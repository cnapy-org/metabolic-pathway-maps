"""Layout engine: load an Escher map, split/merge nodes, place reactions from hand coordinates."""
import json, math, copy
from collections import defaultdict
from .geometry import V, add, nrm, dot

DIRS = {'L': (-1, 0), 'R': (1, 0), 'U': (0, -1), 'D': (0, 1)}

class Layout:
    """Wraps an Escher map ([header, map]).  Coordinates: x right, y DOWN."""
    def __init__(self, path):
        self.path = path
        self.header, self.m = json.load(open(path))
        self.N = self.m['nodes']; self.R = self.m['reactions']; self.T = self.m['text_labels']
        ids = list(self.N)+list(self.R)+list(self.T)+[s for r in self.R.values() for s in r['segments']]
        self._next = max(int(k) for k in ids if str(k).isdigit())+1
        self.zone_of = {}; self.placed = set(); self.zone_nodes = defaultdict(set)
        self.origin = {}; self.flip = {}; self.meta = {}; self.S = {}
        self.exit_owner = {}; self.split_log = []; self.merge_log = []; self.merged_into = {}
        self._topology()

    def newid(self):
        self._next += 1; return str(self._next)

    # ---------------------------------------------------------------- topology
    def _topology(self):
        """per reaction: substrate / product metabolite nodes, midmarker, substrate- and product-side
        multimarker (determined from segment directions, never from ids)"""
        self.top = {}
        for rid, r in self.R.items():
            S, P, mid, mms, mm_sub = [], [], None, set(), None
            for s in r['segments'].values():
                f, t = s['from_node_id'], s['to_node_id']
                for n in (f, t):
                    tp = self.N[n]['node_type']
                    if tp == 'midmarker': mid = n
                    elif tp == 'multimarker': mms.add(n)
                if self.N[f]['node_type'] == 'metabolite': S.append(f); mm_sub = t
                if self.N[t]['node_type'] == 'metabolite': P.append(t)
            if mm_sub is None:                                   # no drawn substrate: pick any
                mm_sub = sorted(mms)[0]
            rest = [x for x in mms if x != mm_sub]
            self.top[rid] = dict(S=S, P=P, mid=mid, mm_sub=mm_sub, mm_prod=rest[0] if rest else mm_sub)

    def uses(self):
        u = defaultdict(set)
        for rid, t in self.top.items():
            for n in t['S']+t['P']: u[n].add(rid)
        return u

    def is_primary(self, n): return bool(self.N[n].get('node_is_primary'))

    def split(self, nid, rid):
        """make a second copy of node `nid` that is used only by reaction `rid`; returns the new id"""
        new = self.newid(); self.N[new] = copy.deepcopy(self.N[nid])
        for s in self.R[rid]['segments'].values():
            if s['from_node_id'] == nid: s['from_node_id'] = new
            if s['to_node_id'] == nid: s['to_node_id'] = new
        self._topology(); self.split_log.append((self.N[nid]['bigg_id'], self.R[rid]['bigg_id'], nid, new))
        return new

    def split_shared_secondaries(self):
        """one secondary (cofactor) node per reaction - avoids cofactor lines spanning the map"""
        for n, rs in sorted(self.uses().items()):
            if self.is_primary(n) or len(rs) < 2: continue
            for rid in sorted(rs, key=int)[1:]: self.split(n, rid)

    def resolve(self, nid, what='node'):
        """node id as written in the spec -> current node id: split names (L.S) are looked up, ids that
        were merged away are followed to the node they were merged into. Unknown ids raise a clear error."""
        n = self.S.get(nid, nid)
        while n in self.merged_into: n = self.merged_into[n]
        if n not in self.N:
            raise KeyError(f'{what}: node {nid!r} does not exist (not a node id of the draft, not a split '
                           f'name in L.S, not merged into an existing node)')
        return n

    def merge(self, keep, *others):
        """re-point all segments of `others` to `keep`, delete `others` (exit labels follow).
        Spec references to a merged-away id keep working through resolve()."""
        keep = self.resolve(keep, 'merge keep')
        for o in others:
            if o in self.merged_into: o = self.resolve(o, 'merge')
            if o not in self.N: raise KeyError(f'merge: node {o} does not exist')
            if o == keep: continue                  # already merged into keep
            if self.N[o]['bigg_id'] != self.N[keep]['bigg_id']:
                raise ValueError(f"merge: {o} is {self.N[o]['bigg_id']}, keep {keep} is {self.N[keep]['bigg_id']}")
            for r in self.R.values():
                for s in r['segments'].values():
                    if s['from_node_id'] == o: s['from_node_id'] = keep
                    if s['to_node_id'] == o: s['to_node_id'] = keep
            for tl, own in list(self.exit_owner.items()):
                if own == o: self.exit_owner[tl] = keep
            del self.N[o]; self.merged_into[o] = keep; self.merge_log.append((self.N[keep]['bigg_id'], keep, o))
        self._topology()

    def drop_duplicate_exits(self):
        """after merges several identical exit labels can belong to one node: keep the one with the
        smallest id, delete the others. Returns [(dropped_id, kept_id, text)]."""
        first, out = {}, []
        for k in sorted(self.exit_owner, key=lambda x: (len(x), x)):
            if k not in self.T: continue
            key = (self.exit_owner[k], self.T[k]['text'])
            if key in first:
                del self.T[k]; out.append((k, first[key], key[1]))
            else: first[key] = k
        for k, _, _ in out: del self.exit_owner[k]
        return out

    def drop_labels(self, ids):
        """delete text labels (exit labels or draft titles) by id"""
        for k in ids:
            if k not in self.T: raise KeyError(f'DROP_LABELS: text label {k!r} does not exist')
            del self.T[k]; self.exit_owner.pop(k, None)

    # ---------------------------------------------------------------- placement
    def gp(self, zone, xy):
        """zone-local -> global coordinates (zone origin + optional vertical mirror)"""
        o = self.origin[zone]; fy = self.flip.get(zone, 1)
        return (o[0]+xy[0], o[1]+fy*xy[1])

    def setp(self, zone, d):
        """set positions of (primary) nodes, coordinates local to `zone`. Ids that no longer exist
        (merged away) are ignored; keys may also be names from L.S (split copies)."""
        for nid, xy in d.items():
            nid = self.S.get(nid, nid)
            if nid not in self.N: continue
            if xy is None or xy[0] is None: raise ValueError(f'setp: node {nid} has no coordinate yet')
            x, y = self.gp(zone, xy); self.N[nid]['x'] = x; self.N[nid]['y'] = y
            self.zone_nodes[zone].add(nid); self.placed.add(nid)

    def P(self, nid): return (self.N[nid]['x'], self.N[nid]['y'])
    def _put(self, nid, xy): self.N[nid]['x'], self.N[nid]['y'] = xy

    def reaction(self, rid, zone, side='R', dist=140, gap=35, step=65, bend=0, shift=(0, 0),
                 axis=None, mid=None, arc=None, ends=None, lab=None, end_cof=False, sides=None, along=0):
        """Place one reaction (markers + its cofactor dots) relative to its already placed primaries.
        side    : 'L','R','U','D' (screen direction), 'out'/'in' (rings) or a vector - cofactor side
        dist    : perpendicular distance of cofactor dots from the axis
        gap     : half distance between the two multimarkers along the axis
        step    : spacing between several cofactors on one side (along the axis)
        bend    : extra shift of cofactors along the axis
        shift   : lane offset vector (parallel reactions between the same two nodes)
        axis/mid: override axis direction / midmarker position (mid is zone-local)
        ends    : (substrate_node, product_node) that define the axis
        arc     : (centre_local, radial_offset) - reaction drawn on a ring (circle through `ends`)
        along   : shift along the arc (map units) - several reactions stacked on one arc lane
        end_cof : cofactors at the free end of the axis when a side has no primary metabolite
        sides   : {cofactor_bigg_id: 'L'|'R'|'U'|'D'} per-cofactor side override
        lab     : provisional label direction (final labels are placed automatically)"""
        if rid not in self.R: raise KeyError(f'reaction {rid} not in map')
        t = self.top[rid]
        ends = tuple(self.resolve(e, f'reaction {rid} ends') for e in ends) if ends else None
        if self.flip.get(zone, 1) < 0:                            # mirrored zone: mirror directions
            sw = {'U': 'D', 'D': 'U'}
            if isinstance(side, str): side = sw.get(side, side)
            if lab: lab = sw.get(lab, lab)
            if axis: axis = (axis[0], -axis[1])
            if isinstance(side, (tuple, list)): side = (side[0], -side[1])
            shift = (shift[0], -shift[1])
            if sides: sides = {k: sw.get(v, v) for k, v in sides.items()}
        self.zone_of[rid] = zone
        pri = lambda l: [n for n in l if self.is_primary(n)]
        sec = lambda l: [n for n in l if not self.is_primary(n)]
        PS, PP = pri(t['S']), pri(t['P'])
        if arc is not None:
            c, off = arc; c = self.gp(zone, c)
            a, b = ends if ends else (PS[0], PP[0])
            A = math.atan2(self.P(a)[1]-c[1], self.P(a)[0]-c[0]); B = math.atan2(self.P(b)[1]-c[1], self.P(b)[0]-c[0])
            d = (B-A+math.pi) % (2*math.pi)-math.pi; sg = 1 if d > 0 else -1
            Rr = math.hypot(*V(c, self.P(a)))+off; am = A+d/2+sg*along/Rr
            M = (c[0]+Rr*math.cos(am), c[1]+Rr*math.sin(am)); u = (-math.sin(am)*sg, math.cos(am)*sg)
            K1 = (c[0]+Rr*math.cos(am-sg*gap/Rr), c[1]+Rr*math.sin(am-sg*gap/Rr))
            K2 = (c[0]+Rr*math.cos(am+sg*gap/Rr), c[1]+Rr*math.sin(am+sg*gap/Rr))
            outv = (math.cos(am), math.sin(am))
            self.meta[rid] = dict(kind='arc', c=c, R=Rr, R0=Rr-off, am=am, sg=sg, gap=gap, ends=(a, b), off=off)
        else:
            mean = lambda l: (sum(self.P(n)[0] for n in l)/len(l), sum(self.P(n)[1] for n in l)/len(l))
            cs = mean(ends[0:1]) if ends else (mean(PS) if PS else None)
            cp = mean(ends[1:2]) if ends else (mean(PP) if PP else None)
            if axis is None and (cs is None or cp is None):
                raise ValueError(f'reaction {rid} ({self.R[rid]["bigg_id"]}): one side has no primary - give axis= and mid=')
            u = nrm(axis) if axis else nrm(V(cs, cp))
            M = self.gp(zone, mid) if mid is not None else ((cs[0]+cp[0])/2, (cs[1]+cp[1])/2)
            M = add(M, shift); K1 = add(M, u, -gap); K2 = add(M, u, gap); outv = None
            self.meta[rid] = dict(kind='lane' if tuple(shift) != (0, 0) else 'line', gap=gap)
        n1 = (-u[1], u[0])
        if isinstance(side, str) and side in DIRS: ref = DIRS[side]
        elif side in ('out', 'in'):
            if outv is None: raise ValueError(f"reaction {rid}: side='{side}' needs arc=")
            ref = outv if side == 'out' else (-outv[0], -outv[1])
        else: ref = side
        n = n1 if dot(n1, ref) >= 0 else (-n1[0], -n1[1])
        self._put(t['mid'], M); self._put(t['mm_sub'], K1); self._put(t['mm_prod'], K2)
        self.zone_nodes[zone].update([t['mid'], t['mm_sub'], t['mm_prod']])
        self.meta[rid].update(u=u, n=n, M=M, K1=K1, K2=K2, zone=zone)
        def place(ns, K, sgn):
            for i, c in enumerate(ns):
                if end_cof and not (PP if sgn > 0 else PS):
                    k = len(ns); p = add(add(K, u, sgn*(dist*0.8)), n, (i-(k-1)/2)*step)
                else:
                    nn = n
                    if sides and self.N[c]['bigg_id'] in sides:
                        r = DIRS[sides[self.N[c]['bigg_id']]]; nn = n if dot(n, r) > 0 else (-n[0], -n[1])
                    p = add(add(K, u, sgn*(bend+i*step)), nn, dist)
                self._put(c, p); self.zone_nodes[zone].add(c); self.placed.add(c)
        place(sec(t['S']), K1, -1); place(sec(t['P']), K2, +1)
        for s in self.R[rid]['segments'].values(): s['b1'] = None; s['b2'] = None
        ld = DIRS[lab] if lab else (-n[0], -n[1])
        self.R[rid]['label_x'] = M[0]+ld[0]*55; self.R[rid]['label_y'] = M[1]+ld[1]*55+8

    # ---------------------------------------------------------------- output
    def save(self, path):
        def clean(o):                       # numpy scalars -> plain floats, coordinates rounded
            if isinstance(o, dict): return {k: (v if k in ('metabolites', 'genes') else clean(v)) for k, v in o.items()}
            if isinstance(o, list): return [clean(v) for v in o]
            if isinstance(o, bool) or o is None or isinstance(o, str): return o
            if isinstance(o, (int, float)) or hasattr(o, 'item'): return round(float(o), 2)
            return o
        json.dump([self.header, clean(self.m)], open(path, 'w'), indent=1, ensure_ascii=False)
