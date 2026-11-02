"""Phase 4 - automatic label placement (candidate search + scoring, no layout algorithm).
Order: reaction labels -> secondary metabolites -> primary metabolites -> exit labels."""
from collections import defaultdict
from .geometry import samples, text_box, FONT, WIDTH_K

DIRS = [(1, 0), (1, -1), (0, -1), (-1, -1), (-1, 0), (-1, 1), (0, 1), (1, 1)]
CELL = 50

class Grid:
    def __init__(s): s.g = defaultdict(list)
    def add(s, x, y, item): s.g[(int(x//CELL), int(y//CELL))].append((x, y, item))
    def query(s, b):
        for i in range(int(b[0]//CELL)-1, int(b[2]//CELL)+2):
            for j in range(int(b[1]//CELL)-1, int(b[3]//CELL)+2):
                yield from s.g.get((i, j), ())

def box_at(cx, cy, d, r, w, h):
    dx, dy = d
    x0 = cx+r if dx > 0 else (cx-r-w if dx < 0 else cx-w/2)
    y0 = cy+r if dy > 0 else (cy-r-h if dy < 0 else cy-h/2)
    if dx and dy: x0 += -dx*r*0.3; y0 += -dy*r*0.3
    return (x0, y0, x0+w, y0+h)

def inter(a, b): return not (a[2] <= b[0] or b[2] <= a[0] or a[3] <= b[1] or b[3] <= a[1])
def circ_hit(b, x, y, r):
    qx = min(max(x, b[0]), b[2]); qy = min(max(y, b[1]), b[3]); return (qx-x)**2+(qy-y)**2 < r*r

class Placer:
    def __init__(s, L, hubs=()):
        s.L = L; s.lines = Grid(); s.nodes = Grid(); s.bgrid = Grid(); s.hubs = set(hubs)
        for x, y, rid in samples(L, step=12): s.lines.add(x, y, rid)
        for nid, n in L.N.items():
            if n['node_type'] == 'metabolite':
                r = (26 if nid in s.hubs else 20) if n.get('node_is_primary') else 10
            elif n['node_type'] == 'midmarker': r = 7
            else: continue
            s.nodes.add(n['x'], n['y'], (nid, r+4))
    def block(s, b, tag): s.bgrid.add((b[0]+b[2])/2, (b[1]+b[3])/2, (b, tag))
    def score(s, b, own=None, pad=4):
        bb = (b[0]-pad, b[1]-pad, b[2]+pad, b[3]+pad)
        ln = sum(1 for x, y, _ in s.lines.query(bb) if bb[0] <= x <= bb[2] and bb[1] <= y <= bb[3])
        nd = sum(1 for x, y, (nid, r) in s.nodes.query(bb) if nid != own and circ_hit(bb, x, y, r))
        lb = sum(1 for x, y, (c, _) in s.bgrid.query((bb[0]-500, bb[1]-80, bb[2]+500, bb[3]+80)) if inter(bb, c))
        return ln*10+nd*120+lb*250
    def best(s, cx, cy, w, h, radii, own=None, prefer=None, dpen=0.08):
        best = None
        for r in radii:
            for d in DIRS:
                b = box_at(cx, cy, d, r, w, h)
                sc = s.score(b, own)+dpen*r+(0 if prefer is None or d in prefer else 6)
                if best is None or sc < best[0]: best = (sc, b)
        return best

def place_all(L, hubs=(), title_ids=()):
    """returns {'rxn':[scores], 'sec':..., 'pri':..., 'exit':...}; score < 10 = clean,
    10..119 = touches line samples, >= 120 = overlaps a node or another label"""
    P = Placer(L, hubs); report = defaultdict(list)
    for k in title_ids: P.block(text_box(L.T[k]['text'], L.T[k]['x'], L.T[k]['y'], 'title'), 'title')
    size = lambda text, kind: (WIDTH_K[kind]*FONT[kind]*len(text), FONT[kind]*1.05, FONT[kind])
    for rid, r in L.R.items():
        m = L.N[L.top[rid]['mid']]; w, h, fs = size(r['bigg_id'], 'rxn')
        sc, b = P.best(m['x'], m['y'], w, h, (30, 46, 66, 92, 125))
        r['label_x'] = b[0]; r['label_y'] = b[3]-0.25*fs; P.block(b, 'rxn'); report['rxn'].append((sc, rid))
    for kind in ('sec', 'pri'):
        for nid, n in L.N.items():
            if n['node_type'] != 'metabolite' or (kind == 'pri') != bool(n.get('node_is_primary')): continue
            text = n['bigg_id'][:-2] if n['bigg_id'][-2:-1] == '_' else n['bigg_id']
            w, h, fs = size(text, kind)
            radii = (14, 22, 32, 46) if kind == 'sec' else (26, 36, 50, 70)
            sc, b = P.best(n['x'], n['y'], w, h, radii, own=nid, prefer=[(1, -1), (1, 0), (-1, -1), (0, -1)])
            n['label_x'] = b[0]; n['label_y'] = b[3]-0.25*fs; P.block(b, kind); report[kind].append((sc, nid))
    for k, o in L.exit_owner.items():
        t = L.T[k]; n = L.N[o]; w, h, fs = size(t['text'], 'exit')
        sc, b = P.best(n['x'], n['y'], w, h, (40, 60, 85, 115, 150), own=o, dpen=0.05)
        t['x'] = b[0]; t['y'] = b[3]-0.25*fs; P.block(b, 'exit'); report['exit'].append((sc, k))
    return report

def node_label_text(n):
    b = n['bigg_id']; return b[:-2] if b[-2:-1] == '_' else b
