"""Phase 3 - Bezier smoothing. Pure geometry from the reaction axes recorded by Layout.reaction.
Only cofactor dots move (along their own reaction axis: 'bend'); primary nodes never move."""
import math
from .geometry import V, add, nrm, dot

DEFAULTS = dict(t=0.6,            # curvature: 0 gentle ... 1 strong (0.6 tested best)
                bend=55,          # cofactor shift along the axis to make room for the curve
                bend_stack=15,    # same for reactions stacked on one lane (rows are close)
                stacked=(),       # reaction ids that share one lane (variant families)
                bend_override={}, # rid -> bend
                special_links={}, # (rid, metabolite_node) -> tangent at the metabolite end
                nested_arc_lanes=True)  # outer ring lanes leave the ring node sooner (stay nested)

def apply_bend(L, cfg):
    for rid, m in L.meta.items():
        t = L.top[rid]
        B = cfg['bend_override'].get(rid, cfg['bend_stack'] if rid in cfg['stacked'] else cfg['bend'])
        u = m['u']
        for n in t['S']:
            if not L.is_primary(n): L._put(n, add(L.P(n), u, -B))
        for n in t['P']:
            if not L.is_primary(n): L._put(n, add(L.P(n), u, +B))

def _rule(M, K, u, t):
    """metabolite M <-> multimarker K: leave M towards the axis, enter K tangentially to the axis"""
    nn = (-u[1], u[0]); MK = V(K, M); d = dot(MK, nn); D = dot(MK, u)
    if abs(d) < 25 or abs(D) < 20: return None
    s = 1 if D > 0 else -1
    g_far = add(M, V(M, K), 0.33)
    g_near = add(K, u, s*min(0.30*abs(D), 0.35*abs(d), 60))
    s_far = add(K, u, D+s*0.2*abs(d))
    s_near = add(K, u, s*min(0.45*abs(D), 0.40*abs(d)))
    mix = lambda a, b: ((1-t)*a[0]+t*b[0], (1-t)*a[1]+t*b[1])
    return mix(g_far, s_far), mix(g_near, s_near)

def _set(s, b1, b2): s['b1'] = {'x': b1[0], 'y': b1[1]}; s['b2'] = {'x': b2[0], 'y': b2[1]}

def apply_curves(L, cfg):
    t = cfg['t']
    special = {(L.S.get(r, r), L.S.get(n, n)): v for (r, n), v in cfg['special_links'].items()}
    for rid, r in L.R.items():
        m = L.meta[rid]; u = m['u']
        for s in r['segments'].values():
            f, to = s['from_node_id'], s['to_node_id']; s['b1'] = s['b2'] = None
            if L.N[f]['node_type'] != 'metabolite' and L.N[to]['node_type'] != 'metabolite': continue
            sub = L.N[f]['node_type'] == 'metabolite'
            Mn, Kn = (f, to) if sub else (to, f); M, K = L.P(Mn), L.P(Kn)
            pri = L.is_primary(Mn)
            if (rid, Mn) in special:                                # long link between clusters
                tm = nrm(special[(rid, Mn)]); Lk = math.dist(M, K)
                sK = 1 if dot(V(K, M), u) > 0 else -1
                cK = add(K, u, sK*0.40*Lk); cM = add(M, tm, 0.40*Lk)
                _set(s, cM, cK) if sub else _set(s, cK, cM); continue
            if m['kind'] == 'arc' and Mn in m['ends']:              # ring: true circular arcs
                c, sg = m['c'], m['sg']
                tan = lambda p: (lambda a: (-math.sin(a)*sg, math.cos(a)*sg))(math.atan2(p[1]-c[1], p[0]-c[0]))
                h = (0.33+0.03*t)*math.dist(M, K); tM, tK = tan(M), tan(K); hM = h
                if cfg['nested_arc_lanes']:
                    off = abs(math.dist(K, c)-math.dist(M, c)); hM = h*150/max(150, off)
                _set(s, add(M, tM, hM), add(K, tK, -h)) if sub else _set(s, add(K, tK, h), add(M, tM, -hM))
                continue
            if m['kind'] == 'lane' and pri:                         # parallel lane: S-shaped bow
                if abs(dot(V(K, M), (-u[1], u[0]))) < 25: continue
                k = (0.32+0.10*t)*math.dist(M, K)
                _set(s, add(M, u, k), add(K, u, -k)) if sub else _set(s, add(K, u, k), add(M, u, -k))
                continue
            cp = _rule(M, K, u, t)
            if cp is None: continue
            far, near = cp
            _set(s, far, near) if sub else _set(s, near, far)

def apply(L, cfg=None):
    c = dict(DEFAULTS); c.update(cfg or {})
    c['stacked'] = {L.S.get(x, x) for x in c['stacked']}
    c['bend_override'] = {L.S.get(k, k): v for k, v in c['bend_override'].items()}
    apply_bend(L, c); apply_curves(L, c)
    return c
