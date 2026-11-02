"""QA checks that replace visual inspection: structural validation against the input map,
curve crossings, near misses (dots too close, foreign lines through dots), label and title collisions."""
import json, math
from collections import defaultdict, Counter
from .geometry import samples, segment_polyline, text_box, is_exit_text, label_kind
from .labels import node_label_text

def validate(path_out, path_in):
    """structural checklist; returns list of (check, ok, detail)"""
    h, m = json.load(open(path_out)); N, R, T, C = m['nodes'], m['reactions'], m['text_labels'], m['canvas']
    src = json.load(open(path_in))[1]
    res = []; ok = lambda name, c, d='': res.append((name, bool(c), d))
    segs = [s for r in R.values() for s in r['segments'].values()]
    ok('segment endpoints exist', all(s['from_node_id'] in N and s['to_node_id'] in N for s in segs))
    ok('coordinates finite', all(math.isfinite(n['x']) and math.isfinite(n['y']) for n in N.values()))
    ok('b1/b2 null or {x,y} pairs', all((s['b1'] is None and s['b2'] is None) or
        (isinstance(s['b1'], dict) and isinstance(s['b2'], dict) and all(math.isfinite(s[b][c]) for b in ('b1', 'b2') for c in 'xy')) for s in segs))
    ok('no temporary keys', not [k for D in (N, R, T) for o in D.values() for k in o if k.startswith('_')])
    diff = [k for k in set(R) | set(src['reactions']) if k not in R or k not in src['reactions'] or
            any(R[k].get(f) != src['reactions'][k].get(f) for f in ('bigg_id', 'name', 'metabolites', 'reversibility', 'genes', 'gene_reaction_rule'))]
    ok('same reaction set (ids, bigg_id, stoichiometry, reversibility, genes)', not diff, diff[:10])
    bad = []
    for rid, r in R.items():
        S = [s['from_node_id'] for s in r['segments'].values() if N[s['from_node_id']]['node_type'] == 'metabolite']
        P = [s['to_node_id'] for s in r['segments'].values() if N[s['to_node_id']]['node_type'] == 'metabolite']
        co = {x['bigg_id']: x['coefficient'] for x in r['metabolites']}
        srcr = src['reactions'][rid]; Ns = src['nodes']
        drawn_src = Counter(Ns[s['from_node_id']].get('bigg_id') for s in srcr['segments'].values() if Ns[s['from_node_id']]['node_type'] == 'metabolite') + \
                    Counter(Ns[s['to_node_id']].get('bigg_id') for s in srcr['segments'].values() if Ns[s['to_node_id']]['node_type'] == 'metabolite')
        drawn = Counter(N[n]['bigg_id'] for n in S+P)
        if not (all(co[N[n]['bigg_id']] < 0 for n in S) and all(co[N[n]['bigg_id']] > 0 for n in P)) or drawn != drawn_src:
            bad.append(r['bigg_id'])
    ok('drawn metabolites identical to input, signs match stoichiometry', not bad, bad[:10])
    mets = [n for n in N.values() if n['node_type'] == 'metabolite']
    ok('every metabolite has position + label', all('label_x' in n and 'label_y' in n for n in mets))
    ok('every reaction has a label position', all('label_x' in r and 'label_y' in r for r in R.values()))
    inside = lambda x, y: C['x'] <= x <= C['x']+C['width'] and C['y'] <= y <= C['y']+C['height']
    ok('all nodes inside canvas', all(inside(n['x'], n['y']) for n in N.values()))
    ok('all text labels inside canvas', all(inside(t['x'], t['y']) for t in T.values()))
    return res

def crossings(L, step=12):
    """crossings between drawn (curved) segments of different reactions, ignoring segments that share
    a node and the immediate surroundings of nodes. Returns Counter of reaction-pair -> count."""
    segs = []
    for rid, r in L.R.items():
        for s in r['segments'].values():
            pl = segment_polyline(L.N, s, step); segs.append((rid, {s['from_node_id'], s['to_node_id']}, pl))
    grid = defaultdict(list); G = 60
    for i, (rid, nodes, pl) in enumerate(segs):
        for j in range(len(pl)-1):
            a, b = pl[j], pl[j+1]
            for gx in range(int(min(a[0], b[0])//G), int(max(a[0], b[0])//G)+1):
                for gy in range(int(min(a[1], b[1])//G), int(max(a[1], b[1])//G)+1):
                    grid[(gx, gy)].append((i, j))
    def cr(p, q, r, s):
        d = lambda a, b, c: (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
        return d(p, q, r)*d(p, q, s) < 0 and d(r, s, p)*d(r, s, q) < 0
    nodes_xy = [(n['x'], n['y']) for n in L.N.values()]
    near_node = lambda p: any(abs(p[0]-x) < 30 and abs(p[1]-y) < 30 for x, y in nodes_xy)
    found = set()
    for cell in grid.values():
        for u in range(len(cell)):
            for v in range(u+1, len(cell)):
                (i, a), (j, b) = cell[u], cell[v]
                if i == j or segs[i][0] == segs[j][0] or segs[i][1] & segs[j][1]: continue
                key = (min(i, j), max(i, j))
                if key in found: continue
                p, q = segs[i][2][a], segs[i][2][a+1]; r, s = segs[j][2][b], segs[j][2][b+1]
                if cr(p, q, r, s) and not near_node(p): found.add(key)
    out = Counter()
    for i, j in found:
        out[tuple(sorted((L.R[segs[i][0]]['bigg_id'], L.R[segs[j][0]]['bigg_id'])))] += 1
    return out

def near_misses(L, dot_dot=45, line_dot=28):
    """metabolite dots closer than dot_dot to each other; lines of a foreign reaction passing closer
    than line_dot to a metabolite dot (not counting reactions that use that node)."""
    mets = {k: n for k, n in L.N.items() if n['node_type'] == 'metabolite'}
    users = defaultdict(set)
    for rid, t in L.top.items():
        for n in t['S']+t['P']: users[n].add(rid)
    dd = []; keys = list(mets)
    grid = defaultdict(list)
    for k in keys: grid[(int(mets[k]['x']//100), int(mets[k]['y']//100))].append(k)
    for k in keys:
        gx, gy = int(mets[k]['x']//100), int(mets[k]['y']//100)
        for i in range(gx-1, gx+2):
            for j in range(gy-1, gy+2):
                for o in grid.get((i, j), ()):
                    if o <= k: continue
                    d = math.dist(L.P(k), L.P(o))
                    if d < dot_dot: dd.append((mets[k]['bigg_id'], k, mets[o]['bigg_id'], o, round(d)))
    ld = []
    pts = defaultdict(list)
    for x, y, rid in samples(L, 8): pts[(int(x//50), int(y//50))].append((x, y, rid))
    for k, n in mets.items():
        r0 = 22 if n.get('node_is_primary') else 12
        gx, gy = int(n['x']//50), int(n['y']//50); worst = {}
        for i in range(gx-1, gx+2):
            for j in range(gy-1, gy+2):
                for x, y, rid in pts.get((i, j), ()):
                    if rid in users[k]: continue
                    d = math.hypot(x-n['x'], y-n['y'])-r0
                    if d < line_dot and d < worst.get(rid, 1e9): worst[rid] = d
        for rid, d in worst.items(): ld.append((n['bigg_id'], k, L.R[rid]['bigg_id'], round(d)))
    return dd, ld

def label_report(L, title_ids, pad=4):
    """labels (all kinds) overlapping other labels, nodes or line samples"""
    from .labels import Grid, inter, circ_hit
    boxes = []
    for rid, r in L.R.items(): boxes.append(('rxn', r['bigg_id'], text_box(r['bigg_id'], r['label_x'], r['label_y'], 'rxn')))
    for nid, n in L.N.items():
        if n['node_type'] == 'metabolite':
            kd = 'pri' if n.get('node_is_primary') else 'sec'
            boxes.append((kd, n['bigg_id'], text_box(node_label_text(n), n['label_x'], n['label_y'], kd)))
    for k, t in L.T.items(): boxes.append((label_kind(t), t['text'], text_box(t['text'], t['x'], t['y'], label_kind(t))))
    over = []
    for i in range(len(boxes)):
        for j in range(i+1, len(boxes)):
            a, b = boxes[i][2], boxes[j][2]
            if inter((a[0]+pad, a[1]+pad, a[2]-pad, a[3]-pad), (b[0]+pad, b[1]+pad, b[2]-pad, b[3]-pad)):
                over.append((boxes[i][0], boxes[i][1], boxes[j][0], boxes[j][1]))
    S = samples(L, 12); lines = []
    big = [b for b in boxes if b[0] == 'title']
    for kd, txt, b in big:
        hits = sum(1 for x, y, _ in S if b[0]-15 <= x <= b[2]+15 and b[1]-15 <= y <= b[3]+15)
        nodes = sum(1 for n in L.N.values() if n['node_type'] == 'metabolite' and b[0]-35 <= n['x'] <= b[2]+35 and b[1]-35 <= n['y'] <= b[3]+35)
        if hits or nodes: lines.append((txt, hits, nodes))
    return over, lines

def report(L, out_json, in_json, label_scores=None):
    """full QA -> dict (and markdown text)"""
    val = validate(out_json, in_json)
    cx = crossings(L); dd, ld = near_misses(L)
    titles = [k for k, t in L.T.items() if not is_exit_text(t['text'])]
    over, tl = label_report(L, titles)
    md = ['# QA report', '', '## Structural validation']
    md += [f"- {'PASS' if ok else 'FAIL'}: {n}" + (f' {d}' if d else '') for n, ok, d in val]
    md += ['', f'## Crossings between different reactions (curves sampled): {sum(cx.values())}']
    md += [f'- {a} × {b}: {c}' for (a, b), c in cx.most_common()]
    md += ['', f'## Metabolite dots closer than 45 units: {len(dd)}'] + [f'- {a} ({ka}) – {b} ({kb}): {d}' for a, ka, b, kb, d in dd]
    md += ['', f'## Foreign lines within 28 units of a dot: {len(ld)}'] + [f'- dot {a} ({k}) ← line of {r}: {d}' for a, k, r, d in sorted(ld, key=lambda x: x[3])]
    md += ['', f'## Overlapping label boxes: {len(over)}'] + [f'- {a}:{b} / {c}:{d}' for a, b, c, d in over]
    md += ['', f'## Titles touching lines/nodes: {len(tl)}'] + [f'- {t}: {h} line pts, {n} nodes' for t, h, n in tl]
    if label_scores:
        md += ['', '## Label placement scores (<10 clean, 10-119 touches a line, >=120 overlaps)']
        for k, v in label_scores.items():
            md.append(f"- {k}: n={len(v)} clean={sum(1 for s, _ in v if s < 10)} line={sum(1 for s, _ in v if 10 <= s < 120)} overlap={sum(1 for s, _ in v if s >= 120)}")
    summary = dict(validation_ok=all(ok for _, ok, _ in val), crossings=sum(cx.values()), dots_too_close=len(dd),
                   lines_near_dots=len(ld), label_overlaps=len(over), titles_on_lines=len(tl))
    md += ['', '## Summary', '```', json.dumps(summary, indent=1), '```']
    return summary, '\n'.join(md)
