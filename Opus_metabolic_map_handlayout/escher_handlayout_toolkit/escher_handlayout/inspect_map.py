"""Phase 0 - inspect a draft Escher map; write a report and a spec skeleton (no coordinates!)."""
import json, math, re
from collections import defaultdict, Counter
from .engine import Layout
from .geometry import text_box, is_exit_text

def exit_owners(L):
    """owner primary node of every exit label ('-> ...').
    1) nearest primary node by distance to the label box; 2) draft generators put exit labels at a few
    fixed vertical offsets (below / above / same row) - offsets that occur >= 3 times are 'patterns';
    re-assign each label to the node that is near AND matches a pattern (bonus 200 units)."""
    pri = {k: v for k, v in L.N.items() if v['node_type'] == 'metabolite' and v.get('node_is_primary')}
    exits = {k: t for k, t in L.T.items() if is_exit_text(t['text'])}
    def bdist(t, p):
        b = text_box(t['text'], t['x'], t['y'], 'exit'); x, y = pri[p]['x'], pri[p]['y']
        return math.hypot(max(b[0]-x, 0, x-b[2]), max(b[1]-y, 0, y-b[3]))
    first = {k: min(pri, key=lambda p: bdist(t, p)) for k, t in exits.items()}
    dyc = Counter(round((t['y']-pri[first[k]]['y'])/4) for k, t in exits.items())
    patterns = {d for d, c in dyc.items() if c >= 3}
    out = {}
    for k, t in exits.items():
        def score(p):
            return bdist(t, p) - (200 if round((t['y']-pri[p]['y'])/4) in patterns else 0)
        out[k] = min(pri, key=score)
    return out

def primary_components(L):
    par = {r: r for r in L.R}
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    by = defaultdict(list)
    for rid, t in L.top.items():
        for n in t['S']+t['P']:
            if L.is_primary(n): by[n].append(rid)
    for l in by.values():
        for x in l[1:]: par[f(x)] = f(l[0])
    g = defaultdict(list)
    for r in L.R: g[f(r)].append(r)
    return sorted(g.values(), key=len, reverse=True)

def all_components(L):
    from .pack import components
    return components(L)

def _bbox(L, rids):
    ns = [n for r in rids for n in L.top[r]['S']+L.top[r]['P']+[L.top[r]['mid']]]
    xs = [L.N[n]['x'] for n in ns]; ys = [L.N[n]['y'] for n in ns]
    return (min(xs), min(ys), max(xs), max(ys))

def propose_subsystems(L):
    """each primary component -> nearest title (titles normally sit above-left of their cluster)"""
    titles = {k: t for k, t in L.T.items() if not is_exit_text(t['text'])}
    res = []
    for c in primary_components(L):
        b = _bbox(L, c); best = None
        for k, t in titles.items():
            tb = text_box(t['text'], t['x'], t['y'], 'title'); cx, cy = (tb[0]+tb[2])/2, tb[3]
            dx = max(b[0]-cx, 0, cx-b[2]); dy = max(b[1]-cy, 0, cy-b[3])
            below = 400 if cy > b[1]+150 else 0                 # title well inside/below the cluster
            sc = math.hypot(dx, dy)+below
            if best is None or sc < best[0]: best = (sc, k, t['text'])
        res.append(dict(reactions=c, title_id=best[1] if best else None, title=best[2] if best else None,
                        distance=round(best[0]) if best else None, bbox=[round(v) for v in b]))
    return res

def eq(L, rid):
    r = L.R[rid]; f = lambda s: ' + '.join(('' if abs(x['coefficient']) == 1 else f"{abs(x['coefficient']):g} ")+x['bigg_id']
                                           for x in r['metabolites'] if (x['coefficient'] < 0) == s)
    return f"{f(True)} {'<=>' if r['reversibility'] else '->'} {f(False)}"

def analyse(path):
    L = Layout(path); rep = {}
    rep['counts'] = dict(reactions=len(L.R), nodes=len(L.N), node_types=dict(Counter(n['node_type'] for n in L.N.values())),
                         text_labels=len(L.T), canvas=L.m['canvas'],
                         primary=sum(1 for n in L.N.values() if n['node_type'] == 'metabolite' and n.get('node_is_primary')),
                         secondary=sum(1 for n in L.N.values() if n['node_type'] == 'metabolite' and not n.get('node_is_primary')))
    rep['titles'] = {k: dict(text=t['text'], x=round(t['x']), y=round(t['y'])) for k, t in L.T.items() if not is_exit_text(t['text'])}
    rep['subsystem_proposal'] = propose_subsystems(L)
    rep['components_all_nodes'] = [len(c) for c in all_components(L)]
    u = L.uses()
    rep['shared_secondary_nodes'] = [dict(node=n, bigg_id=L.N[n]['bigg_id'], reactions=sorted(L.R[r]['bigg_id'] for r in rs))
                                     for n, rs in u.items() if not L.is_primary(n) and len(rs) > 1]
    bc = Counter(r['bigg_id'] for r in L.R.values())
    rep['duplicate_reactions'] = {b: [k for k, r in L.R.items() if r['bigg_id'] == b] for b, c in bc.items() if c > 1}
    fam = defaultdict(list)
    for rid, t in L.top.items():
        fam[frozenset(L.N[n]['bigg_id'] for n in t['S']+t['P'] if L.is_primary(n))].append(rid)
    rep['families_same_primaries'] = [dict(primaries=sorted(k), reactions=v, bigg=[L.R[r]['bigg_id'] for r in v]) for k, v in fam.items() if len(v) > 1]
    dup = defaultdict(list)
    for k, n in L.N.items():
        if n['node_type'] == 'metabolite' and n.get('node_is_primary'): dup[n['bigg_id']].append(k)
    rep['primary_metabolites_with_several_nodes'] = {b: v for b, v in dup.items() if len(v) > 1}
    far = []
    for rid, t in L.top.items():
        m = L.N[t['mid']]
        for n in t['S']+t['P']:
            if L.is_primary(n):
                d = math.hypot(L.N[n]['x']-m['x'], L.N[n]['y']-m['y'])
                if d > 900: far.append(dict(reaction=L.R[rid]['bigg_id'], rid=rid, metabolite=L.N[n]['bigg_id'], node=n, distance=round(d)))
    rep['far_primaries'] = far
    undrawn = {}
    for rid, t in L.top.items():
        drawn = {L.N[n]['bigg_id'] for n in t['S']+t['P']}
        miss = sorted({x['bigg_id'] for x in L.R[rid]['metabolites']}-drawn)
        if miss: undrawn[L.R[rid]['bigg_id']+'#'+rid] = miss
    rep['metabolites_in_stoichiometry_not_drawn'] = undrawn
    rep['exit_label_owner'] = {k: dict(text=L.T[k]['text'], owner=o, bigg_id=L.N[o]['bigg_id']) for k, o in exit_owners(L).items()}
    rep['reactions'] = {rid: dict(bigg_id=r['bigg_id'], equation=eq(L, rid),
                        substrates=[(n, L.N[n]['bigg_id'], L.is_primary(n)) for n in L.top[rid]['S']],
                        products=[(n, L.N[n]['bigg_id'], L.is_primary(n)) for n in L.top[rid]['P']]) for rid, r in L.R.items()}
    return L, rep

def markdown(rep):
    c = rep['counts']; md = ['# Phase 0 report', '',
        f"- reactions {c['reactions']}, nodes {c['nodes']} {c['node_types']}, text labels {c['text_labels']}",
        f"- metabolite nodes: {c['primary']} primary, {c['secondary']} secondary; canvas {({k: round(v) for k, v in c['canvas'].items()})}",
        '', '## Titles'] + [f"- `{k}` {v['text']!r} at ({v['x']}, {v['y']})" for k, v in rep['titles'].items()]
    md += ['', '## Subsystem proposal (primary-node components -> nearest title). CHECK IT.']
    for p in rep['subsystem_proposal']:
        md.append(f"- **{p['title']}** (`{p['title_id']}`, dist {p['distance']}{', LOW CONFIDENCE' if (p['distance'] or 0) > 300 else ''}) n={len(p['reactions'])}: " +
                  ', '.join(f"{rep['reactions'][r]['bigg_id']}#{r}" for r in p['reactions']))
    md += ['', f"## Components over all shared nodes (draft as is): {rep['components_all_nodes']}",
           '', '## Shared secondary nodes (split automatically, one per reaction)']
    md += [f"- {s['bigg_id']} `{s['node']}`: {s['reactions']}" for s in rep['shared_secondary_nodes']]
    md += ['', '## Duplicate reactions (keep; draw as lanes)'] + [f'- {b}: {v}' for b, v in rep['duplicate_reactions'].items()]
    md += ['', '## Families between the same primary metabolites'] + [f"- {f['primaries']}: {list(zip(f['bigg'], f['reactions']))}" for f in rep['families_same_primaries']]
    md += ['', '## Primary metabolites drawn as several nodes (merge candidates)'] + [f'- {b}: {v}' for b, v in rep['primary_metabolites_with_several_nodes'].items()]
    md += ['', '## Primaries > 900 units from their reaction'] + [f"- {f['reaction']}#{f['rid']}: {f['metabolite']} `{f['node']}` {f['distance']}" for f in rep['far_primaries']]
    md += ['', '## In stoichiometry but not drawn (usually protons; leave as is)'] + [f'- {k}: {v}' for k, v in rep['metabolites_in_stoichiometry_not_drawn'].items()]
    md += ['', '## Exit labels -> owner node (row rule)'] + [f"- `{k}` {v['text']!r} -> {v['bigg_id']} `{v['owner']}`" for k, v in rep['exit_label_owner'].items()]
    md += ['', '## Reactions'] + [f"- `{k}` {v['bigg_id']}: {v['equation']}" for k, v in rep['reactions'].items()]
    return '\n'.join(md)

def skeleton(L, rep, input_name):
    """spec skeleton: zones, node ids with draft coordinates in comments, reaction calls.
    Every coordinate is None on purpose - the layout must be chosen by hand."""
    def code(title, taken):
        w = re.sub('[^A-Za-z ]', ' ', title or 'zone').split()
        base = ''.join(x[0] for x in w).upper() if len(w) > 1 else (w[0][:3].upper() if w else 'Z')
        z, i = base, 1
        while z in taken: i += 1; z = f'{base}{i}'
        return z
    zones = {}; used = set(); by_title = {}
    for p in rep['subsystem_proposal']:
        low = p['distance'] is None or p['distance'] > 300
        key = (p['title_id'], low)
        if key not in by_title or low:
            z = code(p['title'], zones); by_title[key] = z
            zones[z] = dict(title=p['title'], title_id=p['title_id'] if not low else None, reactions=[], low=low)
        zones[by_title[key]]['reactions'] += p['reactions']
    out = [f'"""Spec for {input_name} - generated skeleton. Replace EVERY None by hand-chosen coordinates."""',
           'import math', '', 'INPUT = ' + repr(input_name), '',
           '# zone -> title text (None = no title), existing text-label id to reuse, optional colour / anchor',
           'ZONES = {']
    for z, d in zones.items():
        out.append(f"    {z!r}: dict(title={d['title']!r}, title_id={d['title_id']!r}),"+("   # LOW CONFIDENCE: far from any title - check subsystem" if d['low'] else ''))
    out += ['}', '', 'EXIT_OWNER = {}          # automatic owner rule (see report); add {label_id: node_id} to override', '',
            'def prepare(L):', '    """splits of primary nodes: L.S[name] = L.split(node_id, reaction_id)"""', '    pass', '',
            'ORIGIN = {' + ', '.join(f'{z!r}: (None, None)' for z in zones) + '}', 'FLIP = {}', '',
            'MERGES = [   # (keep, [duplicates])  - key metabolites only', ']', '',
            'def layout(L):']
    for z, d in zones.items():
        out.append(f"    # ---------------------------------------------------------------- {z}: {d['title']}")
        nodes = []
        for r in d['reactions']:
            for n in L.top[r]['S']+L.top[r]['P']:
                if L.is_primary(n) and n not in used: used.add(n); nodes.append(n)
        out.append(f'    L.setp({z!r}, {{')
        for n in nodes:
            out.append(f"        {n!r}: None,   # {L.N[n]['bigg_id']}  (draft {round(L.N[n]['x'])}, {round(L.N[n]['y'])})")
        out.append('    })')
        for r in d['reactions']:
            out.append(f"    L.reaction({r!r}, {z!r}, side='R')   # {L.R[r]['bigg_id']}: {eq(L, r)}")
    out += ['', 'NEW_EXITS = []          # (owner node, text) for copies intentionally not merged',
            'EXTRA_TITLES = []       # dict(text=..., anchor=node_id, offset=(dx, dy)) or pos=(x, y)',
            'TITLE_POS = {}          # zone -> (x, y) global position of the zone title',
            'CURVES = dict(t=0.6, stacked=[], bend_override={}, special_links={})',
            'PACK = dict(keep_out=[], auto_ring_keep_out=True)', '']
    return '\n'.join(out)
