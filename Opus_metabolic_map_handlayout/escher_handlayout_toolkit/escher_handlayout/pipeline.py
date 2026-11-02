"""Run a spec end-to-end: splits -> merges -> hand layout -> titles -> curves -> labels -> packing ->
title nudge -> canvas -> JSON -> validation + QA report (+ optional preview PNG)."""
import importlib.util, json, os, sys, time
from collections import Counter
from .engine import Layout
from . import curves, labels, pack, qa, titles as T_
from .inspect_map import exit_owners
from .geometry import is_exit_text

def load_spec(path):
    spec = importlib.util.spec_from_file_location('map_spec', path)
    mod = importlib.util.module_from_spec(spec); sys.modules['map_spec'] = mod; spec.loader.exec_module(mod)
    return mod

def build(spec, input_path):
    L = Layout(input_path)
    L.exit_owner = exit_owners(L); L.exit_owner.update(getattr(spec, 'EXIT_OWNER', None) or {})  # spec overrides auto
    L.split_shared_secondaries()
    if hasattr(spec, 'prepare'): spec.prepare(L)
    for keep, others in getattr(spec, 'MERGES', []):
        L.merge(L.S.get(keep, keep), *[L.S.get(o, o) for o in others])
    bad = [z for z, o in spec.ORIGIN.items() if o is None or None in tuple(o)]
    if bad: raise ValueError(f'ORIGIN not set for zones {bad} - choose zone positions by hand')
    L.origin = dict(spec.ORIGIN); L.flip = dict(getattr(spec, 'FLIP', {}))
    try:
        spec.layout(L); spec.layout(L)      # 2nd pass: every reaction sees final primary positions
    except KeyError as e:
        raise KeyError(f'layout(): unknown zone or id {e} - every zone used in setp/reaction needs an ORIGIN') from None
    missing = [n for n, v in L.N.items() if v['node_type'] == 'metabolite' and n not in L.placed]
    unzoned = [L.R[r]['bigg_id']+'#'+r for r in L.R if r not in L.zone_of]
    if missing or unzoned:
        raise RuntimeError(f'unplaced metabolite nodes {missing[:20]} / reactions without L.reaction call {unzoned[:20]}')
    return L

def place_titles(spec, L):
    """zone titles (reuse existing label ids where given), extra titles; returns {title_id: zone or None},
    {title_id: anchor node}"""
    zone_of_title, anchor_of = {}, {}
    used = set()
    pos = getattr(spec, 'TITLE_POS', {})
    for z, d in spec.ZONES.items():
        if not d.get('title'): continue
        tid = d.get('title_id')
        if tid is None or tid not in L.T: tid = L.newid(); L.T[tid] = {'text': d['title'], 'x': 0, 'y': 0}
        L.T[tid]['text'] = d['title']; used.add(tid); zone_of_title[tid] = z
        if d.get('anchor'):
            n, (dx, dy) = d['anchor']; n = L.S.get(n, n); anchor_of[tid] = n
            L.T[tid]['x'], L.T[tid]['y'] = L.N[n]['x']+dx, L.N[n]['y']+dy
        else:
            if z not in pos: raise KeyError(f'TITLE_POS has no position for zone {z!r}')
            L.T[tid]['x'], L.T[tid]['y'] = pos[z]
    for e in getattr(spec, 'EXTRA_TITLES', []):
        tid = L.newid(); L.T[tid] = {'text': e['text'], 'x': 0, 'y': 0}; used.add(tid); zone_of_title[tid] = e.get('zone')
        if 'anchor' in e:
            n = L.S.get(e['anchor'], e['anchor']); anchor_of[tid] = n; dx, dy = e.get('offset', (0, 0))
            L.T[tid]['x'], L.T[tid]['y'] = L.N[n]['x']+dx, L.N[n]['y']+dy
        else: L.T[tid]['x'], L.T[tid]['y'] = e['pos']
    dropped = [k for k, t in list(L.T.items()) if not is_exit_text(t['text']) and k not in used]
    for k in dropped: del L.T[k]
    return zone_of_title, anchor_of, dropped

def run(spec_path, input_path=None, out_dir='out', render=False, name=None):
    t0 = time.time(); spec = load_spec(spec_path)
    input_path = input_path or os.path.join(os.path.dirname(os.path.abspath(spec_path)), spec.INPUT)
    os.makedirs(out_dir, exist_ok=True)
    name = name or os.path.splitext(os.path.basename(input_path))[0]+'_final'
    log = {}
    L = build(spec, input_path)
    log['splits'] = L.split_log; log['merges'] = L.merge_log
    zone_of_title, anchor_of, dropped = place_titles(spec, L)
    log['dropped_title_labels'] = dropped
    for owner, text in getattr(spec, 'NEW_EXITS', []):
        o = L.S.get(owner, owner); k = L.newid(); L.T[k] = {'text': text, 'x': L.N[o]['x']+30, 'y': L.N[o]['y']+55}; L.exit_owner[k] = o
    cfg = curves.apply(L, getattr(spec, 'CURVES', {}))
    log['curves'] = {k: (sorted(v) if isinstance(v, set) else v) for k, v in cfg.items() if k != 'special_links'}
    hubs = {L.S.get(k, k) for k in getattr(spec, 'HUBS', [m[0] for m in getattr(spec, 'MERGES', [])])}
    title_ids = [k for k, t in L.T.items() if not is_exit_text(t['text'])]
    scores = labels.place_all(L, hubs=hubs, title_ids=title_ids)
    # ---- packing: which text labels travel with which piece
    comps = pack.components(L); comp_of = {r: tuple(sorted(c)) for c in comps for r in c}
    node_comp = {n: comp_of[r] for r, t in L.top.items() for n in t['S']+t['P']}
    main = tuple(sorted(comps[0])); piece_texts = {}
    for k, o in L.exit_owner.items():
        if node_comp[o] != main: piece_texts.setdefault(node_comp[o], []).append(k)
    for tid, z in zone_of_title.items():
        if tid in anchor_of: c = node_comp[anchor_of[tid]]
        elif z is not None:
            c = Counter(comp_of[r] for r, zz in L.zone_of.items() if zz == z).most_common(1)[0][0]
        else: c = main
        if c != main: piece_texts.setdefault(c, []).append(tid)
    moving = {k for v in piece_texts.values() for k in v}
    pk = getattr(spec, 'PACK', {})
    res = pack.pack(L, piece_texts, [k for k in L.T if k not in moving], keep_out=pk.get('keep_out', []),
                    auto_ring_keep_out=pk.get('auto_ring_keep_out', True))
    log['packing'] = res
    log['title_nudges'] = T_.nudge_titles(L)
    log['canvas'] = T_.set_canvas(L, margin=getattr(spec, 'CANVAS_MARGIN', 150))
    out_json = os.path.join(out_dir, name+'.json'); L.save(out_json)
    summary, md = qa.report(L, out_json, input_path, scores); md = md.split('\n')
    a = lambda b: (b[2]-b[0])*(b[3]-b[1])
    b0, b1 = res['before'], res['after']; c = log['canvas']
    md += ['', '## Area', f"- before packing {b0[2]-b0[0]:.0f} x {b0[3]-b0[1]:.0f} = {a(b0)/1e6:.2f} M",
           f"- after packing  {b1[2]-b1[0]:.0f} x {b1[3]-b1[1]:.0f} = {a(b1)/1e6:.2f} M",
           f"- canvas {c['width']:.0f} x {c['height']:.0f} (margin incl.)"]
    md += ['', '## Moves of detached pieces'] + [f"- {m['reactions']}: dx={m['dx']} dy={m['dy']}" for m in res['moves']]
    md += ['', '## Title nudges'] + [f'- {t}: ({dx}, {dy}) {"ok" if ok else "NO CLEAN SPOT"}' for t, dx, dy, ok in log['title_nudges']]
    open(os.path.join(out_dir, 'qa_report.md'), 'w').write('\n'.join(md))
    json.dump(dict(summary=summary, **{k: v for k, v in log.items()}), open(os.path.join(out_dir, 'run_log.json'), 'w'), indent=1, default=str)
    if render:
        from .render import render as R
        colors = {z: d.get('color') for z, d in spec.ZONES.items() if d.get('color')}
        R(L, os.path.join(out_dir, 'overview.png'), scale=getattr(spec, 'RENDER_SCALE', 0.5), hubs=hubs, colors=colors)
    summary['seconds'] = round(time.time()-t0, 1); summary['output'] = out_json
    return L, summary
