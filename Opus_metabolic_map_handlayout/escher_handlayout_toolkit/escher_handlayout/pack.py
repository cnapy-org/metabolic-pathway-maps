"""Phase 5 - pack disconnected pieces around the fixed main component.
Occupancy grids + FFT cross-correlation; only rigid translations of whole pieces."""
import math, numpy as np
from collections import defaultdict
from scipy.signal import fftconvolve
from .geometry import samples, text_box, label_kind
from .labels import node_label_text

CELL = 40; DILATE = 2

def components(L):
    """union-find over reactions sharing any node; largest first"""
    par = {r: r for r in L.R}
    def f(x):
        while par[x] != x: par[x] = par[par[x]]; x = par[x]
        return x
    by = defaultdict(list)
    for rid, r in L.R.items():
        for s in r['segments'].values():
            for n in (s['from_node_id'], s['to_node_id']): by[n].append(rid)
    for l in by.values():
        for x in l[1:]: par[f(x)] = f(l[0])
    g = defaultdict(list)
    for r in L.R: g[f(r)].append(r)
    return sorted(g.values(), key=len, reverse=True)

def piece_nodes(L, rids):
    return {n for rid in rids for s in L.R[rid]['segments'].values() for n in (s['from_node_id'], s['to_node_id'])}

def elements(L, rids, text_ids):
    rs = set(rids)
    pts = [(x, y) for x, y, rid in samples(L, 12, rs)]
    circ, boxes = [], []
    for n in piece_nodes(L, rids):
        v = L.N[n]
        if v['node_type'] == 'metabolite':
            p = bool(v.get('node_is_primary')); circ.append((v['x'], v['y'], 26 if p else 12))
            boxes.append(text_box(node_label_text(v), v['label_x'], v['label_y'], 'pri' if p else 'sec'))
    for rid in rids:
        r = L.R[rid]; boxes.append(text_box(r['bigg_id'], r['label_x'], r['label_y'], 'rxn'))
    for k in text_ids:
        t = L.T[k]; boxes.append(text_box(t['text'], t['x'], t['y'], label_kind(t)))
    return [pts, circ, boxes]

def bbox_of(el):
    pts, circ, boxes = el
    xs = [p[0] for p in pts]+[c[0]-c[2] for c in circ]+[c[0]+c[2] for c in circ]+[b[0] for b in boxes]+[b[2] for b in boxes]
    ys = [p[1] for p in pts]+[c[1]-c[2] for c in circ]+[c[1]+c[2] for c in circ]+[b[1] for b in boxes]+[b[3] for b in boxes]
    return (min(xs), min(ys), max(xs), max(ys))

def raster(el, x0, y0, W, H):
    G = np.zeros((H, W), bool)
    def mark(a, b, c, d):
        i0 = max(0, int((a-x0)//CELL)); i1 = min(W-1, int((c-x0)//CELL))
        j0 = max(0, int((b-y0)//CELL)); j1 = min(H-1, int((d-y0)//CELL))
        if i0 <= i1 and j0 <= j1: G[j0:j1+1, i0:i1+1] = True
    pts, circ, boxes = el
    for x, y in pts: mark(x, y, x, y)
    for c in circ:                         # node dots: bounding square (extra clearance);
        x, y, r = c[:3]                    # keep-out discs (tagged 'disc'): exact disc
        if len(c) > 3: mark_disc(G, x0, y0, x, y, r)
        else: mark(x-r, y-r, x+r, y+r)
    for b in boxes: mark(*b)
    return G

def mark_disc(G, x0, y0, x, y, r):
    """mark the cells that a disc touches (not its bounding square): a cell is occupied when the point
    of the cell rectangle nearest to the centre lies inside the disc"""
    H, W = G.shape
    i0 = max(0, int((x-r-x0)//CELL)); i1 = min(W-1, int((x+r-x0)//CELL))
    j0 = max(0, int((y-r-y0)//CELL)); j1 = min(H-1, int((y+r-y0)//CELL))
    if i0 > i1 or j0 > j1: return
    cx0 = x0+np.arange(i0, i1+1)*CELL; cy0 = y0+np.arange(j0, j1+1)*CELL
    dx = np.clip(x, cx0, cx0+CELL)-x; dy = np.clip(y, cy0, cy0+CELL)-y
    G[j0:j1+1, i0:i1+1] |= dy[:, None]**2+dx[None, :]**2 < r*r

def dilate(G, k):
    out = G.copy()
    for dy in range(-k, k+1):
        for dx in range(-k, k+1): out |= np.roll(np.roll(G, dy, 0), dx, 1)
    return out

def translate(L, rids, text_ids, dx, dy):
    for n in piece_nodes(L, rids):
        v = L.N[n]; v['x'] += dx; v['y'] += dy
        if 'label_x' in v: v['label_x'] += dx; v['label_y'] += dy
    for rid in rids:
        r = L.R[rid]; r['label_x'] += dx; r['label_y'] += dy
        for s in r['segments'].values():
            for b in ('b1', 'b2'):
                if s[b]: s[b]['x'] += dx; s[b]['y'] += dy
    for k in text_ids: L.T[k]['x'] += dx; L.T[k]['y'] += dy

def keep_overlaps(L, comps, piece_texts, main_texts):
    """mode='keep': pairs of pieces whose occupied grid cells overlap (the spec placed them on top
    of each other) -> list of (bigg ids of piece a, bigg ids of piece b), reported in the QA report"""
    els = [elements(L, p, main_texts if i == 0 else piece_texts.get(tuple(sorted(p)), [])) for i, p in enumerate(comps)]
    bbs = [bbox_of(e) for e in els]; out = []
    for i in range(len(comps)):
        for j in range(i+1, len(comps)):
            a, b = bbs[i], bbs[j]
            if a[2] < b[0] or b[2] < a[0] or a[3] < b[1] or b[3] < a[1]: continue
            x0, y0 = min(a[0], b[0]), min(a[1], b[1])
            W = int((max(a[2], b[2])-x0)//CELL)+1; H = int((max(a[3], b[3])-y0)//CELL)+1
            if (raster(els[i], x0, y0, W, H) & raster(els[j], x0, y0, W, H)).any():
                name = lambda p: sorted(L.R[r]['bigg_id'] for r in p)[:3]
                out.append((name(comps[i]), name(comps[j])))
    return out

def ring_keep_out(L, rids):
    """interior discs of all rings (arc reactions) of the main component"""
    seen = {}
    for rid in rids:
        m = L.meta.get(rid)
        if m and m['kind'] == 'arc':
            key = (round(m['c'][0]), round(m['c'][1])); seen[key] = max(seen.get(key, 0), m['R0'])
    return [('circ', (x, y, r)) for (x, y), r in seen.items()]

def pack(L, piece_texts, main_texts, keep_out=(), auto_ring_keep_out=True, mode='pack'):
    """piece_texts: {tuple(sorted(rids)): [text label ids that move with the piece]}
    mode='pack' moves detached pieces to the free spot giving the smallest overall rectangle;
    mode='keep' leaves every piece at its hand-specified position (no moves, no grid snapping)."""
    if mode not in ('pack', 'keep'): raise ValueError(f"PACK mode must be 'pack' or 'keep', not {mode!r}")
    comps = components(L); main = comps[0]; pieces = comps[1:]; log = []
    if mode == 'keep':
        bs = [bbox_of(elements(L, p, main_texts if i == 0 else piece_texts.get(tuple(sorted(p)), [])))
              for i, p in enumerate(comps)]
        u = (min(b[0] for b in bs), min(b[1] for b in bs), max(b[2] for b in bs), max(b[3] for b in bs))
        return dict(before=u, after=u, moves=[], overlaps=keep_overlaps(L, comps, piece_texts, main_texts))
    occ = elements(L, main, main_texts)
    ko = list(keep_out)+(ring_keep_out(L, main) if auto_ring_keep_out else [])
    for kind, data in ko:
        if kind == 'circ': occ[1].append(tuple(data[:3])+('disc',))
        else: occ[2].append(data)
    placed = [occ]
    def union(els):
        bs = [bbox_of(e) for e in els]
        return (min(b[0] for b in bs), min(b[1] for b in bs), max(b[2] for b in bs), max(b[3] for b in bs))
    els = {tuple(sorted(p)): elements(L, p, piece_texts.get(tuple(sorted(p)), [])) for p in pieces}
    before = union([elements(L, main, main_texts)]+list(els.values()))
    area = lambda b: (b[2]-b[0])*(b[3]-b[1])
    pieces.sort(key=lambda p: -area(bbox_of(els[tuple(sorted(p))])))
    for p in pieces:
        key = tuple(sorted(p)); txt = piece_texts.get(key, [])
        el = elements(L, p, txt); pb = bbox_of(el); pw, ph = pb[2]-pb[0], pb[3]-pb[1]
        cur = union(placed)
        x0 = cur[0]-pw-200; y0 = cur[1]-ph-200; x1 = cur[2]+pw+200; y1 = cur[3]+ph+200
        W = int((x1-x0)//CELL)+1; H = int((y1-y0)//CELL)+1
        M = np.zeros((H, W), bool)
        for e in placed: M |= raster(e, x0, y0, W, H)
        M = dilate(M, DILATE)
        F = raster(el, pb[0], pb[1], int(pw//CELL)+2, int(ph//CELL)+2)
        corr = fftconvolve(M.astype(float), F[::-1, ::-1].astype(float), mode='valid')
        best = None
        for j, i in np.argwhere(corr < 0.5):
            nx = float(x0+i*CELL); ny = float(y0+j*CELL)
            ub = (min(cur[0], nx), min(cur[1], ny), max(cur[2], nx+pw), max(cur[3], ny+ph))
            cand = (round(area(ub)/area(cur), 3), math.hypot(nx-pb[0], ny-pb[1]), nx, ny)  # ~equal area -> closest
            if best is None or cand < best: best = cand
        _, _, nx, ny = best; dx, dy = nx-pb[0], ny-pb[1]
        translate(L, p, txt, dx, dy); placed.append(elements(L, p, txt))
        log.append(dict(reactions=sorted(L.R[r]['bigg_id'] for r in p), dx=round(dx), dy=round(dy)))
    return dict(before=before, after=union(placed), moves=log)
