"""title nudging (after curves/packing) and canvas computation"""
import math
from .geometry import samples, text_box, is_exit_text, label_kind, FONT, WIDTH_K
from .labels import node_label_text

def _all_label_boxes(L, skip):
    B = []
    for n in L.N.values():
        if n['node_type'] == 'metabolite':
            kd = 'pri' if n.get('node_is_primary') else 'sec'; B.append(text_box(node_label_text(n), n['label_x'], n['label_y'], kd))
    for r in L.R.values(): B.append(text_box(r['bigg_id'], r['label_x'], r['label_y'], 'rxn'))
    for k, t in L.T.items():
        if k != skip: B.append(text_box(t['text'], t['x'], t['y'], label_kind(t)))
    return B

def nudge_titles(L, max_r=1000):
    """move every title that touches a line, a node or a label to the nearest clean spot (spiral search)"""
    S = samples(L, 12); moved = []
    mets = [(n['x'], n['y']) for n in L.N.values() if n['node_type'] == 'metabolite']
    for k, t in L.T.items():
        if is_exit_text(t['text']): continue
        others = _all_label_boxes(L, k)
        def bad():
            b = text_box(t['text'], t['x'], t['y'], 'title')
            if any(b[0]-15 <= x <= b[2]+15 and b[1]-15 <= y <= b[3]+15 for x, y, _ in S): return True
            if any(b[0]-35 <= x <= b[2]+35 and b[1]-35 <= y <= b[3]+35 for x, y in mets): return True
            return any(not (b[2]+10 < c[0] or c[2] < b[0]-10 or b[3]+10 < c[1] or c[3] < b[1]-10) for c in others)
        if not bad(): continue
        x0, y0 = t['x'], t['y']; done = False
        for r in range(40, max_r, 40):
            for a in range(0, 360, 15):
                t['x'] = x0+r*math.cos(math.radians(a)); t['y'] = y0+r*math.sin(math.radians(a))
                if not bad(): done = True; break
            if done: break
        if not done: t['x'], t['y'] = x0, y0
        moved.append((t['text'], round(t['x']-x0), round(t['y']-y0), done))
    return moved

def set_canvas(L, margin=150):
    xs, ys = [], []
    for n in L.N.values():
        xs += [n['x']-30, n['x']+30]; ys += [n['y']-30, n['y']+30]
        if 'label_x' in n:
            kd = 'pri' if n.get('node_is_primary') else 'sec'
            b = text_box(node_label_text(n), n['label_x'], n['label_y'], kd); xs += [b[0], b[2]]; ys += [b[1], b[3]]
    for r in L.R.values():
        b = text_box(r['bigg_id'], r['label_x'], r['label_y'], 'rxn'); xs += [b[0], b[2]]; ys += [b[1], b[3]]
    for t in L.T.values():
        b = text_box(t['text'], t['x'], t['y'], label_kind(t)); xs += [b[0], b[2]]; ys += [b[1], b[3]]
    L.m['canvas'] = {'x': min(xs)-margin, 'y': min(ys)-margin, 'width': max(xs)-min(xs)+2*margin, 'height': max(ys)-min(ys)+2*margin}
    return L.m['canvas']
