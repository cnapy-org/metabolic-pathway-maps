"""small 2-D helpers shared by all modules"""
import math

def V(a, b): return (b[0]-a[0], b[1]-a[1])
def add(a, b, s=1.0): return (a[0]+s*b[0], a[1]+s*b[1])
def dot(a, b): return a[0]*b[0]+a[1]*b[1]
def nrm(v):
    l = math.hypot(*v)
    return (v[0]/l, v[1]/l) if l > 1e-9 else (0.0, 1.0)

# text metrics (map units). Escher draws text starting at label_x with baseline label_y.
FONT = {'pri': 26, 'sec': 20, 'rxn': 22, 'exit': 22, 'title': 44}
WIDTH_K = {'pri': 0.62, 'sec': 0.62, 'rxn': 0.64, 'exit': 0.62, 'title': 0.64}  # bold -> 0.64

def is_exit_text(text): return text[:1] in '→←'

def text_box(text, x, y, kind):
    """(x0, y0, x1, y1) of a text whose start is x and baseline is y"""
    fs = FONT[kind]
    return (x, y-fs, x+WIDTH_K[kind]*fs*len(text), y+0.25*fs)

def label_kind(T_entry): return 'exit' if is_exit_text(T_entry['text']) else 'title'

def bezier_points(P, n):
    """P = [p0,p1,p2,p3] (cubic) or [p0,p1] (line); n+1 sample points"""
    out = []
    for i in range(n+1):
        t = i/n
        if len(P) == 2:
            out.append((P[0][0]+t*(P[1][0]-P[0][0]), P[0][1]+t*(P[1][1]-P[0][1])))
        else:
            c = ((1-t)**3, 3*t*(1-t)**2, 3*t*t*(1-t), t**3)
            out.append((sum(ci*p[0] for ci, p in zip(c, P)), sum(ci*p[1] for ci, p in zip(c, P))))
    return out

def segment_polyline(N, s, step=12):
    a, b = N[s['from_node_id']], N[s['to_node_id']]
    if s['b1'] is None: P = [(a['x'], a['y']), (b['x'], b['y'])]
    else: P = [(a['x'], a['y']), (s['b1']['x'], s['b1']['y']), (s['b2']['x'], s['b2']['y']), (b['x'], b['y'])]
    n = max(2, int(math.dist(P[0], P[-1])/step)) if len(P) == 4 else max(1, int(math.dist(P[0], P[-1])/step))
    return bezier_points(P, n)

def samples(L, step=12, rids=None):
    """sampled points (x, y, reaction_id) of every drawn segment"""
    out = []
    for rid, r in L.R.items():
        if rids is not None and rid not in rids: continue
        for s in r['segments'].values():
            for x, y in segment_polyline(L.N, s, step): out.append((x, y, rid))
    return out
