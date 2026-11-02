"""Optional matplotlib preview (NOT Escher): coloured lines per zone, arrowheads, labels with a white
halo (the halo exists only in the PNG)."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Circle, PathPatch
from matplotlib.path import Path
from .geometry import V, add, nrm, is_exit_text
from .labels import node_label_text

DEFAULT_COLORS = ['#1f77b4', '#2ca02c', '#d62728', '#9467bd', '#ff7f0e', '#8c564b', '#e377c2',
                  '#7f7f7f', '#bcbd22', '#17becf', '#555555']

def render(L, out, region=None, scale=0.5, hubs=(), colors=None, dpi=100):
    N, R, T = L.N, L.R, L.T
    zones = sorted(set(L.zone_of.values()))
    colors = dict(colors or {})
    for i, z in enumerate(zones): colors.setdefault(z, DEFAULT_COLORS[i % len(DEFAULT_COLORS)])
    if region is None:
        c = L.m.get('canvas'); region = (c['x'], c['y'], c['x']+c['width'], c['y']+c['height'])
    x0, y0, x1, y1 = region
    fig = plt.figure(figsize=((x1-x0)*scale/dpi, (y1-y0)*scale/dpi), dpi=dpi); ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.axis('off'); ax.set_aspect('equal')
    pt = scale*72/dpi; halo = [pe.withStroke(linewidth=24*pt, foreground='white')]
    for rid, r in R.items():
        col = colors.get(L.zone_of.get(rid), '#333')
        for s in r['segments'].values():
            a, b = N[s['from_node_id']], N[s['to_node_id']]
            pts = [(a['x'], a['y'])]+([(s['b1']['x'], s['b1']['y']), (s['b2']['x'], s['b2']['y'])] if s['b1'] else [])+[(b['x'], b['y'])]
            codes = [Path.MOVETO]+([Path.CURVE4]*3 if s['b1'] else [Path.LINETO])
            ax.add_patch(PathPatch(Path(pts, codes), fill=False, ec=col, lw=10*pt, zorder=1, capstyle='round'))
            heads = []
            if b['node_type'] == 'metabolite': heads.append((pts[-2], pts[-1], b))
            if r['reversibility'] and a['node_type'] == 'metabolite': heads.append((pts[1], pts[0], a))
            for p, q, nd in heads:
                d = nrm(V(p, q)); rr = (22 if nd.get('node_is_primary') else 12)+4
                tip = add(q, d, -rr); base = add(tip, d, -30); nn = (-d[1], d[0])
                ax.add_patch(plt.Polygon([tip, add(base, nn, 14), add(base, nn, -14)], color=col, zorder=2))
        m = N[L.top[rid]['mid']]; ax.add_patch(Circle((m['x'], m['y']), 6, color=col, zorder=3))
        ax.text(r['label_x'], r['label_y'], r['bigg_id'], fontsize=22*pt, color='#222', weight='bold', zorder=6, path_effects=halo)
    for nid, n in N.items():
        if n['node_type'] != 'metabolite': continue
        p = bool(n.get('node_is_primary')); rad = 26 if nid in hubs else (20 if p else 10)
        ax.add_patch(Circle((n['x'], n['y']), rad, fc='#f2b54a' if p else '#e6c27a', ec='#a06a00', lw=2*pt, zorder=4))
        ax.text(n['label_x'], n['label_y'], node_label_text(n), fontsize=(26 if p else 20)*pt, color='#000' if p else '#444', zorder=7, path_effects=halo)
    for t in T.values():
        big = not is_exit_text(t['text'])
        ax.text(t['x'], t['y'], t['text'], fontsize=(44 if big else 22)*pt, color='#5b2a86' if big else '#2e7d32',
                weight='bold' if big else 'normal', style='italic' if big else 'normal', zorder=8, path_effects=halo)
    fig.savefig(out, dpi=dpi); plt.close(fig)
