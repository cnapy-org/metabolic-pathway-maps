# Spec reference and cookbook

A spec is a Python file. `run` imports it and uses the names below. Ids are the node / reaction /
text-label ids of the **draft** map (strings). Names stored in `L.S` (split copies) can be used
wherever an id is expected.

## Module-level names

| name | required | meaning |
|---|---|---|
| `INPUT` | no | draft file name relative to the spec (overridden by `--input`) |
| `ZONES` | yes | `{zone: dict(title=str or None, title_id=existing text label id or None, color='#rrggbb', anchor=(node, (dx, dy)))}`. `title_id` reuses (and renames) a draft title label. `anchor` pins the title relative to a node (it travels with a packed piece). Draft titles not referenced anywhere are deleted. |
| `TITLE_POS` | yes* | `{zone: (x, y)}` global title position (*for zones with a title and no `anchor`) |
| `EXTRA_TITLES` | no | `[dict(text=..., anchor=node, offset=(dx, dy), zone=...)]` or `pos=(x, y)` – e.g. titles for detached pieces |
| `EXIT_OWNER` | no | `{exit_label_id: node_id}` overrides of the automatic owner rule |
| `prepare(L)` | no | splits: `L.S['name'] = L.split(node_id, reaction_id)`; early merges allowed |
| `ORIGIN` | yes | `{zone: (x, y)}` global position of each zone's local (0, 0) |
| `FLIP` | no | `{zone: -1}` mirror the zone vertically (side/axis/lab directions are mirrored too) |
| `MERGES` | no | `[(keep, [duplicates])]` – key metabolites only; keepers are drawn larger (hubs) |
| `HUBS` | no | node ids drawn larger in previews (default: merge keepers) |
| `layout(L)` | yes | all `L.setp` / `L.reaction` calls (called twice) |
| `NEW_EXITS` | no | `[(owner_node, '→ text')]` exit labels for copies not merged |
| `CURVES` | no | `dict(t=0.6, bend=55, bend_stack=15, stacked=[rids], bend_override={rid: bend}, special_links={(rid, node): (tx, ty)}, nested_arc_lanes=True)` |
| `PACK` | no | `dict(mode='pack', keep_out=[('box', (x0, y0, x1, y1)), ('circ', (x, y, r))], auto_ring_keep_out=True)`. `mode='keep'` leaves every detached piece exactly where `layout()` put it (no moves, no grid snapping) – use it when you place all pieces by hand; pieces that overlap are listed in the QA report (`piece_overlaps`). |
| `CANVAS_MARGIN`, `RENDER_SCALE` | no | defaults 150, 0.5 |

## `L.setp(zone, {node_id: (x, y)})`
Zone-local coordinates of primary nodes. Every primary node must be placed once. Secondary nodes are
placed by `L.reaction`. Ids that were merged away are ignored.

## `L.reaction(rid, zone, **options)`
Places midmarker, the two multimarkers and the cofactor dots of one reaction. The axis runs from the
mean of the primary substrates to the mean of the primary products (model direction).

| option | default | meaning |
|---|---|---|
| `side` | `'R'` | cofactor side: `'L' 'R' 'U' 'D'` (screen), `'out'`/`'in'` (rings) or a vector |
| `dist` | 140 | perpendicular distance of cofactor dots |
| `gap` | 35 | half distance between the multimarkers (40 on stacked lanes) |
| `step` | 65 | spacing of several cofactors on one side, along the axis |
| `bend` | 0 | extra along-axis shift (Phase 3 adds its own bend, see `CURVES`) |
| `shift` | (0, 0) | lane offset vector (parallel lane) |
| `axis`, `mid` | auto | axis direction / zone-local midmarker position |
| `ends` | auto | `(substrate_node, product_node)` defining axis or ring arc |
| `arc` | None | `(centre_local, radial_offset)` ring reaction (true circular arc) |
| `along` | 0 | shift along the arc (map units) |
| `end_cof` | False | cofactors at the free end when a side has no primary (ATP fans) |
| `sides` | None | `{cofactor_bigg_id: 'L'/'R'/'U'/'D'}` per-cofactor side |
| `lab` | None | provisional label side (final labels are automatic) |

## Cookbook (all from `examples/central_metabolism/spec.py`)

```python
# straight chain with cofactors on the right; reverse reaction as a thin lane on the left
L.setp('GLY', {'480': (0, 0), '482': (0, 300)})
L.reaction('492', 'GLY', side='R')                               # HEX1
L.reaction('1017', 'GLY', side='L', shift=(-120, 0), dist=110)   # G6PP

# fan: one substrate, two products, vertical axis
L.reaction('537', 'GLY', axis=(0, 1), mid=(-70, 1075))           # FBA

# variant family stacked on one straight lane
for rid, y in (('614', 2620), ('1070', 2760), ('1082', 2900), ('1094', 3040)):
    L.reaction(rid, 'GLY', side='R', axis=(0, 1), mid=(0, y), gap=40)
# ... and declare them in CURVES['stacked']

# duplicate reaction as lens lanes
L.reaction('582', 'GLY', shift=(-70, 0)); L.reaction('590', 'GLY', shift=(70, 0))

# ring: put nodes on a circle, draw ring reactions as arcs
L.setp('TCA', {'135': ring(0), '136': ring(1)})                  # ring(k) = R*(cos, sin)(-90+40k)
L.reaction('146', 'TCA', arc=((0, 0), 0), side='out')            # ACONTa

# ring family: one concentric lane per reaction, cofactors in the gaps
for rid, k in (('209', -1), ('221', 0), ('185', 1), ('197', 2), ('1029', 3)):
    L.reaction(rid, 'TCA', arc=((0, 0), 180*k), side='in' if k < 0 else 'out',
               dist={-1: 80, 3: 110}.get(k, 90))

# hourglass whose product is merged into a far hub: pin the waist
L.reaction('944', 'PPP', axis=(-1, 0), mid=(725, 800))           # TALA

# split instead of a long line
def prepare(L): L.S['succ_I'] = L.split('251', '260')           # succinate from ICL, inside the ring

# ATP-type reactions with a single primary: fan around the hub
L.reaction('753', 'ATP', axis=(-0.87, -0.5), mid=(-290, -170), side='L', end_cof=True, step=90)

# long merge link entering a hub horizontally
CURVES = dict(special_links={('944', '493'): (1, 0)})
```

## Reading the QA report
- **Crossings**: sampled Bezier curves of different reactions that intersect away from nodes.
- **Dots closer than 45**: two metabolite circles nearly touching (labels will collide too).
- **Foreign lines within 28 of a dot**: a line of another reaction grazes a node.
- **Label scores**: < 10 clean, 10–119 touches a line, ≥ 120 overlaps a node/label.
