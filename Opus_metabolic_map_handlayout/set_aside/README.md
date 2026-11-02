# set_aside – hand layout

This folder holds the hand-laid-out version of the iML1515 "set aside" draft map. The draft groups
45 reactions from BRENDA pathways that are not drawn in the main split maps. The layout was made with
`escher_handlayout` (`../escher_handlayout_toolkit`), following its `PROMPT.md`.

- Input: `ecmapper/BRENDAmap/draft_map_iML1515/split/set_aside.escher.json`
- Spec (all primary coordinates chosen by hand): `spec.py`
- Final run: `v1/`
  - `set_aside_final.json` is the Escher map.
  - `qa_report.md` and `run_log.json` hold the QA results and the run log.
  - `overview.png` is a matplotlib preview, not an Escher rendering.
- Phase-0 report: `inspect/phase0_report.md`

Reproduce from the repository root:

```bash
PYTHONPATH=Opus_metabolic_map_handlayout/escher_handlayout_toolkit python -m escher_handlayout run Opus_metabolic_map_handlayout/set_aside/spec.py --out Opus_metabolic_map_handlayout/set_aside/v1 --name set_aside_final --render
```

## Counts

| | draft | final |
|---|---|---|
| reactions | 45 | 45 (unchanged) |
| nodes | 269 | 257 |
| metabolite nodes | 134 | 122 |
| text labels | 8 | 10 |

The metabolite count changes as follows:
- 15 primary copies were removed by the merges below.
- 3 nodes were added: the cofactors nad, co2 and nadh were shared by MALDDH and ME1 in the draft and were split, one copy per reaction.

The text labels change as follows:
- 3 piece titles were added.
- 1 duplicate exit label was dropped.

## Subsystems (zones) and their sources

All zone titles come from the draft title labels.

| zone | title | reactions |
|---|---|---|
| MC | methylaspartate cycle (draft title 585) | CS, ACONTa/b, ICDHyr, SUCOAS, MMM, FUM, MDH, SUCDi, FRD2, FRD3 |
| C4 | C4 and CAM-carbon fixation (581) | PPC ×2, PPCK, PYK, MDH ×2, ME1, ME2 ×3, MALDDH, ALATA_L ×2, ASPTA ×2 |
| ACB | acetyl CoA biosynthesis (582) | PDH |
| ATD | aerobic toluene degradation (583) | OP4ENH, HOPNTAL, ACALD, 3OXCOAT |
| BF | butanoate fermentation (584) | ACACT1r (title anchored to aacoa) |
| PHO | photosynthesis (586) | RPI, RPE, RBP4E, X5PL3E, TKT1, TKT2, FBP, FBA, TPI, FBA3 |
| PGK, BF2, BF3 | short piece titles (EXTRA_TITLES) | PGK ('photosynthesis'), ECOAH1 and BUTCT ('butanoate fermentation') |

## Layout

**Methylaspartate cycle**
- Drawn as a TCA-style ring: R = 800, 9 positions every 40°, citrate on top, clockwise.
- The ring is open between 2-oxoglutarate and succinyl-CoA, because the map has no AKGDH.
- Methylmalonyl-CoA (MMM) points outwards into that gap.
- The fumarate ↔ succinate family uses one concentric lane per reaction (FRD3 inner, SUCDi on the ring, FRD2 outer) with the cofactors in the gaps.

**Acetyl-CoA hub (above the ring, feeding CS)**
- The toluene-degradation chain enters it from above.
- 3OXCOAT enters from the right.
- PDH enters from pyruvate on the left.
- ACACT1r leaves it towards acetoacetyl-CoA.

**C4/CAM cycle**
- Drawn as a 1400 × 1400 square to the left of the hub: PEP (top left) → OAA (bottom left) → L-malate (bottom right) → pyruvate (top right) → PEP.
- Pyruvate is the corner shared with PDH.
- The square's edges are drawn as follows:
  - The two PPC copies are stacked on the left edge.
  - PPCK runs as a reverse lane inside the square.
  - The two MDH copies are stacked on the bottom edge.
  - The ME2 copies and ME1 are stacked on the right edge (four rows, 300 apart), with their cofactors inside the square.
  - PYK runs along the top edge.
- The side branches attach at the corners:
  - alanine (two ALATA_L lanes) and D-malate (MALDDH) at pyruvate;
  - aspartate (two ASPTA lanes) at OAA.

**Detached pieces**
- The Calvin-cycle fragment, PGK, ECOAH1 and BUTCT are not connected to the main piece.
- They were laid out by hand and then placed by the automatic packing step.

## Merges

**Merged** (only key metabolites):

| metabolite | kept node | merged copies |
|---|---|---|
| acetyl-CoA | 301 | 183, 223, 249 |
| pyruvate | 32 | 5, 103, 536, 182 |
| PEP | 63 | 45, 480 |
| oxaloacetate | 60 | 46, 131, 481 |
| L-malate | 115 | 1, 87, 532 |

**Not merged:**
- **Methylaspartate cycle:** its own oxaloacetate and malate (nodes 290 and 287) stay separate from the C4 square, to avoid long lines between the ring and the square.
- **HOPNTAL:** its pyruvate is a cofactor (secondary) node and stays one.
- **Cofactors:** none were merged.

## Parameters and size

- Curvature t = 0.6. The stacked lanes use a bend of 15; all other reactions use the default bend of 55.
- Area before packing: 6933 × 4138 = 28.7 M.
- Area after packing: 4303 × 4360 = 18.8 M.
- Canvas: 4603 × 4682 including the 150 margin.
- Keep-out zones: the interior of the C4 square (explicit `PACK` box) and the ring interior (automatic).

## Decisions

1. **One C4/CAM square.** The draft spread the C4 reactions over three pieces, and BRENDA files some of the PPC/ME2 copies under "photosynthesis". All PPC, MDH and ME copies are drawn on one square, using stacked lanes for the duplicates. Their zone membership in the map follows the C4 title.
2. **Acetyl-CoA hub.** The four acetyl-CoA copies are merged into one hub, which connects PDH, the toluene chain, ACACT1r and CS into one main piece. The draft had two identical exit labels '→ Central', one on each of two of these copies. After the merge, one of them (label 579) is removed in `prepare()`; this only touches annotation.
3. **PDH substrate.** The PDH pyruvate is merged into the C4 square's pyruvate, so pyruvate sits between the C4 cycle and acetyl-CoA.
4. **Open ring.** The methylaspartate cycle is drawn as an open ring with the textbook TCA order, because AKGDH is absent.
5. **No biochemistry changed.** Reactions, IDs, stoichiometry, reversibility and genes are unchanged, and all "backwards" directions are kept: SUCOAS succ → succoa, PYK pep → pyr, MDH mal → oaa, PPCK oaa → pep.

## Known weak spots

The final QA report shows 0 for every check:
- validation: all PASS
- crossings: 0
- dots closer than 45: 0
- foreign lines within 28 of a dot: 0
- overlapping labels: 0
- titles on lines: 0

Remaining weak spots:
- **Label scores:** all labels are clean (45 reaction, 82 secondary, 40 primary, 1 exit).
- **Identical titles:** there are three 'butanoate fermentation' titles (ACACT1r, ECOAH1, BUTCT) and two 'photosynthesis' titles (the Calvin fragment and PGK), because those pieces are not connected.
- **Empty interiors:** the C4 square and the ring are deliberately left empty apart from their titles, which makes the canvas larger than strictly necessary.
- **Text widths are estimates:** the QA heuristics use 0.62 × font size (0.64 for bold).

## Testing

The map was validated structurally by the toolkit and its QA checks; it was not loaded in Escher. The
preview `v1/overview.png` was inspected visually.
