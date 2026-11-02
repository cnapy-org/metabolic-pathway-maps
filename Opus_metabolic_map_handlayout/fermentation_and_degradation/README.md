# Fermentation and degradation (iML1515, BRENDA draft): hand layout

Input: `ecmapper/BRENDAmap/draft_map_iML1515/split/fermentation_and_degradation.escher.json`.
The map was laid out with `escher_handlayout_toolkit`: every primary-metabolite coordinate is chosen by hand in `spec.py`.

| file | content |
|---|---|
| `spec.py` | the hand layout, with comments per zone |
| `v3/fermentation_and_degradation_final.json` | the Escher map |
| `v3/qa_report.md`, `v3/run_log.json` | validation, QA, packing moves |
| `v3/overview.png` | matplotlib preview (not Escher) |
| `inspect/` | Phase 0 report of the draft |

Reproduce from this folder:

```bash
PYTHONPATH=../escher_handlayout_toolkit python -m escher_handlayout run spec.py --out v3 --render --name fermentation_and_degradation_final
```

## Counts
| | draft | final |
|---|---|---|
| reactions | 49 | 49 (identical) |
| nodes | 302 (155 metabolite) | 300 (153 metabolite) |
| text labels | 23 (9 titles, 14 exits) | 24 (10 titles, 14 exits) |

The draft has no shared cofactor nodes, so no cofactor copies were added. There are 3 merges (−3 nodes) and 1 split (+1 node).

## Subsystems (zones)
All titles come from the draft labels except one, noted in the table.

| zone | title | reactions |
|---|---|---|
| PF | propionate fermentation | MCITS, MCITD, MICITDr, MCITL2, SUCDi, FRD2, FRD3, MMM, MMCD ×2, PPCSCT, PTA2 |
| PD | 3-phenylpropionate degradation | PPPNDO, DHPPD, 3HPPPNH, HPPPNDO, HKNDDH, OP4ENH, HOPNTAL, ACALD, HKNTDH |
| PDA | phenylacetate degradation (aerobic) | PACCOAL, PACOAT, PACCOAE, REPHACCOAI, OXCOAHDH, OXDHCOAT, DHACOAH, HADPCOADH3, 3OXCOAT |
| AD | allantoin degradation | URIC, ALLTN, ALLTAMH2, UGCIAMH, UGLYCH, URDGLYCD, OXAMTC |
| MD | methylglyoxal degradation | GLYOX3, LGTHL, GLYOX, LDH_D2 |
| LLD | L-lactaldehyde degradation | LCARS, LCARSyi, LCADi |
| OO | octane oxidation | FACOAL80t2pp, FACOAE80 |
| HBZ | **3-hydroxybenzoyl-CoA hydrolysis** (new title) | 3HBZCT |
| CD | creatinine degradation | SARCOX |
| PRD | propanol degradation | ACACCT |

## Merged / not merged
**Merged:**
- **Allantoin:** URIC's product node 160, drawn as a secondary in the draft, is merged into the allantoin chain (169). This connects urate → allantoin.
- **(S)-methylmalonyl-CoA:** 547 is merged into 468, so both MMCD duplicates start from one node.
- **Propionyl-CoA:** 417 (PTA2) is merged into the propionate-fermentation hub 477.

**Split:** propionyl-CoA gets a second copy (`ppcoa_P`) used only by PPCSCT, at the bottom of the propionate block. This avoids a long line from the top hub.

**Not merged:**
- succinyl-CoA in phenylacetate degradation (405) vs. propionate fermentation (467)
- fumarate in 3-phenylpropionate degradation (587) vs. propionate fermentation (489)

They belong to independent pathways, and each keeps its draft exit label.

## Layout
- **Propionate fermentation:** propionyl-CoA hub at the top left.
  - The methylcitrate pathway runs right to 2-methylcitrate and then down: 2mcacn → micit → succinate. Oxaloacetate enters MCITS from above.
  - The sleeping-beauty loop forms the left column: succinyl-CoA → (S)-methylmalonyl-CoA → propionyl-CoA, with MMCD as two lanes. The PPCSCT hourglass closes it at the bottom.
  - Succinate ↔ fumarate is drawn as three nested lanes: SUCDi on top; FRD2 and FRD3 below with their cofactors in the gaps.
  - PTA2 branches off to the left.
- **Linear degradations:** straight chains about 300 apart with cofactors on one side.
  - Phenylacetate degradation is a two-row snake.
  - 3-phenylpropionate degradation has its two entries (pppn, 3hpppn) meeting at dhpppn; hkntd joins op4en from below.
  - Allantoin degradation uses fans for UGLYCH and OXAMTC.
  - Methylglyoxal degradation is a triangle: the direct route GLYOX3 and the glutathione route via lgt__S.
- **Curvature:** t = 0.6.
- **Packing:** every zone except propionate fermentation is a detached piece and is placed by the toolkit's packing step. All titles are anchored to a node of their piece.
- **Area:** 7239 × 3428 (24.8 M) before packing, 2840 × 5278 (15.0 M) after. Final canvas 3294 × 5596.

## Decisions
- **3HBZCT** was proposed for "octane oxidation" with LOW CONFIDENCE (distance 742). It is aromatic, not octane, chemistry, so it gets its own descriptive title.
- **URIC** was proposed as a separate "allantoin degradation" piece. It is attached to the chain by the allantoin merge.
- **Exit labels:** all 14 keep the automatic owners. Some draft targets are odd (octanoate "→ Nucleotides", pyruvate "→ Nucleotides"); they are kept unchanged.
- **ACACCT** stays under the draft title "propanol degradation" as its own piece. PTA2, from the same draft group, is drawn at the propionyl-CoA hub.
- **Fumarate/succinate family:** SUCDi, FRD2 and FRD3 use one lane each rather than a stacked lane. Stacking produced a stack-internal crossing (v1).

## Known weak spots
- **QA (v3):** all structural checks PASS, 0 crossings, 0 dots closer than 45, 0 foreign lines within 28 of a dot, 0 label overlaps, 0 titles on lines; all labels are clean.
- **PPCSCT hourglass:** it spans 700 units and is drawn with strongly bowed curves.
- **propionyl-CoA copy:** the copy at the bottom has no exit label; it is the same metabolite as the hub.
- **Testing:** text widths are estimates. The map was validated structurally by the toolkit and checked in the matplotlib preview, but **not loaded in Escher**.
