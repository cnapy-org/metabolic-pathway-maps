# iML1515 carbohydrates – hand-laid-out Escher map

Input: `ecmapper/BRENDAmap/draft_map_iML1515/split/carbohydrates.escher.json` (BRENDA draft).
Layout: written by hand in `spec.py`, using the `escher_handlayout` toolkit (`../escher_handlayout_toolkit`).

```bash
PYTHONPATH=../escher_handlayout_toolkit python -m escher_handlayout run spec.py --out v3 --render
```

| file | content |
|---|---|
| `v3/carbohydrates.escher_final.json` | final Escher map |
| `v3/qa_report.md`, `v3/run_log.json` | QA report and run log (splits, merges, packing) |
| `v3/overview.png` | matplotlib preview (not Escher) |
| `inspect/phase0_report.md` | Phase 0 report on the draft |

## Counts

| | draft | final |
|---|---|---|
| reactions | 193 | 193 (identical) |
| nodes | 1041 (462 metabolite) | 1040 (461 metabolite) |
| text labels | 51 (15 titles, 36 exits) | 69 (25 titles, 43 exits) |

The toolkit split shared cofactors automatically. On top of that, `prepare()` made 44 splits in total
and there are 45 merges. Most merges remove duplicate drawings of the same metabolite inside one
pathway. The two effects roughly cancel out.

Canvas: **14501 × 8198** (margin included). The bounding box was 14408 × 8390 before packing and
14000 × 8054 after packing. Curvature t = 0.6.

## Layout

- **Main component (fixed):** glucose-1-phosphate hub.
  - **Left:** the maltodextrin column G2…G7. Each level has three stacked reactions:
    - MLTG (hydrolase, up), with a local glucose copy on the left;
    - MLTP (phosphorylase, up), whose line runs to the g1p hub;
    - AMALT (amylomaltase, down), with local maltose donor and glucose copies on the left.
  - **Right:** the glycogen cycle (GLGC → GLCS1 → GLCP).
  - **Below:** g6p via two PGMT lanes, and the G1PP/HEX1 glucose triangle.
  - **Above:** the Leloir pathway. UGLT is drawn with gal1p → g1p on the main line and the
    UDP-glc/UDP-gal pair crossing as co-substrates.
  - **Below g6p:** the trehalose chain down to the glucose hub of the disaccharide block. Around it:
    - lactose/melibiose fans from below;
    - fructose/f1p/f6p on the left;
    - GlcNAc/maltose acetylation and PAI2T on the right.
- **Detached pieces, each with its own title, arranged by the packer:**
  - UDP-GlcNAc + amino sugars, joined by merged f6p/gam6p/UDP-GlcNAc;
  - CMP-KDO;
  - dTDP-L-rhamnose;
  - LPS core glucosylation, plus O16 antigen;
  - D-mannose, drawn as a square with lanes;
  - pentoses: a grid with the D-series on top and the L-series below, plus ribose;
  - D-xylose (Dahms);
  - ketogluconate: a diamond 25dkglcn / 2dhguln / 5dglcn / idonate / gluconate, with variant families
    on stacked lanes;
  - sugar acids: hexuronates → 2ddg6p → EDA → g3p + pyr, and D-galactonate → DDPGALA, which shares
    g3p/pyr;
  - D-glucarate/galactarate;
  - sugar alcohols: glycerol square, galactitol/tagatose column, hexitols → f6p → F6PA, and the
    myristoyl-ACP lipid reactions above glyc3p;
  - glycerophosphodiester degradation;
  - L-rhamnose/L-fucose: two symmetric chains into lald + dhap;
  - sulfoquinovose;
  - sulfated GAG (KDUI);
  - periplasm pieces.

## Decisions

- **Zones:** the draft's 15 titles are reused. Some were renamed or extended, and 10 titles were
  added (`EXTRA_TITLES`).
  - **Moved out of "degradation of pentoses" (LOW CONFIDENCE in the report):**
    - rhamnose/fucose → own block "L-rhamnose / L-fucose degradation";
    - 2DDARAA/GCALDD/PGLYCP → the Dahms d-xylose block.
  - **Moved out of "metabolism of disaccharids" (LOW CONFIDENCE):** the hexuronate reactions
    (GUI1/2, MANAO, MNNH, TAGURr, ALTRH, DDGLK) → sugar acids.
  - **Moved out of "sulfated glycosaminoglycan metabolism" (LOW CONFIDENCE):** EDA/DDGLK → sugar
    acids. KDUI stays alone under the GAG title.
  - **Glycogen block:** XYHDL → pentoses (xylose node).
  - **Disaccharide block:**
    - PFK_2/TGBPA duplicates → sugar alcohols, merged with their twins;
    - FRUK/F1PP → disaccharides (fructose);
    - MANPGH → glucarate block, at glyc__R.
  - **Periplasmic reactions:** grouped into periplasm pieces with "periplasm" titles. Escher labels
    hide the compartment suffix.
- **Merged** (`MERGES`):
  - **g1p hub:** 4 copies, joining glycogen, maltodextrins, Leloir and MLTP.
  - **g6p hub:** 4 copies, joining PGMT, HEX1, TRE6PS and BGLA1.
  - **Within pathways:** duplicate drawings of pentoses, ketogluconates, hexuronates/galactonate,
    sugar alcohols, amino sugars and periplasmic glucose.
- **Intentionally not merged:** these keep their own nodes and, where relevant, an exit label
  (`NEW_EXITS`).
  - glucose of the disaccharide block and of the maltodextrin steps (local copies; they would
    otherwise form a hairball);
  - galactose of LACZ/GALS3 (`→ Leloir pathway`);
  - g1p of the dTDP-rhamnose piece (`← glycogen metabolism (g1p)`);
  - all f6p/dhap/g3p/pyr of the other pieces (`→ Central`);
  - man6p of MANPGH (`→ D-mannose degradation`);
  - the glyc3p hub of the glycerophosphodiesterases (`→ sugar alcohols (glyc3p)`).
- **Local donor copies (no exit label):**
  - UDP-glucose for GLCTR1/GALT1/GLCTR2/GLCTR3/O16GLCT2;
  - maltose donors for AMALT1–4;
  - glucose copies of the MLTG/AMALT/BGLA1 steps;
  - ttdcap for APG3PAT140, ttdca for LPLIPAL2A140.

  These act as co-substrates or side products next to their reaction.
- **Duplicate exit labels removed:** two draft exit labels would have landed twice on the same merged
  node (f6p `2328`, g6p `2317`), so they were deleted. `→ Central, Amino acids` stays on g6p.
- **Packing:** `PACK.keep_out` holds two 10×10 marker boxes. They widen the packing target: without
  them the packer stacked all pieces vertically next to the tall main component, giving 4693 × 18360.
  With them the result is landscape at about the same area.
- **No toolkit change was needed.**

## Known weak spots (from `v3/qa_report.md`)

**Structural validation:** all PASS.

**Crossings: 4.** All are between reactions of one lane family; there are no crossings between
different pathways.

| crossing | explanation |
|---|---|
| DDPGALA × DDPGALA | the two duplicate fans to g3p + pyr cross once |
| TGBPA × TGBPA | the same happens with the duplicate fans to dhap + g3p |
| DKGLCNR2x × DKGLCNR2y | stack-internal crossing on the shared diagonal lane |
| MLTG5 × MLTP3 | stacked on the same maltodextrin level; the MLTP3 line to g1p leaves across the MLTG5 lane |

**Other checks:**
- Dots closer than 45: 0.
- Foreign lines within 28 of a dot: 0.
- Label overlaps: 0.
- Titles on lines: 0.
- Two primary labels touch a line (label score < 120).

**Not covered by QA:**
- Text widths are estimates, and arrowheads are not modelled by QA.
- UGLT intentionally draws its UDP-glc/UDP-gal pair as crossing lines within the reaction itself.

**Validated structurally by the toolkit; not loaded in Escher.**
