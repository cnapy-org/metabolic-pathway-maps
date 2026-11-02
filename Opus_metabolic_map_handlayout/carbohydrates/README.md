# iML1515 carbohydrates – hand-laid-out Escher map

Input: `ecmapper/BRENDAmap/draft_map_iML1515/split/carbohydrates.escher.json` (BRENDA draft).
Layout: written by hand in `spec.py`, using the `escher_handlayout` toolkit (`../escher_handlayout_toolkit`).

```bash
PYTHONPATH=../escher_handlayout_toolkit python -m escher_handlayout run spec.py --out v4 --render
```

| file | content |
|---|---|
| `v4/carbohydrates.escher_final.json` | final Escher map |
| `v4/qa_report.md`, `v4/run_log.json` | QA report and run log (splits, merges) |
| `v4/overview.png` | matplotlib preview (not Escher) |
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

Canvas: **17431 × 8674** (margin included); the hand arrangement spans 17131 × 8370. Packing is off
(`PACK = dict(mode='keep')`), so every piece stays where `spec.py` puts it. Curvature t = 0.6.
v3 used automatic packing (canvas 14501 × 8198); v4 is about 1.27× that area but keeps related
pathways next to each other.

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
- **Detached pieces, each with its own title, placed by hand (see "Global arrangement"):**
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

## Global arrangement (v4, by hand)

The main component sits in the middle (x −2500…1440). Around it:

- **West, top band:** the periplasmic glycerophosphodiesters/cardiolipin piece and the periplasmic
  sugar hydrolases/phosphatases side by side.
- **West, middle band:** below them, the three single periplasmic nucleotide-sugar phosphatases
  (UDPGALPpp, UACGAMPpp, R5PPpp) in a row on the left and cytoplasmic glycerophosphodiester
  degradation on the right.
- **West, far left:** amino sugars + UDP-GlcNAc (f6p/gam6p). Their only links to the rest are exit
  labels to f6p.
- **West, bottom:** sugar alcohols (glycerol / hexitols / galactitol), below the glycerophosphodiester
  block they share glyc3p with.
- **East, row 1 (glycan precursors, next to g1p/Leloir):** LPS core glucosylation (UDP-glucose),
  O16 antigen, CMP-KDO, dTDP-L-rhamnose above the D-mannose block (GDP-mannose), and D-xylose (Dahms).
- **East, row 2 (pentoses and methylpentoses):** pentoses, ribose, L-rhamnose/L-fucose,
  sulfoquinovose and KDUI.
- **East, row 3 (sugar acids):** ketogluconates, then hexuronates/galactonate with
  D-glucarate/galactarate below.

To place the periplasm, ribose and O16 pieces independently, three zones were split off without any
change to the biochemistry: `PGP`, `PNT`, `RIB` and `O16` (from PERI, PEN and LPS). The QA report lists
0 overlapping pieces.

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
  node (f6p `2328`, g6p `2317`), so they are removed with `DROP_LABELS`. `→ Central, Amino acids` stays on g6p.
- **Packing off (v4):** automatic packing stacked all pieces vertically next to the tall main
  component (4693 × 18360). v3 worked around that with two tiny keep-out boxes, but the packer still
  chose where each pathway went. v4 sets `PACK = dict(mode='keep')` and places every piece by hand.
- **No toolkit change was needed.**

## Known weak spots (from `v4/qa_report.md`)

**Structural validation:** all PASS.

**Crossings: 5.** All are between reactions of one lane family; there are no crossings between
different pathways.

| crossing | explanation |
|---|---|
| DDPGALA × DDPGALA | the two duplicate fans to g3p + pyr cross once |
| TGBPA × TGBPA | the same happens with the duplicate fans to dhap + g3p |
| DKGLCNR2x × DKGLCNR2y | stack-internal crossing on the shared diagonal lane |
| MLTG5 × MLTP3 | stacked on the same maltodextrin level; the MLTP3 line to g1p leaves across the MLTG5 lane |
| 2DGULRGx × 2DGULRGy | not a visible crossing: the two stacked rows share one straight lane, and the QA test reports overlapping collinear segments (cross products ≈ 0) depending on floating-point rounding. The geometry is the same as in v3, which by chance did not trigger it |

**Other checks:**
- Dots closer than 45: 0.
- Foreign lines within 28 of a dot: 0.
- Label overlaps: 0.
- Titles on lines: 0.
- All 697 labels are clean (no label touches a line).

**Layout:**
- The left column (amino sugars) leaves empty space below it, and the west middle band is sparse.
  The width is set by the top bands (west: amino sugars + periplasm pieces; east: glycan-precursor
  row), so moving pieces into those gaps would not shrink the canvas.

**Not covered by QA:**
- Text widths are estimates, and arrowheads are not modelled by QA.
- UGLT intentionally draws its UDP-glc/UDP-gal pair as crossing lines within the reaction itself.

**Validated structurally by the toolkit; not loaded in Escher.**
