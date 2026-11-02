# Nucleotides and cofactors (iML1515 split map): hand layout

This is a hand-laid-out Escher map of
`ecmapper/BRENDAmap/draft_map_iML1515/split/nucleotides_and_cofactors.escher.json`, made with
`../escher_handlayout_toolkit` by following its `PROMPT.md`.

**Files**
- `spec.py`: the hand coordinates, commented per zone.
- `v7/nucleotides_and_cofactors.escher_final.json`: the final Escher map.
- `v7/qa_report.md` and `v7/run_log.json`: QA report and run log of the final run.
- `v7/overview.png`: a matplotlib preview. It is not an Escher rendering.
- `inspect/`: the Phase 0 report and the generated skeleton.

**Reproduce**

```bash
PYTHONPATH=../escher_handlayout_toolkit python -m escher_handlayout run spec.py --out v7 --render
```

## Counts

| | draft | final |
|---|---|---|
| reactions | 305 | 305 (identical: ids, bigg_id, stoichiometry, reversibility, genes) |
| nodes | 1825 | 1889 |
| metabolite nodes | 910 | 974 |
| text labels | 51 (18 titles, 33 exit labels) | 54 (21 titles, 33 exit labels) |
| canvas | 12976 × 8246 | 17728 × 10406 |

Nodes go up because of the changes below. Primary-node splits that are later merged back count
twice in the run log, so the totals there do not simply add up.
- **Automatic cofactor splits:** shared cofactor nodes are split into one copy per reaction.
- **My primary-node splits:** listed under Decisions; `run_log.json` lists all 81 splits.
- **Merges:** 17 merges remove nodes again.

**Area:** 17428 × 10083 = 175.7 M map units² (canvas 17728 × 10406 with margin). The draft canvas
was 107 M. v6 reached 91.5 M (canvas 13304 × 7344) with automatic packing, but only by dropping small
pieces into holes of other zones (see "Global arrangement").
**Curvature:** t = 0.6, the default.

## Zones (subsystems)

The subsystem titles come from the draft's text labels.

**Reused titles:** purine, pyrimidine, NAD, 6-hydroxymethyl-dihydropterin diphosphate,
tetrahydrofolate, flavin, heme, vitamin B6, vitamin B1, vitamin K, molybdenum cofactor, biotin, ppGpp,
vitamin B12, ascorbate and lipoate.

**Changed titles:**
- `pantothenate biosynthesis` was renamed **pantothenate and coenzyme A biosynthesis**. The two
  pathways form one chain after merging pantothenate, and the separate "coenzyme A metabolism" title
  was dropped.
- **New titles**, needed because the draft grouped these reactions wrongly by proximity:
  - flavin reductases: siderophore iron and sulfonate reduction. These 20 reactions sat under
    "vitamin B12 metabolism" in the draft.
  - serine and one-carbon (THF) interconversion. SERD_D/L, GHMT2r, MTHFD and MTHFR2 were LOW
    CONFIDENCE "purine" in the draft.
  - acetate kinase / phosphotransacetylase. ACKr and PTAr were LOW CONFIDENCE "purine".
  - periplasmic nucleotidases and phosphatases. These are the 11 `*pp` single reactions plus
    GAL1PPpp, which sat under several unrelated titles in the draft.

**Layout principles**
- **Purine:** the NMP row AMP–IMP–XMP–GMP sits above the nucleoside row (adn, ins, xtsn, gsn), which
  sits above the base row (ade, hxan, xan, gua).
  - Phosphoribosyltransferases run as bypass lanes between the columns.
  - De novo IMP synthesis is a vertical column above IMP.
  - ADP/ATP and the dA nucleotides are on the left; GDP/GTP, cGMP and the dG nucleotides on the right.
  - The ribose-5-phosphate cluster is bottom left.
- **Pyrimidine:** a de novo column ends in UMP. Then come the UMP–UDP–UTP–CTP–CDP row, the CTP/CDP/CMP
  triangle, and a cytidine/uridine/uracil/cytosine salvage rectangle below it. The dU/dT/dC
  nucleotides are above it, with dCTP→dUTP as one long horizontal link.
- **NAD:** a de novo column runs Asp → iminoaspartate (3 stacked oxidases) → quinolinate → NaMN →
  deamido-NAD → NAD. NAD⇄NADP lanes and the transhydrogenases are on the right, the NMN/nicotinamide
  salvage on the left.
- **GTP-derived cofactors share one GTP node:**
  - Above GTP: pterins (2 GTPCI lanes) → 6-hydroxymethyl-dihydropterin-PP → dihydropteroate → DHF →
    THF, plus the pABA branch from chorismate, the preQ0 (queuosine) branch and the tetrahydromonapterin
    branch.
  - Below GTP: riboflavin → FMN → FAD (2 GTPCII2 lanes).
- **Heme:** a straight column from Glu-tRNA to protoheme and heme O, with the siroheme branch on the
  right.
- **Small pathways (B6, B1, K, MoCo, CoA, biotin, ppGpp, B12, ascorbate, lipoate):** chains, snakes or
  hourglasses. Zone titles are anchored to a node of their zone.

## Global arrangement (hand-chosen, since v7)
Every zone and every detached piece is placed by hand (`ORIGIN`). `PACK = dict(mode='keep')` turns
automatic packing off; the QA report confirms that no two pieces overlap. In v6 the packer had placed
49 detached pieces wherever they fit, so e.g. the flavin-reductase rows were scattered over the whole
map and several titles sat far from their pieces.

| column | contents (top → bottom) |
|---|---|
| A | pyrimidine; vitamin B6 + vitamin B1; pantothenate/CoA + lipoate; vitamin K + biotin |
| B | purine, with the periplasmic nucleotidases + acetate kinase/PTA in its empty upper-left corner and serine/one-carbon + ppGpp in its upper-right corner (ppGpp next to GTP/GDP, serine/C1 next to THF in column C); below purine the flavin-reductase block |
| C | GTP-derived cofactors around the shared GTP node (pterin + THF above, riboflavin/FMN/FAD below); molybdenum cofactor below |
| D | NAD; heme; vitamin B12; ascorbate |

The flavin-reductase zone is a compact grid: FMNRx top left, the 6 siderophore rows as 2 rows × 3
columns, the 5 sulfonate monooxygenases as one column on the right, MTRPOX below.

## Decisions

### Merges
These are key metabolites only.

| merged node | what it connects |
|---|---|
| GTP (1 ← 527, 2566) | pterin, THF/queuosine and flavin branches share one GTP |
| ahdt, dhnpt, 6hmhptpp, 25drapp | the pterin/THF/flavin chains (one draft copy each) |
| 4-aminobenzoate | the pABA branch |
| pantothenate | pantothenate → CoA |
| GDP in ppGpp | the ppGpp cycle |
| adocbip | the vitamin B12 chain |
| 5,10-methylene-THF | the C1 piece |
| dscl, scl | the duplicate siroheme branch |
| cytidine, CMP, cytosine | the pyrimidine salvage |

### Kept separate, intentionally
These stay as separate nodes because merging would draw lines thousands of units long across other
pathways:
- GTP in MoCo and ppGpp
- the second purine GTP
- pyruvate (4×), L-alanine, chorismate, dxyl5p, adenine, acetaldehyde
- 2-deoxyribose-1-P (purine vs pyrimidine)
- FMN/riboflavin in the flavin-reductase zone

### Splits of primary nodes, one copy per reaction
- **AMPN:** own ade and r5p, as a short fan left of AMP.
- **5DOAN:** own ade.
- **PUNP4:** own 2dr1p.
- **DCYTD:** own dcyt next to duri.
- **TMDPP:** own 2dr1p.
- **RHAT1 and AADDGT:** own dTDP. These are LPS/ECA reactions.
- **4ABZGLUH:** own glutamate, so it doesn't pull heme's Glu.
- **ALATA_L2/ALATA_D2:** own ala/pyr copies. This avoids crossings with the cysteine desulfurases and
  the vitamer lanes.
- **MoCo:** Cu2+ product copies. The tungsten branch (WCOS, BWCOS) gets its own MPT-AMP copy and
  becomes a separate piece.
- **Flavin-reductase zone:** every reaction gets its own fmnh2/fmn (or rbflvrd/ribflv), fe2 and
  sulfite copy. Each siderophore is one two-lane row (FMNH2 lane above, reduced-riboflavin lane
  below); each sulfonate monooxygenase is one row. In the draft these 20 reactions formed a hairball
  around four shared nodes.

### Other decisions
- **Duplicate reactions** are kept and drawn as lanes: GTPCI×2, GTPCII2×2, NADTRHD×3, UMPK×2,
  SHCHD2×2, UPP3MT×3, GARFT×2.
- **Variant families on one lane:** ASPO6/3/4, DHORD2/5/fum, ATPM/NTP1, CPPPGO2/CPPPGO, PDX5POi/PO2
  and FLVR/FLVRx are stacked on one lane each.
- **Reactions with no primary metabolite on one side** have their cofactors at the free end of the
  axis: UPP3MT×3, MOCDS, AACPS9, CMPN, TDP.
- **Exit labels** use the automatic owner rule; no `NEW_EXITS` were added.
- **GMPR** (GMP → IMP) runs above XMP and enters IMP at a fixed 45° tangent.

## Known weak spots (non-zero QA items, final run)
- **1 crossing, NADTRHD (nad + nadph → nadp) × THD2pp (nadp → nad + nadph).** These are two opposite
  reactions between the same three nodes, so one X is hard to avoid without more copies.
- **4 foreign lines within 28 units of a dot:**
  - GMPR passes a cofactor dot of IMPC (h2o, 11 units).
  - NTPP5 passes DADK's ATP dot (27).
  - NTD5 passes TMDK1's ATP/ADP dots (27); these are neighbouring lanes.
- **Large bulges:** the bypass lanes (ADPT, HXPRT, XPPT, GUAPRT, UPPRT, CYTK2, NTPTP1) and the PYK6
  hourglass draw as large curves. This is readable but bulky.
- **Area:** 175.7 M, 1.9× the packed v6 (91.5 M). Keeping each zone together leaves white space
  (e.g. right of the flavin-reductase block, between purine and the GTP column, under column A), and the
  small zones in the purine corners are the only nesting. A target of ≤ 1.3× v6 was not reachable
  without putting pieces of one zone into holes of another.
- **Estimates:** text widths and the QA heuristics are estimates.
- **No dots closer than 45 units**, no overlapping label boxes and no titles on lines.

## Validation status
The map was validated **structurally** by the toolkit, and every check passes. The checks cover:
- the same reaction set, stoichiometry, reversibility and genes
- the drawn metabolites
- the segment endpoints and Bezier control points
- the canvas

The map was **not loaded in Escher**; only the matplotlib preview was inspected.
