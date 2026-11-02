# Lipids map (iML1515, BRENDA draft): hand layout

- **Input:** `ecmapper/BRENDAmap/draft_map_iML1515/split/lipids.escher.json`
- **Output:** `v5/lipids.escher_final.json`. This folder also holds `spec.py` (the hand coordinates), `v5/qa_report.md`, `v5/run_log.json`, `v5/overview.png` (matplotlib preview, not Escher) and `inspect/` (Phase 0 report).

The map was made with `escher_handlayout` and follows `PROMPT.md`. It is **validated structurally by the toolkit and has not been loaded in Escher.**

Reproduce from the repository root:

```bash
PYTHONPATH=Opus_metabolic_map_handlayout/escher_handlayout_toolkit python -m escher_handlayout run Opus_metabolic_map_handlayout/lipids/spec.py --out Opus_metabolic_map_handlayout/lipids/v5 --render
```

## Counts

| | draft | final |
|---|---|---|
| reactions | 159 | 159 (same ids, stoichiometry, reversibility and genes) |
| nodes | 924 (447 metabolite) | 992 (515 metabolite) |
| text labels | 30 | 35 (8 zone titles, 2 sub-titles, 25 exit labels) |
| canvas | 9556 × 6910 | 12215 × 9378 |

The extra metabolite nodes come from:
- 12 automatic cofactor splits;
- the side-product copies of fatty acids, acetyl-CoA and malonyl-ACP;
- one backbone copy per phospholipid unit.

`run_log.json` lists all 77 splits and 9 merges.

**Area:** the hand arrangement is 11886 × 9046 = 107.5 M. After packing it is 12051 × 9301 = 112.1 M.

## Layout (zones and sources)

Draft titles were reused and renamed where the content differed.

| zone | title | source / decision |
|---|---|---|
| FAS | fatty acid biosynthesis (FAS II) | was 'palmitate biosynthesis' plus the FAS reactions from 'lipid metabolism' |
| BOX | acyl-CoA metabolism / β-oxidation | was 'lipid metabolism' |
| PLS | phosphatidic acid synthesis | was 'arachidonate biosynthesis'; that name did not match its content (PlsB/PlsX/PlsY/PlsC, acyl-ACP synthetases) |
| GLYC | glycerol / glycerol-3-phosphate | new zone |
| PLD | phospholipid degradation and re-acylation | reused title 'mevalonate metabolism' |
| LIPA | lipid A biosynthesis | Raetz pathway + Kdo + EptA + PagP |
| ISO | isoprenoid biosynthesis (MEP pathway) | IPDDI was moved here from 'mevalonate' |
| ECA | ECA / O-antigen precursors | was 'teichoic acid biosynthesis'; E. coli has no teichoic acid, and WecA (ACGAMT) and TDPADGAT belong to ECA |

**Repeated cycles are drawn with one grammar.**
- **FAS snake.** Each elongation round C4…C18 is one row: 3-oxoacyl → 3-hydroxyacyl → trans-2-enoyl → acyl-ACP, with x slots 0/350/700/1150 and rows 600 apart. Rows alternate direction. The condensation (3OAS/KAS) is the vertical step to the next row; its cofactors are on the outer side of the turn.
  - Reduction cofactors are always below the row.
  - The two enoyl reductases (NADH and NADPH, EARx/EARy) are a lens of two lanes.
  - Where the model lacks a step, or a chain intermediate is only a side metabolite (butACP, ocACP, dcaACP, 3hmrsACP, ocdcaACP), the slot stays empty or holds that cofactor dot. The grid is never compressed, so every round sits at the same place in its row.
- **β-oxidation snake** (same grammar, CoA esters). Each row is one chain length: acyl → enoyl → 3-hydroxy → 3-oxo, with slots 550 apart and rows 750 apart. The thiolase is the vertical step.
  - Each acyl-CoA (C18:1, C16:1, C16:0, C6) has the same activation pair on its outer side: periplasmic FA → acyl-CoA (FACOAL…t2pp) above, and acyl-CoA → cytosolic FA (FACOAE…) below.
  - The four btcoa ⇄ b2coa reactions (ACOAD1f ×3, ACOAD1fr) are four lanes with cofactors in the gaps.
- **PA synthesis.** There is one identical column per acyl chain (16:0, 18:1, 16:1):
  - acyl-ACP / acyl-phosphate / free FA on the left, with the ACP ⇄ FA reaction as an outer lane;
  - glycerol-3-P → 1-acyl-G3P with two lanes, G3PAT (from acyl-ACP, above) and APG3PAT (from acyl-P, below).
  - The 16:0 column adds AGPAT160 → PA. The 16:1 column adds 3OAS181 (cis-vaccenate elongation).
- **Phospholipid degradation.** There is one identical unit for each (head group PA/PE/PG × chain 16:0/16:1/18:1).
  - **Periplasm:** a 'Λ' of PLIPA1 (→ 2-acyl) and PLIPA2 (→ 1-acyl), followed by LPLIPAL1 (1-acyl → backbone).
  - **Cytoplasm:** a column of 2AGPxAT re-acylation (PE, PG), then 2-acyl-lyso-PL, then LPLIPAL2 → backbone.

## Decisions

- **Side products.** Free fatty acids released or consumed in the phospholipid-degradation reactions, acetyl-CoA in the thiolases, HXCT and TDPADGAT, malonyl-ACP in 3OAS180/3OAS181, acetate in HXCT, PE and DAG in EptA (PETNT), and palmitate in PagP (LIPAHT) are drawn as small side-product dots. This means one node copy per reaction, with `node_is_primary = false` on those copies. In the draft these were hubs that tied about 45 reactions into a hairball. The biochemistry is unchanged (`validate` PASS).
- **Merges.** These are duplicates of the same metabolite that make a chain continuous:
  - actACP 743 → 336;
  - 3hdecACP 364 → 850 and tdec2eACP 366 → 863, which turns the duplicated 3HAD100 into a lane pair;
  - the side 3ooctdACP of 3OAR180 → 1018, so the C18 round is continuous;
  - glyc3p_p 1273 → 1259;
  - ipdp 705 → 151 and dmpp 706 → 136, so IPDDI joins the MEP pathway.
- **Deliberately not merged** (no long links between zones):
  - palmACP (FAS 954 vs PLS 2), with exit label '← fatty acid biosynthesis' at PLS;
  - hdeACP (PLS 1818 vs lipid A 2025), with exit label '→ lipid A (EDTXS3)';
  - free fatty acid copies in FAS, BOX, PLS, PLD and LIPA;
  - acetyl-CoA: FAS initiation hub 389, ACACT1r copy 715, and side copies elsewhere.
- **Split primaries.**
  - glycerol-3-P: one copy per PLS acyl-chain column, one for the glycerol branch, and one for each cytosolic PA unit;
  - one backbone (G3P / G3PE / G3PG) per phospholipid-degradation unit instead of three-way hubs;
  - Kdo for MOAT3C;
  - FA161ACPHi's hexadecenoate placed in the lipid A zone, next to the second hdeACP.
- **ACACT1r** (2 acetyl-CoA → acetoacetyl-CoA) sits in the 3-oxo slot of the C4 β-oxidation row instead of under 'mevalonate'.
- **Packing.** The packer would move most of the hand-placed pieces into holes of other zones: the first run moved 40 pieces, e.g. the C14–C18 FAS rows by 4900 units. `PACK.keep_out` therefore holds two 1 × 1 boxes at the corners of the hand layout. The reference rectangle is then the whole map and every detached piece stays where it was placed. The trade-off is a larger canvas (112 M instead of the packer's 45 M), but the arrangement stays readable.
- **Other settings:** curvature t = 0.6 (default); bend overrides only for the activation-pair ligases (0) and the 4-lane ACOAD1 family (15).

## QA (final run)

| check | result |
|---|---|
| structural validation | all PASS |
| crossings | 0 |
| metabolite dots closer than 45 | 0 |
| foreign lines within 28 of a dot | 0 |
| overlapping label boxes | 0 |
| titles touching lines or nodes | 0 |
| label scores | all clean (rxn 159, sec 337, pri 178, exit 25) |

## Known weak spots

- **Packing jitter.** Packing quantises translations to its 40-unit grid, so detached pieces move by up to 20 units (TDPADGAT, the largest, by (−5, 18)). Rows of the FAS and β-oxidation snakes that are separate pieces can therefore be offset by a few units from perfect alignment.
- **Empty slots.** The FAS and β-oxidation snakes contain empty slots wherever iML1515 lacks a step, e.g. rows C10/C8 of β-oxidation hold only the thiolases. This is intended (fixed grid) but leaves white space.
- **Fragmented map.** Because of the side-product convention, the map is made of many small pieces. Connections between zones are shown by exit labels only.
- **Estimates.** Text widths and all QA checks are heuristics (see the toolkit README). The map has not been opened in Escher.
