# QA report

## Structural validation
- PASS: segment endpoints exist
- PASS: coordinates finite
- PASS: b1/b2 null or {x,y} pairs
- PASS: no temporary keys
- PASS: same reaction set (ids, bigg_id, stoichiometry, reversibility, genes)
- PASS: drawn metabolites identical to input, signs match stoichiometry
- PASS: every metabolite has position + label
- PASS: every reaction has a label position
- PASS: all nodes inside canvas
- PASS: all text labels inside canvas

## Crossings between different reactions (curves sampled): 4
- DDPGALA × DDPGALA: 1
- DKGLCNR2x × DKGLCNR2y: 1
- MLTG5 × MLTP3: 1
- TGBPA × TGBPA: 1

## Metabolite dots closer than 45 units: 0

## Foreign lines within 28 units of a dot: 0

## Overlapping label boxes: 0

## Titles touching lines/nodes: 0

## Label placement scores (<10 clean, 10-119 touches a line, >=120 overlaps)
- rxn: n=193 clean=193 line=0 overlap=0
- sec: n=247 clean=247 line=0 overlap=0
- pri: n=214 clean=212 line=2 overlap=0
- exit: n=43 clean=43 line=0 overlap=0

## Summary
```
{
 "validation_ok": true,
 "crossings": 4,
 "dots_too_close": 0,
 "lines_near_dots": 0,
 "label_overlaps": 0,
 "titles_on_lines": 0
}
```

## Area
- before packing 14408 x 8390 = 120.88 M
- after packing  14000 x 8054 = 112.76 M
- canvas 14501 x 8198 (margin incl.)

## Moves of detached pieces
- ['ACM6PH', 'ACMAMUT', 'ACNML', 'AGDC', 'AMANAPEr', 'AMANK', 'G1PACT', 'G6PDA', 'GF6PTA', 'PGAMT', 'UACMAMO', 'UAG2E', 'UAGDP']: dx=3 dy=10
- ['AACPS1', 'ACPPAT140', 'AGPAT140', 'ALKP', 'APG3PAT140', 'DHAPT', 'EDTXS4', 'F6PA', 'FA140ACPHi', 'G3PAT140', 'G3PD5', 'G3PD6', 'G3PD7', 'GLTPD', 'GLYCDx', 'GLYK', 'GPDDA4', 'M1PD', 'PFK_2', 'PFK_2', 'SBTPD', 'TGBPA', 'TGBPA']: dx=1093 dy=-280
- ['CLPNS160pp', 'CLPNS161pp', 'CLPNS181pp', 'GPDDA1pp', 'GPDDA2pp', 'GPDDA3pp', 'GPDDA4pp', 'LPLIPAL1A140pp']: dx=775 dy=-2580
- ['ALTRH', 'DDGALK', 'DDGALK', 'DDGLK', 'DDGLK', 'DDPGALA', 'DDPGALA', 'EDA', 'GALCTND', 'GALCTND', 'GALCTND', 'GUI1', 'GUI2', 'MANAO', 'MNNH', 'TAGURr']: dx=9 dy=-16
- ['2DGULRGx', '2DGULRGy', '2DGULRx', '2DGULRy', '5DGLCNR', '5DKGR', 'DKGLCNR1', 'DKGLCNR1', 'DKGLCNR2x', 'DKGLCNR2y', 'GNK', 'GNP', 'IDOND', 'IDOND', 'IDOND2']: dx=97 dy=-70
- ['ARABDI', 'ARABDI', 'ARAI', 'ARAI', 'LYXI', 'RBK_L1', 'RBK_L1', 'RBP4E', 'RBP4E', 'RPE', 'RPE', 'RU5PP', 'X5PL3E', 'X5PL3E', 'XYHDL', 'XYLI1', 'XYLK', 'XYLK', 'XYLK2', 'XYLK2']: dx=-8 dy=-10
- ['FCI', 'FCLK', 'FCLPA', 'LKDRA', 'RHMND', 'RMI', 'RMK', 'RMPA']: dx=19 dy=0
- ['G1PPpp', 'GLCDpp', 'LACZpp', 'TREHpp', 'UDPGPpp']: dx=3657 dy=150
- ['GDPMNH', 'GDPMNP', 'GMAND', 'HEX4', 'MAN1PT2', 'MAN6PI', 'MN6PP', 'PMANM']: dx=475 dy=-170
- ['GALCTD', 'GLCRAL', 'GLCRD', 'GLYCK2', 'HPYRI', 'MANPGH', 'TRSARr']: dx=-2552 dy=-3230
- ['GPDDA1', 'GPDDA2', 'GPDDA3', 'LPLIPAL2A140']: dx=-65 dy=2880
- ['GALT1', 'GLCTR1', 'GLCTR2', 'GLCTR3']: dx=-714 dy=-1690
- ['2DDARAA', 'DXYLTD', 'GCALDD', 'PGLYCP']: dx=9 dy=-20
- ['A5PISO', 'KDOCT2', 'KDOPP', 'KDOPS']: dx=317 dy=-385
- ['3SLAR', '6D6SPA']: dx=-172 dy=10
- ['G1PTT', 'TDPDRE', 'TDPDRR', 'TDPGDH']: dx=-1125 dy=-176
- ['KDUI']: dx=-482 dy=0
- ['O16GLCT2']: dx=-15 dy=-430
- ['R5PP', 'RBK']: dx=16 dy=-606
- ['UACGAMPpp']: dx=627 dy=-350
- ['R5PPpp']: dx=207 dy=650
- ['UDPGALPpp']: dx=127 dy=-630

## Title nudges
- D-xylose degradation (Dahms pathway): (-85, 85) ok
- sulfoquinovose degradation: (-60, -104) ok
- D-glucarate / D-galactarate degradation: (-28, 28) ok
- glycerophosphodiester degradation: (-386, 104) ok