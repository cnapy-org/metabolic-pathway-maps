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

## Crossings between different reactions (curves sampled): 5
- DDPGALA × DDPGALA: 1
- DKGLCNR2x × DKGLCNR2y: 1
- TGBPA × TGBPA: 1
- MLTG5 × MLTP3: 1
- 2DGULRGx × 2DGULRGy: 1

## Metabolite dots closer than 45 units: 0

## Foreign lines within 28 units of a dot: 0

## Overlapping label boxes: 0

## Titles touching lines/nodes: 0

## Label placement scores (<10 clean, 10-119 touches a line, >=120 overlaps)
- rxn: n=193 clean=193 line=0 overlap=0
- sec: n=247 clean=247 line=0 overlap=0
- pri: n=214 clean=214 line=0 overlap=0
- exit: n=43 clean=43 line=0 overlap=0

## Summary
```
{
 "validation_ok": true,
 "crossings": 5,
 "dots_too_close": 0,
 "lines_near_dots": 0,
 "label_overlaps": 0,
 "titles_on_lines": 0
}
```

## Area
- before packing 17131 x 8370 = 143.39 M
- after packing  17131 x 8370 = 143.39 M
- canvas 17431 x 8674 (margin incl.)

## Packing: mode='keep' (pieces at hand positions); overlapping pieces: 0

## Moves of detached pieces

## Title nudges
- D-xylose degradation (Dahms pathway): (-85, 85) ok
- sulfoquinovose degradation: (-60, -104) ok
- D-glucarate / D-galactarate degradation: (-28, 28) ok
- glycerophosphodiester degradation: (-386, 104) ok