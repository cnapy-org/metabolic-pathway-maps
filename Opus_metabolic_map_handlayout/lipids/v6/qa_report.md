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

## Crossings between different reactions (curves sampled): 0

## Metabolite dots closer than 45 units: 0

## Foreign lines within 28 units of a dot: 0

## Overlapping label boxes: 0

## Titles touching lines/nodes: 0

## Label placement scores (<10 clean, 10-119 touches a line, >=120 overlaps)
- rxn: n=159 clean=159 line=0 overlap=0
- sec: n=337 clean=337 line=0 overlap=0
- pri: n=178 clean=178 line=0 overlap=0
- exit: n=25 clean=25 line=0 overlap=0

## Summary
```
{
 "validation_ok": true,
 "crossings": 0,
 "dots_too_close": 0,
 "lines_near_dots": 0,
 "label_overlaps": 0,
 "titles_on_lines": 0
}
```

## Area
- before packing 11886 x 9046 = 107.52 M
- after packing  11886 x 9046 = 107.52 M
- canvas 12186 x 9364 (margin incl.)

## Packing: mode='keep' (pieces at hand positions); overlapping pieces: 0

## Moves of detached pieces

## Title nudges