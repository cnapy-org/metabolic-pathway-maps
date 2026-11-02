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
- rxn: n=45 clean=45 line=0 overlap=0
- sec: n=82 clean=82 line=0 overlap=0
- pri: n=40 clean=40 line=0 overlap=0
- exit: n=1 clean=1 line=0 overlap=0

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
- before packing 6933 x 4138 = 28.68 M
- after packing  5680 x 3260 = 18.51 M
- canvas 5980 x 3568 (margin incl.)

## Moves of detached pieces
- ['FBA', 'FBA3', 'FBP', 'RBP4E', 'RPE', 'RPI', 'TKT1', 'TKT2', 'TPI', 'X5PL3E']: dx=-643 dy=1292
- ['PGK']: dx=-1253 dy=198
- ['BUTCT']: dx=-1005 dy=932
- ['ECOAH1']: dx=-278 dy=2292

## Title nudges