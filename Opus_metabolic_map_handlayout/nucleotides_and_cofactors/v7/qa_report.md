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

## Crossings between different reactions (curves sampled): 1
- NADTRHD × THD2pp: 1

## Metabolite dots closer than 45 units: 0

## Foreign lines within 28 units of a dot: 4
- dot h2o_c (1302) ← line of GMPR: 11
- dot atp_c (1113) ← line of NTPP5: 27
- dot atp_c (2502) ← line of NTD5: 27
- dot adp_c (2503) ← line of NTD5: 27

## Overlapping label boxes: 0

## Titles touching lines/nodes: 0

## Label placement scores (<10 clean, 10-119 touches a line, >=120 overlaps)
- rxn: n=305 clean=305 line=0 overlap=0
- sec: n=606 clean=606 line=0 overlap=0
- pri: n=368 clean=368 line=0 overlap=0
- exit: n=33 clean=33 line=0 overlap=0

## Summary
```
{
 "validation_ok": true,
 "crossings": 1,
 "dots_too_close": 0,
 "lines_near_dots": 4,
 "label_overlaps": 0,
 "titles_on_lines": 0
}
```

## Area
- before packing 17428 x 10083 = 175.72 M
- after packing  17428 x 10083 = 175.72 M
- canvas 17728 x 10406 (margin incl.)

## Packing: mode='keep' (pieces at hand positions); overlapping pieces: 0

## Moves of detached pieces

## Title nudges