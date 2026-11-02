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
- rxn: n=49 clean=49 line=0 overlap=0
- sec: n=95 clean=95 line=0 overlap=0
- pri: n=58 clean=58 line=0 overlap=0
- exit: n=14 clean=14 line=0 overlap=0

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
- before packing 7239 x 3428 = 24.82 M
- after packing  2840 x 5278 = 14.99 M
- canvas 3294 x 5596 (margin incl.)

## Moves of detached pieces
- ['3HPPPNH', 'ACALD', 'DHPPD', 'HKNDDH', 'HKNTDH', 'HOPNTAL', 'HPPPNDO', 'OP4ENH', 'PPPNDO']: dx=-2563 dy=-542
- ['ALLTAMH2', 'ALLTN', 'OXAMTC', 'UGCIAMH', 'UGLYCH', 'URDGLYCD', 'URIC']: dx=-2461 dy=-3704
- ['3OXCOAT', 'DHACOAH', 'HADPCOADH3', 'OXCOAHDH', 'OXDHCOAT', 'PACCOAE', 'PACCOAL', 'PACOAT', 'REPHACCOAI']: dx=-2582 dy=-3686
- ['GLYOX', 'GLYOX3', 'LDH_D2', 'LGTHL']: dx=2720 dy=-106
- ['LCADi', 'LCARS', 'LCARSyi']: dx=2372 dy=-2876
- ['FACOAE80', 'FACOAL80t2pp']: dx=-3171 dy=98
- ['3HBZCT']: dx=-3969 dy=-5202
- ['SARCOX']: dx=-4771 dy=-5928
- ['ACACCT']: dx=3007 dy=-4954

## Title nudges
- 3-phenylpropionate degradation: (-208, 120) ok
- propionate fermentation: (-41, -155) ok