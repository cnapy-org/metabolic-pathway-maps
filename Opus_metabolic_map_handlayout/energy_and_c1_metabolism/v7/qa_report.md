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

## Crossings between different reactions (curves sampled): 2
- GLYCLTDx × GLYCLTDy: 1
- GLYCTO2 × GLYCTO3: 1

## Metabolite dots closer than 45 units: 0

## Foreign lines within 28 units of a dot: 0

## Overlapping label boxes: 0

## Titles touching lines/nodes: 0

## Label placement scores (<10 clean, 10-119 touches a line, >=120 overlaps)
- rxn: n=42 clean=42 line=0 overlap=0
- sec: n=87 clean=87 line=0 overlap=0
- pri: n=53 clean=53 line=0 overlap=0
- exit: n=17 clean=17 line=0 overlap=0

## Summary
```
{
 "validation_ok": true,
 "crossings": 2,
 "dots_too_close": 0,
 "lines_near_dots": 0,
 "label_overlaps": 0,
 "titles_on_lines": 0
}
```

## Area
- before packing 6236 x 4329 = 26.99 M
- after packing  6501 x 4901 = 31.86 M
- canvas 6560 x 4699 (margin incl.)

## Moves of detached pieces
- ['ACACT1r', 'ACCOAC', 'ECOAH1', 'ECOAH1', 'HACD1', 'HACD1']: dx=-2 dy=16
- ['ARGSL', 'ARGSS', 'CBMKr', 'OCBT']: dx=5 dy=-12
- ['CHRPL', 'HBZOPT', 'OCTDPS', 'OHPHM', 'OMPHHXy', 'OPHBDC', 'OPHHX3', 'OPHHXy']: dx=-5 dy=4
- ['DHGLYH', 'GLXCL', 'GLYCLTDx', 'GLYCLTDy', 'GLYCTO2', 'GLYCTO3', 'GLYCTO4']: dx=7 dy=-6
- ['DMQMT', 'OMMBLHX3', 'OMMBLHXy']: dx=6 dy=4
- ['CAT']: dx=16 dy=4
- ['GHBDHx']: dx=171 dy=158
- ['GLYCK']: dx=-4 dy=18
- ['PCNO']: dx=1 dy=-2
- ['MSAR']: dx=-11 dy=16
- ['MMM']: dx=52 dy=9

## Title nudges
- CO2 fixation in Crenarchaeota (3-hydroxypropionate/4-hydroxybutyrate cycle): (-52, -193) ok
- methane metabolism: (-69, 40) ok