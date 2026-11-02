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
- after packing  12051 x 9301 = 112.09 M
- canvas 12215 x 9378 (margin incl.)

## Moves of detached pieces
- ['3HAD140', '3HAD160', '3HAD180', '3OAR160', '3OAR180', '3OAS160', '3OAS180', 'AACPS3', 'EAR140x', 'EAR140y', 'EAR160x', 'EAR160y', 'EAR180x', 'FA160ACPHi']: dx=0 dy=8
- ['EDTXS1', 'EDTXS2', 'EDTXS3', 'FA161ACPHi', 'HEPT1', 'LPADSS', 'MOAT', 'MOAT2', 'TDSK', 'U23GAAT', 'UAGAAT', 'UHGADA', 'USHD']: dx=13 dy=6
- ['ACPPAT160', 'AGPAT160', 'APG3PAT160', 'APH160', 'FA160ACPHi', 'G3PAT160']: dx=-17 dy=-6
- ['CDPMEK', 'DMATT', 'DMPPS', 'DXPRIi', 'DXPS', 'GRTT', 'IPDDI', 'IPDPS', 'MECDPDH5', 'MECDPS', 'MEPCT']: dx=-2 dy=14
- ['3OAS181', 'AACPS4', 'APG3PAT161', 'APH161', 'G3PAT161']: dx=3 dy=-6
- ['3OAR40', 'ACCOAC', 'ACOATA', 'KAS14', 'KAS15', 'MACPD', 'MCOATA', 'OGMEACPS']: dx=15 dy=8
- ['AACPS5', 'APG3PAT181', 'APH181', 'G3PAT181']: dx=-13 dy=-6
- ['ACACT8r', 'ACOAD7f', 'ACOAD7f', 'CTECOAI7', 'FACOAE160', 'FACOAE161', 'FACOAL160t2pp', 'FACOAL161t2pp']: dx=10 dy=13
- ['3HAD60', '3HAD80', '3OAR60', '3OAR80', '3OAS60', '3OAS80', 'EAR60x', 'EAR60y']: dx=-1 dy=8
- ['LPLIPAL1A160pp', 'PLIPA1A160pp', 'PLIPA2A160pp']: dx=-5 dy=-6
- ['ALCD19', 'G3PD2', 'G3PT', 'GAPP']: dx=9 dy=8
- ['3HAD100', '3HAD100', '3OAR100', '3OAS100']: dx=-1 dy=-7
- ['LPLIPAL1E161pp', 'PLIPA1E161pp', 'PLIPA2E161pp']: dx=13 dy=4
- ['LPLIPAL1A181pp', 'PLIPA1A181pp', 'PLIPA2A181pp']: dx=-13 dy=-6
- ['LPLIPAL1G181pp', 'PLIPA1G181pp', 'PLIPA2G181pp']: dx=-13 dy=14
- ['3HAD120', '3OAR120', '3OAS120']: dx=13 dy=8
- ['LPLIPAL1A161pp', 'PLIPA1A161pp', 'PLIPA2A161pp']: dx=13 dy=-6
- ['LPLIPAL1G161pp', 'PLIPA1G161pp', 'PLIPA2G161pp']: dx=13 dy=14
- ['LPLIPAL1G160pp', 'PLIPA1G160pp', 'PLIPA2G160pp']: dx=-7 dy=14
- ['LPLIPAL1E160pp', 'PLIPA1E160pp', 'PLIPA2E160pp']: dx=-7 dy=4
- ['LPLIPAL1E181pp', 'PLIPA1E181pp', 'PLIPA2E181pp']: dx=-7 dy=4
- ['ACGAMT']: dx=7 dy=-12
- ['3OAR121', '3OAS121']: dx=-19 dy=8
- ['MOAT3C']: dx=-7 dy=14
- ['FACOAE181', 'FACOAL181t2pp']: dx=-3 dy=3
- ['PETNT161pp', 'PETNT181pp']: dx=-13 dy=8
- ['2AGPGAT181', 'LPLIPAL2G181']: dx=-13 dy=14
- ['ACACT1r']: dx=9 dy=-6
- ['2AGPEAT181', 'LPLIPAL2E181']: dx=15 dy=4
- ['2AGPGAT160', 'LPLIPAL2G160']: dx=-19 dy=14
- ['2AGPGAT161', 'LPLIPAL2G161']: dx=1 dy=14
- ['2AGPEAT161', 'LPLIPAL2E161']: dx=1 dy=4
- ['2AGPEAT160', 'LPLIPAL2E160']: dx=-19 dy=4
- ['LIPAHT2ex']: dx=-5 dy=8
- ['ACACT4r']: dx=11 dy=4
- ['LIPAHTex']: dx=-13 dy=-2
- ['TDPADGAT']: dx=-5 dy=18
- ['LPLIPAL2A181']: dx=17 dy=-6
- ['LPLIPAL2A161']: dx=1 dy=-6
- ['LPLIPAL2A160']: dx=-19 dy=-6

## Title nudges