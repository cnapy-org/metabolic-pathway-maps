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
- pri: n=368 clean=367 line=1 overlap=0
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
- before packing 18315 x 11490 = 210.43 M
- after packing  13000 x 7040 = 91.52 M
- canvas 13304 x 7344 (margin incl.)

## Moves of detached pieces
- ['ASPCT', 'CSND', 'CTPS2', 'CYTD', 'CYTDH', 'CYTDK2', 'CYTK1', 'CYTK2', 'DCTPD', 'DCYTD', 'DHORD2', 'DHORD5', 'DHORDfum', 'DHORTS', 'DTMPK', 'DURIK1', 'DURIPP', 'DUTPDP', 'NDPK2', 'NDPK3', 'NDPK4', 'NDPK6', 'NDPK7', 'NTD1', 'NTD3', 'NTD4', 'NTD5', 'NTPP3', 'NTPP4', 'NTPP7', 'OMPDC', 'ORPT', 'PYK6', 'PYNP2r', 'RNDR3', 'RNDR4', 'TMDK1', 'TMDPP', 'TMDS', 'UMPK', 'UMPK', 'UPPRT', 'URIH', 'URIK2']: dx=-8507 dy=-766
- ['4ABZGLUH', 'ACP1_FMN', 'ADCL', 'ADCS', 'APRAUR', 'CCGS', 'CDGS', 'CPH4S', 'DB4PS', 'DHFR', 'DHFS', 'DHMPTR', 'DHNPA2r', 'DHNPTE', 'DHPPDA2', 'DHPS2', 'DHPTPE', 'DNMPPA', 'DNTPPA', 'FADRx2', 'FLVR', 'FLVRx', 'FMNAT', 'FTHFD', 'GARFT', 'GTPCI', 'GTPCI', 'GTPCII2', 'GTPCII2', 'HPPK2', 'PMDPHT', 'RBFK', 'RBFSa', 'RBFSb']: dx=-7844 dy=-1308
- ['ASPO3', 'ASPO4', 'ASPO6', 'NADDP', 'NADK', 'NADN', 'NADPPPS', 'NADS1', 'NADTRHD', 'NADTRHD', 'NADTRHD', 'NAMNPP', 'NMNAT', 'NMNDA', 'NMNN', 'NMNt7pp', 'NNAM', 'NNATr', 'NNDPR', 'QULNS', 'THD2pp']: dx=-7298 dy=-284
- ['4HTHRK', 'ALATA_D2', 'ALATA_L2', 'E4PD', 'HYPOE', 'OHPBAT', 'PDX5PO2', 'PDX5POi', 'PDX5PS', 'PDXPP', 'PERD', 'PYAM5PO', 'PYDAMK', 'PYDXK', 'PYDXNK', 'PYDXPP']: dx=-10453 dy=-2048
- ['CPPPGO', 'CPPPGO2', 'FCLT', 'G1SAT', 'GLUTRR', 'GLUTRS', 'HEMEOS', 'HMBS', 'PPBNGS', 'PPPGO', 'PPPGO3', 'UPP3MT', 'UPP3S', 'UPPDC1']: dx=-8476 dy=-178
- ['BMOCOS', 'CPMPS', 'MOCDS', 'MOCOS', 'MPTAT', 'MPTS']: dx=-14325 dy=-4216
- ['DXPS', 'HETZK', 'TDP', 'THZPSN3', 'TMK', 'TMPK', 'TMPPP', 'TYRL']: dx=-20368 dy=-2724
- ['DPCOAK', 'DPR', 'MOHMT', 'PANTS', 'PNTK', 'PPCDC', 'PPNCL2', 'PTPATi', 'TPRDCOAS']: dx=-12202 dy=-8198
- ['AMMQLT8', 'DHNAOT4', 'DHNCOAS', 'DHNCOAT', 'ICHORS_copy1', 'ICHORS_copy2', 'SEPHCHCS', 'SHCHCS3', 'SUCBZL', 'SUCBZS']: dx=-10079 dy=-1624
- ['GHMT2r', 'MTHFD', 'MTHFR2', 'SERD_D', 'SERD_L']: dx=-11655 dy=-9828
- ['ACBIPGT', 'ADOCBIK', 'ADOCBLS', 'NNDMBRT', 'RZ5PP']: dx=-3661 dy=-3344
- ['AMAOTr', 'AOXSr2', 'BTS5', 'DBTS']: dx=-5242 dy=-8078
- ['ADPRDP', 'CMPN', 'PPM', 'YUMPS']: dx=1098 dy=-3508
- ['3KGK', 'ASCBPL', 'DOGULNR', 'KG6PDC']: dx=-8097 dy=-3084
- ['SHCHD2', 'SHCHD2', 'SHCHF', 'UPP3MT', 'UPP3MT']: dx=-9227 dy=-3018
- ['ARBTNR2', 'ARBTNR3']: dx=-229 dy=-6298
- ['CPGNR2', 'CPGNR3']: dx=-2952 dy=-5018
- ['FEENTERR2', 'FEENTERR3']: dx=-5632 dy=-3648
- ['FECRMR2', 'FECRMR3']: dx=3488 dy=-8328
- ['FE3HOXR2', 'FE3HOXR3']: dx=-5552 dy=-8938
- ['FEOXAMR2', 'FEOXAMR3']: dx=-7752 dy=-7280
- ['FMNRx', 'FMNRx2']: dx=-3980 dy=-7511
- ['NTD9pp']: dx=879 dy=-7000
- ['GDPDPK', 'GTPDPDP', 'GTPDPK', 'PPGPPDP']: dx=-3730 dy=-7871
- ['AACPS9', 'OCTNLL']: dx=-7854 dy=-5707
- ['BWCOS', 'WCOS']: dx=-14339 dy=-5350
- ['FDMO']: dx=-6634 dy=-3360
- ['RHAT1']: dx=-2471 dy=250
- ['FDMO2']: dx=1177 dy=-6650
- ['FDMO3']: dx=-2263 dy=-9460
- ['FDMO6']: dx=-9023 dy=-10800
- ['FDMO4']: dx=-3253 dy=-3470
- ['AADDGT']: dx=-4404 dy=-4400
- ['MTRPOX']: dx=-709 dy=-5420
- ['ACKr', 'PTAr']: dx=-8885 dy=-9490
- ['ICYSDS', 'SCYSDS']: dx=366 dy=-6410
- ['5DOAN']: dx=1705 dy=-220
- ['LIPAMPL', 'LIPATPT']: dx=-9824 dy=-9150
- ['NTD5pp']: dx=320 dy=-5720
- ['GAL1PPpp']: dx=-1720 dy=-6160
- ['NTD7pp']: dx=-760 dy=-5160
- ['NTD3pp']: dx=-2640 dy=-5580
- ['NTD10pp']: dx=-1640 dy=-5280
- ['NTD11pp']: dx=-800 dy=-5560
- ['NTD4pp']: dx=2440 dy=-6300
- ['ACP1p']: dx=1840 dy=-6820
- ['NTD1pp']: dx=640 dy=-7060
- ['NTD6pp']: dx=-5080 dy=-6280
- ['NTD8pp']: dx=-480 dy=-7880

## Title nudges