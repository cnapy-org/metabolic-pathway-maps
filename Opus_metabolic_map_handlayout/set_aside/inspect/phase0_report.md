# Phase 0 report

- reactions 45, nodes 269 {'metabolite': 134, 'midmarker': 45, 'multimarker': 90}, text labels 8
- metabolite nodes: 55 primary, 79 secondary; canvas {'x': -93, 'y': -200, 'width': 5214, 'height': 3177}

## Titles
- `581` 'C4 and CAM-carbon fixation' at (3262, 2112)
- `582` 'acetyl CoA biosynthesis' at (2562, 0)
- `583` 'aerobic toluene degradation' at (3985, 0)
- `584` 'butanoate fermentation' at (2760, 1119)
- `585` 'methylaspartate cycle' at (957, 0)
- `586` 'photosynthesis' at (768, 1577)

## Subsystem proposal (primary-node components -> nearest title). CHECK IT.
- **methylaspartate cycle** (`585`, dist 94) n=11: MDH#300, CS#315, FUM#326, FRD2#339, FRD3#351, SUCDi#363, ACONTa#374, SUCOAS#391, ACONTb#402, MMM#411, ICDHyr#426
- **C4 and CAM-carbon fixation** (`581`, dist 94) n=10: ME2#16, ALATA_L#29, ALATA_L#42, PPCK#74, PYK#86, MDH#127, ASPTA#141, ASPTA#153, MALDDH#168, ME1#179
- **photosynthesis** (`586`, dist 94) n=10: RPI#436, RBP4E#460, RPE#468, X5PL3E#477, TKT2#506, FBP#519, TKT1#531, FBA#558, TPI#567, FBA3#578
- **aerobic toluene degradation** (`583`, dist 94) n=4: OP4ENH#209, HOPNTAL#220, ACALD#235, 3OXCOAT#248
- **C4 and CAM-carbon fixation** (`581`, dist 622, LOW CONFIDENCE) n=3: PPC#58, MDH#99, ME2#114
- **acetyl CoA biosynthesis** (`582`, dist 102) n=1: PDH#197
- **butanoate fermentation** (`584`, dist 253) n=1: ACACT1r#260
- **butanoate fermentation** (`584`, dist 326, LOW CONFIDENCE) n=1: ECOAH1#272
- **butanoate fermentation** (`584`, dist 495, LOW CONFIDENCE) n=1: BUTCT#286
- **photosynthesis** (`586`, dist 413, LOW CONFIDENCE) n=1: PGK#450
- **photosynthesis** (`586`, dist 413, LOW CONFIDENCE) n=1: PPC#493
- **photosynthesis** (`586`, dist 655, LOW CONFIDENCE) n=1: ME2#547

## Components over all shared nodes (draft as is): [11, 10, 10, 4, 3, 1, 1, 1, 1, 1, 1, 1]

## Shared secondary nodes (split automatically, one per reaction)
- nad_c `155`: ['MALDDH', 'ME1']
- co2_c `156`: ['MALDDH', 'ME1']
- nadh_c `157`: ['MALDDH', 'ME1']

## Duplicate reactions (keep; draw as lanes)
- ME2: ['16', '114', '547']
- ALATA_L: ['29', '42']
- PPC: ['58', '493']
- MDH: ['99', '127', '300']
- ASPTA: ['141', '153']

## Families between the same primary metabolites
- ['mal__L_c', 'pyr_c']: [('ME2', '16'), ('ME2', '114'), ('ME1', '179'), ('ME2', '547')]
- ['ala__L_c', 'pyr_c']: [('ALATA_L', '29'), ('ALATA_L', '42')]
- ['oaa_c', 'pep_c']: [('PPC', '58'), ('PPCK', '74'), ('PPC', '493')]
- ['mal__L_c', 'oaa_c']: [('MDH', '99'), ('MDH', '127'), ('MDH', '300')]
- ['asp__L_c', 'oaa_c']: [('ASPTA', '141'), ('ASPTA', '153')]
- ['fum_c', 'succ_c']: [('FRD2', '339'), ('FRD3', '351'), ('SUCDi', '363')]

## Primary metabolites drawn as several nodes (merge candidates)
- mal__L_c: ['1', '87', '115', '287', '532']
- pyr_c: ['5', '32', '103', '182', '536']
- pep_c: ['45', '63', '480']
- oaa_c: ['46', '60', '131', '290', '481']
- accoa_c: ['183', '223', '249', '301']

## Primaries > 900 units from their reaction

## In stoichiometry but not drawn (usually protons; leave as is)
- PPC#58: ['h_c']
- PYK#86: ['h_c']
- MDH#99: ['h_c']
- MDH#127: ['h_c']
- ACALD#235: ['h_c']
- MDH#300: ['h_c']
- CS#315: ['h_c']
- PPC#493: ['h_c']

## Exit labels -> owner node (row rule)
- `579` '→ Central' -> accoa_c `249`
- `580` '→ Central' -> accoa_c `301`

## Reactions
- `16` ME2: mal__L_c + nadp_c -> co2_c + nadph_c + pyr_c
- `29` ALATA_L: akg_c + ala__L_c <=> glu__L_c + pyr_c
- `42` ALATA_L: akg_c + ala__L_c <=> glu__L_c + pyr_c
- `58` PPC: co2_c + h2o_c + pep_c -> h_c + oaa_c + pi_c
- `74` PPCK: atp_c + oaa_c -> adp_c + co2_c + pep_c
- `86` PYK: adp_c + h_c + pep_c -> atp_c + pyr_c
- `99` MDH: mal__L_c + nad_c <=> h_c + nadh_c + oaa_c
- `114` ME2: mal__L_c + nadp_c -> co2_c + nadph_c + pyr_c
- `127` MDH: mal__L_c + nad_c <=> h_c + nadh_c + oaa_c
- `141` ASPTA: akg_c + asp__L_c <=> glu__L_c + oaa_c
- `153` ASPTA: akg_c + asp__L_c <=> glu__L_c + oaa_c
- `168` MALDDH: mal__D_c + nad_c -> co2_c + nadh_c + pyr_c
- `179` ME1: mal__L_c + nad_c -> co2_c + nadh_c + pyr_c
- `197` PDH: coa_c + nad_c + pyr_c -> accoa_c + co2_c + nadh_c
- `209` OP4ENH: h2o_c + op4en_c -> 4h2opntn_c
- `220` HOPNTAL: 4h2opntn_c -> acald_c + pyr_c
- `235` ACALD: acald_c + coa_c + nad_c <=> accoa_c + h_c + nadh_c
- `248` 3OXCOAT: coa_c + oxadpcoa_c -> accoa_c + succoa_c
- `260` ACACT1r: 2 accoa_c <=> aacoa_c + coa_c
- `272` ECOAH1: 3hbcoa_c <=> b2coa_c + h2o_c
- `286` BUTCT: accoa_c + but_c -> ac_c + btcoa_c
- `300` MDH: mal__L_c + nad_c <=> h_c + nadh_c + oaa_c
- `315` CS: accoa_c + h2o_c + oaa_c -> cit_c + coa_c + h_c
- `326` FUM: fum_c + h2o_c <=> mal__L_c
- `339` FRD2: fum_c + mql8_c -> mqn8_c + succ_c
- `351` FRD3: 2dmmql8_c + fum_c -> 2dmmq8_c + succ_c
- `363` SUCDi: q8_c + succ_c -> fum_c + q8h2_c
- `374` ACONTa: cit_c <=> acon_C_c + h2o_c
- `391` SUCOAS: atp_c + coa_c + succ_c <=> adp_c + pi_c + succoa_c
- `402` ACONTb: acon_C_c + h2o_c <=> icit_c
- `411` MMM: succoa_c -> mmcoa__S_c
- `426` ICDHyr: icit_c + nadp_c <=> akg_c + co2_c + nadph_c
- `436` RPI: r5p_c <=> ru5p__D_c
- `450` PGK: 3pg_c + atp_c <=> 13dpg_c + adp_c
- `460` RBP4E: ru5p__L_c <=> xu5p__D_c
- `468` RPE: ru5p__D_c <=> xu5p__D_c
- `477` X5PL3E: xu5p__L_c -> ru5p__L_c
- `493` PPC: co2_c + h2o_c + pep_c -> h_c + oaa_c + pi_c
- `506` TKT2: e4p_c + xu5p__D_c <=> f6p_c + g3p_c
- `519` FBP: fdp_c + h2o_c -> f6p_c + pi_c
- `531` TKT1: r5p_c + xu5p__D_c <=> g3p_c + s7p_c
- `547` ME2: mal__L_c + nadp_c -> co2_c + nadph_c + pyr_c
- `558` FBA: fdp_c <=> dhap_c + g3p_c
- `567` TPI: dhap_c <=> g3p_c
- `578` FBA3: s17bp_c <=> dhap_c + e4p_c