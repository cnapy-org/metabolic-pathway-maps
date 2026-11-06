# Phase 0 report

- reactions 49, nodes 302 {'metabolite': 155, 'midmarker': 49, 'multimarker': 98}, text labels 23
- metabolite nodes: 59 primary, 96 secondary; canvas {'x': -475, 'y': -200, 'width': 6056, 'height': 3523}

## Titles
- `663` '3-phenylpropionate degradation' at (2397, 0)
- `664` 'L-lactaldehyde degradation' at (94, 0)
- `665` 'allantoin degradation' at (1630, 2546)
- `666` 'creatinine degradation' at (4765, 0)
- `667` 'methylglyoxal degradation' at (4409, 1508)
- `668` 'octane oxidation' at (4922, 621)
- `669` 'phenylacetate degradation (aerobic)' at (2509, 1273)
- `670` 'propanol degradation' at (4340, 2443)
- `671` 'propionate fermentation' at (618, 1273)

## Subsystem proposal (primary-node components -> nearest title). CHECK IT.
- **propionate fermentation** (`671`, dist 94) n=11: MCITL2#455, MICITDr#466, MMM#476, PPCSCT#488, FRD2#501, FRD3#513, SUCDi#525, MCITD#536, MMCD#546, MMCD#557, MCITS#571
- **3-phenylpropionate degradation** (`663`, dist 94) n=9: 3HPPPNH#18, HPPPNDO#29, HKNDDH#42, OP4ENH#53, HOPNTAL#64, ACALD#79, PPPNDO#95, DHPPD#107, HKNTDH#597
- **phenylacetate degradation (aerobic)** (`669`, dist 94) n=9: REPHACCOAI#330, PACCOAE#347, PACCOAL#364, OXDHCOAT#378, DHACOAH#389, HADPCOADH3#402, 3OXCOAT#415, OXCOAHDH#636, PACOAT#648
- **allantoin degradation** (`665`, dist 94) n=6: ALLTN#180, ALLTAMH2#195, UGCIAMH#208, UGLYCH#223, URDGLYCD#236, OXAMTC#622
- **methylglyoxal degradation** (`667`, dist 94) n=4: GLYOX3#266, LGTHL#277, GLYOX#289, LDH_D2#302
- **L-lactaldehyde degradation** (`664`, dist 49) n=3: LCARS#135, LCADi#150, LCARSyi#609
- **octane oxidation** (`668`, dist 99) n=2: FACOAL80t2pp#320, FACOAE80#584
- **octane oxidation** (`668`, dist 742, LOW CONFIDENCE) n=1: 3HBZCT#121
- **allantoin degradation** (`665`, dist 1363, LOW CONFIDENCE) n=1: URIC#168
- **creatinine degradation** (`666`, dist 98) n=1: SARCOX#254
- **propanol degradation** (`670`, dist 121) n=1: PTA2#429
- **propanol degradation** (`670`, dist 113) n=1: ACACCT#443

## Components over all shared nodes (draft as is): [11, 9, 9, 6, 4, 3, 2, 1, 1, 1, 1, 1]

## Shared secondary nodes (split automatically, one per reaction)

## Duplicate reactions (keep; draw as lanes)
- MMCD: ['546', '557']

## Families between the same primary metabolites
- ['12ppd__S_c', 'lald__L_c']: [('LCARS', '135'), ('LCARSyi', '609')]
- ['pac_c', 'phaccoa_c']: [('PACCOAL', '364'), ('PACOAT', '648')]
- ['fum_c', 'succ_c']: [('FRD2', '501'), ('FRD3', '513'), ('SUCDi', '525')]
- ['mmcoa__S_c', 'ppcoa_c']: [('MMCD', '546'), ('MMCD', '557')]

## Primary metabolites drawn as several nodes (merge candidates)
- succoa_c: ['405', '467']
- ppcoa_c: ['417', '477']
- mmcoa__S_c: ['468', '547']
- fum_c: ['489', '587']

## Primaries > 900 units from their reaction
- MCITL2#455: micit_c `444` 1192

## In stoichiometry but not drawn (usually protons; leave as is)
- 3HPPPNH#18: ['h_c']
- HPPPNDO#29: ['h_c']
- HKNDDH#42: ['h_c']
- ACALD#79: ['h_c']
- PPPNDO#95: ['h_c']
- DHPPD#107: ['h_c']
- LCARS#135: ['h_c']
- LCADi#150: ['h_c']
- ALLTN#180: ['h_c']
- ALLTAMH2#195: ['h_c']
- UGLYCH#223: ['h_c']
- URDGLYCD#236: ['h_c']
- GLYOX3#266: ['h_c']
- GLYOX#289: ['h_c']
- FACOAL80t2pp#320: ['h_c', 'h_p']
- PACCOAE#347: ['h_c']
- HADPCOADH3#402: ['h_c']
- MMCD#546: ['h_c']
- MMCD#557: ['h_c']
- MCITS#571: ['h_c']
- FACOAE80#584: ['h_c']
- HKNTDH#597: ['h_c']
- LCARSyi#609: ['h_c']
- OXCOAHDH#636: ['h_c']
- PACOAT#648: ['h_c']

## Exit labels -> owner node (row rule)
- `649` '→ Carbohydrates' -> lald__L_c `122`
- `650` '→ Nucleotides' -> urate_c `153`
- `651` '→ Amino acids' -> gly_c `241`
- `652` '→ Amino acids' -> mthgxl_c `256`
- `653` '→ Nucleotides' -> pyr_c `290`
- `654` '→ Nucleotides' -> octa_c `574`
- `655` '→ Lipids' -> occoa_c `307`
- `656` '→ Central' -> succoa_c `405`
- `657` '→ Central' -> succ_c `446`
- `658` '→ Central' -> succoa_c `467`
- `659` '→ Central' -> fum_c `489`
- `660` '→ Central' -> oaa_c `559`
- `661` '→ Amino acids' -> fum_c `587`
- `662` '→ Nucleotides' -> cbp_c `611`

## Reactions
- `18` 3HPPPNH: 3hpppn_c + h_c + nadh_c + o2_c -> dhpppn_c + h2o_c + nad_c
- `29` HPPPNDO: dhpppn_c + o2_c -> h_c + hkndd_c
- `42` HKNDDH: h2o_c + hkndd_c -> h_c + op4en_c + succ_c
- `53` OP4ENH: h2o_c + op4en_c -> 4h2opntn_c
- `64` HOPNTAL: 4h2opntn_c -> acald_c + pyr_c
- `79` ACALD: acald_c + coa_c + nad_c <=> accoa_c + h_c + nadh_c
- `95` PPPNDO: h_c + nadh_c + o2_c + pppn_c -> cechddd_c + nad_c
- `107` DHPPD: cechddd_c + nad_c -> dhpppn_c + h_c + nadh_c
- `121` 3HBZCT: 3hbzcoa_c + h2o_c -> 3hbz_c + coa_c
- `135` LCARS: h_c + lald__L_c + nadh_c <=> 12ppd__S_c + nad_c
- `150` LCADi: h2o_c + lald__L_c + nad_c -> 2 h_c + lac__L_c + nadh_c
- `168` URIC: 2 h2o_c + o2_c + urate_c -> alltn_c + co2_c + h2o2_c
- `180` ALLTN: alltn_c + h2o_c -> alltt_c + h_c
- `195` ALLTAMH2: alltt_c + h2o_c + 2 h_c -> co2_c + nh4_c + urdgci_c
- `208` UGCIAMH: h2o_c + urdgci_c -> nh4_c + urdglyc_c
- `223` UGLYCH: h2o_c + 2 h_c + urdglyc_c -> co2_c + glx_c + 2 nh4_c
- `236` URDGLYCD: nad_c + urdglyc_c -> h_c + nadh_c + oxur_c
- `254` SARCOX: h2o_c + o2_c + sarcs_c -> fald_c + gly_c + h2o2_c
- `266` GLYOX3: h2o_c + mthgxl_c -> h_c + lac__D_c
- `277` LGTHL: gthrd_c + mthgxl_c -> lgt__S_c
- `289` GLYOX: h2o_c + lgt__S_c -> gthrd_c + h_c + lac__D_c
- `302` LDH_D2: lac__D_c + q8_c -> pyr_c + q8h2_c
- `320` FACOAL80t2pp: atp_c + coa_c + h_p + octa_p -> amp_c + h_c + occoa_c + ppi_c
- `330` REPHACCOAI: rephaccoa_c <=> 2oxpaccoa_c
- `347` PACCOAE: h_c + nadph_c + o2_c + phaccoa_c -> h2o_c + nadp_c + rephaccoa_c
- `364` PACCOAL: atp_c + coa_c + pac_c -> amp_c + phaccoa_c + ppi_c
- `378` OXDHCOAT: 3oxdhscoa_c + coa_c -> 23dhacoa_c + accoa_c
- `389` DHACOAH: 23dhacoa_c + h2o_c <=> 3hadpcoa_c
- `402` HADPCOADH3: 3hadpcoa_c + nad_c <=> h_c + nadh_c + oxadpcoa_c
- `415` 3OXCOAT: coa_c + oxadpcoa_c -> accoa_c + succoa_c
- `429` PTA2: pi_c + ppcoa_c -> coa_c + ppap_c
- `443` ACACCT: acac_c + accoa_c -> aacoa_c + ac_c
- `455` MCITL2: micit_c <=> pyr_c + succ_c
- `466` MICITDr: 2mcacn_c + h2o_c <=> micit_c
- `476` MMM: succoa_c -> mmcoa__S_c
- `488` PPCSCT: ppcoa_c + succ_c -> ppa_c + succoa_c
- `501` FRD2: fum_c + mql8_c -> mqn8_c + succ_c
- `513` FRD3: 2dmmql8_c + fum_c -> 2dmmq8_c + succ_c
- `525` SUCDi: q8_c + succ_c -> fum_c + q8h2_c
- `536` MCITD: 2mcit_c -> 2mcacn_c + h2o_c
- `546` MMCD: h_c + mmcoa__S_c -> co2_c + ppcoa_c
- `557` MMCD: h_c + mmcoa__S_c -> co2_c + ppcoa_c
- `571` MCITS: h2o_c + oaa_c + ppcoa_c -> 2mcit_c + coa_c + h_c
- `584` FACOAE80: h2o_c + occoa_c -> coa_c + h_c + octa_c
- `597` HKNTDH: h2o_c + hkntd_c -> fum_c + h_c + op4en_c
- `609` LCARSyi: h_c + lald__L_c + nadph_c -> 12ppd__S_c + nadp_c
- `622` OXAMTC: oxur_c + pi_c -> cbp_c + oxam_c
- `636` OXCOAHDH: 2oxpaccoa_c + 2 h2o_c + nadp_c -> 3oxdhscoa_c + 2 h_c + nadph_c
- `648` PACOAT: h2o_c + phaccoa_c -> coa_c + h_c + pac_c