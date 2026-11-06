# Phase 0 report

- reactions 42, nodes 264 {'metabolite': 138, 'midmarker': 42, 'multimarker': 84}, text labels 24
- metabolite nodes: 54 primary, 84 secondary; canvas {'x': -120, 'y': -200, 'width': 5440, 'height': 3168}

## Titles
- `587` 'CO2 fixation in Crenarchaeota' at (2530, 0)
- `588` 'formaldehyde oxidation' at (4204, 1437)
- `589` 'glycolate and glyoxylate degradation' at (2583, 1925)
- `590` 'methane metabolism' at (80, 2085)
- `591` 'reductive acetyl coenzyme A pathway' at (4140, 0)
- `592` 'ubiquinone biosynthesis' at (702, 0)
- `593` 'urea cycle' at (1288, 2085)

## Subsystem proposal (primary-node components -> nearest title). CHECK IT.
- **ubiquinone biosynthesis** (`592`, dist 112) n=8: OCTDPS#281, OPHBDC#293, HBZOPT#305, CHRPL#316, OPHHXy#333, OHPHM#346, OMPHHXy#537, OPHHX3#557
- **glycolate and glyoxylate degradation** (`589`, dist 104) n=7: GLYCLTDy#190, GLYCLTDx#202, GLXCL#213, DHGLYH#447, GLYCTO2#459, GLYCTO3#471, GLYCTO4#483
- **CO2 fixation in Crenarchaeota** (`587`, dist 466, LOW CONFIDENCE) n=6: ECOAH1#50, ECOAH1#59, HACD1#72, HACD1#82, ACACT1r#93, ACCOAC#124
- **formaldehyde oxidation** (`588`, dist 94) n=4: FALGTHLs#136, FALDH2#149, SFGTHi#162, AMPMS2#434
- **reductive acetyl coenzyme A pathway** (`591`, dist 94) n=4: FTHFLi#243, MTHFC#254, MTHFD#267, ULA4NFT#569
- **urea cycle** (`593`, dist 35) n=4: ARGSL#372, OCBT#386, ARGSS#402, CBMKr#417
- **ubiquinone biosynthesis** (`592`, dist 85) n=3: DMQMT#360, OMMBLHX3#504, OMMBLHXy#520
- **CO2 fixation in Crenarchaeota** (`587`, dist 94) n=1: MMM#10
- **CO2 fixation in Crenarchaeota** (`587`, dist 368, LOW CONFIDENCE) n=1: GHBDHx#24
- **CO2 fixation in Crenarchaeota** (`587`, dist 609, LOW CONFIDENCE) n=1: PCNO#38
- **glycolate and glyoxylate degradation** (`589`, dist 946, LOW CONFIDENCE) n=1: MSAR#107
- **glycolate and glyoxylate degradation** (`589`, dist 94) n=1: GLYCK#176
- **methane metabolism** (`590`, dist 122) n=1: CAT#225

## Components over all shared nodes (draft as is): [8, 7, 6, 4, 4, 4, 3, 1, 1, 1, 1, 1, 1]

## Shared secondary nodes (split automatically, one per reaction)
- h2o_c `41`: ['ECOAH1', 'ECOAH1']
- nadh_c `61`: ['HACD1', 'HACD1']
- nad_c `62`: ['HACD1', 'HACD1']

## Duplicate reactions (keep; draw as lanes)
- ECOAH1: ['50', '59']
- HACD1: ['72', '82']

## Families between the same primary metabolites
- ['3hbcoa_c', 'b2coa_c']: [('ECOAH1', '50'), ('ECOAH1', '59')]
- ['3hbcoa_c', 'aacoa_c']: [('HACD1', '72'), ('HACD1', '82')]
- ['glx_c', 'glyclt_c']: [('GLYCLTDy', '190'), ('GLYCLTDx', '202'), ('GLYCTO2', '459'), ('GLYCTO3', '471'), ('GLYCTO4', '483')]
- ['2ohph_c', '2oph_c']: [('OPHHXy', '333'), ('OPHHX3', '557')]
- ['2omhmbl_c', '2ommbl_c']: [('OMMBLHX3', '504'), ('OMMBLHXy', '520')]

## Primary metabolites drawn as several nodes (merge candidates)
- for_c: ['151', '227']

## Primaries > 900 units from their reaction

## In stoichiometry but not drawn (usually protons; leave as is)
- GHBDHx#24: ['h_c']
- PCNO#38: ['h_c']
- HACD1#72: ['h_c']
- HACD1#82: ['h_c']
- MSAR#107: ['h_c']
- ACCOAC#124: ['h_c']
- FALDH2#149: ['h_c']
- SFGTHi#162: ['h_c']
- GLYCK#176: ['h_c']
- GLYCLTDy#190: ['h_c']
- GLYCLTDx#202: ['h_c']
- GLXCL#213: ['h_c']
- MTHFC#254: ['h_c']
- OPHBDC#293: ['h_c']
- OPHHXy#333: ['h_c']
- OHPHM#346: ['h_c']
- DMQMT#360: ['h_c']
- OCBT#386: ['h_c']
- ARGSS#402: ['h_c']
- CBMKr#417: ['h_c']
- AMPMS2#434: ['h_c']
- DHGLYH#447: ['h_c']
- OMMBLHX3#504: ['h_c']
- OMMBLHXy#520: ['h_c']
- OMPHHXy#537: ['h_c']
- OPHHX3#557: ['h_c']
- ULA4NFT#569: ['h_c']

## Exit labels -> owner node (row rule)
- `570` '→ Amino acids' -> malcoa_c `111`
- `571` '→ Nucleotides' -> fald_c `125`
- `572` '→ Central' -> for_c `151`
- `573` '→ Carbohydrates' -> glyc__R_c `164`
- `574` '→ Carbohydrates' -> glyclt_c `179`
- `575` '→ Carbohydrates' -> 2h3oppan_c `203`
- `576` '→ Nucleotides' -> frdp_c `268`
- `577` '→ Amino acids' -> chor_c `306`
- `578` '→ Amino acids' -> argsuc_c `361`
- `579` '→ Amino acids' -> cbp_c `373`
- `580` '→ Amino acids' -> orn_c `374`
- `581` '→ Nucleotides' -> air_c `418`
- `582` '→ Nucleotides' -> dhgly_c `435`
- `583` '→ Amino acids' -> 2ommbl_c `484`
- `584` '→ Amino acids' -> 2ombzl_c `523`
- `585` '→ Amino acids' -> for_c `227`
- `586` '→ Amino acids' -> udpLa4fn_c `559`

## Reactions
- `10` MMM: succoa_c -> mmcoa__S_c
- `24` GHBDHx: h_c + nadh_c + sucsal_c <=> ghb_c + nad_c
- `38` PCNO: nadp_c + ppcoa_c <=> h_c + nadph_c + prpncoa_c
- `50` ECOAH1: 3hbcoa_c <=> b2coa_c + h2o_c
- `59` ECOAH1: 3hbcoa_c <=> b2coa_c + h2o_c
- `72` HACD1: aacoa_c + h_c + nadh_c <=> 3hbcoa_c + nad_c
- `82` HACD1: aacoa_c + h_c + nadh_c <=> 3hbcoa_c + nad_c
- `93` ACACT1r: 2 accoa_c <=> aacoa_c + coa_c
- `107` MSAR: h_c + msa_c + nadph_c -> 3hpp_c + nadp_c
- `124` ACCOAC: accoa_c + atp_c + hco3_c -> adp_c + h_c + malcoa_c + pi_c
- `136` FALGTHLs: fald_c + gthrd_c <=> hmgth_c
- `149` FALDH2: hmgth_c + nad_c <=> Sfglutth_c + h_c + nadh_c
- `162` SFGTHi: Sfglutth_c + h2o_c -> for_c + gthrd_c + h_c
- `176` GLYCK: atp_c + glyc__R_c -> 3pg_c + adp_c + h_c
- `190` GLYCLTDy: glx_c + h_c + nadph_c -> glyclt_c + nadp_c
- `202` GLYCLTDx: glx_c + h_c + nadh_c -> glyclt_c + nad_c
- `213` GLXCL: 2 glx_c + h_c -> 2h3oppan_c + co2_c
- `225` CAT: 2 h2o2_c -> 2 h2o_c + o2_c
- `243` FTHFLi: atp_c + for_c + thf_c -> 10fthf_c + adp_c + pi_c
- `254` MTHFC: h2o_c + methf_c <=> 10fthf_c + h_c
- `267` MTHFD: mlthf_c + nadp_c <=> methf_c + nadph_c
- `281` OCTDPS: frdp_c + 5 ipdp_c -> octdp_c + 5 ppi_c
- `293` OPHBDC: 3ophb_c + h_c -> 2oph_c + co2_c
- `305` HBZOPT: 4hbz_c + octdp_c -> 3ophb_c + ppi_c
- `316` CHRPL: chor_c -> 4hbz_c + pyr_c
- `333` OPHHXy: 2oph_c + h_c + nadph_c + o2_c -> 2ohph_c + h2o_c + nadp_c
- `346` OHPHM: 2ohph_c + amet_c -> 2omph_c + ahcys_c + h_c
- `360` DMQMT: 2omhmbl_c + amet_c -> ahcys_c + h_c + q8h2_c
- `372` ARGSL: argsuc_c <=> arg__L_c + fum_c
- `386` OCBT: cbp_c + orn_c <=> citr__L_c + h_c + pi_c
- `402` ARGSS: asp__L_c + atp_c + citr__L_c -> amp_c + argsuc_c + h_c + ppi_c
- `417` CBMKr: atp_c + co2_c + nh4_c <=> adp_c + cbp_c + 2 h_c
- `434` AMPMS2: air_c + h2o_c + nad_c -> 4ampm_c + 2 for_c + 3 h_c + nadh_c
- `447` DHGLYH: dhgly_c + h2o_c + h_c -> glx_c + nh4_c
- `459` GLYCTO2: glyclt_c + q8_c -> glx_c + q8h2_c
- `471` GLYCTO3: glyclt_c + mqn8_c -> glx_c + mql8_c
- `483` GLYCTO4: 2dmmq8_c + glyclt_c -> 2dmmql8_c + glx_c
- `504` OMMBLHX3: 2ommbl_c + 2 atp_c + 3 h2o_c + nad_c -> 2omhmbl_c + 2 adp_c + 3 h_c + nadh_c + 2 pi_c
- `520` OMMBLHXy: 2ommbl_c + h_c + nadph_c + o2_c -> 2omhmbl_c + h2o_c + nadp_c
- `537` OMPHHXy: 2omph_c + h_c + nadph_c + o2_c -> 2ombzl_c + h2o_c + nadp_c
- `557` OPHHX3: 2oph_c + 2 atp_c + 3 h2o_c + nad_c -> 2ohph_c + 2 adp_c + 3 h_c + nadh_c + 2 pi_c
- `569` ULA4NFT: 10fthf_c + udpLa4n_c -> h_c + thf_c + udpLa4fn_c