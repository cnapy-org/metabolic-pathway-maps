"""Spec for nucleotides_and_cofactors.escher.json - generated skeleton. Replace EVERY None by hand-chosen coordinates."""
import math

INPUT = 'nucleotides_and_cofactors.escher.json'

# zone -> title text (None = no title), existing text-label id to reuse, optional colour / anchor
ZONES = {
    'PM': dict(title='purine metabolism', title_id='4090'),
    'PM2': dict(title='pyrimidine metabolism', title_id='4091'),
    'NM': dict(title='NAD metabolism', title_id='4080'),
    'VBM': dict(title='vitamin B12 metabolism', title_id='4094'),
    'VBM2': dict(title='vitamin B6 metabolism', title_id='4095'),
    'HM': dict(title='heme metabolism', title_id='4085'),
    'FB': dict(title='flavin biosynthesis', title_id='4084'),
    'VKM': dict(title='vitamin K metabolism', title_id='4096'),
    'MCB': dict(title='molybdenum cofactor biosynthesis', title_id='4087'),
    'VBM3': dict(title='vitamin B1 metabolism', title_id='4093'),
    'TM': dict(title='tetrahydrofolate metabolism', title_id='4092'),
    'HDDB': dict(title='6-hydroxymethyl-dihydropterin diphosphate biosynthesis', title_id='4079'),
    'CAM': dict(title='coenzyme A metabolism', title_id='4083'),
    'AM': dict(title='ascorbate metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'BB': dict(title='biotin biosynthesis', title_id='4082'),
    'PB': dict(title='ppGpp biosynthesis', title_id='4089'),
    'PM3': dict(title='purine metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'VBM4': dict(title='vitamin B12 metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'VBM5': dict(title='vitamin B6 metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'PB2': dict(title='pantothenate biosynthesis', title_id='4088'),
    'LB': dict(title='lipoate biosynthesis', title_id='4086'),
    'LB2': dict(title='lipoate biosynthesis', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'PM4': dict(title='purine metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'VBM6': dict(title='vitamin B12 metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'AM2': dict(title='ascorbate metabolism', title_id='4081'),
    'PM5': dict(title='purine metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'PM6': dict(title='purine metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'PM7': dict(title='purine metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'PM8': dict(title='purine metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'FB2': dict(title='flavin biosynthesis', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'TM2': dict(title='tetrahydrofolate metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'VBM7': dict(title='vitamin B12 metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'FB3': dict(title='flavin biosynthesis', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'PM9': dict(title='purine metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'PM10': dict(title='purine metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'PM11': dict(title='purine metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
    'PM12': dict(title='purine metabolism', title_id=None),   # LOW CONFIDENCE: far from any title - check subsystem
}

EXIT_OWNER = {}          # automatic owner rule (see report); add {label_id: node_id} to override

def prepare(L):
    """splits of primary nodes: L.S[name] = L.split(node_id, reaction_id)"""
    pass

ORIGIN = {'PM': (None, None), 'PM2': (None, None), 'NM': (None, None), 'VBM': (None, None), 'VBM2': (None, None), 'HM': (None, None), 'FB': (None, None), 'VKM': (None, None), 'MCB': (None, None), 'VBM3': (None, None), 'TM': (None, None), 'HDDB': (None, None), 'CAM': (None, None), 'AM': (None, None), 'BB': (None, None), 'PB': (None, None), 'PM3': (None, None), 'VBM4': (None, None), 'VBM5': (None, None), 'PB2': (None, None), 'LB': (None, None), 'LB2': (None, None), 'PM4': (None, None), 'VBM6': (None, None), 'AM2': (None, None), 'PM5': (None, None), 'PM6': (None, None), 'PM7': (None, None), 'PM8': (None, None), 'FB2': (None, None), 'TM2': (None, None), 'VBM7': (None, None), 'FB3': (None, None), 'PM9': (None, None), 'PM10': (None, None), 'PM11': (None, None), 'PM12': (None, None)}
FLIP = {}

MERGES = [   # (keep, [duplicates])  - key metabolites only
]

def layout(L):
    # ---------------------------------------------------------------- PM: purine metabolism
    L.setp('PM', {
        '1053': None,   # 2dr1p_c  (draft 2493, 5286)
        '1054': None,   # 2dr5p_c  (draft 2054, 5242)
        '1063': None,   # acald_c  (draft 1626, 5242)
        '1064': None,   # g3p_c  (draft 1626, 5382)
        '1074': None,   # dgtp_c  (draft 282, 5736)
        '1076': None,   # dgsn_c  (draft 1834, 5382)
        '1100': None,   # dad_2_c  (draft 3090, 5286)
        '1114': None,   # damp_c  (draft 3783, 5286)
        '1116': None,   # dadp_c  (draft 4507, 5450)
        '1129': None,   # datp_c  (draft 5457, 5468)
        '1141': None,   # dgmp_c  (draft 1741, 5566)
        '1153': None,   # adp_c  (draft 4491, 5798)
        '1168': None,   # din_c  (draft 2712, 5682)
        '1195': None,   # dgdp_c  (draft 911, 5676)
        '1219': None,   # aicar_c  (draft 4037, 5938)
        '1220': None,   # fprica_c  (draft 2712, 5798)
        '1232': None,   # ap4a_c  (draft 4103, 5798)
        '1243': None,   # atp_c  (draft 5457, 5710)
        '1276': None,   # 25aics_c  (draft 5841, 5877)
        '1287': None,   # gdp_c  (draft 985, 6090)
        '1303': None,   # imp_c  (draft 3107, 6218)
        '1313': None,   # 5aizc_c  (draft 5536, 6157)
        '1330': None,   # gmp_c  (draft 1660, 6090)
        '1347': None,   # dcamp_c  (draft 3652, 5938)
        '1362': None,   # amp_c  (draft 4160, 6334)
        '1385': None,   # gtp_c  (draft 257, 6090)
        '1408': None,   # 5caiz_c  (draft 4956, 6297)
        '1417': None,   # 35cgmp_c  (draft 816, 6162)
        '1441': None,   # xmp_c  (draft 2309, 6242)
        '1500': None,   # air_c  (draft 5457, 6470)
        '1518': None,   # gsn_c  (draft 1251, 6840)
        '1558': None,   # xtsn_c  (draft 2132, 6662)
        '1584': None,   # ins_c  (draft 2810, 6590)
        '1622': None,   # adn_c  (draft 4134, 6730)
        '1662': None,   # fpram_c  (draft 6197, 6793)
        '1676': None,   # ade_c  (draft 4036, 7242)
        '1690': None,   # xan_c  (draft 2120, 7358)
        '1714': None,   # hxan_c  (draft 2810, 7358)
        '1747': None,   # gua_c  (draft 1462, 7212)
        '1860': None,   # fgam_c  (draft 6197, 7115)
        '1894': None,   # r5p_c  (draft 4500, 7242)
        '1907': None,   # gar_c  (draft 4937, 7113)
        '1922': None,   # adprib_c  (draft 4500, 6829)
        '1988': None,   # urate_c  (draft 2120, 7638)
        '2000': None,   # r1p_c  (draft 4356, 7522)
        '2011': None,   # pram_c  (draft 4937, 7451)
        '2028': None,   # prpp_c  (draft 6004, 7485)
        '3353': None,   # csn_c  (draft 4050, 7638)
        '3996': None,   # ura_c  (draft 4844, 7754)
        '3998': None,   # psd5p_c  (draft 4241, 7846)
        '4009': None,   # dad_5_c  (draft 3828, 7730)
        '4011': None,   # 5drib_c  (draft 3417, 7614)
        '3757': None,   # damp_p  (draft 3783, 5354)
        '3759': None,   # dad_2_p  (draft 3090, 5354)
    })
    L.reaction('1062', 'PM', side='R')   # PPM2: 2dr1p_c <=> 2dr5p_c
    L.reaction('1073', 'PM', side='R')   # DRPA: 2dr5p_c -> acald_c + g3p_c
    L.reaction('1087', 'PM', side='R')   # NTPTP1: dgtp_c + h2o_c -> dgsn_c + pppi_c
    L.reaction('1099', 'PM', side='R')   # PUNP4: dgsn_c + pi_c <=> 2dr1p_c + gua_c
    L.reaction('1112', 'PM', side='R')   # PUNP2: dad_2_c + pi_c <=> 2dr1p_c + ade_c
    L.reaction('1126', 'PM', side='R')   # DADK: atp_c + damp_c <=> adp_c + dadp_c
    L.reaction('1139', 'PM', side='R')   # NDPK8: atp_c + dadp_c <=> adp_c + datp_c
    L.reaction('1152', 'PM', side='R')   # NTPP1: dgtp_c + h2o_c -> dgmp_c + h_c + ppi_c
    L.reaction('1167', 'PM', side='R')   # RNDR1: adp_c + trdrd_c -> dadp_c + h2o_c + trdox_c
    L.reaction('1180', 'PM', side='R')   # PUNP6: din_c + pi_c <=> 2dr1p_c + hxan_c
    L.reaction('1192', 'PM', side='R')   # DADA: dad_2_c + h2o_c + h_c -> din_c + nh4_c
    L.reaction('1205', 'PM', side='R')   # DGK1: atp_c + dgmp_c <=> adp_c + dgdp_c
    L.reaction('1217', 'PM', side='R')   # NDPK5: atp_c + dgdp_c <=> adp_c + dgtp_c
    L.reaction('1231', 'PM', side='R')   # AICART: 10fthf_c + aicar_c <=> fprica_c + thf_c
    L.reaction('1242', 'PM', side='R')   # AP4AH: ap4a_c + h2o_c -> 2 adp_c + 2 h_c
    L.reaction('1255', 'PM', side='R')   # ATPM: atp_c + h2o_c -> adp_c + h_c + pi_c
    L.reaction('1265', 'PM', side='R')   # ATPS4rpp: adp_c + 4 h_p + pi_c <=> atp_c + h2o_c + 3 h_c
    L.reaction('1275', 'PM', side='R')   # NTP1: atp_c + h2o_c -> adp_c + h_c + pi_c
    L.reaction('1286', 'PM', side='R')   # ADSL2r: 25aics_c <=> aicar_c + fum_c
    L.reaction('1301', 'PM', side='R')   # RNDR2: gdp_c + trdrd_c -> dgdp_c + h2o_c + trdox_c
    L.reaction('1312', 'PM', side='R')   # IMPC: h2o_c + imp_c <=> fprica_c
    L.reaction('1329', 'PM', side='R')   # PRASCSi: 5aizc_c + asp__L_c + atp_c -> 25aics_c + adp_c + h_c + pi_c
    L.reaction('1344', 'PM', side='R')   # GMPR: gmp_c + 2 h_c + nadph_c -> imp_c + nadp_c + nh4_c
    L.reaction('1361', 'PM', side='R')   # ADSS: asp__L_c + gtp_c + imp_c -> dcamp_c + gdp_c + 2 h_c + pi_c
    L.reaction('1372', 'PM', side='R')   # ADSL1r: dcamp_c <=> amp_c + fum_c
    L.reaction('1382', 'PM', side='R')   # ADK1: amp_c + atp_c <=> 2 adp_c
    L.reaction('1395', 'PM', side='R')   # NDPK1: atp_c + gdp_c <=> adp_c + gtp_c
    L.reaction('1407', 'PM', side='R')   # GK1: atp_c + gmp_c <=> adp_c + gdp_c
    L.reaction('1416', 'PM', side='R')   # AIRC3: 5aizc_c <=> 5caiz_c
    L.reaction('1427', 'PM', side='R')   # GUACYC: gtp_c -> 35cgmp_c + ppi_c
    L.reaction('1437', 'PM', side='R')   # PDE4: 35cgmp_c + h2o_c -> gmp_c + h_c
    L.reaction('1458', 'PM', side='R')   # GMPS2: atp_c + gln__L_c + h2o_c + xmp_c -> amp_c + glu__L_c + gmp_c + 2 h_c + ppi_c
    L.reaction('1472', 'PM', side='R')   # IMPD: h2o_c + imp_c + nad_c -> h_c + nadh_c + xmp_c
    L.reaction('1516', 'PM', side='R')   # AIRC2: air_c + atp_c + hco3_c -> 5caiz_c + adp_c + h_c + pi_c
    L.reaction('1529', 'PM', side='R')   # NTD9: gmp_c + h2o_c -> gsn_c + pi_c
    L.reaction('1555', 'PM', side='R')   # GSNK: atp_c + gsn_c -> adp_c + gmp_c + h_c
    L.reaction('1568', 'PM', side='R')   # NTD10: h2o_c + xmp_c -> pi_c + xtsn_c
    L.reaction('1595', 'PM', side='R')   # INSK: atp_c + ins_c -> adp_c + h_c + imp_c
    L.reaction('1607', 'PM', side='R')   # NTD11: h2o_c + imp_c -> ins_c + pi_c
    L.reaction('1634', 'PM', side='R')   # ADNK1: adn_c + atp_c -> adp_c + amp_c + h_c
    L.reaction('1646', 'PM', side='R')   # NTD7: amp_c + h2o_c -> adn_c + pi_c
    L.reaction('1675', 'PM', side='R')   # PRAIS: atp_c + fpram_c -> adp_c + air_c + 2 h_c + pi_c
    L.reaction('1688', 'PM', side='R')   # ADPT: ade_c + prpp_c -> amp_c + ppi_c
    L.reaction('1701', 'PM', side='R')   # XPPT: prpp_c + xan_c -> ppi_c + xmp_c
    L.reaction('1713', 'PM', side='R')   # ADA: adn_c + h2o_c + h_c -> ins_c + nh4_c
    L.reaction('1724', 'PM', side='R')   # HXPRT: hxan_c + prpp_c -> imp_c + ppi_c
    L.reaction('1759', 'PM', side='R')   # GUAPRT: gua_c + prpp_c -> gmp_c + ppi_c
    L.reaction('1771', 'PM', side='R')   # PUNP1: adn_c + pi_c <=> ade_c + r1p_c
    L.reaction('1798', 'PM', side='R')   # ADNUC: adn_c + h2o_c -> ade_c + rib__D_c
    L.reaction('1810', 'PM', side='R')   # PUNP5: ins_c + pi_c <=> hxan_c + r1p_c
    L.reaction('1822', 'PM', side='R')   # INSH: h2o_c + ins_c -> hxan_c + rib__D_c
    L.reaction('1834', 'PM', side='R')   # PUNP3: gsn_c + pi_c <=> gua_c + r1p_c
    L.reaction('1846', 'PM', side='R')   # PUNP7: pi_c + xtsn_c <=> r1p_c + xan_c
    L.reaction('1858', 'PM', side='R')   # XTSNH: h2o_c + xtsn_c -> rib__D_c + xan_c
    L.reaction('1879', 'PM', side='R')   # PRFGS: atp_c + fgam_c + gln__L_c + h2o_c -> adp_c + fpram_c + glu__L_c + h_c + pi_c
    L.reaction('1904', 'PM', side='R')   # AMPN: amp_c + h2o_c -> ade_c + r5p_c
    L.reaction('1921', 'PM', side='R')   # GART: atp_c + for_c + gar_c -> adp_c + fgam_c + h_c + pi_c
    L.reaction('1934', 'PM', side='R')   # ADPRDP: adprib_c + h2o_c -> amp_c + 2 h_c + r5p_c
    L.reaction('1946', 'PM', side='R')   # ADD: ade_c + h2o_c + h_c -> hxan_c + nh4_c
    L.reaction('1958', 'PM', side='R')   # GUAD: gua_c + h2o_c + h_c -> nh4_c + xan_c
    L.reaction('1972', 'PM', side='R')   # HXAND: h2o_c + hxan_c + nad_c -> h_c + nadh_c + xan_c
    L.reaction('1984', 'PM', side='R')   # GARFT: 10fthf_c + gar_c <=> fgam_c + h_c + thf_c
    L.reaction('1999', 'PM', side='R')   # XAND: h2o_c + nad_c + xan_c -> h_c + nadh_c + urate_c
    L.reaction('2008', 'PM', side='R')   # PPM: r1p_c <=> r5p_c
    L.reaction('2025', 'PM', side='R')   # PRAGSr: atp_c + gly_c + pram_c <=> adp_c + gar_c + h_c + pi_c
    L.reaction('2042', 'PM', side='R')   # GLUPRT: gln__L_c + h2o_c + prpp_c -> glu__L_c + ppi_c + pram_c
    L.reaction('3363', 'PM', side='R')   # CMPN: cmp_c + h2o_c -> csn_c + r5p_c
    L.reaction('3756', 'PM', side='R')   # NTD6: damp_c + h2o_c -> dad_2_c + pi_c
    L.reaction('3782', 'PM', side='R')   # NTD8: dgmp_c + h2o_c -> dgsn_c + pi_c
    L.reaction('3808', 'PM', side='R')   # NTPP5: datp_c + h2o_c -> damp_c + h_c + ppi_c
    L.reaction('4008', 'PM', side='R')   # YUMPS: r5p_c + ura_c <=> h2o_c + psd5p_c
    L.reaction('4021', 'PM', side='R')   # 5DOAN: dad_5_c + h2o_c -> 5drib_c + ade_c
    L.reaction('3770', 'PM', side='R')   # NTD6pp: damp_p + h2o_p -> dad_2_p + pi_p
    # ---------------------------------------------------------------- PM2: pyrimidine metabolism
    L.setp('PM2', {
        '2044': None,   # cbp_c  (draft 759, 105)
        '2045': None,   # cbasp_c  (draft 759, 385)
        '2057': None,   # dhor__S_c  (draft 946, 525)
        '2069': None,   # orot_c  (draft 1138, 897)
        '2105': None,   # orot5p_c  (draft 1311, 1293)
        '2118': None,   # cytd_c  (draft 2713, 1525)
        '2120': None,   # cmp_c  (draft 2563, 1969)
        '2134': None,   # uri_c  (draft 2430, 1245)
        '2146': None,   # ump_c  (draft 1504, 1573)
        '2157': None,   # csn_c  (draft 2981, 1269)
        '2171': None,   # ura_c  (draft 2073, 805)
        '2182': None,   # csn_c  (draft 2215, 245)
        '2233': None,   # cdp_c  (draft 2563, 3003)
        '2246': None,   # udp_c  (draft 1311, 1953)
        '2271': None,   # utp_c  (draft 1928, 2439)
        '2283': None,   # dudp_c  (draft 1446, 2286)
        '2300': None,   # ctp_c  (draft 2647, 2650)
        '2318': None,   # dcyt_c  (draft 2183, 4399)
        '2320': None,   # duri_c  (draft 1985, 3351)
        '2334': None,   # dutp_c  (draft 1777, 2809)
        '2346': None,   # 2dr1p_c  (draft 917, 2801)
        '2359': None,   # cmp_c  (draft 2894, 2650)
        '2373': None,   # dump_c  (draft 1530, 3089)
        '2396': None,   # dctp_c  (draft 2280, 3235)
        '2410': None,   # thymd_c  (draft 790, 3025)
        '2411': None,   # thym_c  (draft 313, 2872)
        '2436': None,   # dtmp_c  (draft 1075, 3795)
        '2448': None,   # cytd_c  (draft 3237, 2650)
        '2475': None,   # dcmp_c  (draft 2330, 4051)
        '2488': None,   # dcdp_c  (draft 2696, 3747)
        '2516': None,   # dtdp_c  (draft 1075, 4331)
        '2541': None,   # dttp_c  (draft 1348, 4172)
        '3904': None,   # pep_c  (draft 1615, 4805)
        '3905': None,   # pyr_c  (draft 1615, 4689)
        '3916': None,   # dtdprmn_c  (draft 418, 4075)
        '3917': None,   # kphphhlipa_c  (draft 456, 4293)
        '3918': None,   # icolipa_c  (draft 456, 4409)
        '4022': None,   # dtdp4aaddg_c  (draft 536, 4525)
        '4023': None,   # unagamu_c  (draft 613, 4641)
        '4024': None,   # unagamuf_c  (draft 754, 4805)
    })
    L.reaction('2056', 'PM2', side='R')   # ASPCT: asp__L_c + cbp_c -> cbasp_c + h_c + pi_c
    L.reaction('2067', 'PM2', side='R')   # DHORTS: dhor__S_c + h2o_c <=> cbasp_c + h_c
    L.reaction('2080', 'PM2', side='R')   # DHORD2: dhor__S_c + q8_c -> orot_c + q8h2_c
    L.reaction('2092', 'PM2', side='R')   # DHORD5: dhor__S_c + mqn8_c -> mql8_c + orot_c
    L.reaction('2104', 'PM2', side='R')   # DHORDfum: dhor__S_c + fum_c -> orot_c + succ_c
    L.reaction('2117', 'PM2', side='R')   # ORPT: orot5p_c + ppi_c <=> orot_c + prpp_c
    L.reaction('2131', 'PM2', side='R')   # CYTDK2: cytd_c + gtp_c -> cmp_c + gdp_c + h_c
    L.reaction('2144', 'PM2', side='R')   # CYTD: cytd_c + h2o_c + h_c -> nh4_c + uri_c
    L.reaction('2155', 'PM2', side='R')   # OMPDC: h_c + orot5p_c -> co2_c + ump_c
    L.reaction('2168', 'PM2', side='R')   # CYTDH: cytd_c + h2o_c -> csn_c + rib__D_c
    L.reaction('2181', 'PM2', side='R')   # PYNP2r: pi_c + uri_c <=> r1p_c + ura_c
    L.reaction('2194', 'PM2', side='R')   # CSND: csn_c + h2o_c + h_c -> nh4_c + ura_c
    L.reaction('2206', 'PM2', side='R')   # URIH: h2o_c + uri_c -> rib__D_c + ura_c
    L.reaction('2218', 'PM2', side='R')   # URIK2: gtp_c + uri_c -> gdp_c + h_c + ump_c
    L.reaction('2230', 'PM2', side='R')   # UPPRT: prpp_c + ura_c -> ppi_c + ump_c
    L.reaction('2243', 'PM2', side='R')   # CYTK1: atp_c + cmp_c <=> adp_c + cdp_c
    L.reaction('2256', 'PM2', side='R')   # UMPK: atp_c + ump_c <=> adp_c + udp_c
    L.reaction('2268', 'PM2', side='R')   # UMPK: atp_c + ump_c <=> adp_c + udp_c
    L.reaction('2281', 'PM2', side='R')   # NDPK2: atp_c + udp_c <=> adp_c + utp_c
    L.reaction('2296', 'PM2', side='R')   # RNDR4: trdrd_c + udp_c -> dudp_c + h2o_c + trdox_c
    L.reaction('2317', 'PM2', side='R')   # CTPS2: atp_c + gln__L_c + h2o_c + utp_c -> adp_c + ctp_c + glu__L_c + 2 h_c + pi_c
    L.reaction('2331', 'PM2', side='R')   # DCYTD: dcyt_c + h2o_c + h_c -> duri_c + nh4_c
    L.reaction('2344', 'PM2', side='R')   # NDPK6: atp_c + dudp_c <=> adp_c + dutp_c
    L.reaction('2357', 'PM2', side='R')   # DURIPP: duri_c + pi_c <=> 2dr1p_c + ura_c
    L.reaction('2370', 'PM2', side='R')   # NTPP4: ctp_c + h2o_c -> cmp_c + h_c + ppi_c
    L.reaction('2383', 'PM2', side='R')   # DURIK1: atp_c + duri_c -> adp_c + dump_c + h_c
    L.reaction('2395', 'PM2', side='R')   # NDPK3: atp_c + cdp_c <=> adp_c + ctp_c
    L.reaction('2408', 'PM2', side='R')   # DCTPD: dctp_c + h2o_c + h_c -> dutp_c + nh4_c
    L.reaction('2421', 'PM2', side='R')   # TMDPP: pi_c + thymd_c <=> 2dr1p_c + thym_c
    L.reaction('2433', 'PM2', side='R')   # DUTPDP: dutp_c + h2o_c -> dump_c + h_c + ppi_c
    L.reaction('2446', 'PM2', side='R')   # TMDS: dump_c + mlthf_c -> dhf_c + dtmp_c
    L.reaction('2459', 'PM2', side='R')   # NTD4: cmp_c + h2o_c -> cytd_c + pi_c
    L.reaction('2486', 'PM2', side='R')   # NTPP3: dctp_c + h2o_c -> dcmp_c + h_c + ppi_c
    L.reaction('2501', 'PM2', side='R')   # RNDR3: cdp_c + trdrd_c -> dcdp_c + h2o_c + trdox_c
    L.reaction('2513', 'PM2', side='R')   # TMDK1: atp_c + thymd_c -> adp_c + dtmp_c + h_c
    L.reaction('2526', 'PM2', side='R')   # DTMPK: atp_c + dtmp_c <=> adp_c + dtdp_c
    L.reaction('2538', 'PM2', side='R')   # NDPK7: atp_c + dcdp_c <=> adp_c + dctp_c
    L.reaction('2551', 'PM2', side='R')   # NDPK4: atp_c + dtdp_c <=> adp_c + dttp_c
    L.reaction('3399', 'PM2', side='R')   # CYTK2: atp_c + dcmp_c <=> adp_c + dcdp_c
    L.reaction('3678', 'PM2', side='R')   # NTD1: dump_c + h2o_c -> duri_c + pi_c
    L.reaction('3704', 'PM2', side='R')   # NTD3: dcmp_c + h2o_c -> dcyt_c + pi_c
    L.reaction('3730', 'PM2', side='R')   # NTD5: dtmp_c + h2o_c -> pi_c + thymd_c
    L.reaction('3820', 'PM2', side='R')   # NTPP7: dttp_c + h2o_c -> dtmp_c + h_c + ppi_c
    L.reaction('3915', 'PM2', side='R')   # PYK6: dtdp_c + h_c + pep_c -> dttp_c + pyr_c
    L.reaction('3928', 'PM2', side='R')   # RHAT1: dtdprmn_c + kphphhlipa_c -> dtdp_c + h_c + icolipa_c
    L.reaction('4034', 'PM2', side='R')   # AADDGT: dtdp4aaddg_c + unagamu_c -> dtdp_c + h_c + unagamuf_c
    # ---------------------------------------------------------------- NM: NAD metabolism
    L.setp('NM', {
        '80': None,   # nad_c  (draft 7215, 366)
        '83': None,   # nadp_c  (draft 6660, 195)
        '95': None,   # nadph_c  (draft 6136, 212)
        '153': None,   # dnad_c  (draft 7945, 207)
        '169': None,   # asp__L_c  (draft 8561, 131)
        '172': None,   # iasp_c  (draft 8316, 595)
        '185': None,   # ncam_c  (draft 6232, 501)
        '198': None,   # nmn_c  (draft 7147, 851)
        '233': None,   # nmn_p  (draft 7147, 919)
        '247': None,   # quln_c  (draft 8316, 922)
        '260': None,   # nac_c  (draft 6231, 1076)
        '273': None,   # nicrnt_c  (draft 7996, 1115)
    })
    L.reaction('93', 'NM', side='R')   # NADTRHD: nad_c + nadph_c -> nadh_c + nadp_c
    L.reaction('106', 'NM', side='R')   # NADTRHD: nad_c + nadph_c -> nadh_c + nadp_c
    L.reaction('118', 'NM', side='R')   # NADPPPS: h2o_c + nadp_c -> nad_c + pi_c
    L.reaction('130', 'NM', side='R')   # NADK: atp_c + nad_c -> adp_c + h_c + nadp_c
    L.reaction('141', 'NM', side='R')   # NADTRHD: nad_c + nadph_c -> nadh_c + nadp_c
    L.reaction('151', 'NM', side='R')   # THD2pp: 2 h_p + nadh_c + nadp_c -> 2 h_c + nad_c + nadph_c
    L.reaction('168', 'NM', side='R')   # NADS1: atp_c + dnad_c + nh4_c -> amp_c + h_c + nad_c + ppi_c
    L.reaction('182', 'NM', side='R')   # ASPO6: asp__L_c + o2_c -> h2o2_c + h_c + iasp_c
    L.reaction('195', 'NM', side='R')   # NADN: h2o_c + nad_c -> adprib_c + h_c + ncam_c
    L.reaction('208', 'NM', side='R')   # NADDP: h2o_c + nad_c -> amp_c + 2 h_c + nmn_c
    L.reaction('220', 'NM', side='R')   # NMNAT: atp_c + h_c + nmn_c -> nad_c + ppi_c
    L.reaction('232', 'NM', side='R')   # NMNN: h2o_c + nmn_c -> h_c + ncam_c + r5p_c
    L.reaction('243', 'NM', side='R')   # NMNt7pp: h2o_c + nmn_p -> h_c + ncam_c + r5p_c
    L.reaction('258', 'NM', side='R')   # QULNS: dhap_c + iasp_c -> 2 h2o_c + pi_c + quln_c
    L.reaction('271', 'NM', side='R')   # NNAM: h2o_c + ncam_c -> nac_c + nh4_c
    L.reaction('284', 'NM', side='R')   # NNATr: atp_c + h_c + nicrnt_c <=> dnad_c + ppi_c
    L.reaction('296', 'NM', side='R')   # NMNDA: h2o_c + nmn_c -> nh4_c + nicrnt_c
    L.reaction('316', 'NM', side='R')   # NAMNPP: atp_c + h2o_c + nac_c + prpp_c -> adp_c + nicrnt_c + pi_c + ppi_c
    L.reaction('330', 'NM', side='R')   # NNDPR: 2 h_c + prpp_c + quln_c -> co2_c + nicrnt_c + ppi_c
    L.reaction('3310', 'NM', side='R')   # ASPO3: asp__L_c + q8_c -> h_c + iasp_c + q8h2_c
    L.reaction('3322', 'NM', side='R')   # ASPO4: asp__L_c + mqn8_c -> h_c + iasp_c + mql8_c
    # ---------------------------------------------------------------- VBM: vitamin B12 metabolism
    L.setp('VBM', {
        '2870': None,   # fmn_c  (draft 9836, 685)
        '2872': None,   # fmnh2_c  (draft 10510, 601)
        '3272': None,   # arbtn_fe3_c  (draft 9321, 453)
        '3273': None,   # arbtn_c  (draft 9075, 569)
        '3274': None,   # fe2_c  (draft 9458, 917)
        '3286': None,   # rbflvrd_c  (draft 9237, 1033)
        '3287': None,   # ribflv_c  (draft 9010, 1149)
        '3364': None,   # cpgn_c  (draft 9332, 569)
        '3365': None,   # cpgn_un_c  (draft 9391, 685)
        '3409': None,   # isetac_c  (draft 10115, 105)
        '3411': None,   # gcald_c  (draft 10367, 105)
        '3413': None,   # so3_c  (draft 10563, 484)
        '3427': None,   # mso3_c  (draft 10901, 573)
        '3429': None,   # fald_c  (draft 11050, 1191)
        '3444': None,   # ethso3_c  (draft 10901, 641)
        '3445': None,   # acald_c  (draft 11050, 1259)
        '3459': None,   # butso3_c  (draft 10901, 709)
        '3460': None,   # btal_c  (draft 11050, 1327)
        '3474': None,   # sulfac_c  (draft 10901, 777)
        '3475': None,   # glx_c  (draft 11050, 1395)
        '3489': None,   # fe3hox_c  (draft 9723, 1265)
        '3490': None,   # fe3hox_un_c  (draft 9723, 1497)
        '3513': None,   # fecrm_c  (draft 9335, 1497)
        '3514': None,   # fecrm_un_c  (draft 9335, 1613)
        '3537': None,   # feenter_c  (draft 9515, 1033)
        '3538': None,   # enter_c  (draft 9550, 1149)
        '3561': None,   # feoxam_c  (draft 9696, 337)
        '3562': None,   # feoxam_un_c  (draft 9321, 221)
        '3650': None,   # Nmtrp_c  (draft 11570, 1613)
        '3654': None,   # trp__L_c  (draft 11740, 1497)
    })
    L.reaction('2883', 'VBM', side='R')   # FMNRx: fmn_c + h_c + nadh_c -> fmnh2_c + nad_c
    L.reaction('2895', 'VBM', side='R')   # FMNRx2: fmn_c + h_c + nadph_c -> fmnh2_c + nadp_c
    L.reaction('3285', 'VBM', side='R')   # ARBTNR2: 2 arbtn_fe3_c + fmnh2_c -> 2 arbtn_c + 2 fe2_c + fmn_c + 2 h_c
    L.reaction('3298', 'VBM', side='R')   # ARBTNR3: 2 arbtn_fe3_c + rbflvrd_c -> 2 arbtn_c + 2 fe2_c + 2 h_c + ribflv_c
    L.reaction('3376', 'VBM', side='R')   # CPGNR2: 2 cpgn_c + fmnh2_c -> 2 cpgn_un_c + 2 fe2_c + fmn_c + 2 h_c
    L.reaction('3387', 'VBM', side='R')   # CPGNR3: 2 cpgn_c + rbflvrd_c -> 2 cpgn_un_c + 2 fe2_c + 2 h_c + ribflv_c
    L.reaction('3426', 'VBM', side='R')   # FDMO: fmnh2_c + isetac_c + o2_c -> fmn_c + gcald_c + h2o_c + h_c + so3_c
    L.reaction('3443', 'VBM', side='R')   # FDMO2: fmnh2_c + mso3_c + o2_c -> fald_c + fmn_c + h2o_c + h_c + so3_c
    L.reaction('3458', 'VBM', side='R')   # FDMO3: ethso3_c + fmnh2_c + o2_c -> acald_c + fmn_c + h2o_c + h_c + so3_c
    L.reaction('3473', 'VBM', side='R')   # FDMO4: butso3_c + fmnh2_c + o2_c -> btal_c + fmn_c + h2o_c + h_c + so3_c
    L.reaction('3488', 'VBM', side='R')   # FDMO6: fmnh2_c + o2_c + sulfac_c -> fmn_c + glx_c + h2o_c + h_c + so3_c
    L.reaction('3501', 'VBM', side='R')   # FE3HOXR2: 2 fe3hox_c + fmnh2_c -> 2 fe2_c + 2 fe3hox_un_c + fmn_c + 2 h_c
    L.reaction('3512', 'VBM', side='R')   # FE3HOXR3: 2 fe3hox_c + rbflvrd_c -> 2 fe2_c + 2 fe3hox_un_c + 2 h_c + ribflv_c
    L.reaction('3525', 'VBM', side='R')   # FECRMR2: 2 fecrm_c + fmnh2_c -> 2 fe2_c + 2 fecrm_un_c + fmn_c + 2 h_c
    L.reaction('3536', 'VBM', side='R')   # FECRMR3: 2 fecrm_c + rbflvrd_c -> 2 fe2_c + 2 fecrm_un_c + 2 h_c + ribflv_c
    L.reaction('3549', 'VBM', side='R')   # FEENTERR2: 2 feenter_c + fmnh2_c -> 2 enter_c + 2 fe2_c + fmn_c + 2 h_c
    L.reaction('3560', 'VBM', side='R')   # FEENTERR3: 2 feenter_c + rbflvrd_c -> 2 enter_c + 2 fe2_c + 2 h_c + ribflv_c
    L.reaction('3573', 'VBM', side='R')   # FEOXAMR2: 2 feoxam_c + fmnh2_c -> 2 fe2_c + 2 feoxam_un_c + fmn_c + 2 h_c
    L.reaction('3584', 'VBM', side='R')   # FEOXAMR3: 2 feoxam_c + rbflvrd_c -> 2 fe2_c + 2 feoxam_un_c + 2 h_c + ribflv_c
    L.reaction('3666', 'VBM', side='R')   # MTRPOX: Nmtrp_c + h2o_c + o2_c -> fald_c + h2o2_c + trp__L_c
    # ---------------------------------------------------------------- VBM2: vitamin B6 metabolism
    L.setp('VBM2', {
        '2963': None,   # e4p_c  (draft 7672, 5697)
        '2966': None,   # 4per_c  (draft 7869, 5982)
        '2981': None,   # ohpb_c  (draft 7869, 6262)
        '2994': None,   # phthr_c  (draft 7755, 6542)
        '3006': None,   # pydam_c  (draft 9285, 7192)
        '3008': None,   # pyam5p_c  (draft 8759, 7286)
        '3020': None,   # pydxn_c  (draft 7122, 6950)
        '3022': None,   # pdx5p_c  (draft 7832, 6928)
        '3034': None,   # pydx_c  (draft 8774, 6494)
        '3036': None,   # pydx5p_c  (draft 8343, 7170)
        '3210': None,   # 4hthr_c  (draft 7190, 6368)
        '3249': None,   # ala__D_c  (draft 8214, 7402)
        '3250': None,   # pyr_c  (draft 8659, 7402)
        '3261': None,   # ala__L_c  (draft 8344, 7402)
        '3624': None,   # cys__L_c  (draft 8232, 7518)
        '3625': None,   # iscs_c  (draft 8497, 7634)
        '3626': None,   # iscssh_c  (draft 8885, 7521)
        '3848': None,   # dxyl5p_c  (draft 7476, 6650)
        '3929': None,   # sufse_c  (draft 7599, 7411)
        '3930': None,   # sufsesh_c  (draft 8000, 7518)
    })
    L.reaction('2978', 'VBM2', side='R')   # E4PD: e4p_c + h2o_c + nad_c <=> 4per_c + 2 h_c + nadh_c
    L.reaction('2991', 'VBM2', side='R')   # PERD: 4per_c + nad_c <=> h_c + nadh_c + ohpb_c
    L.reaction('3004', 'VBM2', side='R')   # OHPBAT: glu__L_c + ohpb_c <=> akg_c + phthr_c
    L.reaction('3018', 'VBM2', side='R')   # PYDAMK: atp_c + pydam_c -> adp_c + h_c + pyam5p_c
    L.reaction('3032', 'VBM2', side='R')   # PYDXNK: atp_c + pydxn_c -> adp_c + h_c + pdx5p_c
    L.reaction('3046', 'VBM2', side='R')   # PYDXK: atp_c + pydx_c -> adp_c + h_c + pydx5p_c
    L.reaction('3062', 'VBM2', side='R')   # PYAM5PO: h2o_c + o2_c + pyam5p_c -> h2o2_c + nh4_c + pydx5p_c
    L.reaction('3074', 'VBM2', side='R')   # PDX5POi: o2_c + pdx5p_c -> h2o2_c + pydx5p_c
    L.reaction('3222', 'VBM2', side='R')   # 4HTHRK: 4hthr_c + atp_c -> adp_c + h_c + phthr_c
    L.reaction('3260', 'VBM2', side='R')   # ALATA_D2: ala__D_c + pydx5p_c -> pyam5p_c + pyr_c
    L.reaction('3271', 'VBM2', side='R')   # ALATA_L2: ala__L_c + pydx5p_c -> pyam5p_c + pyr_c
    L.reaction('3623', 'VBM2', side='R')   # HYPOE: h2o_c + pyam5p_c -> pi_c + pydam_c
    L.reaction('3636', 'VBM2', side='R')   # ICYSDS: cys__L_c + iscs_c -> ala__L_c + iscssh_c
    L.reaction('3847', 'VBM2', side='R')   # PDX5PO2: nad_c + pdx5p_c -> h_c + nadh_c + pydx5p_c
    L.reaction('3867', 'VBM2', side='R')   # PDX5PS: dxyl5p_c + nad_c + phthr_c -> co2_c + 2 h2o_c + h_c + nadh_c + pdx5p_c + pi_c
    L.reaction('3879', 'VBM2', side='R')   # PDXPP: h2o_c + pdx5p_c -> pi_c + pydxn_c
    L.reaction('3903', 'VBM2', side='R')   # PYDXPP: h2o_c + pydx5p_c -> pi_c + pydx_c
    L.reaction('3940', 'VBM2', side='R')   # SCYSDS: cys__L_c + sufse_c -> ala__L_c + sufsesh_c
    # ---------------------------------------------------------------- HM: heme metabolism
    L.setp('HM', {
        '672': None,   # glu__L_c  (draft 10019, 4499)
        '675': None,   # glutrna_c  (draft 10019, 4895)
        '690': None,   # glu1sa_c  (draft 10012, 5291)
        '704': None,   # 5aop_c  (draft 10304, 5259)
        '713': None,   # hmbil_c  (draft 10442, 5883)
        '715': None,   # uppg3_c  (draft 10329, 6175)
        '726': None,   # ppbng_c  (draft 10304, 5639)
        '789': None,   # cpppg3_c  (draft 10847, 6415)
        '814': None,   # pppg9_c  (draft 11146, 6895)
        '841': None,   # ppp9_c  (draft 11023, 7243)
        '854': None,   # pheme_c  (draft 10824, 7518)
        '3197': None,   # 4abzglu_c  (draft 10115, 4151)
        '3199': None,   # 4abz_c  (draft 10231, 4499)
        '3597': None,   # frdp_c  (draft 11330, 7523)
        '3599': None,   # hemeO_c  (draft 10948, 7755)
        '3880': None,   # fum_c  (draft 11437, 6900)
        '3881': None,   # succ_c  (draft 11333, 7243)
    })
    L.reaction('688', 'HM', side='R')   # GLUTRS: atp_c + glu__L_c + trnaglu_c -> amp_c + glutrna_c + ppi_c
    L.reaction('703', 'HM', side='R')   # GLUTRR: glutrna_c + h_c + nadph_c -> glu1sa_c + nadp_c + trnaglu_c
    L.reaction('712', 'HM', side='R')   # G1SAT: glu1sa_c <=> 5aop_c
    L.reaction('724', 'HM', side='R')   # UPP3S: hmbil_c -> h2o_c + uppg3_c
    L.reaction('737', 'HM', side='R')   # HMBS: h2o_c + 4 ppbng_c -> hmbil_c + 4 nh4_c
    L.reaction('747', 'HM', side='R')   # PPBNGS: 2 5aop_c -> 2 h2o_c + h_c + ppbng_c
    L.reaction('760', 'HM', side='R')   # UPP3MT: 2 amet_c + uppg3_c -> 2 ahcys_c + dscl_c + h_c
    L.reaction('798', 'HM', side='R')   # UPPDC1: 4 h_c + uppg3_c -> 4 co2_c + cpppg3_c
    L.reaction('826', 'HM', side='R')   # CPPPGO2: 2 amet_c + cpppg3_c -> 2 co2_c + 2 dad_5_c + 2 met__L_c + pppg9_c
    L.reaction('839', 'HM', side='R')   # CPPPGO: cpppg3_c + 2 h_c + o2_c -> 2 co2_c + 2 h2o_c + pppg9_c
    L.reaction('852', 'HM', side='R')   # PPPGO: 1.5 o2_c + pppg9_c -> 3 h2o_c + ppp9_c
    L.reaction('863', 'HM', side='R')   # FCLT: fe2_c + ppp9_c -> 2 h_c + pheme_c
    L.reaction('3209', 'HM', side='R')   # 4ABZGLUH: 4abzglu_c + h2o_c <=> 4abz_c + glu__L_c
    L.reaction('3611', 'HM', side='R')   # HEMEOS: frdp_c + h2o_c + pheme_c -> hemeO_c + ppi_c
    L.reaction('3891', 'HM', side='R')   # PPPGO3: 3 fum_c + pppg9_c -> ppp9_c + 3 succ_c
    # ---------------------------------------------------------------- FB: flavin biosynthesis
    L.setp('FB', {
        '527': None,   # gtp_c  (draft 4377, 2450)
        '529': None,   # 25drapp_c  (draft 4377, 2730)
        '544': None,   # 5apru_c  (draft 4754, 2870)
        '557': None,   # 5aprbu_c  (draft 5115, 3150)
        '570': None,   # 4r5au_c  (draft 5677, 3172)
        '582': None,   # dmlz_c  (draft 5074, 3420)
        '583': None,   # ribflv_c  (draft 4719, 3454)
        '594': None,   # db4p_c  (draft 5677, 3498)
        '610': None,   # fmn_c  (draft 4371, 3454)
        '623': None,   # rbflvrd_c  (draft 4719, 3855)
        '634': None,   # ru5p__D_c  (draft 6026, 3668)
        '646': None,   # fad_c  (draft 3987, 3512)
        '659': None,   # fadh2_c  (draft 3813, 3729)
    })
    L.reaction('542', 'FB', side='R')   # GTPCII2: gtp_c + 3 h2o_c -> 25drapp_c + for_c + 2 h_c + ppi_c
    L.reaction('555', 'FB', side='R')   # DHPPDA2: 25drapp_c + h2o_c + h_c -> 5apru_c + nh4_c
    L.reaction('568', 'FB', side='R')   # APRAUR: 5apru_c + h_c + nadph_c -> 5aprbu_c + nadp_c
    L.reaction('581', 'FB', side='R')   # PMDPHT: 5aprbu_c + h2o_c -> 4r5au_c + pi_c
    L.reaction('593', 'FB', side='R')   # RBFSb: 2 dmlz_c -> 4r5au_c + ribflv_c
    L.reaction('607', 'FB', side='R')   # RBFSa: 4r5au_c + db4p_c -> dmlz_c + 2 h2o_c + pi_c
    L.reaction('620', 'FB', side='R')   # RBFK: atp_c + ribflv_c -> adp_c + fmn_c + h_c
    L.reaction('633', 'FB', side='R')   # FLVR: h_c + nadph_c + ribflv_c -> nadp_c + rbflvrd_c
    L.reaction('644', 'FB', side='R')   # DB4PS: ru5p__D_c -> db4p_c + for_c + h_c
    L.reaction('657', 'FB', side='R')   # FMNAT: atp_c + fmn_c + h_c -> fad_c + ppi_c
    L.reaction('670', 'FB', side='R')   # FADRx2: fad_c + h_c + nadph_c -> fadh2_c + nadp_c
    L.reaction('3234', 'FB', side='R')   # ACP1_FMN: fmn_c + h2o_c -> pi_c + ribflv_c
    L.reaction('3596', 'FB', side='R')   # FLVRx: h_c + nadh_c + ribflv_c -> nad_c + rbflvrd_c
    # ---------------------------------------------------------------- VKM: vitamin K metabolism
    L.setp('VKM', {
        '3075': None,   # chor_c  (draft 7142, 4526)
        '3076': None,   # ichor_c  (draft 7893, 4526)
        '3094': None,   # 2sephchc_c  (draft 8680, 4526)
        '3106': None,   # 2shchc_c  (draft 8867, 4791)
        '3119': None,   # sucbz_c  (draft 8045, 4769)
        '3122': None,   # sbzcoa_c  (draft 7185, 4769)
        '3145': None,   # 14dhncoa_c  (draft 7618, 5027)
        '3158': None,   # dhna_c  (draft 8545, 5071)
        '3169': None,   # 2dmmql8_c  (draft 9099, 5327)
        '3186': None,   # mql8_c  (draft 7894, 5273)
    })
    L.reaction('3084', 'VKM', side='R')   # ICHORS_copy1: chor_c -> ichor_c
    L.reaction('3092', 'VKM', side='R')   # ICHORS_copy2: chor_c <=> ichor_c
    L.reaction('3105', 'VKM', side='R')   # SEPHCHCS: akg_c + h_c + ichor_c -> 2sephchc_c + co2_c
    L.reaction('3116', 'VKM', side='R')   # SHCHCS3: 2sephchc_c -> 2shchc_c + pyr_c
    L.reaction('3134', 'VKM', side='R')   # SUCBZL: atp_c + coa_c + sucbz_c -> amp_c + ppi_c + sbzcoa_c
    L.reaction('3144', 'VKM', side='R')   # SUCBZS: 2shchc_c -> h2o_c + sucbz_c
    L.reaction('3155', 'VKM', side='R')   # DHNCOAS: h_c + sbzcoa_c -> 14dhncoa_c + h2o_c
    L.reaction('3168', 'VKM', side='R')   # DHNCOAT: 14dhncoa_c + h2o_c -> coa_c + dhna_c + h_c
    L.reaction('3183', 'VKM', side='R')   # DHNAOT4: dhna_c + h_c + octdp_c -> 2dmmql8_c + co2_c + ppi_c
    L.reaction('3196', 'VKM', side='R')   # AMMQLT8: 2dmmql8_c + amet_c -> ahcys_c + h_c + mql8_c
    # ---------------------------------------------------------------- MCB: molybdenum cofactor biosynthesis
    L.setp('MCB', {
        '896': None,   # gtp_c  (draft 8667, 2295)
        '897': None,   # cpmp_c  (draft 9028, 2326)
        '911': None,   # mpt_c  (draft 9028, 2442)
        '912': None,   # mptamp_c  (draft 8949, 2838)
        '926': None,   # moco_c  (draft 9028, 3118)
        '3324': None,   # bmoco_c  (draft 9242, 3398)
        '3325': None,   # cu2_c  (draft 9248, 2606)
        '3337': None,   # wco_c  (draft 9653, 2838)
        '3339': None,   # bwco_c  (draft 9525, 3210)
        '3637': None,   # moadcosh_c  (draft 8849, 1978)
        '3638': None,   # moadcoo_c  (draft 9233, 2094)
        '3981': None,   # tungs_c  (draft 9572, 2600)
    })
    L.reaction('909', 'MCB', side='R')   # CPMPS: gtp_c + h2o_c -> cpmp_c + ppi_c
    L.reaction('923', 'MCB', side='R')   # MPTAT: atp_c + h_c + mpt_c -> mptamp_c + ppi_c
    L.reaction('940', 'MCB', side='R')   # MOCOS: 2 h_c + mobd_c + mptamp_c -> amp_c + cu2_c + h2o_c + moco_c
    L.reaction('953', 'MCB', side='R')   # MOCDS: ctp_c + h_c + moco_c -> mococdp_c + ppi_c
    L.reaction('3336', 'MCB', side='R')   # BMOCOS: moco_c + mptamp_c -> amp_c + bmoco_c + cu2_c
    L.reaction('3350', 'MCB', side='R')   # BWCOS: mptamp_c + wco_c -> amp_c + bwco_c + cu2_c
    L.reaction('3649', 'MCB', side='R')   # MPTS: cpmp_c + cu2_c + 2 moadcosh_c -> 5 h_c + 2 moadcoo_c + mpt_c
    L.reaction('3995', 'MCB', side='R')   # WCOS: 2 h_c + mptamp_c + tungs_c -> amp_c + cu2_c + h2o_c + wco_c
    # ---------------------------------------------------------------- VBM3: vitamin B1 metabolism
    L.setp('VBM3', {
        '2742': None,   # tyr__L_c  (draft 10854, 1983)
        '2745': None,   # dhgly_c  (draft 10793, 2616)
        '2763': None,   # 4mhetz_c  (draft 10421, 2616)
        '2765': None,   # 4mpetz_c  (draft 10657, 3012)
        '2778': None,   # pyr_c  (draft 11463, 1978)
        '2780': None,   # dxyl5p_c  (draft 11428, 2569)
        '2792': None,   # thmmp_c  (draft 10313, 3428)
        '2806': None,   # thmpp_c  (draft 10337, 3781)
        '2818': None,   # thm_c  (draft 10313, 3012)
        '3942': None,   # iscssh_c  (draft 11121, 2616)
        '3947': None,   # iscs_c  (draft 11121, 3080)
    })
    L.reaction('2762', 'VBM3', side='R')   # TYRL: amet_c + nadph_c + tyr__L_c -> 4crsol_c + dad_5_c + dhgly_c + h_c + met__L_c + nadp_c
    L.reaction('2776', 'VBM3', side='R')   # HETZK: 4mhetz_c + atp_c -> 4mpetz_c + adp_c + h_c
    L.reaction('2790', 'VBM3', side='R')   # DXPS: g3p_c + h_c + pyr_c -> co2_c + dxyl5p_c
    L.reaction('2803', 'VBM3', side='R')   # TMPPP: 2mahmp_c + 4mpetz_c + h_c -> ppi_c + thmmp_c
    L.reaction('2816', 'VBM3', side='R')   # TMPK: atp_c + thmmp_c -> adp_c + thmpp_c
    L.reaction('2829', 'VBM3', side='R')   # TMK: atp_c + thm_c -> adp_c + h_c + thmmp_c
    L.reaction('2842', 'VBM3', side='R')   # TDP: h2o_c + thmpp_c -> h_c + pi_c + thmmp_c
    L.reaction('3967', 'VBM3', side='R')   # THZPSN3: atp_c + dhgly_c + dxyl5p_c + h_c + iscssh_c + nadph_c -> 4mpetz_c + amp_c + co2_c + 2 h2o_c + iscs_c + nadp_c + ppi_c
    # ---------------------------------------------------------------- TM: tetrahydrofolate metabolism
    L.setp('TM', {
        '2552': None,   # 10fthf_c  (draft 4939, 1805)
        '2555': None,   # thf_c  (draft 4818, 1383)
        '2580': None,   # dhf_c  (draft 4096, 1869)
        '2615': None,   # dhpt_c  (draft 4090, 1383)
        '2671': None,   # 4abz_c  (draft 4231, 1015)
        '2672': None,   # 6hmhptpp_c  (draft 3842, 1015)
        '2684': None,   # 4adcho_c  (draft 4098, 559)
        '2709': None,   # chor_c  (draft 4098, 105)
        '2566': None,   # gtp_c  (draft 5306, 105)
        '2568': None,   # ahdt_c  (draft 5232, 733)
        '2605': None,   # dhmptp_c  (draft 5081, 1129)
        '2633': None,   # cph4_c  (draft 5549, 1129)
        '2660': None,   # cdg_c  (draft 5549, 1409)
        '2727': None,   # preq0_c  (draft 5549, 1921)
        '2695': None,   # dhmpt_c  (draft 4697, 105)
        '2698': None,   # thmnp_c  (draft 4656, 559)
        '3400': None,   # dhnpt_c  (draft 4940, 105)
    })
    L.reaction('2565', 'TM', side='R')   # GARFT: 10fthf_c + gar_c <=> fgam_c + h_c + thf_c
    L.reaction('2592', 'TM', side='R')   # DHFR: dhf_c + h_c + nadph_c <=> nadp_c + thf_c
    L.reaction('2604', 'TM', side='R')   # FTHFD: 10fthf_c + h2o_c -> for_c + h_c + thf_c
    L.reaction('2630', 'TM', side='R')   # DHFS: atp_c + dhpt_c + glu__L_c -> adp_c + dhf_c + h_c + pi_c
    L.reaction('2683', 'TM', side='R')   # DHPS2: 4abz_c + 6hmhptpp_c -> dhpt_c + ppi_c
    L.reaction('2694', 'TM', side='R')   # ADCL: 4adcho_c -> 4abz_c + h_c + pyr_c
    L.reaction('2721', 'TM', side='R')   # ADCS: chor_c + gln__L_c -> 4adcho_c + glu__L_c
    L.reaction('2579', 'TM', side='R')   # GTPCI: gtp_c + h2o_c -> ahdt_c + for_c + h_c
    L.reaction('2613', 'TM', side='R')   # DHPTPE: ahdt_c <=> dhmptp_c
    L.reaction('2645', 'TM', side='R')   # CPH4S: ahdt_c + h2o_c -> acald_c + cph4_c + h_c + pppi_c
    L.reaction('2670', 'TM', side='R')   # CDGS: cph4_c + h_c -> cdg_c + nh4_c
    L.reaction('2740', 'TM', side='R')   # CCGS: atp_c + cdg_c + nh4_c -> adp_c + h2o_c + h_c + pi_c + preq0_c
    L.reaction('2708', 'TM', side='R')   # DHMPTR: dhmpt_c + h_c + nadph_c -> nadp_c + thmnp_c
    L.reaction('3408', 'TM', side='R')   # DHNPTE: dhnpt_c <=> dhmpt_c
    # ---------------------------------------------------------------- HDDB: 6-hydroxymethyl-dihydropterin diphosphate biosynthesis
    L.setp('HDDB', {
        '1': None,   # gtp_c  (draft 6795, 2552)
        '3': None,   # ahdt_c  (draft 6962, 1620)
        '16': None,   # dhpmp_c  (draft 6962, 2040)
        '29': None,   # dhnpt_c  (draft 7540, 1963)
        '41': None,   # 6hmhpt_c  (draft 8098, 1769)
        '53': None,   # 6hmhptpp_c  (draft 8083, 2152)
        '66': None,   # 25drapp_c  (draft 6795, 2156)
    })
    L.reaction('14', 'HDDB', side='R')   # GTPCI: gtp_c + h2o_c -> ahdt_c + for_c + h_c
    L.reaction('27', 'HDDB', side='R')   # DNTPPA: ahdt_c + h2o_c -> dhpmp_c + h_c + ppi_c
    L.reaction('40', 'HDDB', side='R')   # DNMPPA: dhpmp_c + h2o_c -> dhnpt_c + pi_c
    L.reaction('51', 'HDDB', side='R')   # DHNPA2r: dhnpt_c <=> 6hmhpt_c + gcald_c
    L.reaction('64', 'HDDB', side='R')   # HPPK2: 6hmhpt_c + atp_c -> 6hmhptpp_c + amp_c + h_c
    L.reaction('79', 'HDDB', side='R')   # GTPCII2: gtp_c + 3 h2o_c -> 25drapp_c + for_c + 2 h_c + ppi_c
    # ---------------------------------------------------------------- CAM: coenzyme A metabolism
    L.setp('CAM', {
        '460': None,   # pnto__R_c  (draft 12074, 4156)
        '461': None,   # 4ppan_c  (draft 12074, 4557)
        '475': None,   # 4ppcys_c  (draft 12074, 4953)
        '491': None,   # pan4p_c  (draft 12114, 5233)
        '502': None,   # dpcoa_c  (draft 12114, 5513)
        '516': None,   # coa_c  (draft 11841, 6030)
        '3969': None,   # 2tpr3dpcoa_c  (draft 12216, 6151)
        '3970': None,   # ade_c  (draft 12450, 6267)
    })
    L.reaction('472', 'CAM', side='R')   # PNTK: atp_c + pnto__R_c -> 4ppan_c + adp_c + h_c
    L.reaction('489', 'CAM', side='R')   # PPNCL2: 4ppan_c + ctp_c + cys__L_c -> 4ppcys_c + cmp_c + h_c + ppi_c
    L.reaction('500', 'CAM', side='R')   # PPCDC: 4ppcys_c + h_c -> co2_c + pan4p_c
    L.reaction('513', 'CAM', side='R')   # PTPATi: atp_c + h_c + pan4p_c -> dpcoa_c + ppi_c
    L.reaction('526', 'CAM', side='R')   # DPCOAK: atp_c + dpcoa_c -> adp_c + coa_c + h_c
    L.reaction('3980', 'CAM', side='R')   # TPRDCOAS: atp_c + dpcoa_c -> 2tpr3dpcoa_c + ade_c
    # ---------------------------------------------------------------- AM: ascorbate metabolism
    L.setp('AM', {
        '345': None,   # 3dhguln_c  (draft 7462, 3329)
        '347': None,   # 3dhgulnp_c  (draft 7214, 3477)
        '359': None,   # 23doguln_c  (draft 7980, 3221)
        '373': None,   # xu5p__L_c  (draft 6747, 3453)
        '383': None,   # ascb6p_c  (draft 8200, 3453)
    })
    L.reaction('358', 'AM', side='R')   # 3KGK: 3dhguln_c + atp_c -> 3dhgulnp_c + adp_c + h_c
    L.reaction('371', 'AM', side='R')   # DOGULNR: 23doguln_c + h_c + nadh_c -> 3dhguln_c + nad_c
    L.reaction('382', 'AM', side='R')   # KG6PDC: 3dhgulnp_c + h_c -> co2_c + xu5p__L_c
    L.reaction('393', 'AM', side='R')   # ASCBPL: ascb6p_c + h2o_c -> 3dhgulnp_c + h_c
    # ---------------------------------------------------------------- BB: biotin biosynthesis
    L.setp('BB', {
        '394': None,   # ala__L_c  (draft 3650, 4230)
        '395': None,   # pimACP_c  (draft 3695, 4597)
        '396': None,   # 8aonn_c  (draft 4134, 4341)
        '412': None,   # dann_c  (draft 4838, 4353)
        '426': None,   # dtbt_c  (draft 5507, 4341)
        '441': None,   # btn_c  (draft 6312, 4353)
    })
    L.reaction('409', 'BB', side='R')   # AOXSr2: ala__L_c + pimACP_c -> 8aonn_c + ACP_c + co2_c
    L.reaction('422', 'BB', side='R')   # AMAOTr: 8aonn_c + amet_c <=> amob_c + dann_c
    L.reaction('439', 'BB', side='R')   # DBTS: atp_c + co2_c + dann_c -> adp_c + dtbt_c + 3 h_c + pi_c
    L.reaction('458', 'BB', side='R')   # BTS5: 2fe2s_c + amet_c + dtbt_c -> 2fe1s_c + btn_c + dad_5_c + h_c + met__L_c
    # ---------------------------------------------------------------- PB: ppGpp biosynthesis
    L.setp('PB', {
        '1001': None,   # gtp_c  (draft 11970, 6632)
        '1003': None,   # gdptp_c  (draft 11965, 7144)
        '1015': None,   # gdp_c  (draft 12261, 6772)
        '1017': None,   # ppgpp_c  (draft 12262, 7284)
        '1041': None,   # gdp_c  (draft 12261, 7564)
    })
    L.reaction('1013', 'PB', side='R')   # GTPDPK: atp_c + gtp_c -> amp_c + gdptp_c + h_c
    L.reaction('1027', 'PB', side='R')   # GDPDPK: atp_c + gdp_c -> amp_c + h_c + ppgpp_c
    L.reaction('1039', 'PB', side='R')   # GTPDPDP: gdptp_c + h2o_c -> h_c + pi_c + ppgpp_c
    L.reaction('1052', 'PB', side='R')   # PPGPPDP: h2o_c + ppgpp_c -> gdp_c + ppi_c
    # ---------------------------------------------------------------- PM3: purine metabolism
    L.setp('PM3', {
        '1725': None,   # ser__D_c  (draft 749, 6840)
        '1727': None,   # pyr_c  (draft 238, 6626)
        '1737': None,   # ser__L_c  (draft 749, 6908)
        '1775': None,   # mlthf_c  (draft 265, 7022)
        '1881': None,   # methf_c  (draft 889, 7096)
    })
    L.reaction('1736', 'PM3', side='R')   # SERD_D: ser__D_c -> nh4_c + pyr_c
    L.reaction('1746', 'PM3', side='R')   # SERD_L: ser__L_c -> nh4_c + pyr_c
    L.reaction('1786', 'PM3', side='R')   # GHMT2r: ser__L_c + thf_c <=> gly_c + h2o_c + mlthf_c
    L.reaction('1892', 'PM3', side='R')   # MTHFD: mlthf_c + nadp_c <=> methf_c + nadph_c
    # ---------------------------------------------------------------- VBM4: vitamin B12 metabolism
    L.setp('VBM4', {
        '2910': None,   # adocbip_c  (draft 12302, 809)
        '2912': None,   # agdpcbi_c  (draft 11855, 785)
        '2924': None,   # dmbzid_c  (draft 10699, 949)
        '2926': None,   # 5prdmbz_c  (draft 11474, 949)
        '2940': None,   # rdmbzi_c  (draft 11925, 1089)
        '2951': None,   # adocbl_c  (draft 10699, 833)
    })
    L.reaction('2923', 'VBM4', side='R')   # ACBIPGT: adocbip_c + gtp_c + h_c -> agdpcbi_c + ppi_c
    L.reaction('2937', 'VBM4', side='R')   # NNDMBRT: dmbzid_c + nicrnt_c -> 5prdmbz_c + h_c + nac_c
    L.reaction('2950', 'VBM4', side='R')   # RZ5PP: 5prdmbz_c + h2o_c -> pi_c + rdmbzi_c
    L.reaction('2962', 'VBM4', side='R')   # ADOCBLS: agdpcbi_c + rdmbzi_c -> adocbl_c + gmp_c + h_c
    # ---------------------------------------------------------------- VBM5: vitamin B6 metabolism
    L.setp('VBM5', {
        '761': None,   # dscl_c  (draft 10052, 6884)
        '764': None,   # scl_c  (draft 9977, 7271)
        '800': None,   # sheme_c  (draft 9819, 7580)
    })
    L.reaction('774', 'VBM5', side='R')   # SHCHD2: dscl_c + nad_c -> h_c + nadh_c + scl_c
    L.reaction('787', 'VBM5', side='R')   # UPP3MT: 2 amet_c + uppg3_c -> 2 ahcys_c + dscl_c + h_c
    L.reaction('809', 'VBM5', side='R')   # SHCHF: fe2_c + scl_c -> 3 h_c + sheme_c
    # ---------------------------------------------------------------- PB2: pantothenate biosynthesis
    L.setp('PB2', {
        '954': None,   # 3mob_c  (draft 6787, 3982)
        '957': None,   # 2dhp_c  (draft 7519, 3982)
        '972': None,   # pant__R_c  (draft 7974, 3982)
        '986': None,   # pnto__R_c  (draft 8592, 4021)
    })
    L.reaction('969', 'PB2', side='R')   # MOHMT: 3mob_c + h2o_c + mlthf_c -> 2dhp_c + thf_c
    L.reaction('982', 'PB2', side='R')   # DPR: 2dhp_c + h_c + nadph_c -> nadp_c + pant__R_c
    L.reaction('999', 'PB2', side='R')   # PANTS: ala_B_c + atp_c + pant__R_c -> amp_c + h_c + pnto__R_c + ppi_c
    # ---------------------------------------------------------------- LB: lipoate biosynthesis
    L.setp('LB', {
        '866': None,   # octa_c  (draft 12114, 2215)
        '3823': None,   # octapb_c  (draft 12182, 1978)
    })
    L.reaction('881', 'LB', side='R')   # AACPS9: ACP_c + atp_c + octa_c -> amp_c + ocACP_c + ppi_c
    L.reaction('3835', 'LB', side='R')   # OCTNLL: atp_c + h_c + octa_c -> amp_c + octapb_c + ppi_c
    # ---------------------------------------------------------------- LB2: lipoate biosynthesis
    L.setp('LB2', {
        '883': None,   # lipoate_c  (draft 12038, 3244)
        '884': None,   # lipoamp_c  (draft 12038, 2848)
        '4036': None,   # lipopb_c  (draft 12195, 2616)
    })
    L.reaction('895', 'LB2', side='R')   # LIPATPT: atp_c + lipoate_c -> lipoamp_c + ppi_c
    L.reaction('4045', 'LB2', side='R')   # LIPAMPL: lipoamp_c -> amp_c + lipopb_c
    # ---------------------------------------------------------------- PM4: purine metabolism
    L.setp('PM4', {
        '1473': None,   # ac_c  (draft 1122, 6444)
        '1475': None,   # actp_c  (draft 816, 6444)
        '1487': None,   # accoa_c  (draft 417, 6510)
    })
    L.reaction('1486', 'PM4', side='R')   # ACKr: ac_c + atp_c <=> actp_c + adp_c
    L.reaction('1499', 'PM4', side='R')   # PTAr: accoa_c + pi_c <=> actp_c + coa_c
    # ---------------------------------------------------------------- VBM6: vitamin B12 metabolism
    L.setp('VBM6', {
        '2843': None,   # dscl_c  (draft 11725, 105)
        '2846': None,   # scl_c  (draft 11092, 154)
    })
    L.reaction('2856', 'VBM6', side='R')   # SHCHD2: dscl_c + nad_c -> h_c + nadh_c + scl_c
    L.reaction('2869', 'VBM6', side='R')   # UPP3MT: 2 amet_c + uppg3_c -> 2 ahcys_c + dscl_c + h_c
    # ---------------------------------------------------------------- AM2: ascorbate metabolism
    L.setp('AM2', {
        '331': None,   # gal1p_p  (draft 7463, 2917)
        '333': None,   # gal_p  (draft 7885, 2917)
    })
    L.reaction('344', 'AM2', side='R')   # GAL1PPpp: gal1p_p + h2o_p -> gal_p + pi_p
    # ---------------------------------------------------------------- PM5: purine metabolism
    L.setp('PM5', {
        '1530': None,   # gmp_p  (draft 1660, 6158)
        '1532': None,   # gsn_p  (draft 1251, 6908)
    })
    L.reaction('1543', 'PM5', side='R')   # NTD9pp: gmp_p + h2o_p -> gsn_p + pi_p
    # ---------------------------------------------------------------- PM6: purine metabolism
    L.setp('PM6', {
        '1570': None,   # xmp_p  (draft 2309, 6310)
        '1572': None,   # xtsn_p  (draft 2132, 6730)
    })
    L.reaction('1582', 'PM6', side='R')   # NTD10pp: h2o_p + xmp_p -> pi_p + xtsn_p
    # ---------------------------------------------------------------- PM7: purine metabolism
    L.setp('PM7', {
        '1609': None,   # imp_p  (draft 3107, 6286)
        '1610': None,   # ins_p  (draft 2810, 6658)
    })
    L.reaction('1621', 'PM7', side='R')   # NTD11pp: h2o_p + imp_p -> ins_p + pi_p
    # ---------------------------------------------------------------- PM8: purine metabolism
    L.setp('PM8', {
        '1647': None,   # amp_p  (draft 4160, 6402)
        '1649': None,   # adn_p  (draft 4134, 6798)
    })
    L.reaction('1660', 'PM8', side='R')   # NTD7pp: amp_p + h2o_p -> adn_p + pi_p
    # ---------------------------------------------------------------- FB2: flavin biosynthesis
    L.setp('FB2', {
        '2460': None,   # cmp_p  (draft 2894, 2718)
        '2462': None,   # cytd_p  (draft 3237, 2718)
    })
    L.reaction('2473', 'FB2', side='R')   # NTD4pp: cmp_p + h2o_p -> cytd_p + pi_p
    # ---------------------------------------------------------------- TM2: tetrahydrofolate metabolism
    L.setp('TM2', {
        '2646': None,   # mlthf_c  (draft 5587, 105)
        '2648': None,   # 5mthf_c  (draft 5621, 617)
    })
    L.reaction('2659', 'TM2', side='R')   # MTHFR2: 2 h_c + mlthf_c + nadh_c -> 5mthf_c + nad_c
    # ---------------------------------------------------------------- VBM7: vitamin B12 metabolism
    L.setp('VBM7', {
        '2896': None,   # adocbi_c  (draft 11706, 669)
        '2898': None,   # adocbip_c  (draft 11323, 833)
    })
    L.reaction('2909', 'VBM7', side='R')   # ADOCBIK: adocbi_c + atp_c -> adocbip_c + adp_c + h_c
    # ---------------------------------------------------------------- FB3: flavin biosynthesis
    L.setp('FB3', {
        '3235': None,   # fmn_p  (draft 4371, 3522)
        '3238': None,   # ribflv_p  (draft 4719, 3522)
    })
    L.reaction('3248', 'FB3', side='R')   # ACP1p: fmn_p + h2o_p -> pi_p + ribflv_p
    # ---------------------------------------------------------------- PM9: purine metabolism
    L.setp('PM9', {
        '3679': None,   # dump_p  (draft 1530, 3157)
        '3681': None,   # duri_p  (draft 1985, 3419)
    })
    L.reaction('3692', 'PM9', side='R')   # NTD1pp: dump_p + h2o_p -> duri_p + pi_p
    # ---------------------------------------------------------------- PM10: purine metabolism
    L.setp('PM10', {
        '3705': None,   # dcmp_p  (draft 2330, 4119)
        '3707': None,   # dcyt_p  (draft 2183, 4467)
    })
    L.reaction('3718', 'PM10', side='R')   # NTD3pp: dcmp_p + h2o_p -> dcyt_p + pi_p
    # ---------------------------------------------------------------- PM11: purine metabolism
    L.setp('PM11', {
        '3731': None,   # dtmp_p  (draft 1075, 3863)
        '3734': None,   # thymd_p  (draft 790, 3093)
    })
    L.reaction('3744', 'PM11', side='R')   # NTD5pp: dtmp_p + h2o_p -> pi_p + thymd_p
    # ---------------------------------------------------------------- PM12: purine metabolism
    L.setp('PM12', {
        '3783': None,   # dgmp_p  (draft 1741, 5634)
        '3785': None,   # dgsn_p  (draft 1834, 5450)
    })
    L.reaction('3796', 'PM12', side='R')   # NTD8pp: dgmp_p + h2o_p -> dgsn_p + pi_p

NEW_EXITS = []          # (owner node, text) for copies intentionally not merged
EXTRA_TITLES = []       # dict(text=..., anchor=node_id, offset=(dx, dy)) or pos=(x, y)
TITLE_POS = {}          # zone -> (x, y) global position of the zone title
CURVES = dict(t=0.6, stacked=[], bend_override={}, special_links={})
PACK = dict(keep_out=[], auto_ring_keep_out=True)
