"""SPEC - iML1515 'nucleotides and cofactors' split map (BRENDA draft, 305 reactions) -> hand layout.

Every coordinate of a PRIMARY metabolite below was chosen by hand (zone-local coordinates, x right,
y DOWN, ~300-700 between neighbours). The toolkit only places markers/cofactors, curves, labels and
packs the detached pieces (most cofactor pathways are not connected to each other).

Zones
  PUR  purine metabolism: NMP row AMP-IMP-XMP-GMP, nucleosides below, bases below that,
       de novo IMP synthesis as a column above IMP, ADP/ATP + dA nucleotides left, GDP/GTP + dG right
  PYR  pyrimidine metabolism: de novo column ending in UMP, UMP-UDP-UTP-CTP-CDP row, salvage
       rectangle below, deoxy nucleotides (dU/dT/dC) above
  NAD  NAD metabolism: de novo column asp -> NAD, NAD/NADP lanes, salvage left
  PTE/THF/FLAV  GTP-derived cofactors around one shared GTP node: pterins + THF + preQ0 above,
       riboflavin/FMN/FAD below
  FRED flavin-dependent reductases (siderophore iron reduction, sulfonate monooxygenases); every
       reaction is drawn as its own row with its own fmnh2/fmn/fe2 copies (no hairball)
  HEME, B6, B1, VK, MOCO, COA, BIO, PPG, B12, ASC, LIP, ACE, C1, PERI  (see comments)
"""

INPUT = '../../ecmapper/BRENDAmap/draft_map_iML1515/split/nucleotides_and_cofactors.escher.json'

ZONES = {
 'PUR':  dict(title='purine metabolism', title_id='4090', anchor=('2028', (-750, -250)), color='#1f77b4'),
 'PYR':  dict(title='pyrimidine metabolism', title_id='4091', anchor=('2044', (-900, -150)), color='#2ca02c'),
 'NAD':  dict(title='NAD metabolism', title_id='4080', anchor=('169', (-700, -150)), color='#d62728'),
 'PTE':  dict(title='6-hydroxymethyl-dihydropterin diphosphate biosynthesis', title_id='4079',
              anchor=('1', (-1650, -2300)), color='#17becf'),
 'THF':  dict(title='tetrahydrofolate metabolism', title_id='4092', anchor=('2555', (-750, -150)), color='#17becf'),
 'FLAV': dict(title='flavin biosynthesis', title_id='4084', anchor=('529', (-750, 150)), color='#bcbd22'),
 'FRED': dict(title='flavin reductases: siderophore iron and sulfonate reduction', anchor=('2870', (-100, -250)),
              color='#7f7f7f'),
 'HEME': dict(title='heme metabolism', title_id='4085', anchor=('672', (-400, -150)), color='#8c564b'),
 'B6':   dict(title='vitamin B6 metabolism', title_id='4095', anchor=('2963', (-450, -150)), color='#9467bd'),
 'B1':   dict(title='vitamin B1 metabolism', title_id='4093', anchor=('2742', (-300, -200)), color='#e377c2'),
 'VK':   dict(title='vitamin K metabolism', title_id='4096', anchor=('3075', (-100, -250)), color='#ff7f0e'),
 'MOCO': dict(title='molybdenum cofactor biosynthesis', title_id='4087', anchor=('896', (-500, -150)), color='#7f7f7f'),
 'COA':  dict(title='pantothenate and coenzyme A biosynthesis', title_id='4088', anchor=('954', (-100, -250)),
              color='#8c564b'),
 'BIO':  dict(title='biotin biosynthesis', title_id='4082', anchor=('394', (-100, -250)), color='#2ca02c'),
 'PPG':  dict(title='ppGpp biosynthesis', title_id='4089', anchor=('1001', (-100, -250)), color='#d62728'),
 'B12':  dict(title='vitamin B12 metabolism', title_id='4094', anchor=('2896', (-100, -250)), color='#9467bd'),
 'ASC':  dict(title='ascorbate metabolism', title_id='4081', anchor=('359', (-100, -550)), color='#1f77b4'),
 'LIP':  dict(title='lipoate biosynthesis', title_id='4086', anchor=('866', (-100, -250)), color='#e377c2'),
 'ACE':  dict(title='acetate kinase / phosphotransacetylase', anchor=('1473', (-100, -250)), color='#ff7f0e'),
 'C1':   dict(title='serine and one-carbon (THF) interconversion', anchor=('1725', (-100, -250)), color='#17becf'),
 'PERI': dict(title='periplasmic nucleotidases and phosphatases', anchor=('1530', (-100, -250)), color='#555555'),
}

EXIT_OWNER = {}        # automatic owner rule

# FRED rows: (R2 = fmnh2 reductant, R3 = reduced riboflavin reductant, X_fe3, X_un)
SIDERO = [('3285', '3298', '3272', '3273'), ('3376', '3387', '3364', '3365'), ('3501', '3512', '3489', '3490'),
          ('3525', '3536', '3513', '3514'), ('3549', '3560', '3537', '3538'), ('3573', '3584', '3561', '3562')]
FDMO = [('3426', '3409', '3411'), ('3443', '3427', '3429'), ('3458', '3444', '3445'), ('3473', '3459', '3460'),
        ('3488', '3474', '3475')]          # (rid, sulfonate, aldehyde)

def prepare(L):
    S = L.S
    # purine: AMPN gets its own ade/r5p (short fan left of AMP); 5DOAN its own ade; PUNP4 its own 2dr1p
    S['ade_N'] = L.split('1676', '1904'); S['r5p_A'] = L.split('1894', '1904')
    S['ade_D'] = L.split('1676', '4021')
    S['2dr1p_G'] = L.split('1053', '1099')
    # pyrimidine: DCYTD with its own dcyt next to duri; TMDPP own 2dr1p; dTDP-sugar transferases own dtdp
    S['dcyt_U'] = L.split('2318', '2331'); S['2dr1p_T'] = L.split('2346', '2421')
    S['dtdp_R'] = L.split('2516', '3928'); S['dtdp_A'] = L.split('2516', '4034')
    # folate: 4-aminobenzoyl-glutamate hydrolase must not pull heme's glutamate
    S['glu_F'] = L.split('672', '3209')
    # B6: transaminations ALATA_L2/_D2 get their own ala / pyr copies (cysteine desulfurases keep ala)
    S['ala_A'] = L.split('3261', '3271'); S['pyr_D'] = L.split('3250', '3260')
    # MoCo: cu2 product copies; the tungsten branch gets its own MPT-AMP copy (detached piece)
    S['cu2_B'] = L.split('3325', '3336'); S['cu2_W'] = L.split('3325', '3995'); S['cu2_X'] = L.split('3325', '3350')
    S['mptamp_W'] = L.split('912', '3995'); L.merge(S['mptamp_W'], L.split('912', '3350'))
    # FRED: every reaction gets its own reductant / product copies
    for r2, r3, xf, xu in SIDERO:
        S['fmnh2_'+r2] = L.split('2872', r2); S['fmn_'+r2] = L.split('2870', r2)
        if r2 != '3285': S['fe2_'+r2] = L.split('3274', r2)
        S['fe2_'+r3] = L.split('3274', r3)
        if r3 != '3298':
            S['rbflvrd_'+r3] = L.split('3286', r3); S['ribflv_'+r3] = L.split('3287', r3)
    for rid, x, a in FDMO:
        S['fmnh2_'+rid] = L.split('2872', rid); S['fmn_'+rid] = L.split('2870', rid)
        if rid != '3426': S['so3_'+rid] = L.split('3413', rid)
    S['fald_M'] = L.split('3429', '3666')

MERGES = [('1', ['527', '2566']),        # GTP: pterin, THF (queuosine) and flavin branches share one GTP
          ('3', ['2568']), ('29', ['3400']), ('53', ['2672']), ('529', ['66']),
          ('2671', ['3199']),             # 4-aminobenzoate
          ('986', ['460']),               # pantothenate -> CoA
          ('1015', ['1041']),             # ppGpp cycle
          ('2898', ['2910']),             # adocbip
          ('1775', ['2646']),             # 5,10-methylene-THF
          ('761', ['2843']), ('764', ['2846']),                     # siroheme branch
          ('2118', ['2448']), ('2120', ['2359']), ('2157', ['2182'])]  # cytd, cmp, csn
HUBS = ['1', '1303', '1362', '1330', '2146', '80', '1153', '1243']

ORIGIN = {'PUR': (2700, 4000), 'PYR': (7100, 1000), 'NAD': (12600, 300),
          'PTE': (14700, 3300), 'THF': (14700, 3300), 'FLAV': (14700, 3300),
          'HEME': (17000, 300), 'B6': (8300, 5600), 'FRED': (300, 6900),
          'VK': (10400, 6600), 'MOCO': (12000, 3700), 'B1': (17100, 5000), 'COA': (13600, 7500),
          'BIO': (10400, 7400), 'PPG': (6200, 8000), 'B12': (11300, 8200), 'ASC': (13600, 8500),
          'LIP': (15600, 8400), 'ACE': (17000, 8400), 'C1': (12000, 9300), 'PERI': (6400, 9200)}


def layout(L):
    S = L.S
    R = L.reaction
    # ======================================================================== PUR  purine metabolism
    L.setp('PUR', {
        # NMP row / nucleoside row / base row
        '1362': (0, 0), '1303': (700, 0), '1441': (1400, 0), '1330': (2100, 0),          # amp imp xmp gmp
        '1622': (0, 550), '1584': (700, 550), '1558': (1400, 550), '1518': (2100, 550),  # adn ins xtsn gsn
        '1676': (0, 1100), '1714': (700, 1100), '1690': (1400, 1100), '1747': (2100, 1100),  # ade hxan xan gua
        '1988': (1400, 1550),                                                            # urate
        '1347': (350, -300),                                                             # dcamp
        # de novo IMP synthesis: column above IMP
        '1220': (700, -400), '1219': (700, -780), '1276': (700, -1160), '1313': (700, -1540),
        '1408': (700, -1920), '1500': (700, -2300), '1662': (700, -2680), '1860': (700, -3060),
        '1907': (700, -3460), '2011': (700, -3840), '2028': (700, -4220),
        # adenylates + dA
        '1153': (0, -550), '1243': (0, -1100), '1232': (420, -550),                      # adp atp ap4a
        '1116': (-600, -550), '1129': (-600, -1100), '1114': (-1200, -550),               # dadp datp damp
        '1100': (-1200, 0), '1168': (-1200, 550), '1053': (-1750, 275), '1054': (-2250, 275),
        '1063': (-2650, 75), '1064': (-2650, 475),                                        # acald g3p
        S['ade_N']: (-700, -300), S['r5p_A']: (-700, -50),
        # guanylates + cGMP + dG
        '1287': (2100, -550), '1385': (2100, -1100), '1417': (1650, -800),                # gdp gtp 35cgmp
        '1195': (2700, -550), '1074': (2700, -1100), '1141': (3300, -550), '1076': (3300, 0),
        S['2dr1p_G']: (3300, 550),
        # ribose-5-phosphate cluster (bottom left) and 5'-deoxyadenosine nucleosidase
        '2000': (-1300, 1250), '1894': (-1300, 1750), '1922': (-1800, 1500), '3353': (-800, 1750),
        '3996': (-1300, 2250), '3998': (-1800, 2050),
        '4009': (200, 1700), '4011': (650, 1550), S['ade_D']: (650, 1850)})
    Z = 'PUR'
    # NMP <-> nucleoside lanes, nucleoside -> base lanes
    R('1646', Z, shift=(75, 0), side='R', dist=100); R('1634', Z, shift=(-75, 0), side='L', dist=100)   # NTD7 ADNK1
    R('1607', Z, shift=(75, 0), side='R', dist=100); R('1595', Z, shift=(-75, 0), side='L', dist=100)   # NTD11 INSK
    R('1568', Z, side='R', dist=100)                                                                    # NTD10
    R('1529', Z, shift=(-75, 0), side='L', dist=100); R('1555', Z, shift=(75, 0), side='R', dist=100)   # NTD9 GSNK
    R('1771', Z, shift=(75, 0), side='R', dist=100); R('1798', Z, shift=(-75, 0), side='L', dist=100)   # PUNP1 ADNUC
    R('1810', Z, shift=(75, 0), side='R', dist=100); R('1822', Z, shift=(-75, 0), side='L', dist=100)   # PUNP5 INSH
    R('1846', Z, shift=(75, 0), side='R', dist=100); R('1858', Z, shift=(-75, 0), side='L', dist=100)   # PUNP7 XTSNH
    R('1834', Z, side='R', dist=100)                                                                    # PUNP3
    R('1713', Z, side='U', dist=100)                                                                    # ADA
    R('1946', Z, side='D'); R('1972', Z, side='D'); R('1958', Z, side='D')                              # ADD HXAND GUAD
    R('1999', Z, side='R')                                                                              # XAND
    # phosphoribosyltransferases: bypass lanes between the columns
    R('1688', Z, axis=(0, -1), mid=(-400, 550), side='L', dist=100)       # ADPT  ade -> amp
    R('1724', Z, axis=(0, -1), mid=(1050, 550), side='R', dist=100)       # HXPRT hxan -> imp
    R('1701', Z, axis=(0, -1), mid=(1750, 550), side='R', dist=100)       # XPPT  xan -> xmp
    R('1759', Z, axis=(0, -1), mid=(2450, 550), side='R', dist=100)       # GUAPRT gua -> gmp
    # IMP branch point
    R('1361', Z, side='D', dist=110); R('1372', Z, side='D', dist=100)    # ADSS ADSL1r (via dcamp)
    R('1472', Z, side='U'); R('1458', Z, side='U', step=55)               # IMPD GMPS2
    R('1344', Z, axis=(-1, 0), mid=(1400, -400), side='U', dist=110)      # GMPR gmp -> imp (above xmp)
    # de novo column (cofactors right)
    R('2042', Z, side='R'); R('2025', Z, side='R')                        # GLUPRT PRAGSr
    R('1921', Z, shift=(75, 0), side='R'); R('1984', Z, shift=(-75, 0), side='L')   # GART GARFT
    R('1879', Z, side='R', step=55); R('1675', Z, side='R'); R('1516', Z, side='R')  # PRFGS PRAIS AIRC2
    R('1416', Z, side='R'); R('1329', Z, side='R'); R('1286', Z, side='R')           # AIRC3 PRASCSi ADSL2r
    R('1231', Z, side='R'); R('1312', Z, side='R')                                   # AICART IMPC
    # adenylates
    R('1382', Z, side='L', dist=110)                                      # ADK1 amp -> adp
    for rid, y in (('1255', -720), ('1275', -930)):                       # ATPM, NTP1 stacked
        R(rid, Z, side='R', axis=(0, 1), mid=(0, y), gap=40)
    R('1265', Z, shift=(-110, 0), side='L', dist=100)                     # ATPS4rpp
    R('1242', Z, side='U')                                                # AP4AH
    R('1904', Z, axis=(-1, -0.15), mid=(-330, -130), side='U', dist=100)  # AMPN fan
    R('1167', Z, side='U'); R('1139', Z, side='L'); R('1126', Z, side='U')   # RNDR1 NDPK8 DADK
    R('3808', Z, side='L', dist=150); R('3756', Z, side='L')              # NTPP5 NTD6
    R('1192', Z, side='L'); R('1112', Z, side='U'); R('1180', Z, side='D')  # DADA PUNP2 PUNP6
    R('1062', Z); R('1073', Z)                                            # PPM2 DRPA
    # guanylates
    R('1407', Z, side='R'); R('1395', Z, side='R')                        # GK1 NDPK1
    R('1427', Z, side='L'); R('1437', Z, side='L')                        # GUACYC PDE4
    R('1301', Z, side='U'); R('1217', Z, side='L'); R('1205', Z, side='D')  # RNDR2 NDPK5 DGK1
    R('1152', Z, side='U', dist=110)                                      # NTPP1 dgtp -> dgmp
    R('3782', Z, side='L'); R('1099', Z, side='L')                        # NTD8 PUNP4
    R('1087', Z, axis=(0, 1), mid=(3750, -550), side='R')                 # NTPTP1 dgtp -> dgsn (around dgmp)
    # ribose-5-phosphate cluster, 5DOAN
    R('2008', Z); R('1934', Z); R('4008', Z, side='D')                    # PPM ADPRDP YUMPS
    R('3363', Z, axis=(0, 1), mid=(-1050, 1400), end_cof=True, side='R')  # CMPN (cmp not drawn as primary)
    R('4021', Z, side='L')                                                # 5DOAN

    # ======================================================================== PYR  pyrimidine metabolism
    L.setp('PYR', {
        '2044': (0, 0), '2045': (0, 350), '2057': (0, 750), '2069': (0, 1380), '2105': (0, 1700), '2146': (0, 2050),
        '2246': (650, 2050), '2271': (1300, 2050), '2300': (1950, 2050), '2233': (2600, 2050),
        '2120': (2275, 2600), '2118': (2275, 3150), '2134': (0, 3150), '2171': (0, 3700), '2157': (2275, 3700),
        '2283': (650, 1450), '2334': (650, 900), '2373': (650, 350), '2320': (650, -200), S['dcyt_U']: (650, -750),
        '2346': (150, -450),
        '2436': (1300, 350), '2410': (1300, -200), S['2dr1p_T']: (1100, -650), '2411': (1500, -650),
        '2541': (1950, 350), '2516': (1950, -200), '3904': (2600, -200), '3905': (2600, 350),
        '2488': (3400, 1450), '2396': (3400, 900), '2475': (3400, 350), '2318': (3400, -200),
        '3916': (3100, 2700), '3917': (3100, 3100), S['dtdp_R']: (3700, 2700), '3918': (3700, 3100),
        '4022': (3100, 3550), '4023': (3100, 3950), S['dtdp_A']: (3700, 3550), '4024': (3700, 3950)})
    Z = 'PYR'
    R('2056', Z, side='R'); R('2067', Z, side='R')                        # ASPCT DHORTS
    for rid, y in (('2080', 860), ('2092', 1030), ('2104', 1200)):       # DHORD2/5/fum stacked
        R(rid, Z, side='R', axis=(0, 1), mid=(0, y), gap=40)
    R('2117', Z, side='R'); R('2155', Z, side='R')                        # ORPT OMPDC
    R('2256', Z, shift=(0, -75), side='U', dist=100); R('2268', Z, shift=(0, 75), side='D', dist=100)  # UMPK x2
    R('2281', Z, side='D'); R('2317', Z, side='D', step=55); R('2395', Z, side='U')   # NDPK2 CTPS2 NDPK3
    R('2370', Z, side='L'); R('2243', Z, side='R')                        # NTPP4 CYTK1
    R('2131', Z, shift=(75, 0), side='R', dist=100); R('2459', Z, shift=(-75, 0), side='L', dist=100)  # CYTDK2 NTD4
    R('2144', Z, side='U'); R('2168', Z, side='R'); R('2194', Z, side='D')   # CYTD CYTDH CSND
    R('2181', Z, shift=(75, 0), side='R', dist=100); R('2206', Z, shift=(-75, 0), side='L', dist=100)  # PYNP2r URIH
    R('2218', Z, side='R')                                                # URIK2
    R('2230', Z, axis=(0, -1), mid=(-550, 2600), side='L')                # UPPRT ura -> ump (bypass)
    # deoxy
    R('2296', Z, side='R'); R('2344', Z, side='R'); R('2433', Z, side='R')   # RNDR4 NDPK6 DUTPDP
    R('2408', Z, side='D')                                                # DCTPD dctp -> dutp (long)
    R('2383', Z, shift=(-75, 0), side='L', dist=100); R('3678', Z, shift=(75, 0), side='R', dist=100)  # DURIK1 NTD1
    R('2331', Z, side='L'); R('2357', Z, side='L')                        # DCYTD DURIPP
    R('2446', Z, side='D')                                                # TMDS
    R('2513', Z, shift=(-90, 0), side='L', dist=90); R('3730', Z, shift=(-260, 0), side='L', dist=80)  # TMDK1 NTD5
    R('2421', Z, axis=(0, -1), side='L')                                  # TMDPP
    R('2526', Z, side='U'); R('2551', Z, side='R', dist=110); R('3820', Z, side='D')   # DTMPK NDPK4 NTPP7
    R('3915', Z)                                                          # PYK6 hourglass
    R('2501', Z, side='R'); R('2538', Z, side='L'); R('2486', Z, side='L'); R('3704', Z, side='L')  # RNDR3 NDPK7 NTPP3 NTD3
    R('3399', Z, axis=(0, 1), mid=(3800, 900), side='R')                  # CYTK2 dcmp -> dcdp (bypass)
    R('3928', Z); R('4034', Z)                                            # RHAT1 AADDGT

    # ======================================================================== NAD
    L.setp('NAD', {'169': (0, -150), '172': (0, 600), '247': (0, 900), '273': (0, 1200), '153': (0, 1500),
                   '80': (0, 2000), '83': (700, 2000), '95': (350, 2700),
                   '198': (-700, 2000), '185': (-700, 2700), '233': (-1400, 2700), '260': (-1400, 1200)})
    Z = 'NAD'
    for rid, y in (('182', 50), ('3310', 230), ('3322', 410)):          # ASPO6/3/4 stacked
        R(rid, Z, side='R', axis=(0, 1), mid=(0, y), gap=40)
    R('258', Z, side='R'); R('330', Z, side='R'); R('284', Z, side='R'); R('168', Z, side='L', step=55)
    R('130', Z, shift=(0, -110), side='U', dist=60)                       # NADK
    R('93', Z, shift=(0, -240), side='U', dist=60)                        # NADTRHD
    R('118', Z, shift=(0, -370), side='U', dist=60)                       # NADPPPS
    R('141', Z, mid=(390, 2140))                                          # NADTRHD (nad + nadph -> nadp)
    R('151', Z, mid=(260, 2320))                                          # THD2pp
    R('106', Z, side='R')                                                 # NADTRHD nadph -> nadp
    R('208', Z, shift=(0, -75), side='U', dist=100); R('220', Z, shift=(0, 75), side='D', dist=100)   # NADDP NMNAT
    R('296', Z, side='L'); R('195', Z, side='R'); R('232', Z, side='R'); R('243', Z, side='D')
    R('271', Z, side='L'); R('316', Z, side='U', step=60)                 # NNAM NAMNPP

    # ======================================================================== PTE / THF / FLAV around GTP
    L.setp('PTE', {'1': (0, 0), '3': (0, -450), '16': (0, -800), '29': (0, -1150), '41': (0, -1500),
                   '53': (0, -1850), '2605': (-550, -450), '2633': (550, -650), '2660': (1000, -650),
                   '2727': (1450, -650), '2695': (-550, -1150), '2698': (-1000, -1150)})
    Z = 'PTE'
    R('14', Z, shift=(-75, 0), side='L', dist=100); R('2579', Z, shift=(75, 0), side='R', dist=100)  # GTPCI x2
    R('27', Z, side='L'); R('40', Z, side='L'); R('51', Z, side='L'); R('64', Z, side='L')
    R('2613', Z); R('2645', Z, side='D'); R('2670', Z, side='D'); R('2740', Z, side='D', step=55)
    R('3408', Z); R('2708', Z, side='U')
    L.setp('THF', {'2671': (600, -1850), '2615': (300, -2250), '2684': (1100, -1850), '2709': (1100, -1400),
                   '3197': (600, -1350), S['glu_F']: (900, -1550), '2580': (300, -2650), '2555': (300, -3050),
                   '2552': (900, -3050)})
    Z = 'THF'
    R('2683', Z, side='L'); R('2694', Z, side='U'); R('2721', Z, side='R')            # DHPS2 ADCL ADCS
    R('3209', Z, axis=(0, -1), mid=(620, -1600), side='L')                            # 4ABZGLUH
    R('2630', Z, side='R'); R('2592', Z, side='R')                                    # DHFS DHFR
    R('2565', Z, shift=(0, -75), side='U', dist=100); R('2604', Z, shift=(0, 75), side='D', dist=100)  # GARFT FTHFD
    L.setp('FLAV', {'529': (0, 450), '544': (0, 800), '557': (0, 1150), '570': (0, 1500), '634': (500, 1150),
                    '594': (500, 1500), '582': (250, 1900), '583': (250, 2300), '610': (250, 2800),
                    '646': (250, 3200), '659': (250, 3600), '623': (-450, 2300)})
    Z = 'FLAV'
    R('79', Z, shift=(-75, 0), side='L', dist=100); R('542', Z, shift=(75, 0), side='R', dist=100)  # GTPCII2 x2
    R('555', Z, side='L'); R('568', Z, side='L'); R('581', Z, side='L'); R('644', Z, side='R')
    R('607', Z, side='R'); R('593', Z, side='R')                                      # RBFSa RBFSb
    R('620', Z, shift=(75, 0), side='R', dist=100); R('3234', Z, shift=(-75, 0), side='L', dist=100)  # RBFK ACP1_FMN
    R('657', Z, side='R'); R('670', Z, side='R')                                      # FMNAT FADRx2
    for rid, x in (('633', 90), ('3596', -170)):                                     # FLVR FLVRx stacked
        R(rid, Z, side='U', axis=(-1, 0), mid=(x, 2300), gap=40)

    # ======================================================================== FRED flavin reductases
    pos = {'2870': (0, 0), '2872': (600, 0)}
    for i, (r2, r3, xf, xu) in enumerate(SIDERO):
        x0, y0 = 1800*(i % 3), 1050 + 950*(i // 3)
        pos.update({xf: (x0, y0), xu: (x0+1000, y0),
                    S['fmnh2_'+r2]: (x0, y0-300), S['fmn_'+r2]: (x0+1000, y0-330),
                    S.get('fe2_'+r2, '3274'): (x0+1250, y0-170),
                    S.get('rbflvrd_'+r3, '3286'): (x0, y0+300), S.get('ribflv_'+r3, '3287'): (x0+1000, y0+330),
                    S['fe2_'+r3]: (x0+1250, y0+170)})
    for i, (rid, x, a) in enumerate(FDMO):
        x0, y0 = 5600, 250 + 650*i
        pos.update({x: (x0, y0), S['fmnh2_'+rid]: (x0, y0-250), a: (x0+1000, y0),
                    S['fmn_'+rid]: (x0+1000, y0-250), S.get('so3_'+rid, '3413'): (x0+1250, y0-125)})
    pos.update({'3650': (5600, 3600), S['fald_M']: (6600, 3450), '3654': (6600, 3750)})
    L.setp('FRED', pos)
    Z = 'FRED'
    R('2883', Z, shift=(0, -75), side='U', dist=100); R('2895', Z, shift=(0, 75), side='D', dist=100)  # FMNRx x2
    for r2, r3, xf, xu in SIDERO:
        R(r2, Z, side='U'); R(r3, Z, side='D')
    for rid, x, a in FDMO:
        R(rid, Z, side='D', dist=110)
    R('3666', Z, side='D')                                                            # MTRPOX

    # ======================================================================== HEME
    L.setp('HEME', {'672': (0, 0), '675': (0, 350), '690': (0, 700), '704': (0, 1050), '726': (0, 1400),
                    '713': (0, 1750), '715': (0, 2100), '789': (0, 2450), '814': (0, 3200), '841': (0, 3650),
                    '854': (0, 4000), '3597': (450, 4000), '3599': (225, 4400),
                    '3880': (-500, 3200), '3881': (-500, 3650),
                    '761': (1100, 2100), '764': (1100, 2650), '800': (1100, 3100)})
    Z = 'HEME'
    for rid in ('688', '703', '712', '747', '737', '724', '798'):
        R(rid, Z, side='R')
    for rid, y in (('826', 2650), ('839', 2980)):                        # CPPPGO2 CPPPGO stacked
        R(rid, Z, side='R', axis=(0, 1), mid=(0, y), gap=40)
    R('852', Z, side='R'); R('3891', Z, side='L'); R('863', Z, side='R'); R('3611', Z, side='R')
    R('760', Z, axis=(-1, 0), mid=(-350, 2100), end_cof=True, side='U')   # UPP3MT uppg3 -> (dscl dot)
    R('787', Z, axis=(1, 0), mid=(750, 1950), end_cof=True, side='U')     # UPP3MT -> dscl
    R('2869', Z, axis=(1, 0), mid=(750, 2250), end_cof=True, side='D')    # UPP3MT -> dscl (duplicate)
    R('774', Z, shift=(-75, 0), side='L', dist=100); R('2856', Z, shift=(75, 0), side='R', dist=100)  # SHCHD2 x2
    R('809', Z, side='R')                                                 # SHCHF

    # ======================================================================== B6
    L.setp('B6', {'2963': (0, 0), '2966': (0, 350), '2981': (0, 700), '2994': (0, 1050), '3210': (-450, 1050),
                  '3848': (450, 1050), '3022': (225, 1450), '3036': (925, 1450), '3008': (1625, 1450),
                  '3020': (225, 2050), '3034': (925, 2050), '3006': (1625, 2050),
                  S['ala_A']: (1050, 900), '3250': (1500, 900), '3249': (1050, 1950), S['pyr_D']: (1500, 1950),
                  '3624': (-1150, 350), '3625': (-1150, 100), '3929': (-1150, 600),
                  '3261': (-750, 350), '3626': (-750, 100), '3930': (-750, 600)})
    Z = 'B6'
    R('2978', Z, side='R'); R('2991', Z, side='R'); R('3004', Z, side='R'); R('3222', Z, side='U')
    R('3867', Z, side='L', step=55)                                       # PDX5PS
    for rid, x in (('3074', 475), ('3847', 675)):                        # PDX5POi / PDX5PO2 stacked
        R(rid, Z, side='D', axis=(1, 0), mid=(x, 1450), gap=40, dist=100)
    R('3062', Z, side='U', dist=90)                                       # PYAM5PO
    R('3032', Z, shift=(-75, 0), side='L', dist=100); R('3879', Z, shift=(75, 0), side='R', dist=100)  # PYDXNK PDXPP
    R('3046', Z, side='L', dist=70); R('3903', Z, shift=(-190, 0), side='L', dist=100)  # PYDXK PYDXPP
    R('3018', Z, side='R', dist=70); R('3623', Z, shift=(190, 0), side='R', dist=100)  # PYDAMK HYPOE
    R('3271', Z, axis=(1, 0), mid=(1275, 1180))                           # ALATA_L2 (above)
    R('3260', Z, axis=(1, 0), mid=(1275, 1720))                           # ALATA_D2 (below)
    R('3636', Z); R('3940', Z)                                            # ICYSDS SCYSDS (double hourglass)

    # ======================================================================== B1
    L.setp('B1', {'2742': (0, 0), '2745': (0, 450), '2778': (700, 0), '2780': (700, 450), '3942': (350, 150),
                  '2765': (350, 1150), '3947': (0, 1250), '2763': (800, 1400), '2792': (350, 1650),
                  '2818': (-250, 1650), '2806': (350, 2100)})
    Z = 'B1'
    R('2762', Z, side='L', step=55); R('2790', Z, side='R')               # TYRL DXPS
    R('3967', Z, axis=(0, 1), mid=(350, 780), side='R', step=55)          # THZPSN3
    R('2776', Z, side='D'); R('2803', Z, side='L'); R('2829', Z, side='U'); R('2816', Z, side='R')
    R('2842', Z, axis=(0.6, -1), mid=(120, 1950), end_cof=True, side='L')  # TDP (thmpp not primary)

    # ======================================================================== VK
    L.setp('VK', {'3075': (0, 0), '3076': (450, 0), '3094': (900, 0), '3106': (1350, 0), '3119': (1800, 0),
                  '3122': (1800, 500), '3145': (1350, 500), '3158': (900, 500), '3169': (450, 500), '3186': (0, 500)})
    Z = 'VK'
    R('3084', Z, shift=(0, -60)); R('3092', Z, shift=(0, 60))             # ICHORS x2
    R('3105', Z, side='U'); R('3116', Z, side='U'); R('3144', Z, side='U'); R('3134', Z, side='R')
    R('3155', Z, side='D'); R('3168', Z, side='D'); R('3183', Z, side='D'); R('3196', Z, side='D')

    # ======================================================================== MOCO
    L.setp('MOCO', {'896': (0, 0), '897': (0, 400), '3325': (-400, 250), '3637': (-400, 550), '911': (0, 900),
                    '3638': (-400, 900), '912': (0, 1300), '926': (0, 1750), '3324': (650, 1650), S['cu2_B']: (650, 1400),
                    # tungsten branch (own MPT-AMP copy, detached piece)
                    S['mptamp_W']: (-900, 2300), '3981': (-1350, 2300), '3337': (-900, 2800), S['cu2_W']: (-1350, 2800),
                    '3339': (-350, 2450), S['cu2_X']: (-350, 2750)})
    Z = 'MOCO'
    R('909', Z, side='R'); R('3649', Z, side='L'); R('923', Z, side='R'); R('940', Z, side='L')
    R('953', Z, axis=(0, 1), mid=(0, 2050), end_cof=True, side='R')     # MOCDS
    R('3336', Z, mid=(325, 1525), side='U')                               # BMOCOS
    R('3995', Z, side='L'); R('3350', Z, side='U')                        # WCOS BWCOS

    # ======================================================================== COA
    L.setp('COA', {'954': (0, 0), '957': (400, 0), '972': (800, 0), '986': (1200, 0), '461': (1600, 0),
                   '475': (2000, 0), '491': (2400, 0), '502': (2800, 0), '516': (3200, 0),
                   '3969': (2650, 500), '3970': (3000, 500)})
    Z = 'COA'
    for rid in ('969', '982', '999', '472', '489', '500', '513', '526'):
        R(rid, Z, side='U')
    R('3980', Z, side='L')                                                # TPRDCOAS

    # ======================================================================== BIO, PPG, B12, ASC, LIP, ACE, C1
    L.setp('BIO', {'394': (0, -180), '395': (0, 180), '396': (450, 0), '412': (900, 0), '426': (1350, 0), '441': (1800, 0)})
    R('409', 'BIO', side='D'); R('422', 'BIO', side='U'); R('439', 'BIO', side='U'); R('458', 'BIO', side='U')
    L.setp('PPG', {'1001': (0, 0), '1003': (500, 0), '1017': (500, 500), '1015': (0, 500)})
    R('1013', 'PPG', side='U'); R('1039', 'PPG', side='R')
    R('1027', 'PPG', shift=(0, -75), side='U', dist=90); R('1052', 'PPG', shift=(0, 75), side='D', dist=100)
    L.setp('B12', {'2896': (0, 0), '2898': (450, 0), '2912': (900, 0), '2924': (0, 600), '2926': (450, 600),
                   '2940': (900, 600), '2951': (1350, 300)})
    R('2909', 'B12', side='U'); R('2923', 'B12', side='U'); R('2937', 'B12', side='D'); R('2950', 'B12', side='D')
    R('2962', 'B12', side='R')
    L.setp('ASC', {'359': (0, 0), '345': (450, 0), '347': (900, 0), '373': (1350, 0), '383': (900, -450)})
    R('371', 'ASC', side='D'); R('358', 'ASC', side='D'); R('382', 'ASC', side='D'); R('393', 'ASC', side='R')
    L.setp('LIP', {'866': (0, 0), '3823': (450, 0), '883': (0, 700), '884': (450, 700), '4036': (900, 700)})
    R('3835', 'LIP', side='U'); R('881', 'LIP', axis=(0, 1), mid=(0, 250), end_cof=True, side='L')
    R('895', 'LIP', side='U'); R('4045', 'LIP', side='U')
    L.setp('ACE', {'1473': (0, 0), '1475': (450, 0), '1487': (900, 0)})
    R('1486', 'ACE', side='U'); R('1499', 'ACE', side='U')
    L.setp('C1', {'1737': (0, 0), '1727': (0, -450), '1725': (-450, -450), '1775': (450, 0), '1881': (900, 0),
                  '2648': (450, 450)})
    R('1746', 'C1', side='L'); R('1736', 'C1', side='U'); R('1786', 'C1', side='U'); R('1892', 'C1', side='U')
    R('2659', 'C1', side='R')

    # ======================================================================== PERI periplasmic hydrolases
    PER = [('1543', '1530', '1532'), ('1582', '1570', '1572'), ('1621', '1609', '1610'), ('1660', '1647', '1649'),
           ('2473', '2460', '2462'), ('3248', '3235', '3238'), ('3692', '3679', '3681'), ('3718', '3705', '3707'),
           ('3744', '3731', '3734'), ('3770', '3757', '3759'), ('3796', '3783', '3785'), ('344', '331', '333')]
    for i, (rid, s, p) in enumerate(PER):
        x0, y0 = 600*(i % 4), 700*(i // 4)
        L.setp('PERI', {s: (x0, y0), p: (x0, y0+350)})
        R(rid, 'PERI', side='R', dist=110)


NEW_EXITS = []
EXTRA_TITLES = []
TITLE_POS = {}
CURVES = dict(t=0.6, bend=55, bend_stack=15,
              stacked=['1255', '1275', '2080', '2092', '2104', '182', '3310', '3322', '633', '3596',
                       '826', '839', '3074', '3847'],
              bend_override={'168': 0}, special_links={('1344', '1303'): (1, -1)})
PACK = dict(keep_out=[], auto_ring_keep_out=True)
RENDER_SCALE = 0.25
