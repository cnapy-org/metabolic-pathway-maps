"""SPEC - iML1515 carbohydrates (BRENDA draft, 193 reactions) -> hand-laid-out Escher map.

Every coordinate of a PRIMARY metabolite below was chosen by hand (zone-local, Escher convention:
x right, y DOWN). The toolkit only places markers/cofactors, curves, labels and packs detached pieces.

Main component (fixed): the glucose-1-phosphate hub with glycogen + maltodextrins (left column),
Leloir pathway (above), trehalose / disaccharides / fructose (below g6p).
All other pathways are self-contained pieces (each with its own title) that the packer arranges
around the main component."""

INPUT = '../../ecmapper/BRENDAmap/draft_map_iML1515/split/carbohydrates.escher.json'

ZONES = {
 'GLYC': dict(title='glycogen and maltodextrin metabolism', title_id='2355', color='#1f77b4'),
 'HEX':  dict(title='galactose (Leloir pathway)', title_id='2351', color='#2ca02c'),
 'DIS':  dict(title='metabolism of disaccharides', title_id='2358', color='#9467bd'),
 'RHA':  dict(title='dTDP-L-rhamnose biosynthesis', title_id='2350', anchor=('192', (-60, -230)), color='#8c564b'),
 'LPS':  dict(title=None, color='#7f7f7f'),
 'KDO':  dict(title='CMP-KDO biosynthesis', title_id='2346', anchor=('1', (-60, -230)), color='#e377c2'),
 'UGN':  dict(title='UDP-GlcNAc biosynthesis', title_id='2347', anchor=('90', (-700, -40)), color='#17becf'),
 'AMS':  dict(title='metabolism of amino sugars and derivatives', title_id='2357', anchor=('1293', (-1100, -260)), color='#17becf'),
 'MAN':  dict(title='D-mannose degradation', title_id='2348', anchor=('148', (-60, -260)), color='#bcbd22'),
 'PEN':  dict(title='degradation of pentoses', title_id='2352', anchor=('396', (-60, -300)), color='#ff7f0e'),
 'XYL':  dict(title='D-xylose degradation (Dahms pathway)', title_id='2349', anchor=('179', (-60, -330)), color='#ff7f0e'),
 'KGL':  dict(title='ketogluconate metabolism', title_id='2356', anchor=('1122', (-500, -230)), color='#d62728'),
 'SAC':  dict(title='degradation of sugar acids (hexuronates, D-galactonate)', title_id='2353', anchor=('715', (-60, -330)), color='#d62728'),
 'GLR':  dict(title=None, color='#d62728'),
 'ALC':  dict(title='degradation of sugar alcohols (glycerol, hexitols, galactitol)', title_id='2354', anchor=('855', (-1400, -150)), color='#1f77b4'),
 'GPD':  dict(title=None, color='#1f77b4'),
 'RHM':  dict(title=None, color='#2ca02c'),
 'SQ':   dict(title='sulfoquinovose degradation', title_id='2360', anchor=('1520', (-60, -560)), color='#8c564b'),
 'SGA':  dict(title='sulfated glycosaminoglycan metabolism', title_id='2359', anchor=('1484', (-60, -120)), color='#8c564b'),
 'PERI': dict(title=None, color='#555555'),
}
TITLE_POS = {'GLYC': (-2500, -3700), 'HEX': (-700, -2000), 'DIS': (-1900, 2000)}
EXTRA_TITLES = [
    dict(text='LPS core glucosylation (UDP-glucose)', anchor='2236', offset=(-200, -200), zone='LPS'),
    dict(text='O16 antigen', anchor='2058', offset=(-60, -150), zone='LPS'),
    dict(text='ribose', anchor='430', offset=(-60, -220), zone='PEN'),
    dict(text='D-glucarate / D-galactarate degradation', anchor='782', offset=(-60, -230), zone='GLR'),
    dict(text='glycerophosphodiester degradation', anchor='1845', offset=(-500, -900), zone='GPD'),
    dict(text='L-rhamnose / L-fucose degradation', anchor='385', offset=(-500, -220), zone='RHM'),
    dict(text='periplasm: sugar hydrolases / phosphatases', anchor='1395', offset=(-400, -220), zone='PERI'),
    dict(text='periplasm: glycerophosphodiesters, cardiolipin', anchor='1883', offset=(-300, -220), zone='PERI'),
    dict(text='periplasm', anchor='2136', offset=(-60, -110), zone='PERI'),
    dict(text='periplasm', anchor='2122', offset=(-60, -110), zone='PERI'),
    dict(text='periplasm', anchor='2096', offset=(-60, -110), zone='PERI'),
]
EXIT_OWNER = {}

# ------------------------------------------------------------------ splits (local copies of primaries)
def prepare(L):
    S = L.S
    # glucose: one copy for G1PP + HEX1, otherwise one local copy per maltodextrin reaction
    S['glc_H'] = L.split('1036', '1755'); L.merge(S['glc_H'], L.split('1021', '1033'))
    S['glc_1058'] = L.split('1036', '1058')                 # 1036 stays with AMALT1 1047
    for r in ('2012', '2023'): S['glc_'+r] = L.split('1021', r)   # 1021 stays with AMALT1 1083
    for r in ('1616', '1627', '1638', '1990', '2001', '2175'): S['glc_'+r] = L.split('1085', r)  # 1085 stays MLTG1 1095
    # maltose donor copies for the amylomaltase steps
    for r in ('1047', '1616', '1627', '1638'): S['malt_'+r] = L.split('1034', r)
    # UDP-glucose donor copies for the LPS-core / O16 glycosyltransferases (glucose block keeps 252)
    for r in ('2248', '2285', '2297', '2309', '2070'): S['udpg_'+r] = L.split('252', r)
    # dhap of PAI2T stays local in the disaccharide block
    S['dhap_P'] = L.split('1473', '2082')
    S['glc_B'] = L.split('1345', '1650')         # BGLA1 glucose, local (avoids crossing the trehalose chain)
    # sugar alcohols: glycerophosphodiesterases get their own glyc3p hub (separate piece)
    S['g3p_D'] = L.split('909', '920')
    for r in ('1857', '1869', '1953'): L.merge(S['g3p_D'], L.split('909', r))
    S['ttdca_L'] = L.split('1943', '1953')       # LPLIPAL2A140 product, local
    S['ttdcap_A'] = L.split('2212', '2223')      # APG3PAT140 co-substrate, local (avoids a crossing)
    # duplicate '→ Central' exit labels that would land on the same merged node (f6p 52, g6p 299)
    for k in ('2328', '2317'): del L.T[k]; L.exit_owner.pop(k, None)

ORIGIN = {'GLYC': (0, 0), 'HEX': (0, 0), 'DIS': (0, 0),          # main; detached pieces are packed
          'RHA': (2500, -3000), 'LPS': (2200, 400), 'KDO': (3000, -3000), 'UGN': (-5500, -3000),
          'AMS': (-5500, -3000), 'MAN': (4000, -2500), 'PEN': (2500, -1500), 'XYL': (5500, -3000),
          'KGL': (4500, 1000), 'SAC': (2000, 3200), 'GLR': (4500, 3800), 'ALC': (-4000, 1500),
          'GPD': (-5500, -1000), 'RHM': (6000, -500), 'SQ': (7000, 2500), 'SGA': (7000, 3500), 'PERI': (-3500, -3000)}
FLIP = {}

MERGES = [
    ('262', ['1010', '1098', '1109']),           # g1p hub
    ('299', ['1011', '1023', '1460']),           # g6p hub
    ('1345', ['1358']),                          # glucose of the disaccharide block
    ('1447', ['1698']),                          # fructose
    ('52', ['1305']), ('54', ['1282']), ('90', ['1317']),   # f6p, gam6p, UDP-GlcNAc (amino sugars)
    # pentoses (duplicate drawings of the same metabolites)
    ('396', ['496']), ('397', ['497']), ('406', ['506']), ('407', ['484']), ('455', ['418']),
    ('456', ['515', '884']), ('474', ['532', '896']), ('170', ['570', '882']), ('581', ['895']),
    ('180', ['443']), ('444', ['555']), ('169', ['2165']),          # Dahms pathway
    ('1124', ['1136']),                          # 2-dehydro-L-gulonate
    ('738', ['1494']), ('804', ['1496']), ('1508', ['602']), ('1509', ['603']),   # hexuronates + galactonate
    ('626', ['679']), ('613', ['680']),
    ('936', ['844', '1473']), ('843', ['1409']), ('858', ['1407']), ('845', ['1474']),  # sugar alcohols
    ('1396', ['1372', '1758']), ('1756', ['2151']),   # periplasm
]

def layout(L):
    S = L.S
    # =============================================== MAIN: g1p hub at (0, 0)
    X = -1300                                   # maltodextrin column G2..G7
    G = {'1034': -2925, '1035': -2175, '1060': -1425, '1606': -675, '1617': 75, '1628': 825}
    L.setp('GLYC', {k: (X, y) for k, y in G.items()})
    L.setp('GLYC', {'262': (0, 0), '299': (0, 750), S['glc_H']: (-550, 480),
                    '1096': (1150, 350), '1110': (1150, 1000), '1021': (X-400, -3300)})
    lv = {1: -2550, 2: -1800, 3: -1050, 4: -300, 5: 450}
    def mltg(rid, glc, y):                         # G(n+1) -> G(n) + glucose (copy upper left)
        L.setp('GLYC', {glc: (X-330, y-140)})
        L.reaction(rid, 'GLYC', axis=(0, -1), mid=(X, y), side='R', gap=40)
    def amalt(rid, malt, glc, y):                  # maltose + G(n) -> G(n+1) + glucose (copies left)
        L.setp('GLYC', {malt: (X-330, y-60), glc: (X-330, y+120)})
        L.reaction(rid, 'GLYC', axis=(0, 1), mid=(X, y), side='R', gap=40)
    mltg('1058', S['glc_1058'], lv[1]-110); mltg('1095', '1085', lv[1]+110)          # MLTG1 x2
    L.reaction('1083', 'GLYC', side='R', lab='L')                                      # AMALT1 (malt -> glc)
    mltg('1990', S['glc_1990'], lv[2]-110)                                             # MLTG2
    L.setp('GLYC', {S['malt_1047']: (X-330, lv[2]+50), '1036': (X-330, lv[2]+230)})
    L.reaction('1047', 'GLYC', axis=(0, 1), mid=(X, lv[2]+110), side='R', gap=40)       # AMALT1
    for k, (g, mp, am, mal, glc) in {3: ('2001', '1071', '1616', None, None), 4: ('2012', '2034', '1627', None, None),
                                      5: ('2023', '2045', '1638', None, None)}.items():
        y = lv[k]
        mltg(g, S['glc_'+g], y-200)                                                     # MLTG3-5
        L.reaction(mp, 'GLYC', axis=(0, -1), mid=(X, y), side='R', gap=40)              # MLTP1-3
        amalt(am, S['malt_'+am], S['glc_'+am], y+200)                                   # AMALT2-4
    L.reaction('1019', 'GLYC', shift=(-70, 0), side='L', lab='L')                      # PGMT lanes
    L.reaction('307', 'GLYC', shift=(70, 0), side='R', lab='R')
    L.reaction('1755', 'GLYC', side='L')                                               # G1PP
    L.reaction('1033', 'GLYC', side='D')                                               # HEX1
    L.reaction('1107', 'GLYC', side='U')                                               # GLCP
    L.reaction('1121', 'GLYC', side='D')                                               # GLGC
    L.reaction('1832', 'GLYC', side='R')                                               # GLCS1
    # ----------------------------------------------- Leloir (above g1p)
    L.setp('HEX', {'239': (0, -1300), '241': (0, -650), '252': (600, -250), '253': (600, -950)})
    L.reaction('251', 'HEX', side='L')                                                 # GALKr
    L.reaction('272', 'HEX', axis=(0, 1), mid=(0, -350), lab='L')                      # UGLT (crossed pair)
    L.reaction('298', 'HEX', side='D')                                                 # GALUi
    L.reaction('261', 'HEX', lab='R')                                                  # UDPG4E
    # ----------------------------------------------- trehalose / disaccharides / fructose (below g6p)
    L.setp('DIS', {'1421': (0, 1400), '1384': (0, 2000), '1345': (0, 2600),
                   '1731': (750, 1050), '1733': (520, 620), '1639': (-850, 950), S['glc_B']: (-650, 1400),
                   '1343': (-420, 3200), '1357': (420, 3200), '1344': (0, 3800),
                   '1447': (-950, 2600), '1449': (-950, 3250), '274': (-1600, 2600), '276': (-1600, 3250),
                   '1795': (750, 2200), '1796': (800, 2750), '1954': (1350, 2200), '1955': (1250, 2750),
                   '2072': (750, 1650), S['dhap_P']: (1250, 1800)})
    L.reaction('1472', 'DIS', side='R')                                                # TRE6PS
    L.reaction('1433', 'DIS', axis=(1, 0), mid=(330, 1400), side='D', end_cof=True, step=80)  # TRE6PH
    L.reaction('1445', 'DIS', side='L')                                                # TRE6PP
    L.reaction('1393', 'DIS', side='L')                                                # TREH
    L.reaction('1743', 'DIS', side='R')                                                # FRULYSDG
    L.reaction('1650', 'DIS', side='L')                                                # BGLA1
    L.reaction('1355', 'DIS', side='L')                                                # GALS3
    L.reaction('1368', 'DIS', side='R')                                                # LACZ
    L.reaction('2183', 'DIS', lab='U')                                                 # XYLI2
    L.reaction('1459', 'DIS', shift=(-70, 0), side='L')                                # HEX7
    L.reaction('1730', 'DIS', shift=(70, 0), side='R')                                 # F6PP
    L.reaction('1709', 'DIS', side='U')                                                # F1PP
    L.reaction('286', 'DIS', side='L')                                                 # FRUK
    L.reaction('1807', 'DIS', side='R')                                                # GLCATr
    L.reaction('1966', 'DIS', side='R')                                                # MALTATr
    L.reaction('2082', 'DIS', side='R')                                                # PAI2T
    # =============================================== detached pieces
    L.setp('RHA', {'192': (0, 0), '193': (400, 0), '205': (800, 0), '216': (1200, 0), '226': (1600, 0)})
    L.reaction('204', 'RHA', side='U'); L.reaction('215', 'RHA', side='D')
    L.reaction('224', 'RHA', lab='U'); L.reaction('237', 'RHA', side='U')
    # LPS core chain with local UDP-glucose donors
    L.setp('LPS', {'2236': (0, 0), '2237': (0, 450), '2274': (0, 900), '2286': (0, 1350), '2298': (0, 1800),
                   '2058': (0, 2350), '2059': (550, 2350)})
    for rid, y in (('2248', 225), ('2285', 675), ('2297', 1125), ('2309', 1575)):
        L.setp('LPS', {S['udpg_'+rid]: (-330, y-120)})
        L.reaction(rid, 'LPS', axis=(0, 1), mid=(0, y), side='R', dist=110)
    L.setp('LPS', {S['udpg_2070']: (0, 2700)})
    L.reaction('2070', 'LPS', axis=(1, 0), mid=(260, 2350), side='U', dist=110)        # O16GLCT2
    # CMP-KDO
    L.setp('KDO', {'1': (0, 0), '2': (0, 350), '27': (0, 750), '12': (0, 1150), '13': (0, 1550)})
    for r in ('10', '39', '51', '24'): L.reaction(r, 'KDO', side='R')
    # amino sugars + UDP-GlcNAc
    L.setp('UGN', {'52': (0, 0), '54': (0, 650), '66': (0, 1250), '76': (0, 1850), '90': (0, 2450),
                   '1318': (0, 3050), '1330': (0, 3650), '1581': (-450, 3950), '1583': (0, 4300)})
    L.reaction('65', 'UGN', shift=(-70, 0), side='L', dist=120)                         # GF6PTA
    L.reaction('1316', 'AMS', shift=(70, 0), side='R', dist=120)                        # G6PDA
    L.reaction('74', 'UGN', lab='L')                                                    # PGAMT
    L.reaction('87', 'UGN', side='L'); L.reaction('100', 'UGN', side='L')
    L.reaction('1326', 'UGN', lab='L'); L.reaction('1341', 'UGN', side='L')
    L.reaction('1593', 'UGN', axis=(0, 1), mid=(0, 3975), side='R')                    # ACMAMUT
    L.setp('AMS', {'1271': (650, 650), '1259': (650, 0), '1257': (1250, 0), '1293': (1850, 0),
                   '1568': (1300, 650), '1570': (1050, 1050)})
    L.reaction('1292', 'AMS', side='D'); L.reaction('1279', 'AMS', side='R')
    L.reaction('1270', 'AMS', side='U'); L.reaction('1303', 'AMS', side='U')
    L.reaction('1580', 'AMS', side='U')                                                 # ACM6PH
    # D-mannose (square)
    L.setp('MAN', {'148': (0, 0), '101': (700, 0), '115': (700, 700), '126': (0, 700), '160': (0, 1300),
                   '102': (1350, 0)})
    L.reaction('159', 'MAN', shift=(-70, 0), side='L', dist=110)                        # HEX4
    L.reaction('2057', 'MAN', shift=(70, 0), side='R', dist=110)                        # MN6PP
    L.reaction('1794', 'MAN', side='U')                                                 # GDPMNH
    L.reaction('146', 'MAN', shift=(70, 0), side='R', dist=110)                         # MAN1PT2
    L.reaction('125', 'MAN', shift=(-70, 0), side='L', dist=110)                        # GDPMNP
    L.reaction('134', 'MAN', lab='D'); L.reaction('168', 'MAN', lab='L')
    L.reaction('112', 'MAN', side='U')                                                  # GMAND
    # pentoses grid
    L.setp('PEN', {'396': (0, 0), '397': (450, 0), '465': (900, 0), '169': (1400, 0),
                   '456': (900, 500), '170': (1400, 500),
                   '406': (0, 1000), '407': (450, 1000), '455': (900, 1000), '474': (1400, 1000),
                   '592': (900, 1500), '581': (1400, 1500),
                   '2164': (1400, -650), S['glc_2175']: (1800, -300),
                   '430': (0, 2100), '432': (550, 2100)})
    L.reaction('405', 'PEN', shift=(0, -70), side='U', dist=100); L.reaction('505', 'PEN', shift=(0, 70), side='D', dist=100)
    L.reaction('2120', 'PEN', side='U')                                                 # RU5PP
    L.reaction('473', 'PEN', shift=(-70, 0), side='L', dist=100); L.reaction('531', 'PEN', shift=(70, 0), side='R', dist=100)
    L.reaction('178', 'PEN', side='R')                                                  # XYLI1
    L.reaction('580', 'PEN', shift=(0, -70), side='U', dist=100); L.reaction('894', 'PEN', shift=(0, 70), side='D', dist=100)
    L.reaction('415', 'PEN', shift=(0, -70), side='U', dist=100); L.reaction('514', 'PEN', shift=(0, 70), side='D', dist=100)
    L.reaction('428', 'PEN', shift=(0, -70), side='U', dist=100); L.reaction('495', 'PEN', shift=(0, 70), side='D', dist=100)
    L.reaction('464', 'PEN', shift=(-70, 0), side='L', dist=100); L.reaction('523', 'PEN', shift=(70, 0), side='R', dist=100)
    L.reaction('482', 'PEN', shift=(0, -70), side='U', dist=100); L.reaction('540', 'PEN', shift=(0, 70), side='D', dist=100)
    L.reaction('591', 'PEN', shift=(-70, 0), side='L', dist=100); L.reaction('906', 'PEN', shift=(70, 0), side='R', dist=100)
    L.reaction('600', 'PEN', side='D')                                                  # LYXI
    L.reaction('2175', 'PEN', side='L')                                                 # XYHDL
    L.reaction('442', 'PEN', shift=(0, -70), side='U', dist=100); L.reaction('2094', 'PEN', shift=(0, 70), side='D', dist=100)
    # Dahms pathway
    L.setp('XYL', {'179': (0, 0), '180': (450, 0), '445': (900, -350), '444': (900, 350),
                   '543': (1400, 350), '541': (1400, 850)})
    L.reaction('190', 'XYL', side='D'); L.reaction('454', 'XYL', lab='U')
    L.reaction('569', 'XYL', side='U'); L.reaction('554', 'XYL', side='R')
    # ketogluconate diamond
    L.setp('KGL', {'1122': (0, 0), '1124': (-900, 850), '1172': (900, 850), '1160': (0, 850),
                   '1148': (0, 1700), '1245': (0, 2300)})
    import math
    def unit(a, b): d = math.hypot(b[0]-a[0], b[1]-a[1]); return ((b[0]-a[0])/d, (b[1]-a[1])/d)
    u = unit((0, 0), (-900, 850)); n = (u[1], -u[0])                     # n points up-left (outside)
    L.reaction('1135', 'KGL', shift=(n[0]*70, n[1]*70), side=n, dist=100)
    L.reaction('1146', 'KGL', shift=(-n[0]*70, -n[1]*70), side=(-n[0], -n[1]), dist=100)
    u = unit((0, 0), (900, 850)); n = (-u[1], u[0])                      # outside = up-right
    for rid, t in (('1183', 0.35), ('1696', 0.65)):
        L.reaction(rid, 'KGL', axis=u, mid=(900*t, 850*t), side=(n[0] if n[0] > 0 else -n[0], -abs(n[1])), gap=40)
    for rid, x in (('1170', -620), ('1567', -300)):                      # 2DGULRy/x  2dhguln -> idon
        L.reaction(rid, 'KGL', axis=(1, 0), mid=(x, 850), side='U', gap=40, dist=110)
    for rid, x in (('1195', 600), ('1207', 400), ('1219', 200)):         # IDOND, IDOND2, IDOND
        L.reaction(rid, 'KGL', axis=(-1, 0), mid=(x, 850), side='D', gap=40, dist=90)
    u = unit((-900, 850), (0, 1700))
    for rid, t in (('1159', 0.35), ('1555', 0.65)):                      # 2DGULRGy/x
        L.reaction(rid, 'KGL', axis=u, mid=(-900+900*t, 850+850*t), side=(-1, 1), gap=40)
    u = unit((900, 850), (0, 1700)); n = (u[1], -u[0])
    L.reaction('1231', 'KGL', shift=(n[0]*70, n[1]*70), side=(1, 1), dist=100)         # 5DGLCNR
    L.reaction('1243', 'KGL', shift=(-n[0]*70, -n[1]*70), side=(-1, -1), dist=100)     # 5DKGR
    L.reaction('1256', 'KGL', shift=(-70, 0), side='L', dist=110)                      # GNK
    L.reaction('1844', 'KGL', shift=(70, 0), side='R', dist=110)                       # GNP
    # sugar acids: hexuronates -> 2ddg6p -> g3p + pyr <- galactonate
    L.setp('SAC', {'715': (0, 0), '716': (450, 0), '725': (900, 0), '749': (0, 800), '750': (450, 800),
                   '759': (900, 800), '738': (1350, 400), '804': (1800, 400), '1508': (2250, 100),
                   '1509': (2250, 700), '601': (2750, 400), '613': (3250, 400), '626': (3800, 400)})
    L.reaction('724', 'SAC', lab='U'); L.reaction('737', 'SAC', side='U'); L.reaction('748', 'SAC', side='R')
    L.reaction('758', 'SAC', lab='D'); L.reaction('771', 'SAC', side='D'); L.reaction('781', 'SAC', side='L')
    L.reaction('815', 'SAC', shift=(0, -70), side='U', dist=100); L.reaction('1507', 'SAC', shift=(0, 70), side='D', dist=100)
    L.reaction('1518', 'SAC', lab='U')                                                 # EDA
    L.reaction('612', 'SAC', axis=(-1, 0), mid=(2500, 330), lab='U')                   # DDPGALA x2
    L.reaction('678', 'SAC', axis=(-1, 0), mid=(2500, 470), lab='D')
    L.reaction('625', 'SAC', shift=(0, -70), side='U', dist=100); L.reaction('669', 'SAC', shift=(0, 70), side='D', dist=100)
    L.reaction('636', 'SAC', shift=(0, -160), side='U', dist=100); L.reaction('659', 'SAC', side='U', dist=80, lab='D')
    L.reaction('690', 'SAC', shift=(0, 160), side='D', dist=100)
    # glucarate / galactarate
    L.setp('GLR', {'782': (0, 0), '793': (0, 700), '704': (500, 350), '691': (1050, 350), '1919': (1050, -250),
                   '638': (1600, 350), '1968': (1600, -300), '1969': (2050, -50), '639': (2150, 350)})
    L.reaction('792', 'GLR', side='U'); L.reaction('802', 'GLR', side='D')
    L.reaction('714', 'GLR', side='D'); L.reaction('1927', 'GLR', lab='L')
    L.reaction('703', 'GLR', side='D'); L.reaction('1979', 'GLR', side='L'); L.reaction('650', 'GLR', side='D')
    # sugar alcohols
    L.setp('ALC', {'936': (0, 0), '909': (-900, 0), '986': (0, 600), '973': (-900, 600),
                   '845': (900, 0), '843': (450, -500), '858': (450, -1050), '855': (450, -1600),
                   '818': (450, 1150), '816': (0, 1650), '831': (900, 1650),
                   '1895': (-1500, 450), '1771': (-900, -650), '2201': (-900, -1300),
                   '1770': (-1600, -975), '2212': (-1600, -350), '1943': (-2300, -975),
                   S['ttdcap_A']: (-420, -340), '2261': (-1600, -1650), '2263': (-2400, -1310)})
    for rid, x in (('947', -675), ('959', -450), ('971', -225)):                        # G3PD5/6/7
        L.reaction(rid, 'ALC', axis=(1, 0), mid=(x, 0), side='D', gap=40, dist=110)
    L.reaction('984', 'ALC', side='R'); L.reaction('997', 'ALC', side='D')            # GLYK, GLYCDx
    L.reaction('1009', 'ALC', shift=(70, 0), side='R', dist=110)                      # DHAPT
    L.reaction('1605', 'ALC', shift=(-70, 0), side='L', dist=110)                     # ALKP
    L.reaction('854', 'ALC', shift=(-70, 0), lab='L'); L.reaction('1483', 'ALC', shift=(70, 0), lab='R')  # TGBPA x2
    L.reaction('880', 'ALC', shift=(-70, 0), side='L', dist=110); L.reaction('1419', 'ALC', shift=(70, 0), side='R', dist=110)
    L.reaction('868', 'ALC', side='R')                                                 # GLTPD
    L.reaction('1718', 'ALC', lab='L')                                                 # F6PA
    L.reaction('829', 'ALC', side='L'); L.reaction('842', 'ALC', side='R')            # M1PD, SBTPD
    L.reaction('1906', 'ALC', side='L')                                                # GPDDA4
    L.reaction('1782', 'ALC', shift=(-70, 0), side='L', lab='L')                       # G3PAT140
    L.reaction('2223', 'ALC', shift=(70, 0), side='R', lab='R')                        # APG3PAT140
    L.reaction('2211', 'ALC', axis=(0, -1), mid=(-900, -975), side='R')                # AGPAT140
    L.reaction('2260', 'ALC', side='L')                                                # ACPPAT140
    L.reaction('2199', 'ALC', shift=(0, -70), side='U', dist=110)                      # AACPS1
    L.reaction('2235', 'ALC', shift=(0, 70), side='D', dist=110)                       # FA140ACPHi
    L.reaction('2273', 'ALC', axis=(-1, 0), mid=(-1950, -1310), side='U')              # EDTXS4
    # glycerophosphodiesterases (cytoplasm) around their own glyc3p copy
    L.setp('GPD', {S['g3p_D']: (0, 0), '1858': (-250, -750), '1859': (-650, -250), '1845': (-250, 750),
                   '1847': (-650, 250), '1941': (750, 250), S['ttdca_L']: (350, 650)})
    L.reaction('920', 'GPD', axis=(-0.7, 0.7), mid=(380, -380), side='U', end_cof=True, step=80)   # GPDDA1
    L.reaction('1869', 'GPD', side='R'); L.reaction('1857', 'GPD', side='R')          # GPDDA2, GPDDA3
    L.reaction('1953', 'GPD', side='R')                                                # LPLIPAL2A140
    # L-rhamnose / L-fucose
    L.setp('RHM', {'310': (0, 0), '309': (0, 650), '333': (-650, -900), '321': (-650, -300), '308': (-650, 325),
                   '342': (650, -900), '343': (650, -300), '352': (650, 325),
                   '385': (0, -1300), '374': (0, -700), '375': (-280, -280)})
    L.reaction('341', 'RHM', lab='L'); L.reaction('332', 'RHM', side='L'); L.reaction('319', 'RHM', lab='U')
    L.reaction('351', 'RHM', lab='R'); L.reaction('373', 'RHM', side='R'); L.reaction('361', 'RHM', lab='U')
    L.reaction('395', 'RHM', side='R'); L.reaction('384', 'RHM', lab='R')
    # sulfoquinovose
    L.setp('SQ', {'1520': (0, 0), '1533': (500, 0), '1521': (0, 600)})
    L.reaction('1543', 'SQ', axis=(0, 1), mid=(250, -380), side='R', end_cof=True)
    L.reaction('1532', 'SQ', side='L')
    L.setp('SGA', {'1484': (0, 0), '1485': (0, 450)})
    L.reaction('1493', 'SGA', lab='R')
    # periplasm
    L.setp('PERI', {'1396': (0, 0), '1395': (0, -600), '1810': (600, 0), '1370': (-650, -450), '1371': (-650, 150),
                    '1756': (0, 600), '2150': (0, 1200),
                    '2136': (1400, 600), '2137': (1400, 1100), '2122': (1900, 600), '2123': (1900, 1100),
                    '2096': (2400, 600), '2098': (2400, 1100)})
    L.reaction('1405', 'PERI', side='R'); L.reaction('1822', 'PERI', side='U'); L.reaction('1382', 'PERI', side='L')
    L.reaction('1769', 'PERI', side='R'); L.reaction('2162', 'PERI', side='R')
    L.reaction('2148', 'PERI', side='R'); L.reaction('2134', 'PERI', side='R'); L.reaction('2108', 'PERI', side='R')
    L.setp('PERI', {'923': (0, 3000), '1883': (-250, 2250), '1884': (-650, 2750), '1870': (-250, 3750),
                    '1872': (-650, 3250), '1928': (750, 3250), '1930': (350, 3650),
                    '1907': (-750, 3450), '1653': (-100, 3950)})
    L.reaction('934', 'PERI', axis=(-0.7, 0.7), mid=(380, 2620), side='U', end_cof=True, step=80)   # GPDDA1pp
    L.reaction('1894', 'PERI', side='R'); L.reaction('1882', 'PERI', side='R'); L.reaction('1940', 'PERI', side='R')
    L.reaction('1918', 'PERI', side='L')                                               # GPDDA4pp
    L.setp('PERI', {'1651': (-750, 4450), '1652': (-700, 3950), '1663': (-100, 4650), '1664': (300, 4850),
                    '1674': (550, 4450), '1675': (650, 3950)})
    L.reaction('1662', 'PERI', side='L'); L.reaction('1673', 'PERI', side='R'); L.reaction('1684', 'PERI', side='R')

NEW_EXITS = [('1344', '→ Leloir pathway'), ('192', '← glycogen metabolism (g1p)'), ('1', '← Central'),
             ('g3p_D', '→ sugar alcohols (glyc3p)'), ('1969', '→ D-mannose degradation'),
             ('1533', '→ Central'), ('309', '→ Central'), ('456', '→ Central'), ('dhap_P', '→ Central')]

CURVES = dict(t=0.6, stacked=['1058', '1095', '1990', '1047', '2001', '1071', '1616', '2012', '2034', '1627',
                              '2023', '2045', '1638', '1183', '1696', '1170', '1567', '1195', '1207', '1219',
                              '1159', '1555', '947', '959', '971'],
              bend_override={}, special_links={})
# two tiny keep-out markers widen the packing target (main component is tall) -> landscape canvas
PACK = dict(keep_out=[('box', (-6500, -3400, -6490, -3390)), ('box', (7490, 4300, 7500, 4310))], auto_ring_keep_out=True)
RENDER_SCALE = 0.25
