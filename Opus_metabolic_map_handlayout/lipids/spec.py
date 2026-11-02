"""SPEC - iML1515 lipid metabolism (BRENDA draft split map 'lipids') -> hand-laid-out Escher map.

Every coordinate of a PRIMARY metabolite below was chosen by hand (zone-local coordinates, zone origin
in ORIGIN; Escher convention x right, y DOWN).  Small helper functions only repeat a hand-chosen motif
(rows/slots of an elongation round, one phospholipid-hydrolysis unit, one acyl-chain column) so that
repeated cycles are drawn identically - the same idea as ring(k) in the central-metabolism example.

Zones
  FAS   fatty acid biosynthesis (type II): initiation square + elongation rounds C4 ... C18 drawn as a
        snake (every round = one row  3-oxoacyl -> 3-hydroxyacyl -> trans-2-enoyl -> acyl-ACP; the
        condensation 3OAS/KAS is the vertical step to the next row)
  BOX   acyl-CoA metabolism / beta-oxidation: the same snake grammar for CoA esters (one row per chain
        length: acyl -> enoyl -> 3-hydroxy -> 3-oxo; thiolase = vertical step), activation pair
        (FA_p -> acyl-CoA uptake/ligase, acyl-CoA -> FA_c thioesterase) on the outer side of each acyl-CoA
  PLS   phosphatidic acid synthesis: one identical column per acyl chain (16:0, 18:1, 16:1):
        acyl-ACP / acyl-phosphate / free FA on the left, glycerol-3-P -> 1-acyl-G3P (-> PA) to the right
  GLYC  glycerol / glycerol-3-phosphate side reactions
  PLD   phospholipid degradation + re-acylation: one identical unit per (head group, chain)
  LIPA  lipid A (Raetz) pathway + Kdo, PagP, EptA reactions
  ISO   isoprenoid biosynthesis (MEP pathway)
  ECA   ECA / O-antigen precursor reactions (WecA, dTDP-sugar acetyltransferase)
"""
import math

INPUT = '/home/user/metabolic-pathway-maps/ecmapper/BRENDAmap/draft_map_iML1515/split/lipids.escher.json'

ZONES = {
 'FAS':  dict(title='fatty acid biosynthesis (FAS II)', title_id='2089', color='#1f77b4'),
 'BOX':  dict(title='acyl-CoA metabolism / β-oxidation', title_id='2087', color='#d62728'),
 'PLS':  dict(title='phosphatidic acid synthesis', title_id='2084', color='#2ca02c'),
 'GLYC': dict(title='glycerol / glycerol-3-phosphate', color='#8c564b'),
 'PLD':  dict(title='phospholipid degradation and re-acylation', title_id='2088', color='#9467bd'),
 'LIPA': dict(title='lipid A biosynthesis', title_id='2086', color='#e377c2'),
 'ISO':  dict(title='isoprenoid biosynthesis (MEP pathway)', title_id='2085', color='#17becf'),
 'ECA':  dict(title='ECA / O-antigen precursors', title_id='2090', color='#7f7f7f'),
}

EXIT_OWNER = {}

# ------------------------------------------------------------------ helpers for splits
def copy(L, node, rids):
    """one new primary copy of `node` shared by all reactions `rids`"""
    new = L.split(node, rids[0])
    for r in rids[1:]:
        L.merge(new, L.split(node, r))
    return new

def side(L, node, rids=None):
    """draw `node` as a side product / co-substrate (small dot, like a cofactor) in reactions `rids`
    (default: all its reactions); one copy per reaction"""
    users = sorted(L.uses()[node], key=int)
    rids = users if rids is None else rids
    out = []
    for i, r in enumerate(rids):
        n = node if (i == 0 and set(rids) == set(users)) else L.split(node, r)
        L.N[n]['node_is_primary'] = False; out.append(n)
    return out

# reaction groups of the phospholipid degradation units: (head group, chain) -> ids
PERI = {  # PLIPA1 (-> 2-acyl), PLIPA2 (-> 1-acyl), LPLIPAL1 (1-acyl -> backbone)
 ('A', '160'): ('1477', '1598', '1270'), ('A', '161'): ('1490', '1609', '1283'), ('A', '181'): ('1503', '1620', '1307'),
 ('E', '160'): ('1925', '1937', '2060'), ('E', '161'): ('1515', '1526', '1853'), ('E', '181'): ('1538', '1549', '1864'),
 ('G', '160'): ('1562', '1631', '1333'), ('G', '161'): ('1575', '1642', '1346'), ('G', '181'): ('1587', '1652', '1358')}
CYTO = {  # 2-acyl re-acylation (2AGPxAT, None for PA), LPLIPAL2 (2-acyl -> backbone)
 ('A', '160'): (None, '1370'), ('A', '161'): (None, '1295'), ('A', '181'): (None, '1319'),
 ('E', '160'): ('1966', '1877'), ('E', '161'): ('1669', '1888'), ('E', '181'): ('1686', '1899'),
 ('G', '160'): ('1702', '1384'), ('G', '161'): ('1718', '1397'), ('G', '181'): ('1731', '1409')}

def prepare(L):
    S = L.S
    # ---- FAS: duplicates of chain intermediates -> one node; malonyl-ACP / donors as side substrates
    L.merge('336', '743')                    # acetoacetyl-ACP (KAS product = 3OAR40 substrate)
    L.merge('850', '364'); L.merge('863', '366')   # 3HAD100 duplicate -> lane pair
    sec3o = [n for n in L.top['28']['S'] if L.N[n]['bigg_id'] == '3ooctdACP_c'][0]
    L.N[sec3o]['node_is_primary'] = True; L.merge('1018', sec3o)   # 3OAS180 -> 3OAR180 continuous
    side(L, '1017')                          # malACP of 3OAS180 (as in all other 3OAS rounds)
    S['mal_181'] = side(L, '335', ['1980'])[0]          # malACP of 3OAS181 (PLS zone)
    # ---- acetyl-CoA: 389 stays the FAS initiation hub; thiolases / transferases get a side copy
    for r in ('1059', '1072', '1084', '1215', '1950'):
        side(L, '389', [r])
    side(L, '1205')                          # acetate of HXCT (CoA transferase side product)
    # ---- free fatty acids: own copies per zone
    S['hdcea_box'] = copy(L, '1143', ['1153'])          # FACOAE161 product (BOX)
    S['ocdcea_box'] = copy(L, '1156', ['1166'])         # FACOAE181 product (BOX)
    S['hdcea_lipa'] = copy(L, '1387', ['2049'])         # FA161ACPHi product (LIPA)
    # ---- PLS: one glycerol-3-P per acyl chain; glycerol branch gets its own copy
    S['g3p_181'] = copy(L, '1167', ['1841', '1781'])
    S['g3p_glyc'] = copy(L, '417', ['430', '1191'])
    S['pal_agpat'] = side(L, '2', ['1743'])[0]          # palmitoyl-ACP donor of AGPAT160
    # ---- PLD: fatty acids released / consumed drawn as side products; one backbone per unit
    for rid in ('1370', '1295', '1319'):                # LPLIPAL2A* (cytosol)
        for n in L.top[rid]['P']:
            b = L.N[n]['bigg_id']
            if b in ('hdca_c', 'hdcea_c', 'ocdcea_c'): side(L, n, [rid])
    S['g3p_A160'] = copy(L, '1167', ['1370']); S['g3p_A181'] = copy(L, '1167', ['1319'])
    S['g3p_A161'] = copy(L, '417', ['1295'])
    for n in ('1260', '1323', '1336', '1348', '1374', '1387', '1399'):
        side(L, n)
    for n in ('555', '55'):                  # periplasmic FA also used by FACOAL*t2pp (stays primary)
        side(L, n, [r for r in sorted(L.uses()[n], key=int) if r not in ('570', '70')])
    L.merge('1259', '1273')                  # glyc3p_p
    for (hg, ch), (p1, p2, l1) in PERI.items():
        bb = [n for n in L.top[l1]['P'] if L.is_primary(n)][0]
        if len(L.uses()[bb]) > 1: S[f'bb_p_{hg}{ch}'] = L.split(bb, l1)
        else: S[f'bb_p_{hg}{ch}'] = bb
    for (hg, ch), (at, l2) in CYTO.items():
        if hg == 'A': continue
        bb = [n for n in L.top[l2]['P'] if L.is_primary(n)][0]
        if len(L.uses()[bb]) > 1: S[f'bb_c_{hg}{ch}'] = L.split(bb, l2)
        else: S[f'bb_c_{hg}{ch}'] = bb
    S['pe161_eptA'] = side(L, '1440', ['1452'])[0]; S['pe181_eptA'] = side(L, '1453', ['1464'])[0]
    side(L, '1441'); side(L, '1454')         # diacylglycerol side products of EptA
    side(L, '1230'); side(L, '1243')         # palmitate donors of PagP (LIPAHT*)
    S['ckdo_2'] = copy(L, '271', ['1912'])   # Kdo for MOAT3C (separate piece)
    L.merge('151', '705'); L.merge('136', '706')   # IPP / DMAPP of IPDDI join the MEP pathway

ORIGIN = {'FAS': (900, 1500), 'BOX': (3600, 500), 'PLD': (6400, 500), 'PLS': (450, 7300),
          'GLYC': (9900, 6700), 'LIPA': (6750, 4900), 'ISO': (8700, 5300), 'ECA': (9900, 5000)}
FLIP = {}
MERGES = []

# ------------------------------------------------------------------ FAS / BOX row grammar
HF = 600; SF = (0, 350, 700, 1150)            # FAS: row spacing, slot x (oxo, OH, enoyl, acyl)
def fas(k, s):                                   # row k (C4=0 ... C18=7); odd rows run right->left
    return (SF[s] if k % 2 == 0 else SF[3]-SF[s], k*HF)
HB = 750; SB = (0, 550, 1100, 1650)            # BOX: row spacing, slot x (acyl, enoyl, OH, oxo)
def box(k, s):
    return (SB[s] if k % 2 == 0 else SB[3]-SB[s], k*HB)

def layout(L):
    S = L.S
    # ================================================================ FAS
    L.setp('FAS', {
        '389': (-700, -1300), '377': (0, -1300), '335': (0, -650), '406': (-700, -650),   # accoa malcoa malACP acACP
        '1424': (600, -1100), '1427': (600, -400),                                       # malcoame -> ogmeACP
        '336': fas(0, 0),                                                                # C4: acetoacetyl-ACP
        '728': fas(1, 0), '758': fas(1, 1), '771': fas(1, 2), '782': fas(1, 3),          # C6
        '795': fas(2, 0), '825': fas(2, 1), '827': fas(2, 2),                            # C8
        '810': fas(3, 0), '850': fas(3, 1), '863': fas(3, 2),                            # C10
        '874': fas(4, 0), '889': fas(4, 1), '891': fas(4, 2),                            # C12
        '914': fas(5, 2), '926': fas(5, 3),                                              # C14
        '939': fas(6, 0), '980': fas(6, 1), '968': fas(6, 2), '954': fas(6, 3),          # C16
        '956': (1850, 6*HF),                                                             # palmitate
        '1018': fas(7, 0), '16': fas(7, 1), '30': fas(7, 2),                             # C18
        '444': (-700, 4*HF)})                                                            # 3-oxo-cis-5-dodecenoyl-ACP
    # initiation square
    L.reaction('405', 'FAS', side='U')                     # ACCOAC  accoa -> malcoa
    L.reaction('388', 'FAS', side='R')                     # MCOATA  malcoa -> malACP
    L.reaction('442', 'FAS', side='L')                     # ACOATA  accoa -> acACP
    L.reaction('416', 'FAS', side='U', dist=120)           # MACPD   malACP -> acACP
    L.reaction('1228', 'FAS', ends=('406', '336'), side='L')   # KAS14 acACP + malACP -> actACP
    L.reaction('349', 'FAS', side='R')                     # KAS15   malACP -> actACP
    L.reaction('1438', 'FAS', ends=('1424', '1427'), side='R')  # OGMEACPS (biotin precursor)
    # elongation rounds: rows = 3OAR, 3HAD, EAR (cofactors below, EAR pair as lanes), 3OAS = vertical
    def ear(rx, ry, k):                                    # two enoyl reductases (NADH / NADPH) as lanes
        L.reaction(rx, 'FAS', shift=(0, -75), side='U', dist=120)
        L.reaction(ry, 'FAS', shift=(0, 75), side='D', dist=120)
    L.reaction('756', 'FAS', axis=(1, 0), mid=(175, 0), side='D', end_cof=True)          # 3OAR40 (-> 3haACP)
    L.reaction('742', 'FAS', axis=(0, 1), mid=(SF[3], HF/2), side='R', end_cof=True)     # 3OAS60 (butACP)
    L.reaction('769', 'FAS', side='D'); L.reaction('780', 'FAS', side='D'); ear('793', '1140', 1)
    L.reaction('808', 'FAS', side='L')                                                   # 3OAS80
    L.reaction('848', 'FAS', side='D'); L.reaction('836', 'FAS', side='D')
    L.reaction('824', 'FAS', axis=(0, 1), mid=(SF[3], 2.5*HF), side='R', end_cof=True)   # 3OAS100 (ocACP)
    L.reaction('861', 'FAS', side='D')
    L.reaction('375', 'FAS', shift=(0, -75), side='U', dist=120); L.reaction('872', 'FAS', shift=(0, 75), side='D', dist=120)
    L.reaction('888', 'FAS', axis=(0, 1), mid=(0, 3.5*HF), side='L', end_cof=True)       # 3OAS120 (dcaACP)
    L.reaction('912', 'FAS', side='D'); L.reaction('900', 'FAS', side='D')
    L.reaction('924', 'FAS', axis=(-1, 0), mid=((SF[3]-SF[1]+SF[3]-SF[2])/2, 5*HF), side='D', end_cof=True)  # 3HAD140
    ear('937', '1116', 5)
    L.reaction('952', 'FAS', side='L')                                                   # 3OAS160
    L.reaction('1002', 'FAS', side='D'); L.reaction('990', 'FAS', side='D'); ear('979', '1128', 6)
    L.reaction('966', 'FAS', shift=(0, -75), side='U', dist=110)                         # FA160ACPHi
    L.reaction('1047', 'FAS', shift=(0, 75), side='D', dist=110)                         # AACPS3
    L.reaction('1031', 'FAS', side='R')                                                  # 3OAS180
    L.reaction('28', 'FAS', side='D'); L.reaction('52', 'FAS', side='D')
    L.reaction('42', 'FAS', axis=(-1, 0), mid=(SF[3]-SF[2]-225, 7*HF), side='D', end_cof=True)   # EAR180x
    L.reaction('458', 'FAS', axis=(0, 1), mid=(-700, 3.5*HF), side='L', end_cof=True)    # 3OAS121
    L.reaction('531', 'FAS', axis=(-1, 0), mid=(-875, 4*HF), side='D', end_cof=True)     # 3OAR121

    # ================================================================ BOX
    L.setp('BOX', {
        '57': (0, 0), '55': (-450, -250), S['ocdcea_box']: (-450, 250),                  # odecoa + FA pair
        '557': (SB[2], 0), '555': (SB[2]-450, -250), S['hdcea_box']: (SB[2]-450, 250),     # hdcoa + FA pair
        '1073': box(0, 3),                                                                # 3-oxo C18
        '505': box(1, 0), '506': box(1, 1), '677': (SB[3]+450, HB-250), '694': (SB[3]+450, HB+250),
        '1061': box(2, 3),                                                                # 3-oxo C10
        '1060': box(3, 0), '1048': box(3, 3),                                             # octanoyl, 3-oxo C8
        '608': box(4, 0), '461': box(4, 1), '459': box(4, 2), '471': box(4, 3),           # C6 round
        '647': (-450, 4*HB-250), '664': (-450, 4*HB+250),
        '571': box(5, 0), '485': box(5, 1), '484': box(5, 2), '716': box(5, 3),           # C4 round
        '715': (0, 6*HB)})
    L.reaction('70', 'BOX', side='U', dist=110); L.reaction('1166', 'BOX', side='D', dist=110)   # odecoa pair
    L.reaction('570', 'BOX', side='U', dist=110); L.reaction('1153', 'BOX', side='D', dist=110)  # hdcoa pair
    L.reaction('1104', 'BOX', side='L')                                                  # CTECOAI7
    L.reaction('1084', 'BOX', side='R')                                                  # ACACT8r
    L.reaction('518', 'BOX', shift=(0, -75), side='U', dist=120); L.reaction('552', 'BOX', shift=(0, 75), side='D', dist=120)
    L.reaction('691', 'BOX', side='U', dist=110); L.reaction('704', 'BOX', side='D', dist=110)   # pmtcoa pair
    L.reaction('1072', 'BOX', side='R')                                                  # ACACT4r
    L.reaction('1059', 'BOX', side='L')                                                  # ACACT3r
    L.reaction('661', 'BOX', side='U', dist=110)                                         # FACOAL60t2pp
    u = (-450/math.hypot(450, 250), 250/math.hypot(450, 250)); n = (u[1], -u[0])          # hxcoa -> hxa_c
    L.reaction('674', 'BOX', shift=(n[0]*60, n[1]*60), side='R', dist=100)               # FACOAE60
    L.reaction('1215', 'BOX', shift=(-n[0]*60, -n[1]*60), side='L', dist=100)            # HXCT
    L.reaction('620', 'BOX', shift=(0, -75), side='U', dist=120); L.reaction('632', 'BOX', shift=(0, 75), side='D', dist=120)
    L.reaction('470', 'BOX', shift=(0, -75), side='U', dist=120); L.reaction('540', 'BOX', shift=(0, 75), side='D', dist=120)
    L.reaction('483', 'BOX', side='D')                                                   # HACD2
    L.reaction('644', 'BOX', side='R')                                                   # ACACT2r
    for rid, off, sd, d in (('583', -225, 'U', 110), ('595', -75, 'U', 75), ('607', 75, 'D', 75), ('1096', 225, 'D', 110)):
        L.reaction(rid, 'BOX', shift=(0, off), side=sd, dist=d)                           # ACOAD1f x3, ACOAD1fr
    L.reaction('495', 'BOX', shift=(0, -75), side='U', dist=120); L.reaction('504', 'BOX', shift=(0, 75), side='D', dist=120)
    L.reaction('726', 'BOX', side='L')                                                   # ACACT1r

    # ================================================================ PLS: one column per acyl chain
    def chain(bx, A, P, F, G, LPA, rAP, rPF, rAF, rG3, rAPG):
        L.setp('PLS', {A: (bx, -400), F: (bx, 1100), G: (bx+500, 0), LPA: (bx+1000, 0)})
        if P: L.setp('PLS', {P: (bx, 400)})
        if rAP: L.reaction(rAP, 'PLS', side='R', dist=100)
        L.reaction(rPF, 'PLS', side='R', dist=100)
        L.reaction(rAF, 'PLS', shift=(-300, 0), side='L', dist=110)
        L.reaction(rG3, 'PLS', ends=(G, LPA), shift=(0, -80), side='U', dist=120)
        L.reaction(rAPG, 'PLS', ends=(G, LPA), shift=(0, 80), side='D', dist=120)
    chain(0, '2', '1744', '4', '1167', '1168', '2024', '1793', '14', '1179', '1755')                  # 16:0
    L.setp('PLS', {'1733': (1500, 0)}); L.reaction('1743', 'PLS', ends=('1168', '1733'), side='U')   # AGPAT160
    chain(2200, '1830', '1769', '1156', S['g3p_181'], '1770', None, '1817', '2012', '1841', '1781')  # 18:1
    chain(4400, '1818', '1756', '1143', '417', '1757', None, '1805', '1996', '1829', '1768')          # 16:1
    L.setp('PLS', {'1967': (4400, -900)}); L.reaction('1980', 'PLS', side='R')                       # 3OAS181
    # ================================================================ GLYC
    L.setp('GLYC', {S['g3p_glyc']: (0, 0), '419': (0, -400), '352': (500, 0), '350': (1000, 0), '1192': (1500, 0)})
    L.reaction('430', 'GLYC', side='L'); L.reaction('1191', 'GLYC', side='D')
    L.reaction('363', 'GLYC', side='D'); L.reaction('1204', 'GLYC', side='D')

    # ================================================================ PLD: identical units
    HG = ('A', 'E', 'G'); CH = ('160', '161', '181')
    for i, hg in enumerate(HG):
        for j, ch in enumerate(CH):
            p1, p2, l1 = PERI[(hg, ch)]
            x0, y0 = 1100*j, 1350*i                                  # periplasm unit
            pX = [n for n in L.top[p1]['S'] if L.is_primary(n)][0]
            a2 = [n for n in L.top[p1]['P'] if L.is_primary(n)][0]
            a1 = [n for n in L.top[p2]['P'] if L.is_primary(n)][0]
            L.setp('PLD', {pX: (x0+250, y0), a2: (x0, y0+450), a1: (x0+500, y0+450), S[f'bb_p_{hg}{ch}']: (x0+500, y0+900)})
            L.reaction(p1, 'PLD', side='L', dist=110); L.reaction(p2, 'PLD', side='R', dist=110)
            L.reaction(l1, 'PLD', side='R', dist=110)
            at, l2 = CYTO[(hg, ch)]                                   # cytosol unit (column)
            x0 = 3300 + 700*j
            a2 = [n for n in L.top[l2]['S'] if L.is_primary(n)][0]
            bb = S[f'g3p_{hg}{ch}'] if hg == 'A' else S[f'bb_c_{hg}{ch}']
            L.setp('PLD', {a2: (x0, y0+450), bb: (x0, y0+900)})
            L.reaction(l2, 'PLD', side='R', dist=110)
            if at:
                pc = [n for n in L.top[at]['P'] if L.is_primary(n)][0]
                L.setp('PLD', {pc: (x0, y0)}); L.reaction(at, 'PLD', side='R', dist=110)

    # ================================================================ LIPA
    L.setp('LIPA', {'218': (0, 0), '220': (0, 320), '234': (0, 640), '246': (0, 960), '297': (-450, 1300),
                    '298': (0, 1640), '272': (0, 1960), '271': (-450, 2120), '206': (0, 2280), '191': (0, 2600),
                    '193': (0, 2920), '260': (0, 3240), '2025': (500, 2250), '2027': (700, 2600),
                    S['hdcea_lipa']: (1050, 2250),
                    '1439': (700, 3000), '1442': (1300, 3000), '1244': (700, 3450), '1246': (1300, 3450),
                    '1229': (700, 3800), '1232': (1300, 3800),
                    S['ckdo_2']: (650, 600), '1900': (650, 1000), '1902': (1150, 800)})
    L.reaction('231', 'LIPA', side='R'); L.reaction('244', 'LIPA', side='R'); L.reaction('257', 'LIPA', side='R')
    L.reaction('321', 'LIPA', side='L')                                   # USHD
    L.reaction('309', 'LIPA', ends=('246', '298'), side='R')              # LPADSS
    L.reaction('333', 'LIPA', side='R')                                   # TDSK
    L.reaction('283', 'LIPA', ends=('272', '206'), side='R')              # MOAT
    L.reaction('217', 'LIPA', side='R')                                   # MOAT2
    L.reaction('204', 'LIPA', side='L'); L.reaction('270', 'LIPA', side='L')   # EDTXS1, EDTXS2
    L.reaction('296', 'LIPA', axis=(0, 1), mid=(0, 3420), side='L', end_cof=True)   # HEPT1
    L.reaction('2037', 'LIPA', ends=('191', '2027'), side='D')            # EDTXS3
    L.reaction('2049', 'LIPA', side='U')                                  # FA161ACPHi
    L.reaction('1452', 'LIPA', shift=(0, -70), side='U', dist=110); L.reaction('1464', 'LIPA', shift=(0, 70), side='D', dist=110)
    L.reaction('1256', 'LIPA', side='D'); L.reaction('1242', 'LIPA', side='D')
    L.reaction('1912', 'LIPA', ends=('1900', '1902'), side='D')          # MOAT3C

    # ================================================================ ISO (MEP pathway)
    L.setp('ISO', {'71': (-200, -350), '72': (200, -350), '74': (0, 0), '86': (0, 320), '99': (0, 640),
                   '112': (0, 960), '124': (0, 1280), '135': (0, 1600), '151': (-300, 2000), '136': (300, 2000),
                   '166': (300, 2400), '179': (300, 2800)})
    L.reaction('84', 'ISO', axis=(0, 1), mid=(0, -170), side='R')        # DXPS
    for r in ('97', '110', '123', '134', '1423'): L.reaction(r, 'ISO', side='R')
    L.reaction('164', 'ISO', side='L'); L.reaction('150', 'ISO', side='R')   # IPDPS, DMPPS (fan)
    L.reaction('714', 'ISO', side='D')                                    # IPDDI
    L.reaction('177', 'ISO', side='R'); L.reaction('190', 'ISO', side='R')

    # ================================================================ ECA
    L.setp('ECA', {'1003': (0, 0), '1006': (500, 0), '1938': (0, 450), '1940': (500, 450)})
    L.reaction('1016', 'ECA', side='D'); L.reaction('1950', 'ECA', side='D')

NEW_EXITS = [('2', '← fatty acid biosynthesis'), ('1818', '→ lipid A (EDTXS3)')]
EXTRA_TITLES = [dict(text='periplasm', pos=(6250, 250)), dict(text='cytoplasm', pos=(9600, 250))]
TITLE_POS = {'FAS': (100, -150), 'BOX': (3000, -150), 'PLS': (200, 6400), 'GLYC': (9900, 6000),
             'PLD': (6250, -150), 'LIPA': (6300, 4600), 'ISO': (8400, 4600), 'ECA': (9900, 4600)}
CURVES = dict(t=0.6, stacked=[], special_links={},
              bend_override={'661': 0, '691': 0, '570': 0, '70': 0, '583': 15, '595': 15, '607': 15, '1096': 15})
# The hand arrangement is final: two 1x1 keep-out boxes at the corners of the hand-laid-out extent make
# the packer's reference rectangle the whole map, so every detached piece stays (within one 40-unit
# grid cell) where it was placed instead of being dropped into a hole of another zone.
PACK = dict(keep_out=[('box', (-400, -300, -399, -299)), ('box', (11650, 9000, 11651, 9001))], auto_ring_keep_out=True)
RENDER_SCALE = 0.25
