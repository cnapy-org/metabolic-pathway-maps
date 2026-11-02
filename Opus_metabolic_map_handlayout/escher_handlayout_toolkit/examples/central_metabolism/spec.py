"""SPEC - iML1515 central metabolism (BRENDA draft) -> hand-laid-out Escher map.
Reproduces the delivered v7 layout (one lane per reaction in ring families, cofactors in the gaps).

Every coordinate of a PRIMARY metabolite below was chosen by hand. Zone coordinates are local; the
zone origin (ORIGIN) places the zone on the map. Escher convention: x to the right, y DOWN.
Typical spacing between neighbours of a chain ~300 map units."""
import math

INPUT = 'central_metabolism_escher.json'

# ------------------------------------------------------------------ zones, titles, colours
ZONES = {
 'GLY': dict(title='glycolysis', title_id='1162', color='#1f77b4'),
 'PPP': dict(title='pentose phosphate pathway', title_id='1166', color='#2ca02c'),
 'ED':  dict(title='Entner Doudoroff pathway', title_id='1157', color='#17becf'),
 'ED2': dict(title='ED: glycerate branch', color='#17becf'),
 'PDHB':dict(title='pyruvate → acetyl-CoA', title_id='1164', color='#8c564b'),   # renamed: holds PDH, PFL, POR5
 'LAC': dict(title='lactate fermentation', title_id='1163', color='#e377c2'),
 'ETH': dict(title='ethanol fermentation', title_id='1160', color='#bcbd22'),
 'ACE': dict(title='acetate fermentation', title_id='1158', color='#ff7f0e'),
 'ACX': dict(title='oxalyl-CoA / formate exchange', color='#ff7f0e'),
 'GNG': dict(title='gluconeogenesis', title_id='1161', color='#9467bd'),
 'TCA': dict(title='citric acid cycle', title_id='1159', color='#d62728'),
 'FEX': dict(title='Fe(III)-dicitrate uptake', color='#d62728'),
 'OXN': dict(title=None, color='#7f7f7f'),            # the 3 NADH rows get their own titles (EXTRA_TITLES)
 'ATP': dict(title='oxidative phosphorylation: ATP synthase / ATP use', title_id='1165',
             anchor=('740', (-480, -560)), color='#555555'),
}
TITLE_POS = {'GLY': (-700, -120), 'PPP': (1700, 180), 'ED': (250, -900), 'ED2': (1900, -700),
             'PDHB': (-900, 3800), 'LAC': (900, 2980), 'ETH': (1300, 4120), 'ACE': (1450, 5050),
             'ACX': (2000, 4180), 'GNG': (-1500, 2850), 'TCA': (-400, 6000), 'FEX': (1650, 5480)}
EXTRA_TITLES = [dict(text='NADH dehydrogenase (ubiquinone)', anchor='686', offset=(-40, -260), zone='OXN'),
                dict(text='NADH dehydrogenase (menaquinone)', anchor='672', offset=(-40, -260), zone='OXN'),
                dict(text='NADH dehydrogenase (demethylmenaquinone)', anchor='708', offset=(-40, -260), zone='OXN')]

EXIT_OWNER = {}        # automatic owner rule; add {label_id: node_id} to override single labels.
# note: label 1155 '→ Amino acids' is assigned to trans-aconitate (same-row pattern of the draft);
# the original session used cis-aconitate ('1155': '136').

# ------------------------------------------------------------------ splits (copies of primaries)
def prepare(L):
    S = L.S
    S['accoa_E'] = L.split('343', '369')     # ACALD keeps its own acetyl-CoA copy ...
    L.merge('656', '342'); L.merge('657', '343')   # ... PFL/POR5 join the PDH block (pyr/accoa)
    S['xu5p_K'] = L.split('871', '909')      # xu5p for TKT2 (hourglass chain)
    S['f6p_K']  = L.split('898', '909')      # f6p produced by TKT2 (stays separate)
    S['succ_I'] = L.split('251', '260')      # succinate from ICL (inside the TCA ring)
    S['accoa_M']= L.split('147', '250')      # acetyl-CoA consumed by MALS (inside the ring)
    S['accoa_X']= L.split('78', '956')       # ACOXT block gets its own acetyl-CoA ...
    S['ac_X']   = L.split('63', '956')       # ... and acetate -> becomes a detached piece

# ------------------------------------------------------------------ global arrangement
ORIGIN = {'GLY': (0, 0), 'PPP': (400, 300), 'ED': (1200, 0), 'ED2': (1900, -450),
          'PDHB': (0, 3180), 'LAC': (0, 3180), 'GNG': (0, 3180), 'ETH': (0, 4240), 'ACE': (0, 4240),
          'ACX': (2750, 4300), 'TCA': (0, 6340), 'FEX': (1700, 5600), 'OXN': (-1300, 7950), 'ATP': (600, 8250)}
FLIP = {'ED': -1}       # ED mirrored vertically: hangs ABOVE 6pgc

MERGES = [('480', ['923']), ('482', ['823']), ('493', ['898']), ('528', ['899']), ('592', ['387']),
          ('604', ['656', '616', '414']), ('657', ['147', '78', 'accoa_E']), ('355', ['60']),
          ('149', ['384']), ('807', ['1'])]

def ring(k, R=1050):           # TCA: 9 nodes every 40 deg, clockwise, citrate on top
    a = math.radians(-90+40*k); return (R*math.cos(a), R*math.sin(a))
def polar(deg, r): a = math.radians(deg); return (r*math.cos(a), r*math.sin(a))

def layout(L):
    S = L.S
    # ------------------------------------------------ glycolysis: vertical spine at x = 0
    L.setp('GLY', {'480': (0, 0), '482': (0, 300), '493': (0, 600), '504': (0, 900),   # glc g6p f6p fdp
                   '527': (-320, 1250), '528': (0, 1250),                              # dhap (left), g3p
                   '548': (0, 1550), '561': (0, 1850), '574': (0, 2150), '592': (0, 2450),
                   '604': (0, 3180)})                                                  # pyr (4 stacked PYKs)
    L.reaction('492', 'GLY', side='R')                              # HEX1
    L.reaction('1017', 'GLY', side='L', shift=(-120, 0), dist=110)  # G6PP reverse lane
    L.reaction('501', 'GLY', lab='L')                               # PGI
    L.reaction('514', 'GLY', side='R')                              # PFK
    L.reaction('526', 'GLY', side='L', shift=(-120, 0), dist=110)   # FBP reverse lane
    L.reaction('537', 'GLY', axis=(0, 1), mid=(-70, 1075), lab='L') # FBA fan
    L.reaction('545', 'GLY', lab='D')                               # TPI
    L.reaction('560', 'GLY', side='R')                              # GAPD
    L.reaction('573', 'GLY', side='R')                              # PGK (model direction 3pg -> 13dpg)
    L.reaction('582', 'GLY', shift=(-70, 0), lab='L')               # PGM lane 1
    L.reaction('590', 'GLY', shift=(70, 0), lab='R')                # PGM lane 2
    L.reaction('601', 'GLY', side='R')                              # ENO
    for rid, y in (('614', 2620), ('1070', 2760), ('1082', 2900), ('1094', 3040)):   # PYK, PYK2-4 stacked
        L.reaction(rid, 'GLY', side='R', axis=(0, 1), mid=(0, y), gap=40, lab='L')
    # ------------------------------------------------ pyruvate -> acetyl-CoA: PFL, PDH, POR5 stacked
    L.setp('PDHB', {'656': (0, 0), '657': (0, 1060)})
    for rid, y in (('354', 200), ('671', 490), ('1058', 810)):
        L.reaction(rid, 'PDHB', side='R', axis=(0, 1), mid=(0, y), gap=40, step=60, lab='L')
    # ------------------------------------------------ lactate, right of pyruvate
    L.setp('LAC', {'616': (0, 0), '641': (700, -300), '615': (700, 300)})
    L.reaction('653', 'LAC', side='U')                                           # LDH_D
    for rid, t in (('628', 0.33), ('640', 0.62)):                                # L_LACD2/3 stacked
        L.reaction(rid, 'LAC', side='U', mid=(700*(1-t), 300*(1-t)), gap=40)
    # ------------------------------------------------ gluconeogenesis / anaplerosis, left of spine
    L.setp('GNG', {'414': (0, 0), '387': (0, -700), '454': (-800, 0), '439': (-800, 350), '384': (-650, 1300)})
    L.reaction('411', 'GNG', axis=(0, 1), mid=(-1180, 0), side='R')    # PPC  pep -> oaa, routed far left
    L.reaction('398', 'GNG', axis=(0, -1), mid=(-1320, 0), side='L')   # PPCK
    L.reaction('428', 'GNG', shift=(-140, 0), side='L', dist=150)      # PPS reverse lane of PYK
    for rid, x in (('465', -600), ('478', -330)):                      # ME1, ME2 stacked
        L.reaction(rid, 'GNG', axis=(1, 0), mid=(x, 0), side='U', gap=40, step=60)
    L.reaction('453', 'GNG', side='D')                                 # MALDDH
    L.reaction('438', 'GNG', side='R')                                 # OAADC
    # ------------------------------------------------ TCA ring, R = 1050, centre = zone origin
    L.setp('TCA', {'135': ring(0), '136': ring(1), '163': ring(2), '223': ring(3), '265': ring(4),
                   '251': ring(5), '278': ring(6), '173': ring(7), '149': ring(8),
                   '147': (-250, -1550), '237': (0, -182), S['succ_I']: (500, 150), S['accoa_M']: (-560, 150),
                   '957': polar(245, 720), '1108': polar(-50, 1450), '967': (430, -480), '969': (400, -840)})
    C = (0, 0)
    L.reaction('161', 'TCA', arc=(C, 0), ends=('149', '135'), side='out')     # CS
    L.reaction('966', 'TCA', arc=(C, -150), ends=('135', '149'), side='in')   # CITL inner lane
    L.reaction('146', 'TCA', arc=(C, 0), side='out')                          # ACONTa
    L.reaction('172', 'TCA', arc=(C, 0), side='in')                           # ACONTb
    L.reaction('1116', 'TCA', lab='R')                                        # ACONIs
    L.reaction('236', 'TCA', arc=(C, 0), side='out')                          # ICDHyr
    L.reaction('277', 'TCA', arc=(C, 0), side='out')                          # AKGDH
    L.reaction('304', 'TCA', arc=(C, 0), side='out')                          # SUCOAS
    L.reaction('288', 'TCA', arc=(C, 0), side='in')                           # FUM
    L.reaction('260', 'TCA', ends=('163', '237'), lab='U')                    # ICL (shunt)
    L.reaction('250', 'TCA', ends=('237', '173'), side='U')                   # MALS (shunt)
    L.reaction('981', 'TCA', ends=('967', '135'), side='L')                   # FE3DCITR5
    # ring families: one lane per reaction, 180 apart; cofactors in the gap next to their lane
    GAPW = 180
    for rid, k in (('209', -1), ('221', 0), ('185', 1), ('197', 2), ('1029', 3)):   # MDH2|MDH3|MDH|MDH|MOX
        L.reaction(rid, 'TCA', arc=(C, k*GAPW), side='in' if k < 0 else 'out', dist={-1: 80, 3: 110}.get(k, GAPW/2))
    for rid, k in (('328', -1), ('340', 0), ('316', 1)):                            # FRD3|SUCDi|FRD2
        L.reaction(rid, 'TCA', arc=(C, k*GAPW), side='in' if k < 0 else 'out', dist=110 if k in (-1, 1) else GAPW/2)
    L.setp('FEX', {'982': (0, 0), '983': (0, 300), '984': (400, 150)})
    L.reaction('993', 'FEX', lab='R')                                         # FE3DCITexs
    # ------------------------------------------------ ethanol and acetate, right of acetyl-CoA
    L.setp('ETH', {S['accoa_E']: (0, 0), '355': (650, 0), '370': (1150, 0), '1117': (650, -400)})
    L.reaction('369', 'ETH', side='U'); L.reaction('382', 'ETH', side='U'); L.reaction('1127', 'ETH', side='R')
    L.setp('ACE', {'78': (0, 0), '60': (650, 0), '63': (1150, 500), '94': (500, 700), '119': (1700, 500)})
    L.reaction('75', 'ACE', side='U'); L.reaction('92', 'ACE', side='D'); L.reaction('117', 'ACE', side='L')
    L.reaction('105', 'ACE', side='D'); L.reaction('134', 'ACE', side='D')
    L.setp('ACX', {S['accoa_X']: (0, 0), '945': (0, 300), '994': (0, 600), S['ac_X']: (-400, 0),
                   '946': (-400, 300), '995': (-400, 600), '1096': (-575, 1000), '1097': (-750, 600)})
    L.reaction('956', 'ACX', lab='U'); L.reaction('1005', 'ACX', lab='D'); L.reaction('1107', 'ACX', side='R')
    # ------------------------------------------------ NADH rows and ATP hub
    for i, (q, qh, up, dn) in enumerate([('686', '687', '729', '697'), ('672', '674', '685', '707'), ('708', '709', '739', '719')]):
        L.setp('OXN', {q: (0, 470*i), qh: (650, 470*i)})
        L.reaction(up, 'OXN', shift=(0, -65), side='U', lab='D', dist=100)
        L.reaction(dn, 'OXN', shift=(0, 65), side='D', lab='U', dist=100)
    L.setp('ATP', {'740': (0, 0), '775': (0, 700)})
    L.reaction('753', 'ATP', axis=(-0.87, -0.5), mid=(-290, -170), side='L', end_cof=True, step=90)  # ATPM
    L.reaction('773', 'ATP', axis=(0.87, -0.5), mid=(290, -170), side='R', end_cof=True, step=90)    # NTP1
    L.reaction('763', 'ATP', axis=(0, 1), mid=(0, -330), side='R', end_cof=True, step=90)            # ATPS4rpp
    L.reaction('786', 'ATP', shift=(-260, 0), side='L'); L.reaction('796', 'ATP', side='R', dist=100)
    L.reaction('806', 'ATP', shift=(260, 0), side='R')
    # ------------------------------------------------ pentose phosphate pathway
    L.setp('PPP', {'923': (0, -300), '823': (0, 0), '825': (400, 0), '807': (800, 0), '811': (1250, 0),
                   '861': (1250, 650), '871': (1550, 950), '870': (1850, 950), '888': (2150, 950),
                   '911': (900, 650), '910': (900, 950), '1032': (900, 400), '898': (550, 650), '897': (550, 950),
                   S['xu5p_K']: (550, 1250), S['f6p_K']: (200, 950), '899': (200, 1250)})
    L.reaction('934', 'PPP', shift=(260, 0), side='R', dist=120)      # HEX1 duplicate: 2nd lane
    L.reaction('836', 'PPP', side='D'); L.reaction('860', 'PPP', side='U')           # G6PDH2r, PGL
    L.reaction('822', 'PPP', shift=(0, -70), side='U'); L.reaction('850', 'PPP', shift=(0, 70), side='D', lab='U')  # GND x2
    L.reaction('869', 'PPP', lab='L'); L.reaction('887', 'PPP', lab='R')              # RPI, RPE
    L.reaction('879', 'PPP', lab='U'); L.reaction('896', 'PPP', lab='U')              # RBP4E, X5PL3E
    L.reaction('921', 'PPP', lab='U')                                                 # TKT1 hourglass
    L.reaction('944', 'PPP', axis=(-1, 0), mid=(725, 800), lab='U')                   # TALA (waist pinned)
    L.reaction('909', 'PPP', axis=(-1, 0), mid=(375, 1100), lab='U')                  # TKT2 (waist pinned)
    L.reaction('1042', 'PPP', side='L')                                               # PFK_3
    # ------------------------------------------------ Entner-Doudoroff (zone mirrored: 'D' == up)
    L.setp('ED', {'1': (0, 0), '2': (0, 350), '13': (400, 350), '26': (-250, 750), '27': (250, 750)})
    L.reaction('12', 'ED', side='L'); L.reaction('25', 'ED', side='D'); L.reaction('36', 'ED', lab='L')
    L.setp('ED2', {'48': (0, 0), '37': (350, 0), '38': (700, 0)})
    L.reaction('59', 'ED2', side='U'); L.reaction('46', 'ED2', lab='U')

NEW_EXITS = [('26', '→ Glycolysis'), ('27', '→ Glycolysis'), ('38', '→ Glycolysis'),
             ('f6p_K', '→ Glycolysis'), ('119', '← Glycolysis')]

FAM = ['209', '221', '185', '197', '1029', '328', '340', '316']
CURVES = dict(t=0.6, bend=55, bend_stack=15,
              stacked=['614', '1070', '1082', '1094', '354', '671', '1058', '628', '640', '465', '478'],
              bend_override={'560': 0, '492': 0, '1127': 5, **{r: 30 for r in FAM}},
              special_links={('944', '493'): (1, 0), ('909', '528'): (1, -0.25), ('161', '657'): (-0.25, 1)},
              nested_arc_lanes=True)
PACK = dict(keep_out=[('box', (-1250, 2450, 0, 5300))], auto_ring_keep_out=True)   # + TCA interior (auto)
RENDER_SCALE = 0.5
