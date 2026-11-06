"""SPEC - iML1515 'set aside' draft (BRENDA pathways not drawn in the main maps, 45 reactions)
-> hand-laid-out Escher map.

Every coordinate of a PRIMARY metabolite below was chosen by hand. Zone coordinates are local; the
zone origin (ORIGIN) places the zone on the map. Escher convention: x to the right, y DOWN.

Global arrangement (main connected piece):
  * methylaspartate cycle = TCA-like ring (R = 800, citrate on top, clockwise), open between
    2-oxoglutarate and succinyl-CoA (no AKGDH in this map); succinyl-CoA -> methylmalonyl-CoA outwards.
  * acetyl-CoA hub above the ring (feeds CS).  Into the hub: aerobic toluene degradation chain from
    above, PDH from pyruvate on the left; out of the hub: ACACT1r (butanoate fermentation).
  * C4 / CAM carbon-fixation cycle as a square left of the hub: PEP (top left) -> OAA (bottom left)
    -> malate (bottom right) -> pyruvate (top right) -> PEP.  All PPC / MDH / ME copies of the draft
    (also the ones listed under 'photosynthesis') are drawn on this square as lanes.
Detached pieces (packed automatically): Calvin-cycle fragment ('photosynthesis'), PGK, ECOAH1, BUTCT."""
import math

INPUT = '../../ecmapper/BRENDAmap/draft_map_iML1515/split/set_aside.escher.json'

# ------------------------------------------------------------------ zones, titles, colours
ZONES = {
 'MC':  dict(title='methylaspartate cycle', title_id='585', color='#d62728'),
 'C4':  dict(title='C4 and CAM-carbon fixation', title_id='581', color='#9467bd'),
 'ACB': dict(title='acetyl CoA biosynthesis', title_id='582', color='#8c564b'),
 'ATD': dict(title='aerobic toluene degradation', title_id='583', color='#ff7f0e'),
 'BF':  dict(title='butanoate fermentation', title_id='584', anchor=('250', (-60, -150)), color='#bcbd22'),
 'BF2': dict(title=None, color='#bcbd22'),           # detached pieces, titles in EXTRA_TITLES
 'BF3': dict(title=None, color='#bcbd22'),
 'PHO': dict(title='photosynthesis', title_id='586', anchor=('452', (-350, -620)), color='#2ca02c'),
 'PGK': dict(title=None, color='#2ca02c'),
}
TITLE_POS = {'MC': (-310, 40), 'C4': (-2380, -560), 'ACB': (-1150, -960), 'ATD': (-80, -2050)}
EXTRA_TITLES = [dict(text='butanoate fermentation', anchor='261', offset=(-60, -110), zone='BF2'),
                dict(text='butanoate fermentation', anchor='274', offset=(-60, -110), zone='BF3'),
                dict(text='photosynthesis', anchor='437', offset=(-60, -110), zone='PGK')]

EXIT_OWNER = {}        # automatic owner rule (both '→ Central' labels sit on acetyl-CoA copies)

# acetyl-CoA copies 249 and 301 both carry an identical '→ Central' exit label; after merging them
# into one hub the duplicate label 579 is dropped (annotation only, no biochemistry).
DROP_LABELS = ['579']

# ------------------------------------------------------------------ global arrangement
ORIGIN = {'MC': (0, 0), 'ATD': (0, 0), 'ACB': (0, 0), 'BF': (0, 0), 'C4': (-2700, -1300),
          'PHO': (1900, -2300), 'PGK': (3300, -2300), 'BF2': (1500, -3000), 'BF3': (2300, -3000)}

MERGES = [('301', ['183', '223', '249']),            # acetyl-CoA hub
          ('32', ['5', '103', '536', '182']),        # pyruvate (C4 square + PDH substrate)
          ('63', ['45', '480']),                     # PEP
          ('60', ['46', '131', '481']),              # oxaloacetate
          ('115', ['1', '87', '532'])]               # L-malate

R = 800
def ring(k): a = math.radians(-90+40*k); return (R*math.cos(a), R*math.sin(a))
def polar(deg, r): a = math.radians(deg); return (r*math.cos(a), r*math.sin(a))

def layout(L):
    # ------------------------------------------------ methylaspartate cycle: ring, centre (0, 0)
    L.setp('MC', {'303': ring(0), '364': ring(1), '393': ring(2), '413': ring(3), '379': ring(4),
                  '329': ring(5), '316': ring(6), '287': ring(7), '290': ring(8),
                  '403': polar(50, 1250)})
    C = (0, 0)
    L.reaction('315', 'MC', arc=(C, 0), ends=('290', '303'), side='out')   # CS (+ acetyl-CoA hub)
    L.reaction('374', 'MC', arc=(C, 0), side='out')                        # ACONTa
    L.reaction('402', 'MC', arc=(C, 0), side='in')                         # ACONTb
    L.reaction('426', 'MC', arc=(C, 0), side='out')                        # ICDHyr
    L.reaction('391', 'MC', arc=(C, 0), side='out')                        # SUCOAS (succ -> succoa)
    L.reaction('411', 'MC', side='U')                                      # MMM (outwards)
    L.reaction('326', 'MC', arc=(C, 0), side='in')                         # FUM
    L.reaction('300', 'MC', arc=(C, 0), side='in')                         # MDH (mal -> oaa)
    for rid, k in (('351', -1), ('363', 0), ('339', 1)):                   # FRD3 | SUCDi | FRD2
        L.reaction(rid, 'MC', arc=(C, k*180), side='in' if k < 0 else 'out', dist=110 if k else 90)
    # ------------------------------------------------ acetyl-CoA hub + toluene degradation (global)
    L.setp('ATD', {'301': (-190, -1300), '210': (-190, -1650), '200': (-190, -1950), '199': (-190, -2250),
                   '237': (450, -1300)})
    L.reaction('209', 'ATD', side='L')                                     # OP4ENH
    L.reaction('220', 'ATD', side='L')                                     # HOPNTAL (pyr as cofactor)
    L.reaction('235', 'ATD', side='L')                                     # ACALD
    L.reaction('248', 'ATD', side='D')                                     # 3OXCOAT
    # ------------------------------------------------ butanoate fermentation: ACACT1r out of the hub
    L.setp('BF', {'250': (380, -1800)})
    L.reaction('260', 'BF', side='U')                                      # ACACT1r
    # ------------------------------------------------ PDH: pyruvate (C4 square corner) -> hub
    L.reaction('197', 'ACB', side='D', step=60)
    # ------------------------------------------------ C4 / CAM square, side 1400 (local origin = PEP)
    L.setp('C4', {'63': (0, 0), '60': (0, 1400), '115': (1400, 1400), '32': (1400, 0),
                  '154': (1750, -550), '18': (1050, -550), '129': (-550, 1400)})
    for rid, y in (('58', 450), ('493', 900)):                             # PPC x2 stacked, left edge
        L.reaction(rid, 'C4', side='L', axis=(0, 1), mid=(0, y), gap=40)
    L.reaction('74', 'C4', side='R', shift=(140, 0), dist=110)             # PPCK reverse lane (inside)
    for rid, x in (('99', 350), ('127', 700)):                             # MDH x2 stacked, bottom edge
        L.reaction(rid, 'C4', side='D', axis=(-1, 0), mid=(x, 1400), gap=40)
    for rid, y in (('16', 1150), ('114', 850), ('547', 550), ('179', 250)):   # ME2 x3 + ME1 stacked
        L.reaction(rid, 'C4', side='L', axis=(0, -1), mid=(1400, y), gap=40)
    L.reaction('86', 'C4', side='U')                                       # PYK (top edge)
    L.reaction('168', 'C4', side='R')                                      # MALDDH (mal__D from upper right)
    L.reaction('29', 'C4', shift=(-60, -40), side='L')                     # ALATA_L lane 1
    L.reaction('42', 'C4', shift=(60, 40), side='R')                       # ALATA_L lane 2
    L.reaction('141', 'C4', shift=(0, -70), side='U')                      # ASPTA lane 1
    L.reaction('153', 'C4', shift=(0, 70), side='D')                       # ASPTA lane 2
    # ------------------------------------------------ photosynthesis (Calvin-cycle fragment), detached
    L.setp('PHO', {'452': (0, 0), '427': (-700, 0), '428': (-350, -400), '451': (350, -400),
                   '469': (750, -400), '521': (-350, 450), '495': (650, 0), '507': (650, 450),
                   '548': (650, 900), '559': (1050, 900), '568': (250, 900)})
    L.reaction('436', 'PHO', lab='L')                                      # RPI
    L.reaction('468', 'PHO', lab='R')                                      # RPE
    L.reaction('477', 'PHO', side='U')                                     # X5PL3E
    L.reaction('460', 'PHO', lab='L')                                      # RBP4E
    L.reaction('531', 'PHO', axis=(0, 1), mid=(-350, 225), side='L')       # TKT1 (r5p + xu5p -> s7p)
    L.reaction('506', 'PHO', side='U')                                     # TKT2
    L.reaction('519', 'PHO', side='R')                                     # FBP
    L.reaction('558', 'PHO', side='L')                                     # FBA
    L.reaction('567', 'PHO', side='D')                                     # TPI
    L.reaction('578', 'PHO', side='D')                                     # FBA3
    L.setp('PGK', {'437': (0, 0), '439': (0, 400)})
    L.reaction('450', 'PGK', side='R')                                     # PGK
    L.setp('BF2', {'261': (0, 0), '262': (400, 0)})
    L.reaction('272', 'BF2', side='D')                                     # ECOAH1
    L.setp('BF3', {'274': (0, 0), '276': (500, 0)})
    L.reaction('286', 'BF3', side='D')                                     # BUTCT

NEW_EXITS = []
CURVES = dict(t=0.6, bend=55, bend_stack=15,
              stacked=['58', '493', '99', '127', '16', '114', '547', '179'],
              bend_override={}, special_links={}, nested_arc_lanes=True)
PACK = dict(keep_out=[('box', (-2600, -1200, -1400, 0))], auto_ring_keep_out=True)  # C4 square interior
RENDER_SCALE = 0.5
