"""SPEC - iML1515 energy and C1 metabolism (BRENDA draft, 42 reactions) -> hand-laid-out Escher map.

Every coordinate of a PRIMARY metabolite below was chosen by hand (zone-local coordinates, zone origin
in ORIGIN). Escher convention: x to the right, y DOWN. Typical spacing between chain neighbours ~300.
The map is a set of unconnected pathway fragments; the toolkit packs the detached pieces (Phase 5).
Zone titles are anchored to a node of their piece so that they travel with it."""
import math

INPUT = '../../ecmapper/BRENDAmap/draft_map_iML1515/split/energy_and_c1_metabolism.escher.json'

ZONES = {
 'C1':  dict(title='formaldehyde oxidation', title_id='588', anchor=('125', (-520, -120)), color='#1f77b4'),
 'THF': dict(title='reductive acetyl coenzyme A pathway', title_id='591', anchor=('558', (-60, 230)), color='#17becf'),
 'UB':  dict(title='ubiquinone biosynthesis', title_id='592', anchor=('306', (-150, -170)), color='#2ca02c'),
 'HPB': dict(title='CO2 fixation in Crenarchaeota (3-hydroxypropionate/4-hydroxybutyrate cycle)', title_id='587',
             anchor=('83', (-900, -330)), color='#d62728'),
 'GLX': dict(title='glycolate and glyoxylate degradation', title_id='589', anchor=('179', (-200, -560)), color='#ff7f0e'),
 'UC':  dict(title='urea cycle', title_id='593', anchor=('375', (-100, -170)), color='#9467bd'),
 'MM':  dict(title='methane metabolism', title_id='590', anchor=('214', (-80, -150)), color='#7f7f7f'),
}

EXIT_OWNER = {}

def prepare(L):
    pass

ORIGIN = {'C1': (0, 0), 'THF': (0, 0), 'HPB': (2650, 1000), 'UB': (4200, -100),
          'UC': (1000, 3500), 'GLX': (2500, 2750), 'MM': (4400, 2600)}
FLIP = {}

# formate: SFGTHi/AMPMS2 product = FTHFLi substrate -> links formaldehyde oxidation and the THF C1 pathway
MERGES = [('151', ['227'])]

N_HPB = 13
def ring(k, R=850):            # 3HP/4HB cycle: 13 nodes, clockwise, acetyl-CoA on top
    a = math.radians(-90+360*k/N_HPB); return (R*math.cos(a), R*math.sin(a))
def polar(deg, r): a = math.radians(deg); return (r*math.cos(a), r*math.sin(a))

def layout(L):
    # ------------------------------------------------ formaldehyde oxidation: vertical chain, AMPMS2 fan on the left
    L.setp('C1', {'125': (0, 0), '127': (0, 300), '138': (0, 600), '151': (0, 900),
                  '418': (-550, 600), '421': (-550, 1250)})
    L.reaction('136', 'C1', side='L', dist=120)                    # FALGTHLs
    L.reaction('149', 'C1', side='L', dist=120)                    # FALDH2
    L.reaction('162', 'C1', side='L', dist=120)                    # SFGTHi
    L.reaction('434', 'C1', axis=(0, 1), mid=(-550, 900), side='L')  # AMPMS2: air -> 4ampm + formate
    # ------------------------------------------------ THF C1 pathway (formate -> 10-formyl-THF <- methenyl <- methylene)
    L.setp('THF', {'228': (450, 900), '229': (225, 1250), '245': (225, 1600), '255': (225, 1950),
                   '558': (950, 1250), '559': (950, 900)})
    L.reaction('243', 'THF', axis=(0, 1), mid=(225, 1075), side='L')   # FTHFLi: for + thf -> 10fthf
    L.reaction('569', 'THF', axis=(0, -1), mid=(700, 1075))            # ULA4NFT hourglass (waist pinned)
    L.reaction('254', 'THF', side='R')                                 # MTHFC (written methf -> 10fthf)
    L.reaction('267', 'THF', side='R')                                 # MTHFD (written mlthf -> methf)
    # ------------------------------------------------ ubiquinone biosynthesis: chain at x=250, q8 end at x=1000
    L.setp('UB', {'306': (0, 0), '294': (0, 300), '268': (500, 0), '270': (500, 300),
                  '282': (250, 650), '283': (250, 950), '319': (250, 1450), '335': (250, 1750), '523': (250, 2050),
                  '484': (1000, 2050), '347': (1000, 1550), '350': (1000, 1250)})
    L.reaction('316', 'UB', side='L')                              # CHRPL
    L.reaction('281', 'UB', side='R')                              # OCTDPS
    L.reaction('305', 'UB', side='L')                              # HBZOPT (4hbz + octdp -> 3ophb)
    L.reaction('293', 'UB', side='L')                              # OPHBDC
    L.reaction('333', 'UB', shift=(-130, 0), side='L')             # OPHHXy lane
    L.reaction('557', 'UB', shift=(130, 0), side='R')              # OPHHX3 lane
    L.reaction('346', 'UB', side='L')                              # OHPHM
    L.reaction('537', 'UB', side='L')                              # OMPHHXy
    L.reaction('520', 'UB', shift=(-130, 0), side='L')             # OMMBLHXy lane
    L.reaction('504', 'UB', shift=(130, 0), side='R')              # OMMBLHX3 lane
    L.reaction('360', 'UB', side='R')                              # DMQMT
    # ------------------------------------------------ 3HP/4HB cycle as a ring (missing steps leave gaps)
    ids = ['83', '111', '94', '96', '26', '28', '2', '1', '12', '13', '40', '39', '60']
    L.setp('HPB', {n: ring(k) for k, n in enumerate(ids)})
    L.setp('HPB', {'109': polar(-90+360*0.5/N_HPB, 1350)})        # hco3 outside, next to ACCOAC
    C = (0, 0)
    L.reaction('124', 'HPB', arc=(C, 0), ends=('83', '111'), side='in')   # ACCOAC
    L.reaction('107', 'HPB', arc=(C, 0), side='out')                      # MSAR
    L.reaction('38', 'HPB', arc=(C, 0), side='out')                       # PCNO (written ppcoa -> prpncoa)
    L.reaction('10', 'HPB', arc=(C, 0), side='out')                       # MMM (written succoa -> mmcoa)
    L.reaction('24', 'HPB', arc=(C, 0), side='out')                       # GHBDHx
    L.reaction('93', 'HPB', arc=(C, 0), side='out')                       # ACACT1r (accoa -> aacoa)
    for rid, k in (('50', -1), ('59', 0)):                                # ECOAH1 x2: concentric lanes
        L.reaction(rid, 'HPB', arc=(C, 180*k), side='in' if k < 0 else 'out', dist=80 if k < 0 else 110)
    for rid, k in (('72', -1), ('82', 0)):                                # HACD1 x2
        L.reaction(rid, 'HPB', arc=(C, 180*k), side='in' if k < 0 else 'out', dist=80 if k < 0 else 110)
    # ------------------------------------------------ glycolate / glyoxylate
    L.setp('GLX', {'179': (0, 0), '177': (1000, 0), '203': (1000, 450), '435': (1000, -450),
                   '164': (0, 700), '165': (450, 700)})
    for rid, x in (('459', 330), ('471', 500), ('483', 670)):             # GLYCTO2-4 stacked lane below
        L.reaction(rid, 'GLX', axis=(1, 0), mid=(x, 0), shift=(0, 110), side='D', gap=40)
    for rid, x in (('190', 600), ('202', 400)):                           # GLYCLTDy/x stacked lane above
        L.reaction(rid, 'GLX', axis=(-1, 0), mid=(x, 0), shift=(0, -110), side='U', gap=40)
    L.reaction('213', 'GLX', side='R')                                    # GLXCL
    L.reaction('447', 'GLX', side='R')                                    # DHGLYH
    L.reaction('176', 'GLX', side='D')                                    # GLYCK (detached piece)
    # ------------------------------------------------ urea cycle (arginase not in the map: gap arg -> orn)
    L.setp('UC', {'374': (-450, 0), '375': (0, -450), '361': (450, 0), '362': (0, 450),
                  '373': (-800, -700), '404': (-800, -1100)})
    L.reaction('386', 'UC', arc=((0, 0), 0), ends=('374', '375'), side='in')   # OCBT (cbp from outside)
    L.reaction('402', 'UC', arc=((0, 0), 0), side='out')                       # ARGSS
    L.reaction('372', 'UC', arc=((0, 0), 0), side='out')                       # ARGSL
    L.reaction('417', 'UC', side='L')                                          # CBMKr
    # ------------------------------------------------ catalase
    L.setp('MM', {'214': (0, 0), '215': (400, -150), '216': (400, 150)})
    L.reaction('225', 'MM', axis=(1, 0), mid=(200, 0))

NEW_EXITS = []
EXTRA_TITLES = []
TITLE_POS = {}
CURVES = dict(t=0.6, stacked=['459', '471', '483', '190', '202'], bend_override={}, special_links={})
# Every piece is placed by hand (ORIGIN). Two 1x1 keep-out boxes mark the corners of the intended
# canvas, so the packer finds every piece already at a free spot inside the canvas and leaves it there.
# (A keep-out disc for the 3HP/4HB ring interior is not needed - no piece has to move - and it would
#  collide with the inner cofactor dots because keep-out circles are rasterised as squares.)
PACK = dict(keep_out=[('box', (-900, -700, -899, -699)), ('box', (5600, 4200, 5601, 4201))], auto_ring_keep_out=True)
RENDER_SCALE = 0.5
