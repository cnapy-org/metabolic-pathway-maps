"""SPEC - iML1515 'fermentation and degradation' (BRENDA draft, 49 reactions) -> hand-laid-out Escher map.

Every coordinate of a PRIMARY metabolite below was chosen by hand (zone-local coordinates, zone origin in
ORIGIN). Escher convention: x to the right, y DOWN. Spacing between chain neighbours ~300.
The map consists of ten independent small pathways; only propionate fermentation (with propanol
degradation's PTA2 attached via propionyl-CoA) is larger. Detached pieces are packed automatically."""

INPUT = '../../ecmapper/BRENDAmap/draft_map_iML1515/split/fermentation_and_degradation.escher.json'

# ------------------------------------------------------------------ zones and titles
# All titles are anchored to a node of their piece so that they travel with it during packing.
ZONES = {
 'PF':  dict(title='propionate fermentation', title_id='671', anchor=('477', (-80, -330)), color='#d62728'),
 'PD':  dict(title='3-phenylpropionate degradation', title_id='663', anchor=('82', (-60, -330)), color='#1f77b4'),
 'PDA': dict(title='phenylacetate degradation (aerobic)', title_id='669', anchor=('350', (-60, -330)), color='#2ca02c'),
 'AD':  dict(title='allantoin degradation', title_id='665', anchor=('153', (-60, -330)), color='#9467bd'),
 'MD':  dict(title='methylglyoxal degradation', title_id='667', anchor=('256', (-60, -150)), color='#ff7f0e'),
 'LLD': dict(title='L-lactaldehyde degradation', title_id='664', anchor=('124', (-60, -150)), color='#8c564b'),
 'OO':  dict(title='octane oxidation', title_id='668', anchor=('305', (-60, -150)), color='#e377c2'),
 'HBZ': dict(title='3-hydroxybenzoyl-CoA hydrolysis', anchor=('108', (-60, -230)), color='#7f7f7f'),
 'CD':  dict(title='creatinine degradation', title_id='666', anchor=('239', (-60, -230)), color='#bcbd22'),
 'PRD': dict(title='propanol degradation', title_id='670', anchor=('430', (-60, -230)), color='#17becf'),
}
TITLE_POS = {}
EXTRA_TITLES = []
EXIT_OWNER = {}          # automatic owner rule is correct for all 14 exit labels

# ------------------------------------------------------------------ splits / merges
def prepare(L):
    S = L.S
    S['ppcoa_P'] = L.split('477', '488')     # propionyl-CoA copy for PPCSCT (bottom hourglass),
                                             # avoids a long line from the top hub

MERGES = [('169', ['160']),    # allantoin: URIC product (drawn as secondary) joins the allantoin chain
          ('468', ['547']),    # (S)-methylmalonyl-CoA: both MMCD duplicates start from one node
          ('477', ['417'])]    # propionyl-CoA: PTA2 (propanol degradation) attached to the PF hub

ORIGIN = {'PF': (0, 0), 'PRD': (-1200, 1500), 'PD': (2200, -200), 'PDA': (2200, 900), 'AD': (2200, 2100),
          'MD': (-1200, -100), 'LLD': (-800, -100), 'OO': (4800, 900), 'HBZ': (4800, 2000),
          'CD': (4800, 2500)}
FLIP = {}


def layout(L):
    S = L.S
    # ---------------------------------------------------------------- PF: propionate fermentation
    # propionyl-CoA hub top left; methylcitrate pathway goes right and then down to succinate;
    # the "sleeping beauty" loop succinyl-CoA -> methylmalonyl-CoA -> propionyl-CoA is the left column,
    # closed at the bottom by the PPCSCT hourglass (succ + ppcoa -> succoa + propionate).
    L.setp('PF', {'477': (0, 0), '559': (325, -400),                        # ppcoa, oaa (from above)
                  '526': (700, 0), '456': (700, 300), '444': (700, 600),     # 2mcit, 2mcacn, micit
                  '446': (700, 900), '489': (1300, 900),                     # succ, fum
                  '468': (0, 450), '467': (0, 900),                          # mmcoa, succoa
                  '478': (0, 1400), S['ppcoa_P']: (700, 1400),               # propionate, ppcoa copy
                  '419': (-450, 0)})                                         # propionyl phosphate (PTA2)
    L.reaction('571', 'PF', axis=(1, 0), mid=(360, 0), side='D')            # MCITS (oaa enters from above)
    L.reaction('536', 'PF', side='R')                                        # MCITD
    L.reaction('466', 'PF', side='R')                                        # MICITDr
    L.reaction('455', 'PF', side='L')                                        # MCITL2 (pyr dot inside the loop)
    L.reaction('525', 'PF', shift=(0, -130), side='U', dist=110)             # SUCDi  succ -> fum (upper lane)
    L.reaction('501', 'PF', shift=(0, 130), side='D', dist=85)               # FRD2 fum -> succ (lane 2)
    L.reaction('513', 'PF', shift=(0, 300), side='D', dist=100)              # FRD3 fum -> succ (lane 3, outer)
    L.reaction('476', 'PF', side='L')                                        # MMM  succoa -> mmcoa
    L.reaction('546', 'PF', shift=(-80, 0), side='L', dist=110)              # MMCD lane 1
    L.reaction('557', 'PF', shift=(80, 0), side='R', dist=110)               # MMCD lane 2 (duplicate)
    L.reaction('488', 'PF', axis=(-1, 0), mid=(350, 1150))                   # PPCSCT hourglass
    L.reaction('429', 'PF', side='U')                                        # PTA2 ppcoa -> ppap
    # ---------------------------------------------------------------- PRD: propanol degradation (ACACCT)
    L.setp('PRD', {'430': (0, 0), '432': (450, 0)})
    L.reaction('443', 'PRD', side='D')                                       # ACACCT
    # ---------------------------------------------------------------- PD: 3-phenylpropionate degradation
    # main chain left -> right, two entries (pppn via cechddd, 3hpppn from above) meet at dhpppn;
    # hkntd joins op4en from below (fan op4en + fum)
    L.setp('PD', {'82': (0, 0), '83': (300, 0), '4': (600, 0), '1': (600, -550),
                  '20': (900, 0), '31': (1200, 0), '44': (1500, 0), '54': (1800, 0), '67': (2100, 0),
                  '586': (1200, 500), '587': (1500, 300)})
    L.reaction('95', 'PD', side='U')                                          # PPPNDO
    L.reaction('107', 'PD', side='D')                                         # DHPPD
    L.reaction('18', 'PD', side='R', axis=(0, 1), mid=(600, -330))                                      # 3HPPPNH (vertical)
    L.reaction('29', 'PD', side='D')                                          # HPPPNDO
    L.reaction('42', 'PD', side='U')                                          # HKNDDH (succ dot)
    L.reaction('53', 'PD', side='U')                                          # OP4ENH
    L.reaction('64', 'PD', side='U')                                          # HOPNTAL (pyr dot)
    L.reaction('79', 'PD', side='U')                                          # ACALD
    L.reaction('597', 'PD', axis=(0, -1), mid=(1200, 250), side='L')          # HKNTDH fan
    # ---------------------------------------------------------------- PDA: phenylacetate degradation
    # snake: upper row left -> right, turn down, lower row right -> left ending at succinyl-CoA
    L.setp('PDA', {'350': (0, 0), '333': (450, 0), '321': (800, 0), '322': (1100, 0), '365': (1450, 0),
                   '367': (1450, 450), '380': (1100, 450), '392': (750, 450), '405': (400, 450)})
    L.reaction('364', 'PDA', shift=(0, -75), side='U', dist=110)             # PACCOAL pac -> phaccoa
    L.reaction('648', 'PDA', shift=(0, 75), side='D', dist=110)              # PACOAT  phaccoa -> pac
    L.reaction('347', 'PDA', side='U')                                        # PACCOAE
    L.reaction('330', 'PDA', side='U')                                        # REPHACCOAI
    L.reaction('636', 'PDA', side='U')                                        # OXCOAHDH
    L.reaction('378', 'PDA', side='R')                                        # OXDHCOAT (accoa dot)
    L.reaction('389', 'PDA', side='D')                                        # DHACOAH
    L.reaction('402', 'PDA', side='D')                                        # HADPCOADH3
    L.reaction('415', 'PDA', side='D')                                        # 3OXCOAT
    # ---------------------------------------------------------------- AD: allantoin degradation
    L.setp('AD', {'153': (0, 0), '169': (300, 0), '171': (600, 0), '184': (900, 0), '198': (1200, 0),
                  '210': (1050, -450), '211': (1350, -450),                  # co2, glx (fan up)
                  '226': (1550, 0), '611': (1900, -220), '612': (1900, 220)})   # oxur, cbp, oxam
    for rid in ('168', '180', '195', '208'):                                  # URIC, ALLTN, ALLTAMH2, UGCIAMH
        L.reaction(rid, 'AD', side='D')
    L.reaction('223', 'AD', axis=(0, -1), mid=(1200, -200), side='L')         # UGLYCH fan
    L.reaction('236', 'AD', side='D')                                         # URDGLYCD
    L.reaction('622', 'AD', axis=(1, 0), mid=(1700, 0), side='D', dist=110)   # OXAMTC fan
    # ---------------------------------------------------------------- MD: methylglyoxal degradation
    # direct route (GLYOX3) straight down, glutathione route (LGTHL, GLYOX) via lgt__S on the left
    L.setp('MD', {'256': (0, 0), '268': (-400, 300), '257': (0, 600), '290': (0, 900)})
    L.reaction('266', 'MD', side='R')                                         # GLYOX3
    L.reaction('277', 'MD', side='L')                                         # LGTHL
    L.reaction('289', 'MD', side='L')                                         # GLYOX
    L.reaction('302', 'MD', side='R')                                         # LDH_D2
    # ---------------------------------------------------------------- LLD: L-lactaldehyde degradation
    L.setp('LLD', {'124': (0, 0), '122': (0, 400), '138': (0, 750)})
    L.reaction('135', 'LLD', shift=(-80, 0), side='L', dist=110)              # LCARS   (NADH lane)
    L.reaction('609', 'LLD', shift=(80, 0), side='R', dist=110)               # LCARSyi (NADPH lane)
    L.reaction('150', 'LLD', side='R')                                        # LCADi
    # ---------------------------------------------------------------- OO: octane oxidation
    L.setp('OO', {'305': (0, 0), '307': (0, 400), '574': (450, 400)})
    L.reaction('320', 'OO', side='L', step=70)                                # FACOAL80t2pp (uptake + activation)
    L.reaction('584', 'OO', side='D')                                         # FACOAE80
    # ---------------------------------------------------------------- HBZ / CD: single reactions
    L.setp('HBZ', {'108': (0, 0), '110': (450, 0)})
    L.reaction('121', 'HBZ', side='D')                                        # 3HBZCT
    L.setp('CD', {'239': (0, 0), '241': (450, 0)})
    L.reaction('254', 'CD', side='D')                                         # SARCOX


NEW_EXITS = []
CURVES = dict(t=0.6, stacked=[], bend_override={}, special_links={})
PACK = dict(keep_out=[], auto_ring_keep_out=True)
RENDER_SCALE = 0.5
