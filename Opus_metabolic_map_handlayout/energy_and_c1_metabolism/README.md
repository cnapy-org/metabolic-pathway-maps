# iML1515 energy and C1 metabolism – hand layout

Input: `ecmapper/BRENDAmap/draft_map_iML1515/split/energy_and_c1_metabolism.escher.json` (BRENDA draft).
Made with `escher_handlayout` (`Opus_metabolic_map_handlayout/escher_handlayout_toolkit`). The coordinates of
every primary metabolite and every zone origin were chosen by hand in `spec.py`. Code was used only for
the routine steps: placing markers and cofactors, curves, labels, packing, QA and validation.

```bash
PYTHONPATH=../escher_handlayout_toolkit python -m escher_handlayout run spec.py --out v8 --render
```

## Files
- `v8/energy_and_c1_metabolism.escher_final.json`: the Escher map
- `v8/qa_report.md`, `v8/run_log.json`, `v8/overview.png` (matplotlib preview, not Escher)
- `spec.py`: the hand layout
- `inspect/`: Phase 0 report

## Counts
| | draft | final |
|---|---|---|
| reactions | 42 | 42 (unchanged) |
| nodes | 264 | 266 |
| text labels | 24 (7 titles, 17 exit labels) | 24 |

Node changes:
- +3 cofactor copies. The h2o node is shared by the two ECOAH1 reactions, and the nadh and nad nodes by
  the two HACD1 reactions; each reaction gets its own copy.
- −1 merge (formate).

## Subsystems (zones)
All 7 titles are the draft text labels. One was extended: "CO2 fixation in Crenarchaeota
(3-hydroxypropionate/4-hydroxybutyrate cycle)".

| zone | reactions | notes |
|---|---|---|
| formaldehyde oxidation | FALGTHLs, FALDH2, SFGTHi, AMPMS2 | vertical chain fald → hmgth → S-formylglutathione → formate; AMPMS2 fans in from the left |
| reductive acetyl-CoA pathway | FTHFLi, MTHFC, MTHFD, ULA4NFT | formate + THF → 10-formyl-THF ← methenyl-THF ← methylene-THF; ULA4NFT is an hourglass returning THF |
| ubiquinone biosynthesis | 11 reactions | chorismate + octaprenyl-PP chain; OPHHXy/OPHHX3 and OMMBLHXy/OMMBLHX3 drawn as two lanes each; the 2ommbl → q8h2 piece sits beside the chain (the missing 2ombzl → 2ommbl step leaves a gap) |
| 3HP/4HB cycle | ACCOAC, MSAR, PCNO, MMM, GHBDHx, ECOAH1×2, HACD1×2, ACACT1r | drawn as a ring of 13 nodes (R = 850, clockwise, acetyl-CoA on top); missing steps leave gaps in the ring |
| glycolate / glyoxylate | GLYCLTDx/y, GLYCTO2/3/4, GLXCL, DHGLYH, GLYCK | glycolate ↔ glyoxylate as two stacked lanes, one per direction |
| urea cycle | OCBT, ARGSS, ARGSL, CBMKr | ring orn → citr → argsuc → arg (arginase is not in the map, so arg → orn is a gap); carbamoyl-P from the upper left |
| methane metabolism | CAT | draft title kept |

## Decisions
- **MSAR moved to the 3HP/4HB cycle.** The draft proposed it for glycolate/glyoxylate (LOW CONFIDENCE).
  Malonate-semialdehyde → 3-hydroxypropionate is a step of that cycle.
- **The 3HP/4HB cycle is drawn as one ring.** In the map it is 1 connected piece (accoa → aacoa →
  3hbcoa → b2coa, plus ACCOAC) and 4 single-reaction pieces (MSAR, PCNO, MMM, GHBDHx).
  - All are placed on one circle in cycle order: accoa, malcoa, msa, 3hpp, ppcoa, prpncoa, mmcoa,
    succoa, sucsal, ghb, b2coa, 3hbcoa, aacoa.
  - Reactions are not flipped: PCNO (ppcoa → prpncoa) and MMM (succoa → mmcoa) keep their model
    direction.
  - The duplicates ECOAH1×2 and HACD1×2 use concentric lanes: one inner lane with cofactors inside, and
    the ring lane with cofactors outside.
- **Merges:** only formate (`227` → `151`). It joins formaldehyde oxidation with the THF C1 pathway.
  - No other key metabolite occurs twice as a primary node.
  - Cofactor-type primaries were not merged. For example, q8h2 is a primary in DMQMT but a cofactor in
    GLYCTO2.
- **Whole-map arrangement:**
  - Top row: formaldehyde + THF (main component, left), 3HP/4HB ring (middle), ubiquinone (right).
  - Bottom row: urea cycle, glycolate/glyoxylate (+ GLYCK), catalase.
- **Packing is off** (`PACK = dict(mode='keep')`, since v8).
  - Every detached piece got a hand-chosen position and stays exactly there; the QA report lists no
    overlapping pieces.
  - Automatic packing scattered the single-reaction pieces of the 3HP/4HB cycle all over the map.
    v7 avoided that with two 1×1 keep-out boxes at the canvas corners, which still left grid-snap
    offsets of up to 25 units.
- **Zone titles are anchored to a node of their piece.** v8 moves the reductive acetyl-CoA title below
  the THF column; in v7 it ran into the 3HP/4HB ring.
- **Curvature:** t = 0.6 (default).
- **Exit labels:** 17, all assigned by the automatic owner rule.

## Area
- Hand arrangement: 6236 × 4329 (packing off, so unchanged).
- Final canvas: **6554 × 4727** (margin 150 included).
- The map is less compact than an automatically packed one. Earlier versions with automatic packing
  reached 3880 × 3623, but they broke up the 3HP/4HB ring and placed pieces at random.

## QA (final run v8)
| check | result |
|---|---|
| structural validation | all PASS |
| crossings | 2, both inside a stack |
| dots closer than 45 | 0 |
| foreign lines within 28 of a dot | 0 |
| overlapping labels / titles on lines | 0 / 0 |
| exit labels touching a line | 0 of 17 |

## Known weak spots
- **GLYCLTDx × GLYCLTDy and GLYCTO2 × GLYCTO3 crossings.** These are stacked rows on one lane. Each row
  bows separately from glycolate/glyoxylate, so the bows of neighbouring rows cross. PROMPT.md accepts
  this.
- **The ring interior and the space left of the urea cycle are empty.** That costs area; it was kept for
  readability.
- **The "methane metabolism" title is kept from the draft.** It covers only catalase (CAT).

Validated structurally by the toolkit (`validate`: all PASS). **Not loaded in Escher.**

## Toolkit observations
- **Packing had no off-switch** (fixed: the toolkit now has `PACK = dict(mode='keep')`, used since v8).
- **Keep-out circles are rasterised as their bounding square.** A ring-interior keep-out disc therefore
  collides with inner cofactor dots near the diagonals. The same is true of `auto_ring_keep_out`.
