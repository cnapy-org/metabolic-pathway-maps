# escher_handlayout

Toolkit for turning an auto-generated draft **Escher** map into a textbook-style map whose layout is
**specified by hand**: you (or an AI agent) write the coordinates of the primary metabolites in a small
Python *spec* file; the toolkit does the mechanical rest:
- marker and cofactor placement relative to your coordinates
- Bezier curves
- label placement
- packing of detached pieces
- QA checks
- structural validation
- valid Escher JSON output

There is deliberately no automatic graph-layout algorithm in it.

Developed while laying out the iML1515 central-metabolism draft; that map is included as a worked,
regression-tested example (`examples/central_metabolism`).

## Contents

```
PROMPT.md            autonomous prompt for an AI agent (no user interaction) – start here for agents
SPEC_REFERENCE.md    spec format, L.reaction() options, cookbook, how to read the QA report
escher_handlayout/   the package
  engine.py          Layout: load map, split/merge nodes, setp(), reaction() (lines, lanes, stacks, arcs)
  inspect_map.py     Phase 0 report + spec skeleton generator, exit-label owner rule
  curves.py          Phase 3 Bezier rules (cofactor rule, S-lanes, circular arcs, special links)
  labels.py          Phase 4 label placement
  pack.py            Phase 5 packing of detached pieces (FFT occupancy search, ring keep-out)
  titles.py          title nudging, canvas
  qa.py              validation against the draft + crossings / near misses / label collisions
  render.py          optional matplotlib preview (not Escher)
  pipeline.py        runs a spec end-to-end
  __main__.py        command line interface
examples/central_metabolism/
  central_metabolism_escher.json   draft input
  spec.py                          hand layout (v7)
  reference_v7.json                expected output
tests/test_example.py              regression test
```

## Install

Python ≥ 3.9 with numpy, scipy, matplotlib.

```bash
pip install -e .            # provides the command  escher-handlayout
# or without installing:   PYTHONPATH=<this dir> python -m escher_handlayout ...
python tests/test_example.py    # should print: all tests passed
```

## Quick start: reproduce the example

```bash
escher-handlayout run examples/central_metabolism/spec.py --out out_example --render
```

This writes to `out_example/`:
- `central_metabolism_escher_final.json`, the Escher map
- `qa_report.md`, the validation and QA findings plus areas, packing moves and title nudges
- `run_log.json`, which records splits, merges, curve settings and packing
- `overview.png`, written only with `--render`

## Workflow for a new map

1. **Inspect**

   ```bash
   escher-handlayout inspect my_draft.json --out work/inspect --skeleton work/spec.py
   ```

   `work/inspect/phase0_report.md` lists:
   - subsystem proposals; entries marked LOW CONFIDENCE need a decision
   - shared cofactors, which are split automatically
   - duplicate reactions and reaction families
   - merge candidates
   - far-away primaries
   - exit-label owners
   - all reaction equations

   `work/spec.py` is a skeleton. It contains one zone per proposed subsystem, every primary node id
   (draft coordinates appear in comments) and one `L.reaction` call per reaction. **All coordinates
   are `None`** on purpose.
2. **Decide zones**: merge, rename or move reactions between zones as the biochemistry demands.
3. **Write the hand layout**: fill in `ORIGIN`, the `setp` coordinates, the reaction options and the
   splits/merges. Use the principles in `PROMPT.md` §4 and the cookbook in `SPEC_REFERENCE.md`; the
   example spec shows every pattern.
4. **Run and iterate**

   ```bash
   escher-handlayout run work/spec.py --input my_draft.json --out work/v1 [--render]
   ```

   Read `qa_report.md` and fix the spec until the acceptance targets in `PROMPT.md` §5 are met.
   Use a new `--out` folder per version.
5. **Check in Escher.** The toolkit validates structurally only; load the final JSON in Escher once
   before publishing.

To let an AI agent do steps 1–4 autonomously, give it `PROMPT.md`, this toolkit and the draft map.

## Other commands

```bash
escher-handlayout validate final.json draft.json      # structural checklist against the draft
escher-handlayout render final.json out.png --scale 0.2 [--region x0 y0 x1 y1]
```

## What the pipeline does (`run`)

1. **Load and split cofactors:** load the draft and give every reaction its own cofactor nodes.
2. **Splits and merges:** apply the spec's `prepare()` (splits), then `MERGES`.
3. **Hand layout:** call `layout()` twice, so every reaction sees final primary positions.
   Unplaced nodes or reactions abort the run with a message.
4. **Titles and exit labels:**
   - Zone titles reuse the draft labels given as `title_id`, or are created.
   - `EXTRA_TITLES` are added.
   - Unused draft titles are removed.
   - `NEW_EXITS` are added.
5. **Curves:** t = 0.6 by default. Cofactor dots move along their own reaction axis only.
6. **Labels:** reaction labels first, then cofactor, primary and exit labels. Each label picks the
   best of 8 directions × 4–5 radii against line samples, nodes, titles and earlier labels.
7. **Packing:** every piece not connected to the main map is translated to the free position that
   gives the smallest overall rectangle (closest to its old position among near-equal ones).
   Ring interiors are kept free. Titles anchored to a piece travel with it.
8. **Titles, canvas and output:** nudge titles that touch something, compute the canvas (+150 margin),
   write the JSON, validate and write the QA report.

## Limitations

- The preview renderer is not Escher. Text widths are estimates (0.62 × font size, 0.64 for bold),
  so keep clearance generous.
- QA checks are heuristics on sampled curves. Arrowheads are not modelled.
- The pipeline keeps the reaction set, stoichiometry, reversibility and genes byte-identical
  (`validate` checks this). Metabolites that the draft does not draw (usually protons) stay undrawn.
