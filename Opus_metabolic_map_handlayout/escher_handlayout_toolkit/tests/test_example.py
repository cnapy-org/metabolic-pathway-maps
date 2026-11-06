"""Regression test: the example spec reproduces the reference map (v7) and passes validation.
Run:  python -m pytest tests   (or:  python tests/test_example.py)"""
import json, math, os, sys, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from escher_handlayout.pipeline import run

EX = os.path.join(ROOT, 'examples', 'central_metabolism')

def test_example_reproduces_reference():
    out = tempfile.mkdtemp()
    L, s = run(os.path.join(EX, 'spec.py'), None, out)
    assert s['validation_ok']
    assert s['label_overlaps'] == 0 and s['titles_on_lines'] == 0
    got = json.load(open(s['output']))[1]; ref = json.load(open(os.path.join(EX, 'reference_v7.json')))[1]
    assert set(got['nodes']) == set(ref['nodes'])
    worst = max(math.dist((got['nodes'][k]['x'], got['nodes'][k]['y']), (ref['nodes'][k]['x'], ref['nodes'][k]['y'])) for k in ref['nodes'])
    assert worst < 0.05, worst
    for rid, r in ref['reactions'].items():
        for sid, seg in r['segments'].items():
            g = got['reactions'][rid]['segments'][sid]
            assert (g['b1'] is None) == (seg['b1'] is None)

def test_inspect_and_skeleton():
    from escher_handlayout.inspect_map import analyse, skeleton
    import ast
    L, rep = analyse(os.path.join(EX, 'central_metabolism_escher.json'))
    assert rep['counts']['reactions'] == 91 and len(rep['exit_label_owner']) == 29
    ast.parse(skeleton(L, rep, 'central_metabolism_escher.json'))

def test_pack_mode_keep():
    """PACK mode='keep': detached pieces stay where layout() put them (no moves)"""
    out = tempfile.mkdtemp(); spec = os.path.join(out, 'spec_keep.py'); ex = os.path.join(EX, 'spec.py')
    open(spec, 'w').write(f"exec(open({ex!r}).read())\nINPUT = {os.path.join(EX, 'central_metabolism_escher.json')!r}\n"
                          "PACK = dict(mode='keep')\n")
    L, s = run(spec, None, out)
    log = json.load(open(os.path.join(out, 'run_log.json')))
    assert s['validation_ok'] and log['packing']['moves'] == [] and 'piece_overlaps' in s

def _wrapped(extra):
    out = tempfile.mkdtemp(); spec = os.path.join(out, 'spec_x.py'); ex = os.path.join(EX, 'spec.py')
    open(spec, 'w').write(f"exec(open({ex!r}).read())\nINPUT = {os.path.join(EX, 'central_metabolism_escher.json')!r}\n" + extra)
    return spec, out

def test_merged_away_anchor_and_unknown_id():
    """a merged-away id ('923' -> glucose '480') is followed; an unknown id gives a clear error"""
    spec, out = _wrapped("EXTRA_TITLES = EXTRA_TITLES + [dict(text='glucose here', anchor='923', offset=(0, -200))]\n")
    L, s = run(spec, None, out)
    assert s['validation_ok'] and any(t['text'] == 'glucose here' for t in L.T.values())
    spec, out = _wrapped("EXTRA_TITLES = [dict(text='x', anchor='no_such_node')]\n")
    try: run(spec, None, out); assert False, 'expected KeyError'
    except KeyError as e: assert 'no_such_node' in str(e) and 'does not exist' in str(e)

def test_duplicate_exits_and_drop_labels():
    """identical exit labels on one node after a merge are reduced to one; DROP_LABELS removes labels"""
    from escher_handlayout.engine import Layout
    L = Layout(os.path.join(EX, 'central_metabolism_escher.json'))
    k1, k2 = [k for k, t in L.T.items() if t['text'][:1] == '→'][:2]
    L.T[k2]['text'] = L.T[k1]['text']; L.exit_owner = {k1: '480', k2: '923'}
    L.merge('480', '923')
    dropped = L.drop_duplicate_exits()
    assert len(dropped) == 1 and len(L.exit_owner) == 1
    gone, kept, _ = dropped[0]
    assert {gone, kept} == {k1, k2} and kept in L.T and gone not in L.T
    spec, out = _wrapped(f"DROP_LABELS = ['1155']\n")
    L, s = run(spec, None, out)
    assert '1155' not in L.T and s['validation_ok']

def test_keep_out_disc_is_round():
    """a keep-out disc occupies the cells it touches, not its bounding square"""
    import numpy as np
    from escher_handlayout.pack import raster, CELL
    G = raster([[], [(400, 400, 300, 'disc')], []], 0, 0, 20, 20)
    assert G[10, 10] and not G[3, 3] and not G[16, 16]          # centre in, corners of the square out
    S = raster([[], [(400, 400, 300)], []], 0, 0, 20, 20)       # node dot: bounding square (unchanged)
    assert S[3, 3]

if __name__ == '__main__':
    test_example_reproduces_reference(); test_inspect_and_skeleton(); test_pack_mode_keep(); test_merged_away_anchor_and_unknown_id(); test_duplicate_exits_and_drop_labels(); test_keep_out_disc_is_round(); print('all tests passed')
