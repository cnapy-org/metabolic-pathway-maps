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

if __name__ == '__main__':
    test_example_reproduces_reference(); test_inspect_and_skeleton(); test_pack_mode_keep(); print('all tests passed')
