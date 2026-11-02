"""CLI:
  python -m escher_handlayout inspect DRAFT.json --out DIR [--skeleton spec.py]
  python -m escher_handlayout run SPEC.py [--input DRAFT.json] --out DIR [--render] [--name NAME]
  python -m escher_handlayout validate FINAL.json DRAFT.json
  python -m escher_handlayout render FINAL.json OUT.png [--region x0 y0 x1 y1] [--scale S]"""
import argparse, json, os, sys

def main(argv=None):
    ap = argparse.ArgumentParser(prog='escher_handlayout'); sub = ap.add_subparsers(dest='cmd', required=True)
    a = sub.add_parser('inspect'); a.add_argument('draft'); a.add_argument('--out', default='inspect'); a.add_argument('--skeleton')
    b = sub.add_parser('run'); b.add_argument('spec'); b.add_argument('--input'); b.add_argument('--out', default='out')
    b.add_argument('--render', action='store_true'); b.add_argument('--name')
    c = sub.add_parser('validate'); c.add_argument('final'); c.add_argument('draft')
    d = sub.add_parser('render'); d.add_argument('map'); d.add_argument('png'); d.add_argument('--region', nargs=4, type=float); d.add_argument('--scale', type=float, default=0.5)
    args = ap.parse_args(argv)
    if args.cmd == 'inspect':
        from .inspect_map import analyse, markdown, skeleton
        os.makedirs(args.out, exist_ok=True); L, rep = analyse(args.draft)
        open(os.path.join(args.out, 'phase0_report.md'), 'w').write(markdown(rep))
        json.dump(rep, open(os.path.join(args.out, 'phase0_report.json'), 'w'), indent=1, default=str)
        if args.skeleton: open(args.skeleton, 'w').write(skeleton(L, rep, os.path.basename(args.draft)))
        print('report written to', args.out, '(skeleton: %s)' % args.skeleton if args.skeleton else '')
    elif args.cmd == 'run':
        from .pipeline import run
        try:
            L, s = run(args.spec, args.input, args.out, args.render, args.name)
        except (ValueError, KeyError, RuntimeError) as e:
            print('SPEC ERROR:', e, file=sys.stderr); sys.exit(2)
        print(json.dumps(s, indent=1))
        sys.exit(0 if s['validation_ok'] else 1)
    elif args.cmd == 'validate':
        from .qa import validate
        res = validate(args.final, args.draft)
        for n, ok, d in res: print('PASS' if ok else 'FAIL', n, d or '')
        sys.exit(0 if all(ok for _, ok, _ in res) else 1)
    elif args.cmd == 'render':
        from .engine import Layout
        from .render import render
        L = Layout(args.map)
        for rid in L.R: L.zone_of[rid] = 'map'
        render(L, args.png, region=args.region, scale=args.scale)

if __name__ == '__main__': main()
