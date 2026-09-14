import os, re
cDir = os.path.dirname(os.path.abspath(os.sys.argv[0]))

svgdir = os.path.join(cDir, 'iJR904_Output')

font_size = re.compile('font-size=\"[\.|0-9]*\"')
new_font_size = 'font-size=\"1.0\"'

for f in os.listdir(svgdir):
    if not f.endswith('.svg'):
        continue
    print((os.path.join(svgdir, f)))
    F = file(os.path.join(svgdir, f), 'r')
    outfile = ''
    for l in F:
        rex = re.findall(font_size, l)
        if len(rex) > 0:
            for fz in rex:
                l = l.replace(fz, new_font_size)
        outfile += l
    F.close()
    F = file(os.path.join(svgdir, f), 'w')
    F.write(outfile)
    F.close()
