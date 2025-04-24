"""Prints average viral percentage."""

import re
from collections import defaultdict
from glob import glob

from matplotlib import pyplot as plt
from myplot import ctx

from violin import violin


def read_file(f: str):
    rgx = re.compile('\\S+')
    root, vir = None, None
    for line in open(f):
        parts = rgx.findall(line)
        if parts[3] == 'R':
            root = int(parts[1])
        if parts[3] == 'D' and parts[5] == 'Viruses':
            vir = int(parts[1])
        if root and vir:
            break
    return vir / root


def main():
    print('Reading files')
    files = glob('../../../Data/ww2-kraken/*.krk.txt')
    rgx = re.compile('_INF_|_SOL_')
    names = {'_INF_': 'Influent', '_SOL_': 'Solid'}
    d = defaultdict(list)
    for f in files:
        g = rgx.findall(f)
        if not g:
            continue
        d[names[g[0]]].append(read_file(f))

    print([(k, len(v)) for k, v in d.items()])

    with ctx('virperc', sizeratio=0.75):
        violin(d)
        plt.ylabel('Viral read ratio')
    with ctx('virperc-log', sizeratio=0.75):
        violin(d)
        plt.ylabel('Viral read ratio')
        plt.yscale('log')


if __name__ == '__main__':
    main()
