"""Displays read count distribution."""

import json

from matplotlib import pyplot as plt
from myplot import ctx

from violin import violin


def read_counts():
    lines = open('readcounts.txt')
    items = (x.split(' ') for x in lines)
    items = ((a, int(int(b) / 4)) for a, b in items)
    items = list(items)
    assert all(len(x) == 2 for x in items)
    return items


def derive(xx: list) -> list:
    return [1 / (x2 - x1) for x1, x2 in zip(xx, xx[1:])]


rc = read_counts()
json.dump({a: b for a, b in rc}, open('readcounts.json', 'w'), indent=2)

inf = [x for x in rc if '_INF_' in x[0]]
sol = [x for x in rc if '_SOL_' in x[0]]
oth = [x for x in rc if '_INF_' not in x[0] and '_SOL_' not in x[0]]
assert len(inf) > 0
assert len(sol) > 0
assert len(oth) > 0
assert len(set(inf) | set(sol) | set(oth)) == len(rc)
assert len(inf) + len(sol) + len(oth) == len(rc)

vinf = sorted(x[1] / 1000000 for x in inf)
vsol = sorted(x[1] / 1000000 for x in sol)
voth = sorted(x[1] / 1000000 for x in oth)


with ctx('readcounts', sizeratio=0.75):
    d = {'Influent': vinf, 'Solid': vsol, 'Other': voth}
    violin(d)
    plt.ylabel('Read count (millions)')


plt.style.use('bmh')
with ctx('readcounts_cdf'):
    yy = [i / (len(vinf) - 1) for i in range(len(vinf))]
    plt.plot(vinf, yy, label='Influent')
    yy = [i / (len(vsol) - 1) for i in range(len(vsol))]
    plt.plot(vsol, yy, label='Solid')
    yy = [i / (len(voth) - 1) for i in range(len(voth))]
    plt.plot(voth, yy, label='Other')

    # plt.plot(vinf[:-1], derive(vinf), label='Influent')
    # plt.plot(vsol[:-1], derive(vsol), label='Solid')
    # plt.plot(voth[:-1], derive(voth), label='Other')

    plt.title('Read count CDF')
    plt.xlabel('Read count (millions)')
    plt.ylabel('Ratio')
    plt.legend()
