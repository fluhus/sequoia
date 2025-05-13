"""Plots rarefaction curves."""

import json

from matplotlib import pyplot as plt
from myplot import ctx

from samplenaming import LUNA_GROUPS, fix_name


def main():
    data: dict[str, list[list[int]]] = json.load(open('rrf.json'))
    data = {
        fix_name(k): v for k, v in data.items() if not k.startswith('Undetermined')
    }
    print(len(data), 'samples')

    plt.style.use('bmh')
    with ctx('rrf', dpi=400, sizeratio=1):
        for g in LUNA_GROUPS:
            first = True
            for k, (x, y) in data.items():
                if not k.startswith(g):
                    continue
                if first:
                    plt.plot(x, y, label=LUNA_GROUPS[g], lw=0.5)
                    c = plt.gca().lines[-1].get_color()
                    plt.plot(x[-1], y[-1], '*', c=c)
                    first = False
                else:
                    c = plt.gca().lines[-1].get_color()
                    plt.plot(x, y, c=c, lw=0.5)
                    plt.plot(x[-1], y[-1], '*', c=c)
        plt.xlabel('Sequencing depth (reads)')
        plt.ylabel('Number of OTUs')
        # plt.xscale('log')

        for line in plt.legend().get_lines():
            line.set_linewidth(3.0)


if __name__ == '__main__':
    main()
