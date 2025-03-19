"""rRNA-related plots."""

from glob import glob

import numpy as np
from matplotlib import pyplot as plt
from myplot import ctx
from scipy.stats import mannwhitneyu

from config import DATA_DIR
from samplenaming import LUNA_GROUPS, fix_name

NREADS_FILES = f'{DATA_DIR}/ww-greengenes/*.nreads'


def load_nreads_single(f) -> tuple[int, int]:
    lines = list(open(f))
    # assert len(lines) == 2, f'found {len(lines)} lines, want 2'
    if len(lines) != 2:
        return None
    return int(lines[0]), int(lines[1])


def load_nreads() -> dict[str, tuple[int, int]]:
    return {
        fix_name(f): nr for f in glob(NREADS_FILES) if (nr := load_nreads_single(f))
    }


def plot_bars_single(d: dict[str, float], ymax=None):
    if not d:
        return
    plt.bar(list(d.keys()), list(d.values()))
    plt.xticks(rotation=45, ha='right')
    plt.ylabel('Ratio of reads mapped to GreenGenes')
    if ymax:
        plt.ylim(0, ymax)


def plot_bars(data: dict[str, tuple[int, int]]):
    data = {k: v[1] / v[0] for k, v in data.items()}
    mx = max(data.values()) * 1.1
    with ctx('rrna', sizeratio=2):
        for i, g in enumerate(LUNA_GROUPS):
            plt.subplot(221 + i)
            gdata = {k: v for k, v in data.items() if k.startswith(g)}
            plot_bars_single(gdata, ymax=mx)


def violin(d: dict[str, list]):
    plt.violinplot(list(d.values()), showmedians=True, showextrema=False)
    # plt.violinplot(
    #     list(d.values()),
    #     quantiles=[[0.25, 0.5, 0.75]] * len(d),
    #     showextrema=False,
    # )
    plt.xticks(list(range(1, len(d) + 1)), d.keys())
    for i, yy in enumerate(d.values()):
        xx = np.array([float(i + 1)] * len(yy))
        xx += np.random.normal(scale=0.03, size=len(xx))
        plt.plot(xx, yy, 'ko', markersize=2)


def plot_violins(data: dict[str, tuple[int, int]]):
    data = {k: v[1] / v[0] for k, v in data.items()}
    d = {}
    for g in LUNA_GROUPS:
        gdata = {k: v for k, v in data.items() if k.startswith(g)}
        if not gdata:
            continue
        d[LUNA_GROUPS[g]] = list(gdata.values())

    p1 = mannwhitneyu(d['v1 solid'], d['v1 influent']).pvalue
    p2 = mannwhitneyu(d['v2 solid'], d['v2 influent']).pvalue
    print(f'Mann-Whitney 1 p={p1:.2f}')
    print(f'Mann-Whitney 2 p={p2:.2f}')

    mx = max(x for v in d.values() for x in v) * 1.1

    # plt.style.use('bmh')
    with ctx('rrna_violin', sizeratio=0.75, dpi=400):
        violin(d)
        plt.ylabel('Ratio of reads mapped to GreenGenes')

        # Draw p-value bars.
        plt.plot([1, 1, 2, 2], [mx, mx * 1.05, mx * 1.05, mx], 'k')
        plt.text(1.5, mx * 1.07, f'p={p1:.2f}', horizontalalignment='center')
        plt.plot([3, 3, 4, 4], [mx, mx * 1.05, mx * 1.05, mx], 'k')
        plt.text(3.5, mx * 1.07, f'p={p2:.2f}', horizontalalignment='center')

        # Extend Y-limit because it doesn't automatically take the text in
        # consideration.
        ylim = list(plt.ylim())
        ylim[1] *= 1.05
        plt.ylim(ylim)


def main():
    data = load_nreads()
    # data = {k: v for k, v in list(data.items())[:9]}
    # plt.style.use('bmh')
    # plot_bars(data)
    plot_violins(data)


if __name__ == '__main__':
    main()
