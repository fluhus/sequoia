"""Abundance plots."""

import colorsys
import itertools
import re
from typing import Iterable

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from myplot import ctx

from abundance import AbundancePaths, load_data
from config import DATA_DIR, WS_DATA_DIR
from samplenaming import LUNA_GROUPS, fix_name2

# Run parameters.
ALT_TOP = None
WITH_CMAP = True

# Constant data.
REFSEQ_FILE = f'{DATA_DIR}/refseq/viral.1.genomic.json'
KRAKEN_SPECIES_FILE = f'{WS_DATA_DIR}/species_krk.json'


def colors(name: str, n=None, reverse=False):
    cmap = plt.get_cmap(name)
    if n is None:
        return (cmap(i) for i in itertools.count())
    if reverse:
        return (cmap(1 - i / (n - 1)) for i in range(n))
    return (cmap(i / (n - 1)) for i in range(n))


def to_pastel(c):
    h, ll, s = colorsys.rgb_to_hls(*c[:3])
    s = 0.6
    ll = 0.65  # 1 - ((1 - l) / 2)
    return colorsys.hls_to_rgb(h, ll, s)


def nsubplots(n) -> np.ndarray:
    a = int(np.ceil(n**0.5))
    b = int(np.ceil(n / a))
    return np.array([a, b])


def sorted_by_time(a) -> list:
    no_prefix = re.compile('.*?_(.*)')
    return sorted(a, key=lambda x: no_prefix.findall(x)[0])


def repl_strings(s: str, d: dict) -> str:
    for k, v in d.items():
        s = s.replace(k, v)
    return s


def fix_locations(s: Iterable[str]) -> list[str]:
    d = {
        '_LB_': 'Site 1, ',
        '_Tur_': 'Site 2, ',
        '_Wod_': 'Site 3, ',
    }
    return [repl_strings(x, d) for x in s]


def fix_dates(s: Iterable[str]) -> list[str]:
    rgx = re.compile('(\\d\\d)(\\d\\d)(\\d\\d)')
    return [rgx.sub('\\1-\\2-\\3', x) for x in s]


def jason_plot(df: pd.DataFrame, groups=None):
    df = df.loc[sorted(df.index.tolist())]
    # df.index = fix_locations(df.index.tolist())
    # df.index = fix_dates(df.index.tolist())
    sums = df.sum().sort_values(ascending=False)
    tops = sums.index.tolist()[:10]
    if ALT_TOP:
        if type(ALT_TOP) is int:
            tops = sums.index.tolist()[:ALT_TOP]
        elif type(ALT_TOP) is list:
            tops = ALT_TOP
        else:
            raise TypeError(f'bad type: {type(ALT_TOP)=}')
    # print('Tops:', tops)
    if 'Other' in tops:
        tops.remove('Other')
    plt.style.use('ggplot')

    if not groups:
        groups = ('',)

    nsp = nsubplots(len(groups))
    sp = nsp[1] * 100 + nsp[0] * 10 + 1

    with ctx('jason', dpi=500, sizeratio=(nsp * 0.5).tolist()):
        for i, g in enumerate(groups):
            c = colors('plasma', len(tops))
            # c = colors('tab10')
            # c = (to_pastel(x) for x in c)
            plt.subplot(sp + i)
            df1 = df.loc[[x for x in df.index.tolist() if x.startswith(g)]]
            bottom = np.zeros(len(df1))
            xx = df1.index.tolist()
            xx = [fix_name2(x) for x in xx]
            for s in tops:
                yy = df1[s].values
                if WITH_CMAP:
                    plt.bar(xx, yy, 0.95, bottom=bottom, label=s, color=next(c))
                else:
                    plt.bar(xx, yy, 0.95, bottom=bottom, label=s)
                bottom += yy
            plt.bar(
                xx,
                1 - bottom,
                0.95,
                bottom=bottom,
                label='Other',
                color='lightgrey',
            )
            plt.xticks(rotation=45, ha='right')
            plt.ylabel('Relative abundance')

        legend_handles = plt.subplot(sp + 0).get_legend_handles_labels()

    with ctx('jason_legend', dpi=300):
        plt.legend(*legend_handles, loc='center')
        plt.gca().axis('off')


def main():
    df = load_data(AbundancePaths.BY_NAME, remove_spike=True)
    jason_plot(df, LUNA_GROUPS)


if __name__ == '__main__':
    main()
