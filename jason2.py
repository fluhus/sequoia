"""Abundance plots."""

import colorsys
import itertools
import json
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from myplot import ctx
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.manifold import MDS

from abundance import AbundancePaths, load_data2
from confidence_ellipse import confidence_ellipse

# Run parameters.
ALT_TOP = None
WITH_CMAP = True


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


def df_top(df: pd.DataFrame) -> list[str]:
    sums = df.sum().sort_values(ascending=False)
    tops = sums.index.tolist()[:10]
    if ALT_TOP:
        if type(ALT_TOP) is int:
            tops = sums.index.tolist()[:ALT_TOP]
        elif type(ALT_TOP) is list:
            tops = ALT_TOP
        else:
            raise TypeError(f'bad type: {type(ALT_TOP)=}')
    if 'Other' in tops:
        tops.remove('Other')
    return tops


def jason_plot(df: pd.DataFrame):
    df = df.loc[sorted(df.index.tolist())]
    # df.index = fix_locations(df.index.tolist())
    # df.index = fix_dates(df.index.tolist())
    # tops = df_top(df)
    plt.style.use('ggplot')

    # Divide DF into groups.
    rgx = re.compile('^([^_]+)_([^_]+)_(INF|SOL)_')
    groups = defaultdict(list)
    for x in df.index.tolist():
        m = rgx.findall(x)
        if not m:
            continue
        m = m[0]  # First result.
        key = (m[1], m[2])
        val = (x, m[0])
        groups[key].append(val)
    # print(groups.keys())
    # return

    # legend_handles = None
    type_titles = {'INF': 'influent', 'SOL': 'solid'}

    for grp in groups:
        with ctx(f'abundance_{grp[0]}_{grp[1]}', dpi=500, sizeratio=[1.5, 0.75]):
            plt.subplot(1, 2, 1)
            gdf = df.loc[[x[0] for x in groups[grp]]]
            assert len(gdf) > 0
            tops = df_top(gdf)
            c = colors('plasma', len(tops))
            # c = (to_pastel(x) for x in c)
            bottom = np.zeros(len(gdf))
            xx = [x[1] for x in groups[grp]]
            for s in tops:
                yy = gdf[s].values
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
            plt.title(f'{grp[0]}, {type_titles[grp[1]]}')

            legend_handles = plt.gca().get_legend_handles_labels()

            # with ctx('abundance_legend', dpi=300):
            plt.subplot(1, 2, 2)
            plt.legend(*legend_handles, loc='center')
            plt.gca().axis('off')


def nmds_plot(df: pd.DataFrame):
    d = squareform(pdist(df.values, 'braycurtis'))
    print('Distances:', d.shape)

    pca = PCA(2).fit_transform(d)

    mds_model = MDS(2, metric=False, dissimilarity='precomputed')
    mds = mds_model.fit_transform(d, init=pca)
    stress = mds_model.stress_
    print('Stress:', stress)
    del mds_model
    suffix = ''  # '-global' if GLOBAL_MODE else '-viral'

    rgx = re.compile('_(INF|SOL)_')
    rgx = re.compile('_(ESP|LOB|MER|MOD|TUR|UCD|WIN)_')

    def groupf(x):
        m = rgx.findall(x)
        if not m:
            return None
        return m[0]

    groups = group_indices(df.index.tolist(), groupf)
    print(list(groups))

    plt.style.use('bmh')
    with ctx(f'mds{suffix}', sizeratio=0.75, dpi=400):
        for g, ii in groups.items():
            sc = plt.scatter(mds[ii, 0], mds[ii, 1], label=g, alpha=0.5)
            color = sc.get_facecolors()[0]
            confidence_ellipse(
                mds[ii, 0],
                mds[ii, 1],
                plt.gca(),
                2,
                facecolor=color,
                zorder=0,
                alpha=0.1,
            )
        plt.legend()
        prefix = ''  # 'Global ' if GLOBAL_MODE else 'Viral '
        plt.title(f'{prefix}NMDS (stress={stress:.2f})')


def group_indices(a, f) -> dict:
    d = defaultdict(list)
    for i, x in enumerate(a):
        fx = f(x)
        if fx is None:
            continue
        d[fx].append(i)
    return d


def filter_by_read_count(df: pd.DataFrame, min_count: int) -> pd.DataFrame:
    rc = json.load(open('readcounts.json'))
    loc = [rc[x] >= min_count for x in df.index.tolist()]
    print('Small samples removed:', len(loc) - sum(loc))
    return df.loc[loc]


def main():
    df = load_data2(
        AbundancePaths.GENUS2,
        remove_spike=True,
    )
    df = filter_by_read_count(df, 5000000)
    jason_plot(df)
    # nmds_plot(df)


if __name__ == '__main__':
    main()
