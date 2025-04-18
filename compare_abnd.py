"""Calculates PCA and MDS on abundances."""

import re

import pandas as pd
from matplotlib import pyplot as plt
from myplot import ctx
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix, permanova
from sklearn.decomposition import PCA
from sklearn.manifold import MDS

from abundance import AbundancePaths, load_data
from confidence_ellipse import confidence_ellipse
from config import WS_DATA_DIR
from samplenaming import LUNA_GROUPS

GLOBAL_MODE = False
"""If true, compare abundances of ALL species."""

EXPORT_TO_R = False
"""If true, create CSVs for running permanova in R."""


def enumerator():
    e = {}

    def f(x):
        if x not in e:
            e[x] = len(e)
        return e[x]

    return f


def main():
    if GLOBAL_MODE:
        print('Running in GLOBAL MODE')
        df = load_data(AbundancePaths.ALL_GEN, remove_spike=True)
    else:
        df = load_data(AbundancePaths.BY_GENUS, remove_spike=True)
    # print(df)

    d = squareform(pdist(df.values, 'braycurtis'))
    print('Distances:', d.shape)

    pca = PCA(2).fit_transform(d)
    plt.style.use('bmh')
    with ctx('pca', sizeratio=1.5):
        for g in LUNA_GROUPS:
            ii = [i for i, x in enumerate(df.index.tolist()) if x.startswith(g)]
            if not ii:
                raise Exception(f'group {g!r} not found')
            plt.scatter(pca[ii, 0], pca[ii, 1], label=LUNA_GROUPS[g])
        plt.legend()
        plt.title('PCA')

    mds_model = MDS(2, metric=False, dissimilarity='precomputed')
    mds = mds_model.fit_transform(d, init=pca)
    print('Stress:', mds_model.stress_)
    suffix = '-global' if GLOBAL_MODE else '-viral'
    with ctx(f'mds{suffix}', sizeratio=0.75, dpi=400):
        for g in ('1.Euro', '1.Inh', '2.Euro', '2.Inh'):
            ii = [i for i, x in enumerate(df.index.tolist()) if x.startswith(g)]
            if not ii:
                raise Exception(f'group {g!r} not found')
            sc = plt.scatter(mds[ii, 0], mds[ii, 1], label=LUNA_GROUPS[g], alpha=0.5)
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
        prefix = 'Global ' if GLOBAL_MODE else 'Viral '
        plt.title(f'{prefix}NMDS (stress={mds_model.stress_:.2f})')

    for pat in ('^(1\\.|2\\.)', '(Euro|Inh)', '(_LB_|_Tur_|_Wod_)'):
        grouper = re.compile(pat)
        prm_grouping = []
        prm_indexes = []
        enm = enumerator()
        for i, x in enumerate(df.index.tolist()):
            m = grouper.findall(x)
            if not m:
                continue
            m = m[0]
            prm_grouping.append(enm(m))
            prm_indexes.append(i)

        # print(prm_indexes)
        # print(prm_grouping)

        prm_d = d[prm_indexes][:, prm_indexes].copy()
        prmnv = permanova(
            DistanceMatrix(prm_d),
            prm_grouping,
            permutations=100000,
        )
        print('Permanova on', grouper.pattern, ':')

        prmnv = prmnv.tolist()
        n_samples = float(prmnv[2])
        n_groups = float(prmnv[3])
        f = float(prmnv[4])
        pval = float(prmnv[5])

        df1 = n_groups - 1
        df2 = n_samples - n_groups
        r2 = (f * df1) / (f * df1 + df2)
        print(f'{r2=:.2f} {pval=:.2g}')

    if EXPORT_TO_R:
        print('Exporting!')
        open(WS_DATA_DIR + '/prmnv_dist.csv', 'w').write(
            '\n'.join(','.join(str(x) for x in line) for line in d)
        )
        rgx = re.compile('^(\\d)\\.(Inh|Euro)_([^_]+)')
        meta = [rgx.findall(x)[0] for x in df.index.tolist()]
        meta = pd.DataFrame(meta, columns=['batch', 'part', 'loc'], index=df.index)
        meta.to_csv(WS_DATA_DIR + '/prmnv_meta.csv', index=False)


if __name__ == '__main__':
    main()
