"""Calculates PCA and MDS on abundances."""

import re

from matplotlib import pyplot as plt
from myplot import ctx
from scipy.spatial.distance import pdist, squareform
from skbio.stats.distance import DistanceMatrix, permanova
from sklearn.decomposition import PCA
from sklearn.manifold import MDS

from abundance import AbundancePaths, load_data
from confidence_ellipse import confidence_ellipse
from samplenaming import LUNA_GROUPS


def enumerator():
    e = {}

    def f(x):
        if x not in e:
            e[x] = len(e)
        return e[x]

    return f


def main():
    df = load_data(AbundancePaths.BY_TID, remove_spike=True)
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
    with ctx('mds', sizeratio=0.75, dpi=400):
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
        plt.title(f'NMDS (stress={mds_model.stress_:.2f})')

    for pat in ['^(1\\.|2\\.)', '(Euro|Inh)', '(_LB_|_Tur_|_Wod_)']:
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


if __name__ == '__main__':
    main()
