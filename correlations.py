"""Plots correlations between corresponding samples."""

import pandas as pd
from matplotlib import pyplot as plt
from myplot import ctx
from numpy import corrcoef

from abundance import AbundancePaths, load_data
from samplenaming import sample_group
from violin import violin


def sample_names(df: pd.DataFrame):
    return {x.removeprefix(sample_group(x) + '_') for x in df.index.tolist()}


def corr(x, y):
    return corrcoef(x, y)[0][1]


def main():
    df = load_data(AbundancePaths.BY_GENUS, remove_spike=True)
    print(df.shape)

    corrs_v1 = []
    corrs_v2 = []
    corrs_sol = []
    corrs_inf = []

    for name in sample_names(df):
        i1 = df.loc['1.Inh_' + name].values
        i2 = df.loc['2.Inh_' + name].values
        s1 = df.loc['1.Euro_' + name].values
        s2 = df.loc['2.Euro_' + name].values

        corrs_v1.append(corr(i1, s1))
        corrs_v2.append(corr(i2, s2))
        corrs_sol.append(corr(s1, s2))
        corrs_inf.append(corr(i1, i2))

    # plt.style.use('bmh')
    with ctx('correlations', sizeratio=0.6, dpi=500):
        violin(
            {
                'solid-influent, v1': corrs_v1,
                'solid-influent, v2': corrs_v2,
                'v1-v2, solid': corrs_sol,
                'v1-v2, influent': corrs_inf,
            },
            rotate_xticks=20,
        )
        plt.ylabel('Pearson correlation')


if __name__ == '__main__':
    main()
