"""Counts viral reads in our data."""

from matplotlib import pyplot as plt
from myplot import ctx

from abundance import AbundancePaths, load_data
from samplenaming import LUNA_GROUPS, sample_group


def main():
    df = load_data(AbundancePaths.VIR_OR, remove_spike=False)
    vir = (1 - df['Other'].values) * 100
    print('with spike: {:.2g}% +- {:.2g}%'.format(vir.mean(), vir.std()))

    m_with_spike = {}
    for g, gdf in df.groupby(sample_group):
        vir = (1 - gdf['Other'].values) * 100
        print('  {: <8} {:.2g}% +- {:.2g}%'.format(g, vir.mean(), vir.std()))
        m_with_spike[g] = vir.tolist()

    df = load_data(AbundancePaths.VIR_OR, remove_spike='Other')
    vir = (1 - df['Other'].values) * 100
    print('without spike: {:.2g}% +- {:.2g}%'.format(vir.mean(), vir.std()))

    m_no_spike = {}
    for g, gdf in df.groupby(sample_group):
        vir = (1 - gdf['Other'].values) * 100
        print('  {: <8} {:.2g}% +- {:.2g}%'.format(g, vir.mean(), vir.std()))
        m_no_spike[g] = vir.tolist()

    plt.style.use('bmh')
    size = [0.9, 0.3]
    dpi = 600
    with ctx('vreads', dpi=dpi, sizeratio=size):
        plt.boxplot(
            list(m_with_spike.values()),
            showfliers=False,
            whis=0,
            vert=False,
            widths=0.75,
        )
        plt.boxplot(
            list(m_no_spike.values()),
            showfliers=False,
            whis=0,
            vert=False,
            widths=0.75,
        )

        for i, g in enumerate(LUNA_GROUPS):
            data = m_with_spike[g]
            plt.plot(data, [i + 1] * len(data), '^', alpha=0.5)
        plt.gca().set_prop_cycle(None)
        for i, g in enumerate(LUNA_GROUPS):
            data = m_no_spike[g]
            plt.plot(data, [i + 1] * len(data), 'o', alpha=0.5)

        plt.yticks(list(range(1, len(m_with_spike) + 1)), LUNA_GROUPS.values())
        plt.xlabel('% reads assigned to viruses per sample')
        plt.gca().invert_yaxis()
        plt.xscale('log')


if __name__ == '__main__':
    main()
