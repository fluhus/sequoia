import numpy as np
from matplotlib import pyplot as plt
from myplot import ctx
from scipy.stats import mannwhitneyu

from abundance import AbundancePaths, load_data


def main():
    SIG_THRESHOLD = 0.05
    MODE = 'vsp'  # vsp or state
    WITH_TITLE = False

    match MODE:
        case 'vsp':
            title = 'VSP2 vs VSP1'
            groups = [list(range(18, 36)), list(range(18))]
        case 'state1':
            title = 'Solids vs Influent (VSP1)'
            groups = [
                list(range(9)),
                list(range(9, 18)),
            ]
        case 'state2':
            title = 'Solids vs Influent (VSP2)'
            groups = [
                list(range(18, 27)),
                list(range(27, 36)),
            ]
        case _:
            raise ValueError(f'bad mode: {MODE!r}')

    df = load_data(AbundancePaths.BY_NAME, remove_spike=True)
    species = df.columns.tolist()
    tests = []

    for spc, sdata in zip(species, df.values.T):
        a, b = sdata[groups[0]], sdata[groups[1]]

        # Make sure the species appears in both groups.
        if len(np.nonzero(a)[0]) < 1 or len(np.nonzero(b)[0]) < 1:
            continue

        mw = mannwhitneyu(a, b)
        effect_size = mw.statistic / (18 * 18)
        tests.append([spc, float(effect_size), float(mw.pvalue)])
    print(len(tests), 'tests done')

    tests = [x for x in tests if x[2] <= SIG_THRESHOLD]
    print(len(tests), 'significant before correction')

    # BH correction.
    tests = sorted(tests, key=lambda x: x[2])
    for i, spc in enumerate(tests):
        spc[2] *= len(species) / (i + 1)
    tests = [x for x in tests if x[2] <= SIG_THRESHOLD]
    print(len(tests), 'significant after correction')
    # print([f'{x[0]}:{x[1]:.2f}' for x in tests])

    if not tests:
        return

    cm = plt.get_cmap('plasma')
    plt.style.use('bmh')
    with ctx('mannwhit-' + MODE, sizeratio=[0.15 * len(tests), 1.5]):
        tests = sorted(tests, key=lambda x: -x[1])
        c = [cm(x[1]) for x in tests]
        plt.bar([x[0] for x in tests], [x[1] for x in tests], color=c)
        plt.xticks(rotation=45, ha='right')
        plt.ylabel('Mann-Whitney effect size')
        if WITH_TITLE:
            plt.title(title)


if __name__ == '__main__':
    main()
