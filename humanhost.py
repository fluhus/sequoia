"""Find human-host viruses in our data."""

import json

from matplotlib import pyplot as plt
from myplot import ctx

from abundance import AbundancePaths, load_data
from samplenaming import LUNA_GROUPS

HUHO_TID_FILE = '../data/human_host.tid.txt'


def main_plot():
    df = load_data(AbundancePaths.BY_TID, remove_spike=False)
    print(df.shape)
    # huho = load_human_host()
    huho = {x.strip() for x in open(HUHO_TID_FILE)}
    df = df.loc[:, [x in huho for x in df]]
    if df.shape[1] == 0:
        raise RuntimeError('df left with 0 columns after retaining human-host :(')
    print(df.shape)
    # print(df.columns.tolist())
    # return

    df[df > 0] = 1
    counts = df.sum(axis=1).tolist()
    idx = df.index.tolist()
    gcounts = [
        [x for x, l in zip(counts, idx) if l.startswith(g)] for g in LUNA_GROUPS
    ]

    plt.style.use('bmh')
    with ctx('humanhost', dpi=600, sizeratio=[0.9, 0.3]):
        plt.boxplot(gcounts, showfliers=False, whis=0, vert=False, widths=0.75)
        for i, g in enumerate(gcounts):
            plt.plot(g, [i + 1] * len(g), 'o', alpha=0.5)
        plt.yticks(list(range(1, len(gcounts) + 1)), list(LUNA_GROUPS.values()))
        plt.xlabel('Human-targeting viral species per sample')
        plt.gca().invert_yaxis()


def main_print():
    df = load_data(AbundancePaths.tids, remove_spike=False)
    print(df.shape)
    huho = {x.strip() for x in open(HUHO_TID_FILE)}
    df = df.loc[:, [x in huho for x in df]]
    if df.shape[1] == 0:
        raise RuntimeError('df left with 0 columns after retaining human-host :(')
    print(df.shape)
    a = df.columns.tolist()
    print(a)

    del df, huho

    krk = json.load(open('../data/krk_std.json'))
    for x in a:
        xx = krk[x]
        print(','.join([xx['Name']] + xx['Accs']))


if __name__ == '__main__':
    main_print()
