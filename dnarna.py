import json

import numpy as np
from matplotlib import pyplot as plt
from myplot import ctx

from abundance import AbundancePaths, load_data, sort_columns
from samplenaming import LUNA_GROUPS, fix_name2

PHYLUM_TO_MOLTYPE_FILE = '../data/vir_phylum_moltype.json'


def main():
    df = load_data(AbundancePaths.BY_PHY, remove_spike=True, phylum_mode=True)
    idx = df.index.tolist()
    cols = df.columns.tolist()

    mol_types = json.load(open(PHYLUM_TO_MOLTYPE_FILE))
    df.columns = [mol_types.get(c, 'Unknown') for c in cols]
    df = df.groupby(df.columns, axis=1).sum()
    df = sort_columns(df)

    # plt.style.use('bmh')
    with ctx('dnarna', dpi=500):
        for ig, g in enumerate(LUNA_GROUPS):
            plt.subplot(221 + ig)
            irow = [i for i, x in enumerate(idx) if x.startswith(g)]
            xx = list(range(len(irow)))
            bottom = np.zeros(len(xx))

            for mol in df.columns.tolist():
                vals = df.iloc[irow][mol].values
                plt.bar(xx, vals, bottom=bottom, label=mol)
                bottom += vals

            plt.xticks(
                list(range(len(irow))),
                [fix_name2(idx[i]) for i in irow],
                rotation=45,
                ha='right',
            )
            # plt.legend(loc='lower right') if ig == 1 else ...
            plt.ylabel('Relative abundance')

        legend_handles = plt.subplot(221).get_legend_handles_labels()

    with ctx('dnarna_legend'):
        plt.legend(*legend_handles, loc='center')
        plt.gca().axis('off')


if __name__ == '__main__':
    main()
