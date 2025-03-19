"""Exports abundance data for analysis with external tools."""

import re

import pandas as pd
from scipy.spatial.distance import pdist, squareform

from abundance import AbundancePaths, load_data


def export():
    df = load_data(AbundancePaths.BY_NAME)
    df.to_csv('data.csv')


def export_permanova():
    df = load_data(AbundancePaths.BY_TID, remove_spike=True)
    d = pd.DataFrame(squareform(pdist(df.values, 'braycurtis')))
    print(d.shape)
    # print(d)
    d.to_csv('dist.csv', index=False, header=False)

    idx = df.index.tolist()
    batch = [int(x[0]) for x in idx]
    # print(batch)

    rgx = re.compile('^\\d+\\.([^_]+)')
    part = [rgx.findall(x)[0] for x in idx]
    # print(part)

    rgx = re.compile('^\\d+\\.[^_]+_([^_]+)')
    loc = [rgx.findall(x)[0] for x in idx]
    # print(loc)

    pd.DataFrame(
        {
            'batch': batch,
            'part': part,
            'loc': loc,
        }
    ).to_csv('meta.csv', index=False)


if __name__ == '__main__':
    export()
