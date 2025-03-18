"""K-mer analysis plotting logic."""

import json
import re

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from myplot import ctx
from sklearn.decomposition import PCA

from samplenaming import fix_name


def load_data():
    data = (json.loads(x) for x in open('kmers.json'))
    names = next(data)
    data = {x[0]: np.array(x[1], dtype='int8') for x in data}
    return pd.DataFrame(data, index=names)


def main():
    print('Loading data')
    df = load_data()
    print(df.shape)

    idx = df.index.tolist()
    idx = [fix_name(x) for x in idx]
    group_re = re.compile(r'^([12]\.(Euro|Inh))')
    groups = [group_re.findall(x)[0][0] for x in idx]
    igroups = {g: i for i, g in enumerate(set(groups))}
    cmap = plt.get_cmap('Dark2')
    colors = [cmap(igroups[g]) for g in groups]
    print(groups)

    print('PCAing')
    pca = PCA(2).fit_transform(df.values)
    print(pca.shape)

    with ctx('kmers-pca'):
        plt.scatter(pca[:, 0], pca[:, 1], c=colors)
        for i in range(len(df)):
            plt.text(pca[i, 0], pca[i, 1], idx[i], fontsize=6)


if __name__ == '__main__':
    main()
