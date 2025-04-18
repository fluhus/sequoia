"""Compares our results with DRAGEN's."""

import json
from glob import glob
from os.path import basename

import numpy as np
from matplotlib import pyplot as plt
from myplot import ctx

from abundance import AbundancePaths, load_data
from config import DATA_DIR, WS_DATA_DIR
from samplenaming import fix_name2, sample_group

V1_FILES = f'{DATA_DIR}/ww-dragen/v1/*.json'
V2_FILES = f'{DATA_DIR}/ww-dragen/v2/*.json'
KRK_FILE = f'{WS_DATA_DIR}/krk_std2.json'
MIN_ABUNDANCE = 0.00001

VSP_V1 = True


def load_dragen(files_glob: str):
    d = {}
    for f in glob(files_glob):
        data = json.load(open(f))
        b = basename(f).removesuffix('.json')
        d[b] = data
    return d


def to_dragen_name(s: str) -> str:
    g = sample_group(s)
    n = fix_name2(s)
    return f'{g[2:]}-{n[:-1]}-{n[-1:]}'


def remove_spike(d: dict) -> dict:
    spike = {'Betacoronavirus', 'Mastadenovirus'}
    d = {k: v for k, v in d.items() if k not in spike}
    s = sum(d.values())
    d = {k: v / s for k, v in d.items()}
    return d


def load_phage_genera() -> dict:
    d: dict = json.load(open(KRK_FILE))
    phages = []
    for k, v in d.items():
        if v['Level'][0] not in {'S', 'G'}:
            continue
        if 'phage' not in v['Name'].lower():
            continue
        phages.append(k)
    genera = set()
    for p in phages:
        while d[p]['Level'] != 'G':
            if d[p]['Level'][0] not in {'S', 'G'}:
                break
            p = d[p]['ParentTID']
        else:
            genera.add(d[p]['Name'])
    return genera


def main():
    phages = load_phage_genera()
    drg = load_dragen(V1_FILES if VSP_V1 else V2_FILES)
    krk = load_data(AbundancePaths.BY_GENUS, remove_spike=False)
    krk = krk[:18] if VSP_V1 else krk[18:]
    krk.index = [to_dragen_name(x) for x in krk.index.tolist()]
    # print(df)
    # return

    assert set(krk.index.tolist()) == set(drg.keys())

    for k in drg:
        dkrk = krk.loc[k]
        dkrk = {
            x: v for x, v in zip(dkrk.index.tolist(), dkrk.values.tolist()) if v > 0
        }
        ddrg = drg[k]

        dkrk = remove_spike(dkrk)
        ddrg = remove_spike(ddrg)
        n1, n2 = len(dkrk), len(ddrg)
        dkrk = {k: v for k, v in dkrk.items() if k not in phages}
        ddrg = {k: v for k, v in ddrg.items() if k not in phages}
        print('Phages removed:', n1 - len(dkrk), n2 - len(ddrg))

        species = sorted(set(ddrg) | set(dkrk))
        # print(set(ddrg) - set(dkrk))
        v1 = [np.log10(dkrk.get(s, 0) + MIN_ABUNDANCE) for s in species]
        v2 = [np.log10(ddrg.get(s, 0) + MIN_ABUNDANCE) for s in species]
        corr = float(np.corrcoef(v1, v2)[0][1])
        # print(f'{k}: {corr:.2f}')

        plt.style.use('bmh')
        namev = 'v1' if VSP_V1 else 'v2'
        with ctx(f'dragen-{namev}-{k}', sizeratio=0.75):
            plt.scatter(v1, v2)
            plt.xlabel('Kraken abundance')
            plt.ylabel('Dragen abundance')
            plt.title(f'{k}  :  $\\rho$={corr:.2f}')
            plt.xticks(
                [-5, -4, -3, -2, -1, 0],
                ['0', '0.0001', '0.001', '0.01', '0.1', '1'],
            )
            plt.yticks(
                [-5, -4, -3, -2, -1, 0],
                ['0', '0.0001', '0.001', '0.01', '0.1', '1'],
            )

        # break


if __name__ == '__main__':
    # load_phage_genera()
    main()
