"""Prints top level statistics about the abundances."""

import csv
import json
import sys

import numpy as np
from matplotlib import pyplot as plt
from myplot import ctx
from scipy.stats import entropy

from abundance import AbundancePaths, load_data
from config import WS_DATA_DIR
from samplenaming import LUNA_GROUPS, fix_name

PRINT_FOUND_VSP2_SPECIES = False

PHYLUM_TO_MOLTYPE_FILE = f'{WS_DATA_DIR}/vir_phylum_moltype.json'


def load_vsp2_taxids() -> dict[str, str]:
    """Returns a dict from VSP2 tax IDs to VSP2 name."""
    data = json.load(open(f'{WS_DATA_DIR}/vsp2_graph.json'))
    return {tid: k for k, x in data.items() for tid in x.get('TIDs', [])}


def load_phyla() -> dict[str, str]:
    krk = json.load(open(f'{WS_DATA_DIR}/krk_std2.json'))
    phyla = {}
    for tid, data in krk.items():
        if not data['Level'].startswith('S'):
            continue
        ptid = tid
        while krk[ptid]['Level'] != 'P':
            if krk[ptid]['Level'] == 'R':
                break
            ptid = krk[ptid]['ParentTID']
        else:
            assert krk[ptid]['Level'] == 'P'
            phyla[tid] = krk[ptid]['Name']
    return phyla


def load_vsp2_phyla() -> dict[str, str]:
    tids = load_vsp2_taxids()
    phyla = load_phyla()
    return {name: phyla[tid] for tid, name in tids.items() if tid in phyla}


def main():
    assert len(sys.argv) == 2, 'Usage: diversity.py remove_spike?(0/1)'
    rmspike = {'0': False, '1': True}[sys.argv[1]]
    print(f'{rmspike=}')
    df = load_data(AbundancePaths.BY_TID, remove_spike=rmspike)
    idx: list[str] = df.index.tolist()

    ents = [entropy(x) for x in df.values]
    nz = [sum(x > 0) for x in df.values]
    assert len(ents) == len(df)
    assert len(nz) == len(df)

    nreads: dict[str, int] = {}
    for x in csv.DictReader(open('samples-ww.csv')):
        if x['File'].startswith('Undetermined'):
            continue
        nreads[fix_name(x['File'])] = int(x['Number of reads (raw)'])

    # STD number of reads per VSP.
    # print(
    #     pd.DataFrame({'sample': list(nreads), 'reads': list(nreads.values())})
    #     .set_index('sample')
    #     .sort_index()
    #     .groupby(lambda x: str(x)[0])
    #     .std()
    # )
    # return

    # Read count ratio euro and inh.
    # a = (
    #     pd.DataFrame({'sample': list(nreads), 'reads': list(nreads.values())})
    #     .set_index('sample')
    #     .sort_index()
    # )
    # euro = a.iloc[list(range(9)) + list(range(18, 27))]
    # inh = a.iloc[list(range(9, 18)) + list(range(27, 36))]
    # print(euro)
    # print(inh)
    # assert [x[6:] for x in euro.index.tolist()] == [
    #     x[5:] for x in inh.index.tolist()
    # ], 'indexes don\'t match'
    # ratio = (inh.values / euro.values - 1) * 100
    # print(ratio)
    # print(ratio.mean(), ratio.std())
    # return

    groups = list(LUNA_GROUPS)
    labels = list(LUNA_GROUPS.values())
    # groups = ('',)

    greads = []
    gents = []
    gspecies = []
    for g in groups:
        greads.append([nreads[k] for k in idx if k.startswith(g)])
        gents.append([v for k, v in zip(idx, ents) if k.startswith(g)])
        gspecies.append([v for k, v in zip(idx, nz) if k.startswith(g)])

    print('Species counts:')
    for g, x in zip(groups, gspecies):
        print(
            f'- {g} ({min(x)}-{max(x)}) '
            f'Avg={np.mean(x):.1f}+-{np.std(x):.1f} '
            f'Med={np.median(x)}'
        )

    plt.style.use('bmh')

    size = [0.9, 0.3]
    dpi = 600
    with ctx('reads', dpi=dpi, sizeratio=size):
        plt.boxplot(greads, showfliers=False, whis=0, vert=False, widths=0.75)
        for i, g in enumerate(greads):
            plt.plot(g, [i + 1] * len(g), 'o', alpha=0.5)
        plt.yticks(list(range(1, len(groups) + 1)), labels)
        plt.xlabel('Reads per sample')
        plt.gca().invert_yaxis()

    with ctx('ents', dpi=dpi, sizeratio=size):
        plt.boxplot(gents, showfliers=False, whis=0, vert=False, widths=0.75)
        for i, g in enumerate(gents):
            plt.plot(g, [i + 1] * len(g), 'o', alpha=0.5)
        plt.yticks(list(range(1, len(groups) + 1)), labels)
        plt.xlabel('Viral Shannon-diversity per sample')
        plt.gca().invert_yaxis()

    with ctx('species', dpi=dpi, sizeratio=size):
        plt.boxplot(gspecies, showfliers=False, whis=0, vert=False, widths=0.75)
        for i, g in enumerate(gspecies):
            plt.plot(g, [i + 1] * len(g), 'o', alpha=0.5)
        plt.yticks(list(range(1, len(groups) + 1)), labels)
        plt.xlabel('Species per sample')
        plt.gca().invert_yaxis()

    with ctx('reads-species', sizeratio=0.75):
        for i in range(len(groups)):
            plt.plot(greads[i], gspecies[i], 'o', alpha=0.5, label=labels[i])
        plt.legend()
        plt.xlabel('Number of reads')
        plt.ylabel('Number of species')
        plt.title('Reads vs Species')
        plt.xscale('log')
        plt.yscale('log')

    with ctx('reads-ents', sizeratio=0.75):
        for i in range(len(groups)):
            plt.plot(greads[i], gents[i], 'o', alpha=0.5, label=labels[i])
        plt.legend()
        plt.xlabel('Number of reads')
        plt.ylabel('Shannon diversity')
        plt.title('Reads vs Diversity')
        plt.xscale('log')

    # top_species = df.columns.tolist()[:10]
    vsp2 = load_vsp2_taxids()
    print(len(vsp2), 'vsp2 taxids')
    print('before', df.shape, '', end='')
    df = df[[col for col in df if col in vsp2]]
    print('after1', df.shape, '', end='')
    df.columns = [vsp2[c] for c in df.columns]
    df = df.T.groupby(df.columns).sum().T
    print('after2', df.shape)
    del vsp2

    # vsp2_species = df.columns.tolist()
    # json.dump(
    #     list(set(top_species + vsp2_species)),
    #     open('krkcov/taxids.json', 'w'),
    # )
    # json.dump(
    #     {'top': top_species, 'vsp2': vsp2_species},
    #     open('cov_species.json', 'w'),
    # )

    vsp2_totals = df.sum(1)
    # print(vsp2_totals)
    # print(df)
    print(
        f'v1 enriched total abundance: {vsp2_totals[:18].mean():.2f} '
        f'+- {vsp2_totals[:18].std():.2f}'
    )
    print(
        f'v2 enriched total abundance: {vsp2_totals[18:].mean():.2f} '
        f'+- {vsp2_totals[18:].std():.2f}'
    )

    nz = [sum(x > 0) for x in df.values]
    assert len(nz) == len(df)
    gspecies = []
    for g in groups:
        gspecies.append([v for k, v in zip(idx, nz) if k.startswith(g)])

    with ctx('species-vsp2', dpi=dpi, sizeratio=size):
        # plt.boxplot(gspecies, showfliers=False, whis=0, vert=False, widths=0.75)
        for i, g in enumerate(gspecies):
            plt.plot(g, [i + 1] * len(g), 'o', alpha=0.5)
        plt.yticks(list(range(1, len(groups) + 1)), labels)
        plt.xlabel(
            'VSP v2 enriched viruses per sample'
            + (' (without spike)' if rmspike else ' (with spike)')
        )
        plt.gca().invert_yaxis()

    if PRINT_FOUND_VSP2_SPECIES:
        mol_types = json.load(open(PHYLUM_TO_MOLTYPE_FILE))
        phyla = load_vsp2_phyla()
        print()
        print('Species in the table:')
        for c in df:
            print(f'{c},{mol_types.get(phyla.get(c, ''), 'Unknown')}')


if __name__ == '__main__':
    main()
