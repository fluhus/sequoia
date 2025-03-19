"""Main abundance plotting logic."""

import colorsys
import csv
import itertools
import json
import re
import sys
from collections import defaultdict, namedtuple
from typing import Iterable

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from myplot import ctx
from scipy.spatial.distance import pdist, squareform
from scipy.stats import entropy, mannwhitneyu
from skbio.stats.distance import DistanceMatrix, permanova
from sklearn.decomposition import PCA
from sklearn.manifold import MDS

from abundance import AbundancePaths, load_data, sort_columns
from confidence_ellipse import confidence_ellipse
from samplenaming import LUNA_GROUPS, fix_name, fix_name2, sample_group

# Run parameters.
BATCH = 'kraken-viror'
REMOVE_SPIKE = False
ALT_TOP = None
WITH_CMAP = True

Params = namedtuple('Params', ('raw', 'tax', 'fix_names'))

# Constant data.
DATA_DIR = '/dfs7/whitesonlab/alavon/Data'
TAX_FILE = f'{DATA_DIR}/ncbi/sequence.tax.json'
REFSEQ_FILE = f'{DATA_DIR}/refseq/viral.1.genomic.json'
KRAKEN_SPECIES_FILE = '../data/species_krk.json'
PHYLUM_TO_MOLTYPE_FILE = '../data/vir_phylum_moltype.json'
HUHO_TID_FILE = '../data/human_host.tid.txt'


def load_human_host():
    refseq = json.load(open(REFSEQ_FILE))
    tax = {v['TaxID'] for v in refseq if v.get('Host', '') == 'Homo sapiens'}
    if len(tax) == 0:
        raise ValueError('found no human-host species')
    return tax


def load_refseq() -> dict[str, set[str]]:
    """Returns a dict from NC code to a set of tax IDs."""
    refseq = json.load(open(REFSEQ_FILE))
    d = defaultdict(set)
    for x in refseq:
        for acc in x['Accs']:
            d[acc].add(x['TaxID'])
    return d


def load_kraken_species() -> dict[str, set[str]]:
    krk = json.load(open(KRAKEN_SPECIES_FILE))
    d = defaultdict(set)
    for tid, accs in krk.items():
        for acc in accs:
            d[acc].add(tid)
    return d


def load_vsp2_taxids() -> set[str]:
    """Returns VSP2 tax IDs."""
    data = json.load(open('../data/vsp2_graph.json'))
    return {tid for x in data.values() for tid in x.get('TIDs', [])}


def colors(name: str, n=None, reverse=False):
    cmap = plt.get_cmap(name)
    if n is None:
        return (cmap(i) for i in itertools.count())
    if reverse:
        return (cmap(1 - i / (n - 1)) for i in range(n))
    return (cmap(i / (n - 1)) for i in range(n))


def to_pastel(c):
    h, ll, s = colorsys.rgb_to_hls(*c[:3])
    s = 0.6
    ll = 0.65  # 1 - ((1 - l) / 2)
    return colorsys.hls_to_rgb(h, ll, s)


def nsubplots(n) -> np.ndarray:
    a = int(np.ceil(n**0.5))
    b = int(np.ceil(n / a))
    return np.array([a, b])


def sorted_by_time(a) -> list:
    no_prefix = re.compile('.*?_(.*)')
    return sorted(a, key=lambda x: no_prefix.findall(x)[0])


def repl_strings(s: str, d: dict) -> str:
    for k, v in d.items():
        s = s.replace(k, v)
    return s


def fix_locations(s: Iterable[str]) -> list[str]:
    d = {
        '_LB_': 'Site 1, ',
        '_Tur_': 'Site 2, ',
        '_Wod_': 'Site 3, ',
    }
    return [repl_strings(x, d) for x in s]


def fix_dates(s: Iterable[str]) -> list[str]:
    rgx = re.compile('(\\d\\d)(\\d\\d)(\\d\\d)')
    return [rgx.sub('\\1-\\2-\\3', x) for x in s]


def enumerator():
    e = {}

    def f(x):
        if x not in e:
            e[x] = len(e)
        return e[x]

    return f


def jason_plot(df: pd.DataFrame, groups=None):
    df = df.loc[sorted(df.index.tolist())]
    # df.index = fix_locations(df.index.tolist())
    # df.index = fix_dates(df.index.tolist())
    sums = df.sum().sort_values(ascending=False)
    tops = sums.index.tolist()[:10]
    if ALT_TOP:
        if type(ALT_TOP) is int:
            tops = sums.index.tolist()[:ALT_TOP]
        elif type(ALT_TOP) is list:
            tops = ALT_TOP
        else:
            raise TypeError(f'bad type: {type(ALT_TOP)=}')
    # print('Tops:', tops)
    if 'Other' in tops:
        tops.remove('Other')
    plt.style.use('ggplot')

    if not groups:
        groups = ('',)

    nsp = nsubplots(len(groups))
    sp = nsp[1] * 100 + nsp[0] * 10 + 1

    with ctx('jason', dpi=500, sizeratio=(nsp * 0.5).tolist()):
        for i, g in enumerate(groups):
            c = colors('plasma', len(tops))
            # c = colors('tab10')
            # c = (to_pastel(x) for x in c)
            plt.subplot(sp + i)
            df1 = df.loc[[x for x in df.index.tolist() if x.startswith(g)]]
            bottom = np.zeros(len(df1))
            xx = df1.index.tolist()
            xx = [fix_name2(x) for x in xx]
            for s in tops:
                yy = df1[s].values
                if WITH_CMAP:
                    plt.bar(xx, yy, 0.95, bottom=bottom, label=s, color=next(c))
                else:
                    plt.bar(xx, yy, 0.95, bottom=bottom, label=s)
                bottom += yy
            plt.bar(
                xx,
                1 - bottom,
                0.95,
                bottom=bottom,
                label='Other',
                color='lightgrey',
            )
            plt.xticks(rotation=45, ha='right')
            plt.ylabel('Relative abundance')

        legend_handles = plt.subplot(sp + 0).get_legend_handles_labels()

    with ctx('jason_legend', dpi=300):
        plt.legend(*legend_handles, loc='center')
        plt.gca().axis('off')


def compare_abnd():
    df = load_data(AbundancePaths.tids, remove_spike=True)
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


def diversity():
    assert len(sys.argv) == 3, 'Usage: diversity remove_spike?(0/1)'
    rmspike = {'0': False, '1': True}[sys.argv[2]]
    print(f'{rmspike=}')
    df = load_data(AbundancePaths.tids, remove_spike=rmspike)
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
        plt.xlabel('Shannon-diversity per sample')
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
    print('after', df.shape)
    # vsp2_species = df.columns.tolist()
    # json.dump(
    #     list(set(top_species + vsp2_species)),
    #     open('krkcov/taxids.json', 'w'),
    # )
    # json.dump(
    #     {'top': top_species, 'vsp2': vsp2_species},
    #     open('cov_species.json', 'w'),
    # )

    del vsp2

    if False:
        # Print VSP names.
        krk = json.load(open('../data/krk_std2.json'))
        for sp in sorted([krk[x]['Name'] for x in df.columns.tolist()]):
            print(sp)
        return

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
            'VSP v2 enriched species per sample'
            + (' (without spike)' if rmspike else ' (with spike)')
        )
        plt.gca().invert_yaxis()


def plot_stuff():
    df = load_data(AbundancePaths.names, remove_spike=True)
    jason_plot(df, LUNA_GROUPS)


def dnarna():
    df = load_data(AbundancePaths.phy, remove_spike=True, phylum_mode=True)
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


def mannwhit():
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

    df = load_data(AbundancePaths.names, remove_spike=True)
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


def export():
    df = load_data(AbundancePaths.names)
    df.to_csv('data.csv')


def export_permanova():
    df = load_data(AbundancePaths.tids, remove_spike=True)
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


def plot_human_host():
    df = load_data(AbundancePaths.tids, remove_spike=False)
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


def print_human_host():
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


def rarefy():
    data: dict[str, list[list[int]]] = json.load(open('rrf.json'))
    data = {
        fix_name(k): v for k, v in data.items() if not k.startswith('Undetermined')
    }
    print(len(data), 'samples')

    plt.style.use('bmh')
    with ctx('rrf', dpi=400, sizeratio=1):
        for g in LUNA_GROUPS:
            first = True
            for k, (x, y) in data.items():
                if not k.startswith(g):
                    continue
                if first:
                    plt.plot(x, y, label=LUNA_GROUPS[g], lw=0.5)
                    c = plt.gca().lines[-1].get_color()
                    plt.plot(x[-1], y[-1], '*', c=c)
                    first = False
                else:
                    c = plt.gca().lines[-1].get_color()
                    plt.plot(x, y, c=c, lw=0.5)
                    plt.plot(x[-1], y[-1], '*', c=c)
        plt.xlabel('Sequencing depth (reads)')
        plt.ylabel('Number of OTUs')
        # plt.xscale('log')
        plt.legend()


def count_viral_reads():
    df = load_data(AbundancePaths.viror, remove_spike=False)
    vir = (1 - df['Other'].values) * 100
    print('with spike: {:.2g}% +- {:.2g}%'.format(vir.mean(), vir.std()))

    m_with_spike = {}
    for g, gdf in df.groupby(sample_group):
        vir = (1 - gdf['Other'].values) * 100
        print('  {: <8} {:.2g}% +- {:.2g}%'.format(g, vir.mean(), vir.std()))
        m_with_spike[g] = vir.tolist()

    df = load_data(AbundancePaths.viror, remove_spike='Other')
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


def main():
    actions = {
        'cvreads': count_viral_reads,
        'dnarna': dnarna,
        'diversity': diversity,
        'jason-vir': plot_stuff,
        'rarefy': rarefy,
        'huho': plot_human_host,
        'huho2': print_human_host,
    }
    args = sys.argv
    if len(args) == 1 or args[1] not in actions:
        print('Available actions:', ', '.join(sorted(actions)))
        exit(1)
    actions[args[1]]()


if __name__ == '__main__':
    # plot_stuff()
    # compare_abnd()
    # diversity()
    # dnarna()
    # rrna()
    # mannwhit()
    # lg()
    # export()
    # export_permanova()
    # plot_human_host()
    # rarefy()
    # count_viral_reads()
    # json.dump(list(load_vsp2_taxids().keys()), open('vsp2_taxids.json', 'w'))
    main()
