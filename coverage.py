"""Creates coverage plots."""

import json
import re
import sys
from collections import defaultdict
from glob import glob
from os import path
from typing import Iterable

import numpy as np
from matplotlib import pyplot as plt
from myplot import ctx
from pympler.asizeof import asizeof
from scipy.stats import kendalltau

from samplenaming import LUNA_GROUPS

WITH_TAU = True
WITH_COVERAGE = True
WITH_LOG_Y = False

SAMPLES_FILE = 'samples.txt'
KRK_FILE = '../data/krk_viral.json'
COV_SPECIES = json.load(open('../data/cov_species.json'))
RAW_COV_FILE_NAME_RE = re.compile(r'^([^\.]+)\.([^\.]+)\.cov\.json$')
COV_ANNOT = {
    (False, False): '',
    (True, False): ' (top)',
    (False, True): ' (vsp)',
    (True, True): ' (top,vsp)',
}


def load_data(glb: str) -> Iterable[tuple[str, dict[str, np.ndarray]]]:
    """Yields pairs of (species, counts by group)."""
    return (
        (
            extract_species(f),
            json.load(open(f), object_hook=json_hook),
        )
        for f in glob(glb)
    )


def load_data_nz(glb: str) -> Iterable[tuple[str, dict[str, float]]]:
    """Yields pairs of (species, nonzero map)."""
    return (
        (
            extract_species(f),
            json.load(open(f)),
        )
        for f in glob(glb)
    )


def json_hook(d: dict):
    for k in d:
        if d[k]:
            d[k] = np.array(d[k])
    return d


def extract_species(f: str) -> str:
    return (
        path.basename(f)
        .removesuffix('.json')
        .removesuffix('.cov')
        .removesuffix('.nz')
    )


def sum_covs(covs: Iterable[np.ndarray]) -> np.ndarray:
    s = mysum(covs)
    if s is None:
        return None
    return s  # / s.max()


def mysum(a):
    s = None
    for x in a:
        if x is None:
            continue
        if s is None:
            s = x
        else:
            s += x
    return s


def sum_species(d: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    dd = {}
    for g in LUNA_GROUPS:
        dd[g] = sum_covs((v for k, v in d.items() if k.startswith(g)))
    return dd


def plot_species(d: dict[str, list[float]]):
    ax1 = plt.gca()
    ax2 = plt.twinx()
    for g, gg in LUNA_GROUPS.items():
        yy = d[g]
        if yy is None:
            yy = []
        try:
            xx = list(range(len(yy)))
        except TypeError:
            raise TypeError(f'cannot len type {type(yy)}')
        if g.startswith('1'):
            ax1.plot(xx, yy, label=gg, alpha=0.5, linewidth=0.75)
            ax2.plot([], [], label=gg, alpha=0.5, linewidth=0.75)
        if g.startswith('2'):
            ax2.plot(xx, yy, label=gg, alpha=0.5, linewidth=0.75)
            ax1.plot([], [], label=gg, alpha=0.5, linewidth=0.75)
    ax1.set_ylabel('VSP v1')
    ax2.set_ylabel('VSP v2')
    if WITH_LOG_Y:
        ax1.set_yscale('log')
        ax2.set_yscale('log')
    # plt.legend()


def load_names() -> dict[str, str]:
    raw = json.load(open(KRK_FILE))
    return {x: raw[x]['Name'] for x in raw}


def tau_wrapper(a, b):
    if a is None or b is None:
        return None
    nz = [i for i in range(len(a)) if a[i] != 0 or b[i] != 0]
    a = a[nz]  # [a[i] for i in nz]
    b = b[nz]  # [b[i] for i in nz]
    return kendalltau(a, b).statistic  # mytau.tauf(a, b)


def violin(d: dict[str, list], rotate_xticks=None):
    plt.violinplot(list(d.values()), showmedians=True, showextrema=True)
    plt.xticks(list(range(1, len(d) + 1)), d.keys())
    if rotate_xticks:
        plt.xticks(rotation=rotate_xticks)
    for i, yy in enumerate(d.values()):
        xx = np.array([float(i) + 1] * len(yy))
        xx += np.random.normal(scale=0.03, size=len(xx))
        plt.plot(xx, yy, 'ko', markersize=1)


def main():
    if len(sys.argv) not in {2, 3}:
        print('Usage:', sys.argv[0], 'input_dir [title]')
        exit(1)

    print('Starting')
    if len(sys.argv) == 2:
        in_dir, title = sys.argv[1], '{name}'
    elif len(sys.argv) == 3:
        in_dir, title = sys.argv[1:]
    in_glob = path.join(in_dir, '*.covs.json')
    in_glob_nz = path.join(in_dir, '*.nz.json')
    names = load_names()

    t1s, t2s, t12s = [], [], []
    legend = None
    for key, kdata in load_data(in_glob):
        with ctx(f'cov-{key}', dpi=300, sizeratio=0.66):
            print(f'Size: {asizeof(kdata)/(2**20):.1f}mb')
            plot_species(kdata)
            if not legend:
                legend = plt.gca().get_legend_handles_labels()

            t1, t2, t12 = None, None, None
            if WITH_TAU:
                t1 = tau_wrapper(kdata['1.Inh'], kdata['1.Euro'])
                t2 = tau_wrapper(kdata['2.Inh'], kdata['2.Euro'])
                t12 = tau_wrapper(
                    sum_covs([kdata['1.Inh'], kdata['1.Euro']]),
                    sum_covs([kdata['2.Inh'], kdata['2.Euro']]),
                )
            if t1:
                t1s.append(t1)
            if t2:
                t2s.append(t2)
            if t12:
                t12s.append(t12)

            suf = COV_ANNOT[(key in COV_SPECIES['top'], key in COV_SPECIES['vsp2'])]
            tsuf1 = f' $\\tau_1$={t1:.2f}' if t1 else ''
            tsuf2 = f' $\\tau_2$={t2:.2f}' if t2 else ''
            tsuf12 = f' $\\tau_{{12}}$={t12:.2f}' if t12 else ''
            tsuf = (tsuf1 + tsuf2 + tsuf12).strip()
            if tsuf:
                tsuf = '\n' + tsuf
            plt.title(title.format(name=names.get(key, key)) + suf + tsuf)

    with ctx('cov-legend', dpi=400):
        plt.legend(*legend, loc='center')
        plt.gca().axis('off')
    del legend

    # plt.style.use('bmh')

    if WITH_TAU:
        print(
            'Tau medians: {:.2f} {:.2f} {:.2f}'.format(
                np.median(t1s), np.median(t2s), np.median(t12s)
            )
        )
        with ctx('cov-tau', dpi=500, sizeratio=[0.5, 0.5]):
            tdata = {
                'VSP v1\nsolid vs.\ninfluent': t1s,
                'VSP v2\nsolid vs.\ninfluent': t2s,
                'VSP v1\nvs. v2': t12s,
            }
            violin({k: tdata[k] for k in tdata if tdata[k]})
            plt.ylabel('Kendall\'s tau')

    if WITH_COVERAGE:
        pp = defaultdict(list)
        group_re = re.compile('|'.join(re.escape(x) for x in LUNA_GROUPS))
        for _, spc_nz in load_data_nz(in_glob_nz):
            for smpl, p in spc_nz.items():
                g = group_re.findall(smpl)[0]
                pp[g].append(p)
        pp = {x: pp[x] for x in sorted(pp)}  # Sort keys
        with ctx('cov-prc', dpi=500, sizeratio=0.5):
            violin(pp)
            plt.ylabel('Percent coverage')


main()
