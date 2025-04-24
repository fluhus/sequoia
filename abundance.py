import json
from glob import glob
from os.path import basename

import numpy as np
import pandas as pd

from config import DATA_DIR
from samplenaming import fix_name


class AbundancePaths:
    BY_NAME = f'{DATA_DIR}/ww-kraken/*.vir.json'
    BY_TID = f'{DATA_DIR}/ww-kraken/*.tid.json'
    BY_GENUS = f'{DATA_DIR}/ww-kraken/*.gen.json'
    BY_PHY = f'{DATA_DIR}/ww-kraken/*.phy.json'
    VIR_OR = f'{DATA_DIR}/ww-kraken/*.viror.json'
    ALL_SP = f'{DATA_DIR}/ww-kraken/*.allsp.json'
    ALL_GEN = f'{DATA_DIR}/ww-kraken/*.allgen.json'
    MINDY = f'{DATA_DIR}/ww-mindy/*.json'

    SPECIES2 = f'{DATA_DIR}/ww2-kraken/*.vir.json'
    GENUS2 = f'{DATA_DIR}/ww2-kraken/*.gen.json'


SPIKE_TAXA = [
    'NC_003045',  # Bovine covid
    'NC_006213',  # Human coronavirus OC43
    'NC_001405',  # Human mastadenovirus C
    'NC_001454',  # Human mastadenovirus F
    'NC_001401',  # Adeno-associated virus 2
    # Tax IDs in Kraken
    '694003',
    '129951',
    '130309',
    '2955291',
    '1986019',
    '130310',
    # COVID spike in bundy
    'Bovine coronavirus',
    'Human coronavirus OC43',
    # COVID spike in kraken
    'Betacoronavirus 1',
    'Rabbit coronavirus HKU14',
    # COVID spike in CZID
    'Betacoronavirus',
    # Adenovirus
    'Human adenovirus 1',
    'Human mastadenovirus A',
    'Human mastadenovirus B',
    'Human mastadenovirus C',
    'Human mastadenovirus D',
    'Human mastadenovirus E',
    'Human mastadenovirus F',
    'Human mastadenovirus G',
    'Simian mastadenovirus B',
    'Simian mastadenovirus C',
    'Simian mastadenovirus F',
    'adeno-associated virus 2',
    'Mastadenovirus',
    # Flu
    'Alphainfluenzavirus influenzae',
]


def load_data(
    glb: str,
    remove_spike: bool | str,
    exp_files=36,
    normalize=True,
    phylum_mode=False,
) -> pd.DataFrame:
    files = glob(glb)
    files = [f for f in files if 'Undetermined' not in f]
    if exp_files:
        assert len(files) == exp_files, f'expected {exp_files}, got {len(files)}'
    print(len(files), 'files')
    jsons = (json.load(open(f)) for f in files)

    df = (
        pd.DataFrame(jsons)
        .assign(name=[fix_name(f) for f in files])
        .set_index('name')
        .replace(np.nan, 0)
        .sort_index()
    )
    df = sort_columns(df)

    if phylum_mode:
        # Split columns into species and phylum.
        # Using species for spike removal, then assigning phylum again.
        spc, phy = [], []
        for c in df.columns.tolist():
            s = c.split(',')
            assert len(s) == 2, f'{c!r}: {len(s)} parts'
            spc.append(s[0])
            phy.append(s[1])
        df.columns = spc

    if remove_spike:
        keep = [True] * df.shape[1]
        print('Removing spike:')
        for x in SPIKE_TAXA:
            if x in df:
                if type(remove_spike) is str:
                    df[remove_spike] += df[x]
                keep[df.columns.get_loc(x)] = False
                # del df[x]
                print('-', x)
        df = df.loc[:, keep]

    if phylum_mode:
        # Sum by phylum.
        phy = [x for x, k in zip(phy, keep) if k]
        df.columns = phy
        df = df.groupby(df.columns, axis=1).sum()

    if normalize:  # Normalize to 1
        df = df.div(df.sum(axis=1), axis=0)

    return df


def load_data2(
    glb: str,
    remove_spike: bool | str,
    exp_files=281,
    normalize=True,
    phylum_mode=False,
) -> pd.DataFrame:
    files = glob(glb)
    if exp_files:
        assert len(files) == exp_files, f'expected {exp_files}, got {len(files)}'
    print(len(files), 'files')
    jsons = (json.load(open(f)) for f in files)

    df = (
        pd.DataFrame(jsons)
        .assign(name=[clean_file_name(f) for f in files])
        .set_index('name')
        .replace(np.nan, 0)
        .sort_index()
    )
    df = sort_columns(df)

    # if phylum_mode:
    #     # Split columns into species and phylum.
    #     # Using species for spike removal, then assigning phylum again.
    #     spc, phy = [], []
    #     for c in df.columns.tolist():
    #         s = c.split(',')
    #         assert len(s) == 2, f'{c!r}: {len(s)} parts'
    #         spc.append(s[0])
    #         phy.append(s[1])
    #     df.columns = spc

    if remove_spike:
        keep = [True] * df.shape[1]
        print('Removing spike:')
        for x in SPIKE_TAXA:
            if x in df:
                if type(remove_spike) is str:
                    df[remove_spike] += df[x]
                keep[df.columns.get_loc(x)] = False
                # del df[x]
                print('-', x)
        df = df.loc[:, keep]

    # Remove empty samples.
    nbefore = len(df)
    df = df[df.sum(axis=1) != 0]
    nafter = len(df)
    print('Empty samples removed:', nbefore - nafter)

    # if phylum_mode:
    #     # Sum by phylum.
    #     phy = [x for x, k in zip(phy, keep) if k]
    #     df.columns = phy
    #     df = df.groupby(df.columns, axis=1).sum()

    if normalize:  # Normalize to 1
        df = df.div(df.sum(axis=1), axis=0)

    return df


def sort_columns(df: pd.DataFrame) -> pd.DataFrame:
    column_sums = df.sum(axis=0)
    column_sums = column_sums.sort_values(ascending=False)
    return df.reindex(columns=column_sums.index)


def clean_file_name(f: str) -> str:
    f = basename(f)
    for s in ('.json', '.vir', '.gen'):
        f = f.removesuffix(s)
    return f
