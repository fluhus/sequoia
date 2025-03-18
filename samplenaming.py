"""Helps with processing WW sample names."""

import re
from os.path import basename

SAMPLE_NAME_MAPPING = {
    "A1": "Euro_Tur_111622",
    "A2": "Inh_Tur_111622",
    "B1": "Euro_Tur_113022",
    "B2": "Inh_Tur_113022",
    "C1": "Euro_Tur_040722",
    "C2": "Inh_Tur_040722",
    "D1": "Euro_Wod_111622",
    "D2": "Inh_Wod_111622",
    "E1": "Euro_Wod_113022",
    "H1": "Inh_Wod_113022",
    "F1": "Euro_Wod_041522",
    "G1": "Inh_Wod_041522",
    "H2": "Euro_LB_111622",
    "E2": "Inh_LB_111622",
    "A3": "Euro_LB_113022",
    "F2": "Inh_LB_113022",
    "B3": "Euro_LB_041422",
    "G2": "Inh_LB_041422",
}
SAMPLE_NAME_MAPPING2 = {
    "Euro_Tur_040722": "Tur1",
    "Inh_Tur_040722": "Tur1",
    "Euro_Tur_111622": "Tur2",
    "Inh_Tur_111622": "Tur2",
    "Euro_Tur_113022": "Tur3",
    "Inh_Tur_113022": "Tur3",
    "Euro_Wod_041522": "Wod1",
    "Inh_Wod_041522": "Wod1",
    "Euro_Wod_111622": "Wod2",
    "Inh_Wod_111622": "Wod2",
    "Euro_Wod_113022": "Wod3",
    "Inh_Wod_113022": "Wod3",
    "Euro_LB_041422": "LB1",
    "Inh_LB_041422": "LB1",
    "Euro_LB_111622": "LB2",
    "Inh_LB_111622": "LB2",
    "Euro_LB_113022": "LB3",
    "Inh_LB_113022": "LB3",
}
LUNA_GROUPS = {
    '1.Euro': 'v1 solid',
    '1.Inh': 'v1 influent',
    '2.Euro': 'v2 solid',
    '2.Inh': 'v2 influent',
}

_sample_name_re = re.compile(r'((Euro|Inh)_([a-zA-Z]+)_([0-9]+))')

_sample_group_re = re.compile(
    '^(' + '|'.join([re.escape(g) for g in LUNA_GROUPS] + ['']) + ')'
)


def fix_name(s: str) -> str:
    s = basename(s)
    if s.startswith('Euro_') or s.startswith('Inh_'):
        return '1.' + _sample_name_re.findall(s)[0][0]
    return '2.' + SAMPLE_NAME_MAPPING[s[:2]]


def fix_name2(s: str) -> str:
    return SAMPLE_NAME_MAPPING2[s[2:]]


def sample_group(s: str) -> str:
    return _sample_group_re.findall(s)[0]
