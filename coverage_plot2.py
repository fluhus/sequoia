"""Coverage plots for second batch results."""

import json
import sys
from array import array
from glob import glob

import numpy as np
from matplotlib import pyplot as plt
from myplot import ctx


def load_data(dir: str):
    files = glob(dir + '/*.cov.json')
    print(len(files), 'files')

    d = {}
    for f in files:
        j = json.load(open(f))
        for k, v in j.items():
            if k not in d:
                d[k] = array('Q', v)
                continue
            dk = d[k]
            if len(dk) < len(v):
                dif = len(v) - len(dk)
                dk.extend([0] * dif)
            for i, x in enumerate(v):
                dk[i] += x

    return d


def main():
    data = load_data(sys.argv[1])
    print([(k, len(v)) for k, v in data.items()])
    for k, v in data.items():
        with ctx(k):
            plt.plot(np.arange(len(v)), v)
            plt.title(k)


if __name__ == '__main__':
    main()
