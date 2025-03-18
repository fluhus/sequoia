"""Prints average viral percentage."""

import sys

import numpy as np

in_file = sys.argv[1]
data = [float(x.strip()) for x in open(in_file)]

print(data)
mean, std = np.mean(data), np.std(data)
print(f'{mean:.1f} +- {std:.1f}')
