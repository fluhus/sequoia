import numpy as np
from matplotlib import pyplot as plt


def violin(d: dict[str, list], rotate_xticks=None):
    plt.violinplot(list(d.values()), showmedians=True, showextrema=True)
    plt.xticks(list(range(1, len(d) + 1)), d.keys())
    if rotate_xticks:
        plt.xticks(rotation=rotate_xticks)
    for i, yy in enumerate(d.values()):
        xx = np.array([float(i) + 1] * len(yy))
        xx += np.random.normal(scale=0.03, size=len(xx))
        plt.plot(xx, yy, 'ko', markersize=1)
