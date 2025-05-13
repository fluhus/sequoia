import numpy as np
from matplotlib import pyplot as plt


def violin(d: dict[str, list], separate_colors=False, rotate_xticks=None):
    if separate_colors:
        for i, v in enumerate(d.values()):
            plt.violinplot(v, [i + 1], showmedians=True, showextrema=True)
    else:
        plt.violinplot(list(d.values()), showmedians=True, showextrema=True)
    plt.xticks(list(range(1, len(d) + 1)), d.keys())
    if rotate_xticks:
        plt.xticks(rotation=rotate_xticks, ha='right')
    for i, yy in enumerate(d.values()):
        xx = np.array([float(i) + 1] * len(yy))
        xx += np.random.normal(scale=0.03, size=len(xx))
        plt.plot(xx, yy, 'ko', markersize=1)
