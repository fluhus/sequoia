import json
from glob import glob

from matplotlib import pyplot as plt

PREFIX = '../../data/'
SUFFIX = '_R1_001.fastq.gz'

COUNT_READS = False


def load_nreads():
    data = json.load(open('lines.json'))
    data = {
        k[len(PREFIX) : -len(SUFFIX)]: v // 4 for k, v in data.items() if '_R1_' in k
    }
    euro = [v for k, v in data.items() if k.startswith('Euro')]
    inh = [v for k, v in data.items() if k.startswith('Inh')]
    return euro, inh


def load_nphages():
    files = glob('*.json')
    data = {f: len(json.load(open(f))) for f in files}
    data = {k[:-5]: v for k, v in data.items()}  # Remove .json
    euro = [v for k, v in data.items() if k.startswith('Euro')]
    inh = [v for k, v in data.items() if k.startswith('Inh')]
    return euro, inh


euro, inh = load_nreads() if COUNT_READS else load_nphages()
print(euro, inh)


plt.style.use('bmh')
plt.figure(dpi=200, figsize=(4, 1.5))
plt.boxplot([euro, inh], showfliers=False, whis=0, vert=False, widths=0.7)
plt.plot(euro, [1] * len(euro), 'go')
plt.plot(inh, [2] * len(inh), 'o')
plt.yticks([1, 2], ['Solids', 'Influent'])
plt.xlabel('Reads per sample' if COUNT_READS else 'Viral species detected')
# plt.show()
plt.tight_layout()
fname = 'reads' if COUNT_READS else 'species'
plt.savefig(fname)
print(f'./{fname}.png')
