import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import os

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import auxplt; importlib.reload(auxplt)
import auxsys; importlib.reload(auxsys)

time_dir = 'prop/'

runs = sys.argv[1]

parts = ['wrcont', 'como', 'etl', 'steal', 'cycle']

ct = np.zeros((len(parts), len(runs.split())))

j = 0

for run in runs.split():

    for i in range(len(parts)):

        path = paths.it0h + time_dir + run + '/'

        time = np.loadtxt(path + 'cpu_time.'  + parts[i])

        ct[i, j] = np.mean(time)

    j += 1

plt.close('all')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

fig.suptitle('Total Time: NLTE 30 - 32 min; NLTE 11 - 21 min', y = 1.02, fontsize = 30)

ax.scatter(np.arange(len(parts)), ct[:, 0], s = 100, label = 'NLTE 30')
ax.scatter(np.arange(len(parts)), ct[:, 1], s = 100, label = 'NLTE 11')

ax.set_ylabel('Mean CPU time, [s]', fontsize = 30)

ax.yaxis.grid(True)

ax.set_ylim(0, 51)
ax.set_xlim(-0.5, 4.5)

ax.yaxis.set_major_locator(MultipleLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))

ax.xaxis.set_major_locator(MultipleLocator(1))

parts.insert(0, '0.0')

parts.append('')

ax.set_xticklabels(parts, minor = False, rotation = 45)

ax.tick_params(axis = 'both', labelsize = 20)

leg = ax.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 30})

auxplt.savepdf('ct')
