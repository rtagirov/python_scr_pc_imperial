import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

import importlib
import sys
import os

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)
import phys;   importlib.reload(phys)

ga = np.loadtxt(paths.atlruns + 'test/lev.out', usecols = [0])
gn = np.loadtxt(paths.atlruns + 'test/lev.out', usecols = [1])
ea = np.loadtxt(paths.atlruns + 'test/lev.out', usecols = [2])
en = np.loadtxt(paths.atlruns + 'test/lev.out', usecols = [3])

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12, 10))

fig.tight_layout()

fig.suptitle('Statistical weights and energies of the first 48 terms of FeI: ATLAS vs. NIST Atomic Database', y = 1.02)

ax[0].scatter(np.arange(1, 49), ga, color = 'r',                            label = 'ATLAS')
ax[0].scatter(np.arange(1, 49), gn, color = 'k', s = 10.4, marker = (5, 2), label = 'NIST')

ax[0].set_xlim(1, 48)

ax[0].yaxis.grid(True)

ax[0].set_ylabel('Statistical weights', fontsize = 12.5)

leg = ax[0].legend(framealpha = 1, loc = 'best', handletextpad = 1, prop = {'size': 20.5})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

ax[1].scatter(np.arange(1, 49), (ea - en) * 100.0 / ea, color = 'k')

for j in range(2):

    ax[j].xaxis.set_minor_locator(AutoMinorLocator(5))
    ax[j].xaxis.set_major_locator(MultipleLocator(5))

    for i in range(1, 49):

        ax[j].axvline(x = i, linewidth = 0.3)

ax[1].set_xlim( 1.0, 48)
ax[1].set_ylim(-2.5, 20)

ax[1].yaxis.grid(True)

ax[1].set_xlabel('Term Number',                                              fontsize = 12.5)
ax[1].set_ylabel(r'$(E_{ATLAS} - E_{NIST}) \times 100\ /\ E_{ATLAS}, [\%]$', fontsize = 12.5)

auxplt.savepdf('nist_comp')
