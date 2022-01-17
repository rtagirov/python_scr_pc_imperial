import numpy as np
import matplotlib.pyplot as plt

import itertools

import importlib

import sys
import paths;  importlib.reload(paths)
import pltaux; importlib.reload(pltaux)

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

h_r  = np.loadtxt(paths.it0h + 'murmean_r/ATM_STR', skiprows = 2, usecols = [1])
T_r  = np.loadtxt(paths.it0h + 'murmean_r/ATM_STR', skiprows = 2, usecols = [2])
n_r  = np.loadtxt(paths.it0h + 'murmean_r/ATM_STR', skiprows = 2, usecols = [3])
Tg_r = np.loadtxt(paths.it0h + 'murmean_r/ATM_STR', skiprows = 2, usecols = [6])
ng_r = np.loadtxt(paths.it0h + 'murmean_r/ATM_STR', skiprows = 2, usecols = [7])

h_f  = np.loadtxt(paths.it0h + 'murmean_f/ATM_STR', skiprows = 2, usecols = [1])
T_f  = np.loadtxt(paths.it0h + 'murmean_f/ATM_STR', skiprows = 2, usecols = [2])
n_f  = np.loadtxt(paths.it0h + 'murmean_f/ATM_STR', skiprows = 2, usecols = [3])
Tg_f = np.loadtxt(paths.it0h + 'murmean_f/ATM_STR', skiprows = 2, usecols = [6])
ng_f = np.loadtxt(paths.it0h + 'murmean_f/ATM_STR', skiprows = 2, usecols = [7])

fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (15, 10))

fig.tight_layout(h_pad = 1.5, w_pad = 2.5)

ax[0, 0].plot(h_f, T_f,              label = '256 points')
ax[0, 0].plot(h_r, T_r, color = 'r', label = '86   points')

ax[0, 1].plot(h_f, n_f)
ax[0, 1].plot(h_r, n_r, color = 'r')

ax[1, 0].plot(h_f, Tg_f)
ax[1, 0].plot(h_r, Tg_r, color = 'r')

ax[1, 1].plot(h_f, ng_f)
ax[1, 1].plot(h_r, ng_r, color = 'r')

ax[0, 1].set_yscale('log')

for i, j in itertools.product(range(2), range(2)):

    ax[i, j].set_xlabel('Height, [km]')

    ax[i, j].set_xlim(0, 950)

    ax[i, j].grid(True)

    ax[i, j].xaxis.set_minor_locator(AutoMinorLocator(4))

    if i == 0 and j == 0: ax[i, j].yaxis.set_minor_locator(AutoMinorLocator(4))
    if i == 1 and j == 0: ax[i, j].yaxis.set_minor_locator(AutoMinorLocator(5))
    if i == 1 and j == 1: ax[i, j].yaxis.set_minor_locator(AutoMinorLocator(4))
        
#    ax[i, j].yaxis.set_major_locator(MultipleLocator(0.1))

ax[0, 0].set_ylabel(r'$T$, [K]')
ax[0, 1].set_ylabel(r'$N$, [cm$^{-3}$]')
ax[1, 0].set_ylabel(r'$\Delta T / \Delta h$, [K / km]')
ax[1, 1].set_ylabel(r'$\max{(\Delta r)\Delta N} / \max{(\Delta N)}\Delta r$')

leg = ax[0, 0].legend(frameon = True, framealpha = 1, loc = 'best', prop={'size': 23})

for handle in leg.legendHandles: handle.set_linewidth(5.0)

pltaux.savepdf('mean_muram_atm')
