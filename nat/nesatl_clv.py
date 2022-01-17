import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

from scipy import interpolate

import importlib
import sys
import os

import paths;  importlib.reload(paths)
import pltaux; importlib.reload(pltaux)
import sysaux; importlib.reload(sysaux)
import phys;   importlib.reload(phys)

def clv_band(wvl, clv, band):

    idx_b = np.where((wvl >= band[0]) & (wvl <= band[1]))

    frq_b = phys.c / (wvl[idx_b] * 1.0e-7)

    its = np.zeros(len(clv[0, :]))

    for i in range(len(clv[0, :])):

        its[i] = np.trapz(clv[idx_b[0], i], frq_b)

    return its / its[0]

def ip(mu):

    p = np.sqrt(1.0 - mu**2.0)

    return ['%.3f' % z for z in p]
    
mu = np.array([1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.05])

clv_n = np.loadtxt(paths.it0f + 'murmean_r_atl_abd/CLV')
clv_a = np.loadtxt(paths.out  + 'atl_clv_86.dat')

wvl_n = clv_n[:, 0] / 10.0
wvl_a = clv_a[:, 0]

bands = [[117.00000, 242.00000], \
         [242.00000, 360.00000], \
         [410.00000, 750.00000], \
         [1500.0000, 3000.0000], \
         [3000.0000, 3500.0000], \
         [3500.0000, 4000.0000], \
         [10020.000, 160001.02]]

#col = ['k', 'r', 'b', 'g', 'o', 'c', 'm']
col = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6']

#sysaux.clean_dir(paths.figdir)

pltaux.figpar(fontsize = 20)

fig, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize = (12, 10))

fig.tight_layout()

idx_b = 0

for band in bands:

    its_n = clv_band(wvl_n, clv_n[:, 1 : 12], band)
    its_a = clv_band(wvl_a, clv_a[:, 1 : 12], band)

    if idx_b != 6: l = str(int(band[0])) + r' nm $\rightarrow$ ' + str(int(band[1])) + ' nm'
    if idx_b == 6: l = str(int(band[0] / 1e3)) + r' $\mu$m $\rightarrow$ ' + str(int(band[1] / 1e3)) + ' $\mu$m'

    ax1.plot(mu, its_n,       color = col[idx_b], label = l)
    ax1.plot(mu, its_a, '--', color = col[idx_b])

    idx_b = idx_b + 1

ax1.xaxis.set_major_locator(MultipleLocator(0.1))
ax1.yaxis.set_major_locator(MultipleLocator(0.1))

ax1.set_xlabel(r'$\mu = \cos\theta$')
ax1.set_ylabel('$I(\mu) / I(\mu = 1)$')

ax1.grid(True)

ax1.set_xlim(1.0, 0.05)
ax1.set_ylim(0.0, 1.00)

ax2 = ax1.twiny()

ax2.set_xlim(ax1.get_xlim())

ax2.set_xticks(mu[0 : len(mu) - 1])

ax2.set_xticklabels(ip(mu[0 : len(mu) - 1]))

ax2.set_xlabel(r'$p = R / R_\odot$', labelpad = 10.5)

leg = ax1.legend(frameon = True, framealpha = 1, loc = 'best', prop={'size': 23})

for handle in leg.legendHandles: handle.set_linewidth(5.0)

pltaux.savepdf('nesatl_clv')
