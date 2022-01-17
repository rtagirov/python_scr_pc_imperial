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

mu = ['1.00', '0.90', '0.80', '0.70', '0.60', '0.50', '0.40', '0.30', '0.20', '0.10', '0.05']

idxs_mu = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])

ip = np.zeros(len(mu))

for i in range(len(mu)): ip[i] = np.sqrt(1.0 - float(mu[i])**2.0)

clv_n = np.loadtxt(paths.it0f + 'murmean_r_atl_abd/CLV_UVI')
clv_a = np.loadtxt(paths.out  + 'atl_clv_86.dat')

wvl_n = clv_n[:, 0] / 10.0
wvl_a = clv_a[:, 0]; wvl = wvl_a

idx_n = np.where((wvl > min(wvl_n)) & (wvl < max(wvl_n)))

wvl = wvl[idx_n]

#sysaux.clean_dir(paths.figdir)

col = ['g', 'r', 'g', 'r', 'r', 'r', 'r']

pltaux.figpar(fontsize = 20)

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12, 10))

fig.tight_layout()

p = str(ip[0])[0 : 6]

its_a = clv_a[idx_n[0], 1]

its_n = interpolate.interp1d(wvl_n, clv_n[:, 1])(wvl)

ax[0].plot(wvl, its_n / its_a, color = 'k', label = '$\mu =$ ' + mu[0] + ', $p =$ ' + p + '000')
ax[1].plot(wvl, its_n / its_a, color = 'k')

idx_r = np.where((idxs_mu != 1) & \
                 (idxs_mu != 2) & \
                 (idxs_mu != 3) & \
                 (idxs_mu != 4) & \
                 (idxs_mu != 5) & \
                 (idxs_mu != 6) & \
                 (idxs_mu != 7) & \
                 (idxs_mu != 8) & \
                 (idxs_mu != 9))

idx_c = 0

for idx_mu in idxs_mu[idx_r]:

    p = str(ip[idx_mu - 1])[0 : 6]

    its_a = clv_a[idx_n[0], idx_mu]

    its_n = interpolate.interp1d(wvl_n, clv_n[:, idx_mu])(wvl)

    if idx_c == 0: lw = 2.0
    if idx_c == 1: lw = 0.5

    ax[0].plot(wvl, its_n / its_a, color = col[idx_c], linewidth = lw, label = '$\mu =$ ' + mu[idx_mu - 1] + ', $p =$ ' + p)
    ax[1].plot(wvl, its_n / its_a, color = col[idx_c], linewidth = lw)

    idx_c = idx_c + 1

ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[1].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[0].yaxis.set_major_locator(MultipleLocator(0.4))
ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[0].set_ylabel('NESSY / ATLAS')
ax[1].set_ylabel('NESSY / ATLAS')

ax[0].set_xlabel('Wavelength, [nm]')
ax[1].set_xlabel('Wavelength, [nm]')

ax[0].grid(True)
ax[1].grid(True)

#ax.set_xlim(100, 4001)
ax[0].set_xlim(120, 500)
ax[1].set_xlim(500, 4001)
ax[0].set_ylim(0.0, 3.2)
ax[1].set_ylim(0.85, 1.1)

leg = ax[0].legend(frameon = True, framealpha = 1, loc = 1, prop={'size': 30})

for handle in leg.legendHandles: handle.set_linewidth(5.0)

pltaux.savepdf('nesatl_uvi')
