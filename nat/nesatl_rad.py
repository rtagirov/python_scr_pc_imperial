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

idxs_mu = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

ip = np.zeros(len(mu))

for i in range(len(mu)): ip[i] = np.sqrt(1.0 - float(mu[i])**2.0)

clv_n = np.loadtxt(paths.it0f + 'murmean_r_atl_abd/CLV_RAD')
clv_a = np.loadtxt(paths.out  + 'atl_clv_86.dat')

wvl_n = clv_n[:, 0] * 1e4 / 1e8
wvl_a = clv_a[:, 0] * 1e4 / 1e7; wvl = wvl_a

idx_n = np.where((wvl > min(wvl_n)) & (wvl < max(wvl_n)))

wvl = wvl[idx_n]

#sysaux.clean_dir(paths.figdir)

pltaux.figpar(fontsize = 20)

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12, 10))

fig.tight_layout()

p = str(ip[0])[0 : 6]

its_a = clv_a[idx_n[0], 1]

its_n = interpolate.interp1d(wvl_n, clv_n[:, 1])(wvl)

ax.plot(wvl, its_n / its_a, color = 'k', label = '$\mu =$ ' + mu[0] + ', $p =$ ' + p + '000')

for idx_mu in idxs_mu:

    if idx_mu != 1:

        p = str(ip[idx_mu - 1])[0 : 6]

        if len(p) == 3: p = p + '000'

        its_a = clv_a[idx_n[0], idx_mu]

        its_n = interpolate.interp1d(wvl_n, clv_n[:, idx_mu])(wvl)

        ax.plot(wvl, its_n / its_a, label = '$\mu =$ ' + mu[idx_mu - 1] + ', $p =$ ' + p)

ax.xaxis.set_major_locator(MultipleLocator(10))
ax.yaxis.set_major_locator(MultipleLocator(0.002))

ax.set_xlabel('Wavelength, [$\mu$m]')
ax.set_ylabel('NESSY / ATLAS')

ax.grid(True)

ax.set_xlim(7, 160)
#ax.set_ylim(1.002, 1.022)
ax.set_ylim(0.998, 1.020)

leg = ax.legend(frameon = True, framealpha = 1, loc = 2, prop={'size': 15})

for handle in leg.legendHandles: handle.set_linewidth(3.0)

pltaux.savepdf('nesatl_rad')
