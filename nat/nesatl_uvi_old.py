import numpy as np
import pylab as pl
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import os

import paths;  importlib.reload(paths)
import pltaux; importlib.reload(pltaux)
import sysaux; importlib.reload(sysaux)

mu =    ['1.0', '0.9', '0.8', '0.7', '0.6', '0.5', '0.4', '0.3', '0.2', '0.1', '0.05']

mu_fn = ['1',   '09',  '08',  '07',  '06',  '05',  '04',  '03',  '02',  '01',  '005']

idxs_mu = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]

clv_n = np.loadtxt(paths.it0f + 'murmean/CLV_UVI')
clv_a = np.loadtxt(paths.out  + 'atl_clv.dat')

wvl_n = clv_n[:, 0] / 10.0; wvl = wvl_n
wvl_a = clv_a[:, 0]

idx_n = np.where((wvl_a > min(wvl_n)) & (wvl_a < max(wvl_n)))

sysaux.clean_dir(paths.figdir)

pltaux.figpar()

for idx_mu in idxs_mu:

    its_n = clv_n[:, idx_mu]
    its_a = np.interp(wvl_n, wvl_a[idx_n], clv_a[idx_n[0], idx_mu])

    figname = mu_fn[idx_mu - 1]

    fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    ax[0].set_title('$\mu =$ ' + mu[idx_mu - 1])

    ax[0].plot(wvl, its_n, color = 'b',                  label = 'NESSY')
    ax[0].plot(wvl, its_a, color = 'r', linewidth = 0.5, label = 'ATLAS')

    ax[0].xaxis.grid(True)
    ax[0].yaxis.grid(True)

    ax[0].set_ylabel(r'Intensity, [erg / (cm$^2$ $\times$ Hz $\times$ s $\times$ sterad)]')

    ax[0].set_yscale('log')

    ax[0].set_xlim(185, 4000)

    ax[1].plot(wvl, its_n / its_a)

    ax[1].set_xlabel('Wavelength, [nm]')
    ax[1].set_ylabel('NESSY / ATLAS')

    ax[1].xaxis.grid(True)
    ax[1].yaxis.grid(True)

    ax[1].set_xlim(185, 4000)

    ax[1].set_ylim(0.1, 1.3)

    leg = ax[0].legend(framealpha = 0, loc = 4, handlelength = 0, handletextpad=0, prop={'size': 30})

    for handle in leg.legendHandles: handle.set_visible(False)

    tn, ta = leg.get_texts()

    pl.setp(tn, color = 'b')
    pl.setp(ta, color = 'r')

    pltaux.savepdf(paths.figdir, figname)

os.system('pdftk ' + paths.figdir + '*' + ' output ' + paths.figdir + 'allmus.pdf')
