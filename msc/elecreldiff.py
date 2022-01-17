import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as patches

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import os

from tqdm import tqdm

import paths;   importlib.reload(paths)
import spec;    importlib.reload(spec)
import nessy;   importlib.reload(nessy)
import sysaux;  importlib.reload(sysaux)
import pltaux;  importlib.reload(pltaux)
import auxfunc; importlib.reload(auxfunc)

#h = np.loadtxt(paths.it1h + 'spec/nlte0/ATM_MOD', usecols = [0])
T = np.loadtxt(paths.it1h + 'spec/nlte0/ATM_MOD', usecols = [1])

if len(sys.argv) == 2 and sys.argv[1] == 'r':

    for j in range(0, 31):

        if j == 2: continue

        name = 'nlte' + str(j)

        if j == 0:

            f1, rne_lte, f2 = nessy.read_popnum(paths.it1h + 'spec/nlte0')

            rne = np.zeros((len(rne_lte), 31))

            rne[:, 0] = rne_lte

        if j != 0: rne[:, j] = np.loadtxt(paths.it1h + 'spec/' + name + '/' + paths.lev + 'ELECTR', skiprows = 2, usecols = [3])

sysaux.clean_dir(paths.figdir + 'elecreldiff', mode = 'noverbose')

l, an = nessy.get_elems(paths.it1h + 'spec/nlte0/DATOM_FULL'); l.insert(0, 'lte')

idx = np.arange(1, len(rne[:, 0]) + 1)

names_pdf = ''

for j in tqdm(range(3, 30), ncols = auxfunc.term_width(), desc = 'Plotting'):

    plt.close('all')

    names_pdf = names_pdf + str(j) + '.pdf '

    pltaux.figpar(fontsize = 20)

    fig, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    fig.suptitle('NLTE up to ' + l[j], y = 1.01)

    ax1.plot(idx, (rne[:, j] - rne[:, 30]) * 100.0 / rne[:, 30], color = 'k')

    ax1.plot(idx, np.zeros(len(rne[:, 0])), '--', color = 'k', linewidth = 0.5)

    ax1.set_xlim(len(rne[:, 0]), 1)
    ax1.set_ylim(-2, 10)

    ax1.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator(10))

    ax1.yaxis.set_major_locator(ticker.MultipleLocator(1))
    ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))

    ax1.add_patch(patches.Rectangle((48, -0.5), -8, 2.5))

    ax2 = ax1.twinx()

    ax2.plot(idx, T, color = 'g', label = 'Temperature')

#    ax1.set_yscale('log')
    ax2.set_yscale('log')

    ax2.set_xlim(82, 1)

    ax1.set_xlabel('Depth index', labelpad = 15)
    ax1.set_ylabel(r'$(n_\mathrm{e}^\mathrm{' + l[j] + '} - n_\mathrm{e}^\mathrm{All})\ /\ n_\mathrm{e}^\mathrm{All}$, [%]')
    ax2.set_ylabel('Temperature, [K]')

#    leg1 = ax1.legend(framealpha = 1, loc = 'upper right', handletextpad = 1, prop = {'size': 25})
    leg2 = ax2.legend(framealpha = 1, loc = 'lower right', handletextpad = 1, prop = {'size': 25})

    left, bottom, width, height = [0.5, 0.7, 0.3, 0.2]

    axi = fig.add_axes([left, bottom, width, height])

    axi.plot(idx, (rne[:, j] - rne[:, 30]) * 100.0 / rne[:, 30], 'k')
    axi.plot(idx, np.zeros(len(rne[:, 0])), '--', color = 'k', linewidth = 0.5)

    axi.set_xlim(48, 40)
    axi.set_ylim(-0.5, 2)

    axi.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    axi.xaxis.set_major_locator(ticker.MultipleLocator(5))

    axi.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))

    pltaux.savepdf('elecreldiff/' + str(j))

os.chdir(paths.figdir + 'elecreldiff/')

os.system('pdftk ' + names_pdf + 'output ' + 'elecreldiff.pdf')

os.chdir(paths.pydir)
