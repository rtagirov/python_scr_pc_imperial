import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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

sysaux.clean_dir(paths.figdir + 'nlte_edccomp', mode = 'noverbose')

l, an = nessy.get_elems(paths.it1h + 'spec/nlte0/DATOM_FULL'); l.insert(0, 'lte')

names_pdf = ''

for j in tqdm(range(3, 30), ncols = auxfunc.term_width(), desc = 'Plotting'):

    plt.close('all')

    names_pdf = names_pdf + str(j) + '.pdf '

    fig, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    pltaux.figpar(fontsize = 20)

    ax1.plot(rne[:, 30] / rne[:, 0], color = 'r', linewidth = 4, label = 'All')
    ax1.plot(rne[:, 1]  / rne[:, 0], color = 'c', linewidth = 4, label = 'only H')

    ax1.plot(np.ones(len(rne[:, 0])) * 1.0, '--', color = 'k', linewidth = 0.5)

    ax1.set_xlim(len(rne[:, 0]), 1)
    ax1.set_ylim(1e-1, 2)

#    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

#    I[:, j] = np.load(paths.npz + 'I' + str(j) + '.npz')['I'] 

    ax1.plot(rne[:, j] / rne[:, 0], color = 'k', label = 'up to ' + l[j])

    ax2 = ax1.twinx()

    ax2.plot(T, color = 'g', label = 'Temperature')

    ax1.set_yscale('log')
    ax2.set_yscale('log')

    ax1.set_xlabel('Depth index', labelpad = 15)
    ax1.set_ylabel(r'$n_\mathrm{e}^\mathrm{NLTE} / n_\mathrm{e}^\mathrm{LTE}$')
    ax2.set_ylabel('Temperature, [K]')

    leg1 = ax1.legend(framealpha = 1, loc = 'upper left',  handletextpad = 1, prop = {'size': 25})
    leg2 = ax2.legend(framealpha = 1, loc = 'lower right', handletextpad = 1, prop = {'size': 25})

    pltaux.savepdf('nlte_edccomp/' + str(j))

os.chdir(paths.figdir + 'nlte_edccomp/')

os.system('pdftk ' + names_pdf + 'output ' + 'nlte_edccomp.pdf')

os.chdir(paths.pydir)
