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

h = np.loadtxt(paths.it1h + 'spec/nlte0/ATM_MOD', usecols = [0])
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

sysaux.clean_dir(paths.figdir + 'nlte_edccomp_height', mode = 'noverbose')

l, an = nessy.get_elems(paths.it1h + 'spec/nlte0/DATOM_FULL'); l.insert(0, 'lte')

names_pdf = ''

for j in tqdm(range(3, 30), ncols = auxfunc.term_width(), desc = 'Plotting'):

    plt.close('all')

    names_pdf = names_pdf + str(j) + '.pdf '

    fig, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    pltaux.figpar(fontsize = 20)

    ax1.plot(h, (rne[:, 30] - rne[:, 0]) / rne[:, 0], color = 'r', linewidth = 4, label = 'All')
    ax1.plot(h, (rne[:, 1]  - rne[:, 0]) / rne[:, 0], color = 'c', linewidth = 4, label = 'only H')

    ax1.plot(np.zeros(len(rne[:, 0])), '--', color = 'k', linewidth = 0.5)

    ax1.plot(h, (rne[:, j] - rne[:, 0]) / rne[:, 0], color = 'k', label = 'up to ' + l[j])

    ax2 = ax1.twinx()

    ax2.plot(h, T, color = 'g', label = 'Temperature')

    ax2.set_yscale('log')

    ax1.set_xlim(0,    1100)
    ax2.set_ylim(4e+3, 1e+4)

    ax1.set_xlabel('Height, [km]', labelpad = 15)
    ax1.set_ylabel(r'$\left(n_\mathrm{e}^\mathrm{NLTE} - n_\mathrm{e}^\mathrm{LTE}\right) / n_\mathrm{e}^\mathrm{LTE}$')
    ax2.set_ylabel('Temperature, [K]')

    leg1 = ax1.legend(framealpha = 1, loc = 'lower left',  handletextpad = 1, prop = {'size': 25})
    leg2 = ax2.legend(framealpha = 1, loc = 'lower right', handletextpad = 1, prop = {'size': 25})

    pltaux.savepdf('nlte_edccomp_height/' + str(j))

os.chdir(paths.figdir + 'nlte_edccomp_height/')

os.system('pdftk ' + names_pdf + 'output ' + 'nlte_edccomp_height.pdf')

os.chdir(paths.pydir)
