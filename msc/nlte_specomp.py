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

if len(sys.argv) == 2 and sys.argv[1] == 'r':

    for j in range(0, 31):

        if j == 2: continue

        name = 'nlte' + str(j)

        w, I =  nessy.read_spec(paths.it1f + 'spec/' + name, wvl1 = 1005, wvl2 = 10000)

        w = w / 10.0

        ws, Is = spec.running_mean(w, I, 1.0)

        np.savez(paths.npz + 'I' + str(j), I = Is)

    np.savez(paths.npz + 'w', w = ws)

sysaux.clean_dir(paths.figdir + 'nlte_specomp', mode = 'noverbose')

l, an = nessy.get_elems(paths.it1f + 'spec/nlte0/DATOM_FULL'); l.insert(0, 'lte')

w = np.load(paths.npz + 'w.npz')['w']

I = np.zeros((len(w), 31))

I[:, 0] =  np.load(paths.npz + 'I0.npz')['I']
I[:, 1] =  np.load(paths.npz + 'I1.npz')['I']
I[:, 30] = np.load(paths.npz + 'I30.npz')['I']

names_pdf = ''

for j in tqdm(range(3, 30), ncols = auxfunc.term_width(), desc = 'Plotting'):

    plt.close('all')

    names_pdf = names_pdf + str(j) + '.pdf '

    pltaux.figpar(fontsize = 20)

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    ax.set_xlabel('Wavelength, [nm]', labelpad = 15)
    ax.set_ylabel(r'$F_\mathrm{NLTE} / F_\mathrm{LTE}$')

    ax.plot(w, I[:, 30] / I[:, 0], color = 'r', linewidth = 4, label = 'All')
    ax.plot(w, I[:, 1] / I[:, 0],  color = 'c', linewidth = 4, label = 'only H')

    ax.plot(np.ones(len(I[:, 0])) * 1.0, '--', color = 'k', linewidth = 0.5)

#    ax.set_yscale('log')

    ax.set_xlim(100, 400)
#    ax.set_xlim(300, 1000)

#    ax.set_ylim(7e-3, 2)
#    ax.set_ylim(0.8, 1.1)
    ax.set_ylim(0.0, 1.25)

    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    I[:, j] = np.load(paths.npz + 'I' + str(j) + '.npz')['I']

    ax.plot(w, I[:, j] / I[:, 0], color = 'k', label = 'up to ' + l[j])

    leg = ax.legend(framealpha = 1, loc = 'best', handletextpad = 1, prop = {'size': 25})

    pltaux.savepdf('nlte_specomp/' + str(j))

#    sys.exit()

os.chdir(paths.figdir + 'nlte_specomp/')

os.system('pdftk ' + names_pdf + 'output ' + 'nlte_specomp.pdf')

os.chdir(paths.pydir)
