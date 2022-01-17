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

sysaux.clean_dir(paths.figdir + 'nlte_specomp_abs', mode = 'noverbose')

l, an = nessy.get_elems(paths.it1f + 'spec/nlte0/DATOM_FULL'); l.insert(0, 'lte')

w = np.load(paths.npz + 'w.npz')['w']

I = np.zeros((len(w), 31))

I[:, 0] =  np.load(paths.npz + 'I0.npz')['I']
I[:, 1] =  np.load(paths.npz + 'I1.npz')['I']
I[:, 30] = np.load(paths.npz + 'I30.npz')['I']

names_pdf = ''

label = ''

for j in tqdm(range(1, 30), ncols = auxfunc.term_width(), desc = 'Plotting'):

    if j == 2:  continue
    if j == 3:  continue
    if j == 4:  continue
    if j == 5:  continue
    if j == 7:  continue
    if j == 9:  continue
    if j == 10: continue
    if j == 11: continue
    if j == 15: continue
    if j == 17: continue
    if j == 18: continue
    if j == 19: continue
    if j == 21: continue
    if j == 22: continue
    if j == 23: continue
    if j == 24: continue
    if j == 25: continue
    if j == 27: continue
    if j == 29: continue
    if j == 30: continue

    plt.close('all')

    names_pdf = names_pdf + str(j) + '.pdf '

    pltaux.figpar(fontsize = 20)

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    ax.set_xlabel('Wavelength, [nm]', labelpad = 15)
#    ax.set_ylabel(r'$F_\mathrm{NLTE} / F_\mathrm{LTE}$')
    ax.set_ylabel(r'Flux, $[\mathrm{W} / \mathrm{m}^2 / \mathrm{nm}]$')

#    ax.plot(w, I[:, 30] / I[:, 0], color = 'r', linewidth = 4, label = 'All')
    ax.plot(w, I[:, 0],  color = 'c', linewidth = 2, label = 'LTE (All)')
    ax.plot(w, I[:, 30], color = 'r', linewidth = 2, label = 'NLTE (All)')

#    ax.plot(np.ones(len(I[:, 0])) * 1.0, '--', color = 'k', linewidth = 0.5)

#    ax.set_yscale('log')

    ax.set_xlim(100, 400)
#    ax.set_xlim(300, 1000)

#    ax.set_ylim(7e-3, 2)
#    ax.set_ylim(0.8, 1.1)
#    ax.set_ylim(0.0, 1.25)

    ax.set_yscale('log')

    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    I[:, j] = np.load(paths.npz + 'I' + str(j) + '.npz')['I']

    if j == 1: label = label + l[j]

    if j != 1: label = label + '+' + l[j]

    ax.plot(w, I[:, j], color = 'k', label = 'NLTE (' + label + ')')

    leg = ax.legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 25})

    pltaux.savepdf('nlte_specomp_abs/' + str(j))

#    sys.exit()

os.chdir(paths.figdir + 'nlte_specomp_abs/')

os.system('pdftk ' + names_pdf + 'output ' + 'nlte_specomp.pdf')

os.chdir(paths.pydir)
