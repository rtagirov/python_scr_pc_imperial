import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import os

from tqdm import tqdm

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;   importlib.reload(paths)
import spec;    importlib.reload(spec)
import nessy;   importlib.reload(nessy)
import auxsys;  importlib.reload(auxsys)
import auxplt;  importlib.reload(auxplt)
import auxfunc; importlib.reload(auxfunc)

if len(sys.argv) == 2 and sys.argv[1] == 'r':

    for j in range(0, 30):

        if j == 2: continue

        name = str(j)

        w, I =  nessy.read_spec(paths.it1f + 'test_def/' + name, wvl1 = 1005, wvl2 = 10000)

        w = w / 10.0

        ws, Is = spec.mean_within_delta(w, I, 1.0)

        np.savez(paths.npz + 'Eon' + str(j), I = Is)

    np.savez(paths.npz + 'w', w = ws)

auxsys.clean_dir(paths.figdir + 'E_old_new', mode = 'noverbose')

l, an = nessy.get_elems(paths.it1f + 'test_def/0/datom.inp')

l.insert(0, 'old')

w = np.load(paths.npz + 'w.npz')['w']

I = np.zeros((len(w), 31))

I[:, 0] =  np.load(paths.npz + 'Eon0.npz')['I']

names_pdf = ''

for j in tqdm(range(1, 30), ncols = auxfunc.term_width(), desc = 'Plotting'):

    if j == 2: continue

    plt.close('all')

    names_pdf = names_pdf + str(j) + '.pdf '

    auxplt.figpar(fontsize = 20)

    fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    fig.subplots_adjust(hspace = 0.05)

    ax[0].set_ylabel(r'Flux, [W / m$^2$ / nm]')
    ax[1].set_xlabel('Wavelength, [nm]', labelpad = 15)
    ax[1].set_ylabel(r'$\mathrm{Flux}_\mathrm{new}\ /\ \mathrm{Flux}_\mathrm{old}$')

    ax[1].plot(np.ones(len(I[:, 0])) * 1.0, '--', color = 'k', linewidth = 0.5)

    ax[0].set_xlim(100, 800)
    ax[1].set_xlim(100, 800)

    ax[0].set_ylim(-0.1, 2.4)
    ax[1].set_ylim(0.15, 1.125)

    ax[0].xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax[1].xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax[0].yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    ax[1].yaxis.set_minor_locator(ticker.AutoMinorLocator(4))

    ax[0].xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax[1].xaxis.set_major_locator(ticker.MultipleLocator(50))

    ax[0].tick_params(labelbottom = 'off')

    I[:, j] = np.load(paths.npz + 'Eon' + str(j) + '.npz')['I']

    ax[0].plot(w, I[:, 0], color = 'k', label = 'OLD')
    ax[0].plot(w, I[:, j], color = 'r', label = 'NEW', alpha = 0.5)

    ax[1].plot(w, I[:, j] / I[:, 0], color = 'k', label = 'up to ' + l[j])

    leg0 = ax[0].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 25})
    leg1 = ax[1].legend(framealpha = 0, loc = 4, handletextpad = 0, prop = {'size': 25})

    for handle in leg0.legendHandles: handle.set_linewidth(3.0)
    for handle in leg1.legendHandles: handle.set_visible(False)

    auxplt.savepdf('E_old_new/' + str(j))

os.chdir(paths.figdir + 'E_old_new/')

os.system('pdftk ' + names_pdf + 'output ' + 'E_old_new.pdf')

os.chdir(paths.pydir + '/msc')
