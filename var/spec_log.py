import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

import importlib
import scipy.io
import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;   importlib.reload(paths)
import spec;    importlib.reload(spec)
import nessy;   importlib.reload(nessy)
import auxsys;  importlib.reload(auxsys)
import auxplt;  importlib.reload(auxplt)
import varfunc; importlib.reload(varfunc)

w = scipy.io.readsav(paths.satnlt + 'satire_nlte.sav')['wl']

qnf = varfunc.spec('nlt', 'Q_fal')
qnk = varfunc.spec('nlt', 'Q_kur')
fnf = varfunc.spec('nlt', 'F_fal')
pnk = varfunc.spec('nlt', 'P_kur')
unk = varfunc.spec('nlt', 'U_kur')

qlk = varfunc.spec('lte', 'qs_1.dat')
flk = varfunc.spec('lte', 'fac_5.dat')
plk = varfunc.spec('lte', 'penumbra_1.dat')
ulk = varfunc.spec('lte', 'umbra_1.dat')

plt.close('all')

fig, ax = plt.subplots(nrows = 4, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

ax[0].semilogy(w, qnf, color = 'g', label = 'NESSY - FAL', linewidth = 3)
ax[0].semilogy(w, qnk, color = 'r', label = 'NESSY - KUR')
ax[0].semilogy(w, qlk, color = 'b', label = 'ATLAS - KUR')

ax[1].semilogy(w, fnf, color = 'r', label = 'NESSY - FAL')
ax[1].semilogy(w, flk, color = 'b', label = 'ATLAS - KUR')

ax[2].semilogy(w, pnk, color = 'r', label = 'NESSY - KUR')
ax[2].semilogy(w, plk, color = 'b', label = 'ATLAS - KUR')

ax[3].semilogy(w, unk, color = 'r', label = 'NESSY - KUR')
ax[3].semilogy(w, ulk, color = 'b', label = 'ATLAS - KUR')

ax[3].set_xlabel('Wavelength, [nm]')

tlist = ['Quiet Sun', 'Facula', 'Penumbra', 'Umbra']

props = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.5)

ax[0].set_ylim(1e-13, 1e+1)
ax[1].set_ylim(1e-7, 1e+1)
ax[2].set_ylim(1e-13, 1e+1)
ax[3].set_ylim(1e-18, 1e+1)

for i in range(len(ax)):

    ax[i].xaxis.set_minor_locator(AutoMinorLocator(5))

    ax[i].set_xlim(100, 300)

    ax[i].set_ylabel(r'Flux, [W / m$^2$ / nm]')

    ax[i].text(0.02, 0.93, tlist[i], transform = ax[i].transAxes, fontsize = 14, verticalalignment = 'top', bbox = props)

    leg = ax[i].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 15.5})

    for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/inp_spec_log')
