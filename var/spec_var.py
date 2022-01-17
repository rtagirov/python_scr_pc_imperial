import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

import scipy.io

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;   importlib.reload(paths)
import spec;    importlib.reload(spec)
import nessy;   importlib.reload(nessy)
import auxsys;  importlib.reload(auxsys)
import auxplt;  importlib.reload(auxplt)
import varfunc; importlib.reload(varfunc)

var_nlt = scipy.io.readsav(paths.satnlt + 'satire_nlte.sav')
var_lte = scipy.io.readsav(paths.satlte + 'satire_lte.sav')

w = var_nlt['wl']

t = time(var_nlt['date'])

tsi_nlt = var_nlt['tsi']
tsi_lte = var_lte['tsi']

svar_nlt = svar(var_nlt['ssi'], t)
svar_lte = svar(var_lte['ssi'], t)

np.savez(paths.npz + 'svar', nlt = svar_nlt, lte = svar_lte)

svar = np.load(paths.npz + 'svar.npz')

svar_nlt = svar['nlt']
svar_lte = svar['lte']

rdsv = svar_nlt / svar_lte

rdtsi = (tsi_nlt / tsi_lte - 1.0) * 100.0

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

ax[0].plot(w, np.ones(len(w)), color = 'k', linestyle = '--')

ax[0].semilogy(w, abs(rdsv), color = 'r')
ax[1].plot(t, rdtsi, color = 'k')

ax[0].set_xlabel('Wavelength, [nm]')
ax[1].set_xlabel('Year')

ax[0].set_xlim(100, 500)
ax[0].set_ylim(1e-1, 1e+5)

ax[1].ticklabel_format(useOffset = False)

ax[1].xaxis.set_minor_locator(AutoMinorLocator(10))

tlist = ['Spectral Variability', 'TSI']

xtext = [0.82, 0.93]

ax[0].set_ylabel('NLTE / LTE')
ax[1].set_ylabel('(NLTE - LTE) / LTE, [%]')

for i in range(len(ax)):

    ax[i].text(xtext[i], 0.93, tlist[i], transform = ax[i].transAxes, fontsize = 14, verticalalignment = 'top', bbox = props)

auxplt.savepdf('satire_nlte_svar_tsi')
