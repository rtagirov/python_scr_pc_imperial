import matplotlib.pyplot as plt
import numpy as np

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

var_lte = scipy.io.readsav(paths.satlte + 'satire_lte.sav')

w = var_lte['wl']

t = varfunc.time(var_lte['date'])

wr, tr, ssi_r = varfunc.readssi(paths.satrnm + 'ssi.dat')

#tsi_nlt = var_nlt['tsi']
#tsi_lte = var_lte['tsi']

bvar1_rnm = varfunc.bvar(ssi_r, wr, 175.0, 242.0)
bvar2_rnm = varfunc.bvar(ssi_r, wr, 242.0, 360.0)
bvar3_rnm = varfunc.bvar(ssi_r, wr, 410.0, 750.0)

bvar1_lte = varfunc.bvar(var_lte['ssi'], w, 175.0, 242.0)
bvar2_lte = varfunc.bvar(var_lte['ssi'], w, 242.0, 360.0)
bvar3_lte = varfunc.bvar(var_lte['ssi'], w, 410.0, 750.0)

idx = np.where((tr >= min(t)) & (tr <= max(t)))

t_rnm = tr[idx]

bvar1_r = bvar1_rnm[idx]
bvar2_r = bvar2_rnm[idx]
bvar3_r = bvar3_rnm[idx]

np.savez(paths.npz + 'bvar_rnm', t = t_rnm, bv1 = bvar1_r,\
                                            bv2 = bvar2_r,\
                                            bv3 = bvar3_r)

np.savez(paths.npz + 'bvar_lte', bv1 = bvar1_lte,\
                                 bv2 = bvar2_lte,\
                                 bv3 = bvar3_lte)

bv_rnm = np.load(paths.npz + 'bvar_rnm.npz')
bv_lte = np.load(paths.npz + 'bvar_lte.npz')

t_rnm = bv_rnm['t']

bvar1_rnm = bv_rnm['bv1']
bvar2_rnm = bv_rnm['bv2']
bvar3_rnm = bv_rnm['bv3']

bvar1_lte = bv_lte['bv1']
bvar2_lte = bv_lte['bv2']
bvar3_lte = bv_lte['bv3']

fig, ax = plt.subplots(nrows = 4, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

fig.suptitle('Absolute variability values', y = 1.01)

ax[0].plot(t_rnm, bvar1_rnm, color = 'r')
ax[0].plot(t, bvar1_lte, color = 'k')

ax[1].plot(t_rnm, bvar2_rnm, color = 'r')
ax[1].plot(t, bvar2_lte, color = 'k')

ax[2].plot(t_rnm, bvar3_rnm, color = 'r')
#ax[2].plot(t, bvar3_lte, color = 'k')

#ax[3].plot(t, tsi_nlt, color = 'r')
#ax[3].plot(t, tsi_lte, color = 'k')

ax[3].set_xlabel('Year')

tlist = ['175 nm - 242 nm', '242 nm - 360 nm', '410 nm - 750 nm', 'TSI']

props = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.5)

for i in range(len(ax)):

    ax[i].ticklabel_format(useOffset = False)
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))

    ax[i].set_ylabel(r'Flux, [W / m$^2$]')

    ax[i].text(0.02, 0.93, tlist[i], transform = ax[i].transAxes, fontsize = 14, verticalalignment = 'top', bbox = props)

auxplt.savepdf('var/satire_rnm_var_abs')
