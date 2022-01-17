import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

import importlib
import scipy.io
import varfunc
import sys

if not '../aux/' in sys.path:

    sys.path.append('../aux/')

import paths
import auxplt

importlib.reload(varfunc)
importlib.reload(paths)

#wr, tr, ssi_r = varfunc.readssi(paths.satrnm + 'ssi_rnm.dat')

#np.savez(paths.npz + 'ssi_r', wr = wr, tr = tr, ssi_r = ssi_r)

ssir = np.load(paths.npz + 'ssi_r.npz')

wr = ssir['wr']
tr = ssir['tr']

tt, tsi_r, tsi_l, tsi_u, s = np.loadtxt(paths.satrnm + 'tsi_rnm.dat', skiprows = 24, unpack = True)

var_nlt = scipy.io.readsav(paths.inp + 'satire_nlte.sav')
var_lte = scipy.io.readsav(paths.inp + 'satire_lte.sav')

tsi_nlt = var_nlt['tsi']
tsi_lte = var_lte['tsi']

w = var_nlt['wl']

t = varfunc.time(var_nlt['date'])

#bvar0_nlt = varfunc.bvar(var_nlt['ssi'], w, 115.0, 175.0)
#bvar1_nlt = varfunc.bvar(var_nlt['ssi'], w, 175.0, 242.0)
#bvar2_nlt = varfunc.bvar(var_nlt['ssi'], w, 242.0, 360.0)
#bvar3_nlt = varfunc.bvar(var_nlt['ssi'], w, 410.0, 750.0)

#bvar0_lte = varfunc.bvar(var_lte['ssi'], w, 115.0, 175.0)
#bvar1_lte = varfunc.bvar(var_lte['ssi'], w, 175.0, 242.0)
#bvar2_lte = varfunc.bvar(var_lte['ssi'], w, 242.0, 360.0)
#bvar3_lte = varfunc.bvar(var_lte['ssi'], w, 410.0, 750.0)

#bvar0_rnm = varfunc.bvar(ssir['ssi_r'], wr, 115.0, 175.0)
#bvar1_rnm = varfunc.bvar(ssir['ssi_r'], wr, 175.0, 242.0)
#bvar2_rnm = varfunc.bvar(ssir['ssi_r'], wr, 242.0, 360.0)
#bvar3_rnm = varfunc.bvar(ssir['ssi_r'], wr, 410.0, 750.0)

#np.savez(paths.npz + 'bvar_nlt', bv0 = bvar0_nlt,\
#                                 bv1 = bvar1_nlt,\
#                                 bv2 = bvar2_nlt,\
#                                 bv3 = bvar3_nlt)

#np.savez(paths.npz + 'bvar_lte', bv0 = bvar0_lte,\
#                                 bv1 = bvar1_lte,\
#                                 bv2 = bvar2_lte,\
#                                 bv3 = bvar3_lte)

#np.savez(paths.npz + 'bvar_rnm', bv0 = bvar0_rnm,\
#                                 bv1 = bvar1_rnm,\
#                                 bv2 = bvar2_rnm,\
#                                 bv3 = bvar3_rnm)

bv_nlt = np.load(paths.npz + 'bvar_nlt.npz')
bv_lte = np.load(paths.npz + 'bvar_lte.npz')
bv_rnm = np.load(paths.npz + 'bvar_rnm.npz')

bvar0_nlt = bv_nlt['bv0']
bvar1_nlt = bv_nlt['bv1']
bvar2_nlt = bv_nlt['bv2']
bvar3_nlt = bv_nlt['bv3']

bvar0_lte = bv_lte['bv0']
bvar1_lte = bv_lte['bv1']
bvar2_lte = bv_lte['bv2']
bvar3_lte = bv_lte['bv3']

bvar0_rnm = bv_rnm['bv0']
bvar1_rnm = bv_rnm['bv1']
bvar2_rnm = bv_rnm['bv2']
bvar3_rnm = bv_rnm['bv3']

#rd0 = (bvar0_nlt / bvar0_lte - 1.0) * 100.0
#rd1 = (bvar1_nlt / bvar1_lte - 1.0) * 100.0
#rd2 = (bvar2_nlt / bvar2_lte - 1.0) * 100.0
#rd3 = (bvar3_nlt / bvar3_lte - 1.0) * 100.0

lr0 = bvar0_lte / bvar0_rnm
lr1 = bvar1_lte / bvar1_rnm
lr2 = bvar2_lte / bvar2_rnm
lr3 = bvar3_lte / bvar3_rnm

nr0 = bvar0_nlt / bvar0_rnm
nr1 = bvar1_nlt / bvar1_rnm
nr2 = bvar2_nlt / bvar2_rnm
nr3 = bvar3_nlt / bvar3_rnm

plt.close('all')

fig, ax = plt.subplots(nrows = 4, ncols = 1, figsize = (6.0, 6.0))

fig.tight_layout()

plt.subplots_adjust(hspace = 0.0)

lw = 1.0

ax[0].plot(t, bvar0_rnm, color = 'b', linewidth = lw)
ax[0].plot(t, bvar0_lte, color = 'g', linewidth = lw)
ax[0].plot(t, bvar0_nlt, color = 'r', linewidth = lw)

ax[1].plot(t, bvar1_rnm, color = 'b', linewidth = lw)
ax[1].plot(t, bvar1_lte, color = 'g', linewidth = lw)
ax[1].plot(t, bvar1_nlt, color = 'r', linewidth = lw)

ax[2].plot(t, bvar2_rnm, color = 'b', linewidth = lw)
ax[2].plot(t, bvar2_lte, color = 'g', linewidth = lw)
ax[2].plot(t, bvar2_nlt, color = 'r', linewidth = lw)

ax[3].plot(t, bvar3_rnm, color = 'b', linewidth = lw)
ax[3].plot(t, bvar3_lte, color = 'g', linewidth = lw)
ax[3].plot(t, bvar3_nlt, color = 'r', linewidth = lw)

ax[3].set_xlabel('Year', fontsize = 12)

ax[0].tick_params(labelbottom = 'off')
ax[1].tick_params(labelbottom = 'off')
ax[2].tick_params(labelbottom = 'off')

tlist = ['115 nm - 175 nm', '175 nm - 242 nm', '242 nm - 360 nm', '410 nm - 750 nm']

props1 = dict(boxstyle = 'round', facecolor = 'green', alpha = 0.5)
props2 = dict(boxstyle = 'round', facecolor = 'blue',  alpha = 0.5)
props3 = dict(boxstyle = 'round', facecolor = 'red',   alpha = 0.5)

for i in range(len(ax)):

    ax[i].ticklabel_format(useOffset = False)

    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))

    ax[i].set_ylabel(r'Flux, [W / m$^2$]')

    ax[i].set_xlim(min(t), max(t))

    axr = ax[i].twinx()

    axr.set_ylabel(tlist[i])

    axr.tick_params(axis = 'y', which = 'both', right = False, labelright = False)

#    ax[i].yaxis.set_minor_locator(AutoMinorLocator(4))

ax[0].text(2012,   0.025, 'LTE',                          fontsize = 10, bbox = props1)
ax[0].text(2013.3, 0.025, '$\mathrm{LTE}_\mathrm{CORR}$', fontsize = 10, bbox = props2)
ax[0].text(2015,   0.025, 'NLTE',                         fontsize = 10, bbox = props3)

#ax[1].text(2010.55, 0.89955, tlist[1], fontsize = 10, bbox = props)
#ax[2].text(2010.55, 1.00694, tlist[2], fontsize = 10, bbox = props)
#ax[2].text(2015.10, 1.00810, tlist[3], fontsize = 10, bbox = props)

auxplt.savepdf('var/bvar')
