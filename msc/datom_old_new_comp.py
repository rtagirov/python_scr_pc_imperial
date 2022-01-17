import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

from itertools import product

import importlib
import math
import sys
import scipy.io

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import phys;   importlib.reload(phys)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)

prefix0 = paths.it0f
prefix1 = paths.it1f

#wn, qfo_lte = nessy.read_spec(prefix0 + 'VAR/Q_fal_old', wvl1 = 1005, wvl2 = 16000)
#wn, qfn_lte = nessy.read_spec(prefix0 + 'VAR/Q_fal_new', wvl1 = 1005, wvl2 = 16000)
#wn, qko_lte = nessy.read_spec(prefix0 + 'VAR/Q_kur_old', wvl1 = 1005, wvl2 = 16000)
#wn, qkn_lte = nessy.read_spec(prefix0 + 'VAR/Q_kur_new', wvl1 = 1005, wvl2 = 16000)
#wn, ffo_lte = nessy.read_spec(prefix0 + 'VAR/F_fal_old', wvl1 = 1005, wvl2 = 16000)
#wn, ffn_lte = nessy.read_spec(prefix0 + 'VAR/F_fal_new', wvl1 = 1005, wvl2 = 16000)
#wn, fko_lte = nessy.read_spec(prefix0 + 'VAR/F_kur_old', wvl1 = 1005, wvl2 = 16000)
#wn, fkn_lte = nessy.read_spec(prefix0 + 'VAR/F_kur_new', wvl1 = 1005, wvl2 = 16000)
#wn, qfo_nlt = nessy.read_spec(prefix1 + 'VAR/Q_fal_old', wvl1 = 1005, wvl2 = 16000)
#wn, qfn_nlt = nessy.read_spec(prefix1 + 'VAR/Q_fal_new', wvl1 = 1005, wvl2 = 16000)
#wn, qko_nlt = nessy.read_spec(prefix1 + 'VAR/Q_kur_old', wvl1 = 1005, wvl2 = 16000)
#wn, qkn_nlt = nessy.read_spec(prefix1 + 'VAR/Q_kur_new', wvl1 = 1005, wvl2 = 16000)
#wn, ffo_nlt = nessy.read_spec(prefix1 + 'VAR/F_fal_old', wvl1 = 1005, wvl2 = 16000)
#wn, ffn_nlt = nessy.read_spec(prefix1 + 'VAR/F_fal_new', wvl1 = 1005, wvl2 = 16000)
#wn, fko_nlt = nessy.read_spec(prefix1 + 'VAR/F_kur_old', wvl1 = 1005, wvl2 = 16000)
#wn, fkn_nlt = nessy.read_spec(prefix1 + 'VAR/F_kur_new', wvl1 = 1005, wvl2 = 16000)
#
#wn = wn / 10.0
#
#wns, qfo_lte = spec.mean_within_delta(wn, qfo_lte, 10)
#wns, qfn_lte = spec.mean_within_delta(wn, qfn_lte, 10)
#wns, qko_lte = spec.mean_within_delta(wn, qko_lte, 10)
#wns, qkn_lte = spec.mean_within_delta(wn, qkn_lte, 10)
#wns, ffo_lte = spec.mean_within_delta(wn, ffo_lte, 10)
#wns, ffn_lte = spec.mean_within_delta(wn, ffn_lte, 10)
#wns, fko_lte = spec.mean_within_delta(wn, fko_lte, 10)
#wns, fkn_lte = spec.mean_within_delta(wn, fkn_lte, 10)
#wns, qfo_nlt = spec.mean_within_delta(wn, qfo_nlt, 10)
#wns, qfn_nlt = spec.mean_within_delta(wn, qfn_nlt, 10)
#wns, qko_nlt = spec.mean_within_delta(wn, qko_nlt, 10)
#wns, qkn_nlt = spec.mean_within_delta(wn, qkn_nlt, 10)
#wns, ffo_nlt = spec.mean_within_delta(wn, ffo_nlt, 10)
#wns, ffn_nlt = spec.mean_within_delta(wn, ffn_nlt, 10)
#wns, fko_nlt = spec.mean_within_delta(wn, fko_nlt, 10)
#wns, fkn_nlt = spec.mean_within_delta(wn, fkn_nlt, 10)
#
#np.savez(paths.npz + 'datom_new_old_comp', w = wns,
#                                           qfo_lte = qfo_lte,\
#                                           qfn_lte = qfn_lte,\
#                                           qko_lte = qko_lte,\
#                                           qkn_lte = qkn_lte,\
#                                           ffo_lte = ffo_lte,\
#                                           ffn_lte = ffn_lte,\
#                                           fko_lte = fko_lte,\
#                                           fkn_lte = fkn_lte,\
#                                           qfo_nlt = qfo_nlt,\
#                                           qfn_nlt = qfn_nlt,\
#                                           qko_nlt = qko_nlt,\
#                                           qkn_nlt = qkn_nlt,\
#                                           ffo_nlt = ffo_nlt,\
#                                           ffn_nlt = ffn_nlt,\
#                                           fko_nlt = fko_nlt,\
#                                           fkn_nlt = fkn_nlt)

contr = np.load(paths.npz + 'datom_new_old_comp.npz')

w = contr['w']

qfo_lte = contr['qfo_lte']
qfn_lte = contr['qfn_lte']
qko_lte = contr['qko_lte']
qkn_lte = contr['qkn_lte']
ffo_lte = contr['ffo_lte']
ffn_lte = contr['ffn_lte']
fko_lte = contr['fko_lte']
fkn_lte = contr['fkn_lte']
qfo_nlt = contr['qfo_nlt']
qfn_nlt = contr['qfn_nlt']
qko_nlt = contr['qko_nlt']
qkn_nlt = contr['qkn_nlt']
ffo_nlt = contr['ffo_nlt']
ffn_nlt = contr['ffn_nlt']
fko_nlt = contr['fko_nlt']
fkn_nlt = contr['fkn_nlt']

plt.close('all')

#fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (8.00, 8.00))
#fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (8.00, 4.00))
fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (18.00, 6.03))

auxplt.figpar(3, 3, 20)

fig.tight_layout()

lw = 1.0

#ax[0, 0].title.set_text('FAL99-C, old datom')

#ax[0, 0].plot(w, qfo_lte, color = 'k', linewidth = lw, linestyle = '--', label = 'LTE')
#ax[0, 0].plot(w, qfo_nlt, color = 'k', linewidth = lw,                   label = 'NLTE')

ax.title.set_text('Kurucz-QS, old datom')

ax.plot(w, qko_lte, color = 'k', linewidth = lw, linestyle = '--', label = 'LTE')
ax.plot(w, qko_nlt, color = 'k', linewidth = lw,                   label = 'NLTE')

#ax[1, 0].title.set_text('FAL99-C, new datom')

#ax[1, 0].plot(w, qfn_lte, color = 'k', linewidth = lw, linestyle = '--', label = 'LTE')
#ax[1, 0].plot(w, qfn_nlt, color = 'k', linewidth = lw,                   label = 'NLTE')

#ax[1, 1].title.set_text('Kurucz-QS, new datom')

#ax[1, 1].plot(w, qkn_lte, color = 'k', linewidth = lw, linestyle = '--', label = 'LTE')
#ax[1, 1].plot(w, qkn_nlt, color = 'k', linewidth = lw,                   label = 'NLTE')

ax.set_xlim(170, 800)
ax.set_ylim(0.84, 2.2)
ax.set_ylim(bottom = 0.0)

ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))

ax.set_ylabel(r'Flux, [W / m$^2$ / nm]', fontsize = 20)
#ax.set_ylabel(r'Flux, [W / m$^2$ / nm]', fontsize = 20)

ax.set_xlabel('Wavelength, [nm]', fontsize = 20)
#ax.set_xlabel('Wavelength, [nm]', fontsize = 20)

leg0 = ax.legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 15})

auxplt.savepdf('var/datom_old_new_comp_qs')

fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (8.00, 8.00))

auxplt.figpar(3, 3, 20)

fig.tight_layout()

lw = 1.0

ax[0, 0].title.set_text('FAL99-P, old datom')

ax[0, 0].plot(w, ffo_lte, color = 'r', linewidth = lw, linestyle = '--', label = 'LTE')
ax[0, 0].plot(w, ffo_nlt, color = 'r', linewidth = lw,                   label = 'NLTE')

ax[0, 1].title.set_text('Kurucz-F, old datom')

ax[0, 1].plot(w, fko_lte, color = 'r', linewidth = lw, linestyle = '--')
ax[0, 1].plot(w, fko_nlt, color = 'r', linewidth = lw,                 )

ax[1, 0].title.set_text('FAL99-P, new datom')

ax[1, 0].plot(w, ffn_lte, color = 'r', linewidth = lw, linestyle = '--')
ax[1, 0].plot(w, ffn_nlt, color = 'r', linewidth = lw,                 )

ax[1, 1].title.set_text('Kurucz-F, new datom')

ax[1, 1].plot(w, fkn_lte, color = 'r', linewidth = lw, linestyle = '--')
ax[1, 1].plot(w, fkn_nlt, color = 'r', linewidth = lw,                 )

for i, j in product(range(2), range(2)):

    ax[i, j].set_xlim(170, 800)
    ax[i, j].set_ylim(0.84, 2.4)
    ax[i, j].set_ylim(bottom = 0.0)

    ax[i, j].xaxis.set_major_locator(MultipleLocator(200))
    ax[i, j].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[0, 0].set_ylabel(r'Flux, [W / m$^2$ / nm]', fontsize = 20)
ax[1, 0].set_ylabel(r'Flux, [W / m$^2$ / nm]', fontsize = 20)

ax[1, 0].set_xlabel('Wavelength, [nm]', fontsize = 20)
ax[1, 1].set_xlabel('Wavelength, [nm]', fontsize = 20)

leg0 = ax[0, 0].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 15})

auxplt.savepdf('var/datom_old_new_comp_fac')
