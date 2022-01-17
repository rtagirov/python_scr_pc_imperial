import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

from scipy import interpolate 

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

#wr, tr, ssi_r = varfunc.readssi(paths.satrnm + 'ssi_rnm.dat')

#np.savez(paths.npz + 'ssi_r', wr = wr, tr = tr, ssi_r = ssi_r)

ssir = np.load(paths.npz + 'ssi_r.npz')

wr = ssir['wr']
tr = ssir['tr']

tt, tsi_r, tsi_l, tsi_u, s = np.loadtxt(paths.satrnm + 'tsi_rnm.dat', skiprows = 24, unpack = True)

var_nlt = scipy.io.readsav(paths.inp + 'satire_nlte.sav')
var_lte = scipy.io.readsav(paths.inp + 'satire_lte.sav')

wff = var_nlt['wl']

w = np.zeros(len(wff))

for i in range(len(wff)):

    w[i] = wff[i, 0]

t = varfunc.time(var_nlt['date'])

tsi_nlt = var_nlt['tsi']
tsi_lte = var_lte['tsi']

t1 = 2010.42
t2 = 2011.25
#t3 = 2014.02; t4 = 2014.85
t3 = 2014.925; dt = t2 - t1; t4 = t3 + dt

mvar_nlt = varfunc.mmmvar(var_nlt['ssi'], t, t1, t2, t3, t4)
mvar_lte = varfunc.mmmvar(var_lte['ssi'], t, t1, t2, t3, t4)

mvar_sir = varfunc.mmmvar(ssir['ssi_r'], tr, t1, t2, t3, t4)

#mvar_sir_i = interpolate.interp1d(wr, mvar_sir)(w)
mvar_sir_i = np.interp(w, wr, mvar_sir)

mvar_rat_1 = mvar_nlt / mvar_lte
mvar_rat_2 = mvar_sir_i / mvar_lte

idx1 = np.where(mvar_nlt < 0.0)[0]
idx2 = np.where(mvar_lte < 0.0)[0]
idx3 = np.where(mvar_rat_1 < 0.0)[0]
idx4 = np.where(mvar_rat_2 < 0.0)[0]

w1 = w[idx1]
w2 = w[idx2]
w3 = w[idx3]
w4 = w[idx4]

mvar_nlt_neg = mvar_nlt[idx1]
mvar_lte_neg = mvar_lte[idx2]

mvar_rat_1_neg = mvar_rat_1[idx3]
mvar_rat_2_neg = mvar_rat_2[idx4]

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12, 6))

bbox = dict(boxstyle = 'round', ec = (1.0, 0.5, 0.5), fc = (1.0, 0.8, 0.8),)

fig.tight_layout()

auxplt.figpar(3, 3, 15)

plt.subplots_adjust(hspace = 0.3)

#ax[0].fill_between(np.array([t1, t2]), 1358.5, 1362.5, facecolor = 'gray', linewidth = 0, alpha = 0.5)
#ax[0].fill_between(np.array([t3, t4]), 1358.5, 1362.5, facecolor = 'gray', linewidth = 0, alpha = 0.5)

#ax[0].axvline(x = t1, linestyle = '--', color = 'k')
ax[0].axvline(x = t2, linestyle = '--', color = 'k')
ax[0].axvline(x = t3, linestyle = '--', color = 'k')
ax[0].axvline(x = t4, linestyle = '--', color = 'k')

ax[0].fill_between(varfunc.time(tt), tsi_l + 0.15, tsi_u + 0.15, color = 'gray', label = 'RENORM', alpha = 0.5)

ax[0].plot(t, tsi_lte, color = 'k', linewidth = 1.0, label = 'LTE')
ax[0].plot(t, tsi_nlt, color = 'r', linewidth = 0.5, label = 'NLTE')

ax[0].text(2010.725, 1359, 'min', bbox = bbox)
ax[0].text(2015.225, 1359, 'max', bbox = bbox)

ax[1].semilogy(w[np.where(w >= 150.0)], mvar_lte[np.where(w >= 150.0)], color = 'k', label = 'ATLAS9, LTE, U99')
#ax[1].semilogy(w, mvar_sir_i, color = 'b', label = r'$\mathrm{LTE}_\mathrm{corr}$')
ax[1].semilogy(w, mvar_sir_i, color = 'gray', label = 'ATLAS9, LTE, U99, emp. corr.')

ax[1].semilogy(w, mvar_nlt, color = 'r', label = 'NESSY, NLTE, FAL99')

ax[1].scatter(w1, -mvar_nlt_neg, color = 'r', s = 10.5)
ax[1].scatter(w2, -mvar_lte_neg, color = 'k', s = 10.5)

ax[0].set_xlabel('Year')
ax[1].set_xlabel('Wavelength, [nm]')

ax[1].set_xlim(100, 500)
ax[1].set_ylim(5e-4, 2e+2)

ax[0].set_ylim(1358.5, 1362.5)

ax[0].set_xlim(min(t), max(t))

ax[0].ticklabel_format(useOffset = False)

ax[0].xaxis.set_minor_locator(AutoMinorLocator(4))
ax[0].yaxis.set_minor_locator(AutoMinorLocator(4))
ax[1].xaxis.set_minor_locator(AutoMinorLocator(5))

ax[1].yaxis.set_major_locator(LogLocator(10, numticks = 7))
ax[1].yaxis.set_minor_locator(LogLocator(base = 10.0, subs = (2, 3, 4, 5, 6, 7, 8, 9), numticks = 8))

leg1 = ax[1].legend(framealpha = 1, loc = 3, handletextpad = 1, prop = {'size': 17.0})

for obj in leg1.legendHandles: obj.set_linewidth(3.0)

ax[0].set_ylabel(r'TSI, $\left[\mathrm{W} / \mathrm{m}^2\right]$')
ax[1].set_ylabel(r'$\Delta \mathrm{SSI}, \mathrm{[\%]}$')

auxplt.savepdf('var/tsi_var_spec')
