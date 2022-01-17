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

var_nlt = scipy.io.readsav(paths.inp + 'satire_nlte.sav')
#var_nss = scipy.io.readsav(paths.inp + 'satire_nlte_ss_fd.sav')
var_nss = scipy.io.readsav(paths.inp + 'satire_nlte_ss.sav')

wff = var_nlt['wl']

w = np.zeros(len(wff))

for i in range(len(wff)):

    w[i] = wff[i, 0]

t = varfunc.time(var_nlt['date'])

tsi_nlt = var_nlt['tsi']
tsi_nss = var_nss['tsi']

t1 = 2010.42
t2 = 2011.25
#t3 = 2014.02; t4 = 2014.85
t3 = 2014.925; dt = t2 - t1; t4 = t3 + dt

mvar_nlt = varfunc.mmmvar(var_nlt['ssi'], t, t1, t2, t3, t4)
mvar_nss = varfunc.mmmvar(var_nss['ssi'], t, t1, t2, t3, t4)

mvar_rat = mvar_nlt / mvar_nss

idx1 = np.where(mvar_nlt < 0.0)[0]
idx2 = np.where(mvar_nss < 0.0)[0]
idx3 = np.where(mvar_rat < 0.0)[0]

w1 = w[idx1]
w2 = w[idx2]
w3 = w[idx3]

mvar_nlt_neg = mvar_nlt[idx1]
mvar_nss_neg = mvar_nss[idx2]

mvar_rat_neg = mvar_rat[idx3]

plt.close('all')

fig, ax = plt.subplots(nrows = 3, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

ax[1].plot(w, np.ones(len(w)), color = 'k', linestyle = '--')

ax[0].semilogy(w, mvar_nss, color = 'g', label = 'RE-SCALED')
ax[0].semilogy(w, mvar_nlt, color = 'r', label = 'NLTE', linewidth = 0.5)

ax[0].scatter(w1, -mvar_nlt_neg, color = 'r', s = 10.5)
ax[0].scatter(w2, -mvar_nss_neg, color = 'k', s = 10.5)

#ax[1].semilogy(w, mvar_rat, color = 'r')
ax[1].plot(w, mvar_rat, color = 'r')

ax[1].scatter(w3, -mvar_rat_neg, color = 'r', s = 10.5)

ax[2].fill_between(np.array([t1, t2]), 1358.5, 1362.5, facecolor = 'wheat', linewidth = 0)
ax[2].fill_between(np.array([t3, t4]), 1358.5, 1362.5, facecolor = 'wheat', linewidth = 0)

ax[2].plot(t, tsi_nlt, color = 'r', linewidth = 0.5, label = 'NLTE')
ax[2].plot(t, tsi_nss, color = 'g', linewidth = 0.5, label = 'RE-SCALED')

ax[1].set_xlabel('Wavelength, [nm]')
ax[2].set_xlabel('Year')

ax[0].set_xlim(100,  1000)
ax[1].set_xlim(100,  1000)
ax[0].set_ylim(6e-4, 1e+2)
#ax[1].set_ylim(6e-2, 1e+1)
ax[1].set_ylim(0.75, 1.5)

ax[2].set_ylim(1358.5, 1362.5)

ax[2].set_xlim(min(t), max(t))

ax[2].ticklabel_format(useOffset = False)

ax[0].xaxis.set_minor_locator(AutoMinorLocator(10))
ax[1].xaxis.set_minor_locator(AutoMinorLocator(10))
ax[2].xaxis.set_minor_locator(AutoMinorLocator(10))

leg0 = ax[0].legend(framealpha = 1, loc = 3, handletextpad = 1, prop = {'size': 15.5})
#leg1 = ax[1].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 20.5})
leg2 = ax[2].legend(framealpha = 1, loc = 3, handletextpad = 1, prop = {'size': 15.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
#for obj in leg1.legendHandles: obj.set_linewidth(3.0)
for obj in leg2.legendHandles: obj.set_linewidth(3.0)

ax[0].set_ylabel('Flux difference, (MAX - MIN) / AVERAGE, [%]')
ax[1].set_ylabel('Ratio')
ax[2].set_ylabel(r'TSI, [W / m$^2$]')

auxplt.savepdf('var/satire_mmmvar_scaling')
