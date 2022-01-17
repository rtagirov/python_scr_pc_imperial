import matplotlib.pyplot as plt
import numpy as np

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

import itertools
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

wvl, Iwl_max, Iw, Iwu_max = np.loadtxt(paths.inp + 'whi_unc_max.dat', unpack = True)
wvl, Iwl_min, Iw, Iwu_min = np.loadtxt(paths.inp + 'whi_unc_min.dat', unpack = True)

waf, Faf = varfunc.spec('lte', 'qs_1.dat')

nw_conv = scipy.io.readsav('/mnt/SSD/sim/idl/src/var/neswhiconv.sav')
aw_conv = scipy.io.readsav('/mnt/SSD/sim/idl/src/var/atlwhiconv.sav')

Fn = nw_conv['inf']
Fa = aw_conv['Fa']

#idx_wac = np.where((waf > min(wvl)) & (waf < max(wvl)))

#wac = waf[idx_wac]

#Fac = Faf[idx_wac]

#Fa = np.interp(wvl, wac, Fac)

lfs = 30

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (20.0, 15.0))

fig.tight_layout()

plt.subplots_adjust(wspace = 0.20, hspace = 0.1)

ax[0, 0].fill_between(wvl, Iwl_max, Iwu_max, color = 'grey', label = 'SIRS WHI', alpha = 0.3)
ax[0, 0].fill_between(wvl, Iwl_min, Iwu_min, color = 'grey', alpha = 1.0)

ax[0, 0].plot(wvl, Fa, color = 'b', label = 'ATLAS')
ax[0, 0].plot(wvl, Fn, color = 'r', label = 'NESSY')

ax[0, 0].set_xlim(170, 400)
ax[0, 0].set_ylim(0.0, 1.7)

ax[0, 1].fill_between(wvl, Iwl_max, Iwu_max, color = 'grey', alpha = 0.3)
ax[0, 1].fill_between(wvl, Iwl_min, Iwu_min, color = 'grey', alpha = 1.0)

ax[0, 1].plot(wvl, Fa, color = 'b')
ax[0, 1].plot(wvl, Fn, color = 'r')

ax[0, 1].set_xlim(400, 1000)
ax[0, 1].set_ylim(0.65, 2.2)

#ax[1, 0].fill_between(wvl, -7, 7, color = 'grey', alpha = 0.3)
#ax[1, 0].fill_between(wvl, -3.5, 3.5, color = 'grey', alpha = 1.0)

ax[1, 0].plot(wvl, (Fa - Iw) * 100.0 / Iw, color = 'b')
ax[1, 0].plot(wvl, (Fn - Iw) * 100.0 / Iw, color = 'r')

ax[1, 0].set_xlim(170, 400)
ax[1, 0].set_ylim(-75, 250)

#ax[1, 1].fill_between(wvl, , (Iwu_max - Iw) * 100.0 / Iw, color = 'grey', alpha = 0.3)
#ax[1, 1].fill_between(wvl, , (Iwu_min - Iw) * 100.0 / Iw, color = 'grey', alpha = 1.0)

ax[1, 1].plot(wvl, (Fa - Iw) * 100.0 / Iw, color = 'b')
ax[1, 1].plot(wvl, (Fn - Iw) * 100.0 / Iw, color = 'r')

ax[1, 1].set_xlim(400, 1000)
ax[1, 1].set_ylim(-10, 6)

for i, j in itertools.product(range(2), range(2)):

    ax[i, j].tick_params(axis = 'both', labelsize = 19)
    ax[i, j].tick_params(axis = 'both', labelsize = 19)

leg = ax[0, 0].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 25.5})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

ax[0, 0].set_ylabel(r'Flux, [W / m$^2$ / nm]', fontsize = lfs)

ax[1, 0].set_ylabel('Deviation, [%]', fontsize = lfs)
ax[1, 0].set_xlabel('Wavelengh, [nm]', fontsize = lfs)
ax[1, 1].set_xlabel('Wavelengh, [nm]', fontsize = lfs)

auxplt.savepdf('var/whi_comp')
