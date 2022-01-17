import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)

wvl, flu1 = nessy.read_spec(paths.it1f + 'atlodf/old/kur/base/',            wvl1 = 1005, wvl2 = 10000)
wvl, flu2 = nessy.read_spec(paths.it1f + 'atlodf/new/kur/nopre/',           wvl1 = 1005, wvl2 = 10000)
wvl, flu3 = nessy.read_spec(paths.it1f + 'atlodf/new/kur/nopre_noH_nopif/', wvl1 = 1005, wvl2 = 10000)

wvl = wvl / 10.0

wvls, flus1 = spec.mean_within_delta(wvl, flu1, 1.0)
wvls, flus2 = spec.mean_within_delta(wvl, flu2, 1.0)
wvls, flus3 = spec.mean_within_delta(wvl, flu3, 1.0)

outname = 'spec_comp_kurucz_QS_atl_vs_fioss_opac'

np.savez(paths.npz + outname, w = wvls, f1 = flus1, f2 = flus2, f3 = flus3)

specs = np.load(paths.npz + outname + '.npz')

w =  specs['w']
f1 = specs['f1']
f2 = specs['f2']
f3 = specs['f3']

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6.0, 6.0))

fig.tight_layout()

for i in range(len(ax)):

    ax[i].set_xlim(100, 500)

    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))
    ax[i].xaxis.set_major_locator(MultipleLocator(100))

#ax[0].set_ylim(1e-4, 1e-1)
#ax[1].set_ylim(-30, 10)

#ax[0].set_xlim(100, 500)

ax[0].plot(w, f1, label = 'FIOSS',                       color = 'k')
ax[0].plot(w, f2, label = 'ATLAS (no PS)',               color = 'b', linewidth = 0.5)
ax[0].plot(w, f3, label = 'ATLAS (no PS, no H, no PIF)', color = 'r', linewidth = 0.5)

ax[1].plot(w, f1 / f2, color = 'b')
ax[1].plot(w, f1 / f3, color = 'r')

ax[1].axhline(y = 1.0, color = 'k', linestyle = '--', linewidth = 0.5)

ax[1].set_xlabel('Wavelength, [nm]')
ax[0].set_ylabel('Flux, [W / m$^2$ / nm]')
ax[1].set_ylabel('Ratio')

ax[0].set_yscale('log')

leg = ax[0].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 10.5})

for obj in leg.legendHandles:

    obj.set_linewidth(3.0)

auxplt.savepdf('spec_comp_kurucz_QS_atl_vs_fioss_opac')
