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

wvl, flu1 = nessy.read_spec(paths.it1f + 'runtime/def',                                        wvl1 = 1005, wvl2 = 10000)
wvl, flu2 = nessy.read_spec(paths.it1f + 'runtime/def_atlodf_cut_off',                         wvl1 = 1005, wvl2 = 10000)
wvl, flu3 = nessy.read_spec(paths.it1f + 'runtime/nlte11_fgrid_red_5_20_eps_2',                wvl1 = 1005, wvl2 = 10000)
wvl, flu4 = nessy.read_spec(paths.it1f + 'runtime/nlte11_fgrid_red_5_20_eps_2_atlodf_cut_off', wvl1 = 1005, wvl2 = 10000)

wvl = wvl / 10.0

wvls, flus1 = spec.mean_within_delta(wvl, flu1, 1.0)
wvls, flus2 = spec.mean_within_delta(wvl, flu2, 1.0)
wvls, flus3 = spec.mean_within_delta(wvl, flu3, 1.0)
wvls, flus4 = spec.mean_within_delta(wvl, flu4, 1.0)

outname = 'runtime_specomp_atlodf_cut_off'

np.savez(paths.npz + outname, w = wvls, f1 = flus1,
                                        f2 = flus2,
                                        f3 = flus3,
                                        f4 = flus4)

specs = np.load(paths.npz + outname + '.npz')

w =  specs['w']
f1 = specs['f1']
f2 = specs['f2']
f3 = specs['f3']
f4 = specs['f4']

plt.close('all')

fig, ax = plt.subplots(nrows = 3, ncols = 1, figsize = (6.0, 6.0))

fig.tight_layout()

fig.suptitle('Best combination + ATLAS ODF impact', y = 1.01)

ax[0].plot(w, f1 / f3, color = 'k', label = 'def / best_combination (???)')
ax[1].plot(w, f1 / f2, color = 'k', label = 'def / def_atlodf')
ax[2].plot(w, f3 / f4, color = 'k', label = 'best_combination / best_combination_atlodf')

for j in [0, 1, 2]:

    ax[j].set_xlim(100, 1000)

    ax[j].xaxis.set_minor_locator(AutoMinorLocator(10))
    ax[j].xaxis.set_major_locator(MultipleLocator(100))

    ax[j].set_ylabel('Ratio')

#ax.set_ylim(0.95, 1.08)

    ax[j].axhline(y = 1.0, color = 'k', linestyle = '--', linewidth = 0.5)

ax[2].set_xlabel('Wavelength, [nm]')
#ax[0].set_ylabel('Flux, [W / m$^2$ / nm]')

#ax[0].set_yscale('log')

leg0 = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size':7.5})
leg1 = ax[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size':7.5})
leg2 = ax[2].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size':7.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg1.legendHandles: obj.set_linewidth(3.0)
for obj in leg2.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('runtime_spec_comp_atlodf_cut_off')
