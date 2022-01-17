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

if len(sys.argv) == 1 or len(sys.argv) == 2:

    sysaux.abort('Two directory names have to be provided.')

name1 = sys.argv[1]
name2 = sys.argv[2]

odf_it = '0'

if len(sys.argv) > 3: odf_it = sys.argv[3]

if odf_it == '0': prefix = paths.it0f
if odf_it == '1': prefix = paths.it1f
if odf_it == '2': prefix = paths.it2f
if odf_it == '3': prefix = paths.it3f

wvl, flu1 = nessy.read_spec(prefix + name1, wvl1 = 1005, wvl2 = 10000)
wvl, flu2 = nessy.read_spec(prefix + name2, wvl1 = 1005, wvl2 = 10000)

wvl = wvl / 10.0

wvls, flus1 = spec.mean_within_delta(wvl, flu1, 1.0)
wvls, flus2 = spec.mean_within_delta(wvl, flu2, 1.0)

outname = 'specomp'

np.savez(paths.npz + outname, w = wvls, f1 = flus1, f2 = flus2)

specs = np.load(paths.npz + outname + '.npz')

w =  specs['w']
f1 = specs['f1']
f2 = specs['f2']

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6.0, 6.0))

fig.tight_layout()

for i in range(len(ax)):

    ax[i].set_xlim(100, 1000)

    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))
    ax[i].xaxis.set_major_locator(MultipleLocator(100))

#ax[0].set_ylim(1e-4, 1e-1)
#ax[1].set_ylim(-30, 10)

ax[0].set_xlim(100, 500)

ax[0].plot(w, f1, label = 'FIOSS', color = 'k')
ax[0].plot(w, f2, label = 'ATLAS', color = 'r', linewidth = 0.5)

ax[1].plot(w, f1 / f2, color = 'purple')

ax[1].axhline(y = 1.0, color = 'k', linestyle = '--', linewidth = 0.5)

ax[1].set_xlabel('Wavelength, [nm]')
ax[0].set_ylabel('Flux, [W / m$^2$ / nm]')
ax[1].set_ylabel('FIOSS / ATLAS')

ax[0].set_yscale('log')

leg = ax[0].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 10.5})

for obj in leg.legendHandles:

    obj.set_linewidth(3.0)

auxplt.savepdf('specomp_mem_test')
