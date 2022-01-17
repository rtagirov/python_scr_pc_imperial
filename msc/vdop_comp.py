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

wvl0, chi0 = np.loadtxt('600_620v0.dat', unpack = True)
wvl1, chi1 = np.loadtxt('600_620v1.dat', unpack = True)
wvl2, chi2 = np.loadtxt('600_620v2.dat', unpack = True)

wvl0 = wvl0[ : len(wvl0) - 1]

#chi0 = chi0[ : len(wvl0) - 1]
#chi1 = chi1[ : len(wvl0) - 1]
#chi2 = chi2[ : len(wvl0) - 1]

delta02 = 0.2
delta1 =  1.0
delta10 = 10

wvls02 = np.linspace(600.1, 619.9, int(20.0 / delta02))
wvls1  = np.linspace(600.5, 619.5, int(20.0 / delta1))
wvls10 = np.linspace(605.0, 615.0, int(20.0 / delta10))

wvls = wvls1
delta = delta1

chi0s = spec.sort_within_delta_over_grid(chi0, wvl0, delta, wvls)
chi1s = spec.sort_within_delta_over_grid(chi1, wvl0, delta, wvls)
chi2s = spec.sort_within_delta_over_grid(chi2, wvl0, delta, wvls)

chi0m = spec.mean_within_delta_over_grid(chi0s, wvl0, delta10, wvls10)
chi1m = spec.mean_within_delta_over_grid(chi1s, wvl0, delta10, wvls10)
chi2m = spec.mean_within_delta_over_grid(chi2s, wvl0, delta10, wvls10)

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6.0, 6.0))

ax[0].plot(wvl0, chi0s, color = 'k', label = 'vdop = 0')
ax[0].plot(wvl0, chi1s, color = 'r', label = 'vdop = 1 km / s', linewidth = 0.5)
ax[0].plot(wvl0, chi2s, color = 'g', label = 'vdop = 2 km / s', linewidth = 0.5)

ax[1].step(wvls10 + delta10 / 2.0, chi0m, color = 'k', label = 'vdop = 0')
ax[1].step(wvls10 + delta10 / 2.0, chi1m, color = 'r', label = 'vdop = 1 km / s', linewidth = 0.5)
ax[1].step(wvls10 + delta10 / 2.0, chi2m, color = 'g', label = 'vdop = 2 km / s', linewidth = 0.5)

ax[0].set_xlim(600, 620)
ax[1].set_xlim(600, 620)
#ax[0].set_xlim(600, 605)
#ax[1].set_xlim(600, 605)

ax[1].set_xlabel('Wavelength, nm')
ax[0].set_ylabel(r'Opacity, cm$^2$ / g')
ax[1].set_ylabel(r'Opacity, cm$^2$ / g')

ax[0].set_yscale('log')
ax[1].set_yscale('log')

ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[1].xaxis.set_minor_locator(AutoMinorLocator(5))
#ax[0].xaxis.set_major_locator(MultipleLocator(2))
#ax[1].xaxis.set_major_locator(MultipleLocator(2))
#ax[0].xaxis.set_minor_locator(AutoMinorLocator(2))
#ax[1].xaxis.set_minor_locator(AutoMinorLocator(2))
ax[0].xaxis.set_major_locator(MultipleLocator(5))
ax[1].xaxis.set_major_locator(MultipleLocator(5))
ax[0].xaxis.set_minor_locator(AutoMinorLocator(5))
ax[1].xaxis.set_minor_locator(AutoMinorLocator(5))

for i, elem in enumerate(wvls):

    ax[0].axvline(x = elem, color = 'b', linewidth = 0.5)

for i, elem in enumerate(wvls10):

    ax[1].axvline(x = elem, color = 'b', linewidth = 0.5)

leg = ax[0].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 7.5}, bbox_to_anchor=(0, 1.12))

for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('vdop_comp')
