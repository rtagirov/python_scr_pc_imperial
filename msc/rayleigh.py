import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib

import pylab as plb

#import sys

import paths;      importlib.reload(paths)
import spec;       importlib.reload(spec)
import nessy_spec; importlib.reload(nessy_spec)

#if len(sys.argv) == 1 or len(sys.argv) == 2:

#   print('\nTwo directory names have to be provided. Abort.'); sys.exit()

#name1 = sys.argv[1]
#name2 = sys.argv[2]

#odf_it = '0'

#if len(sys.argv) > 3: odf_it = sys.argv[3]

#if odf_it == '0': prefix = paths.it0f
#if odf_it == '1': prefix = paths.it1f
#if odf_it == '2': prefix = paths.it2f
#if odf_it == '3': prefix = paths.it3f

#wvl_o, flu_o = nessy_spec.read(paths.it3f + 'OLD_SPEC',    wvl1 = 905, wvl2 = 1700)
#wvl_n, flu_n = nessy_spec.read(paths.it3f + 'no_rayleigh', wvl1 = 705, wvl2 = 1700)
#wvl_r, flu_r = nessy_spec.read(paths.it3f + 'rayleigh',    wvl1 = 705, wvl2 = 1700)

#wvl_o = wvl_o / 10.0
#wvl_n = wvl_n / 10.0
#wvl_r = wvl_r / 10.0

#wvls_o, flus_o = spec.running_mean(wvl_o, flu_o, 0.1)
#wvls_n, flus_n = spec.running_mean(wvl_n, flu_n, 0.1)
#wvls_r, flus_r = spec.running_mean(wvl_r, flu_r, 0.1)

outname = 'rayleigh'

#np.savez(paths.npz + outname, wo = wvls_o, \
#                              wn = wvls_n, \
#                              wr = wvls_r, \
#                              fo = flus_o, \
#                              fr = flus_r, \
#                              fn = flus_n)

specs = np.load(paths.npz + outname + '.npz')

wo = specs['wo']
wn = specs['wn']
wr = specs['wr']
fo = specs['fo']
fn = specs['fn']
fr = specs['fr']

ww = np.loadtxt('/mnt/SSD/sim/idl/output/whi_uv_unc.dat', usecols = [0])
fl = np.loadtxt('/mnt/SSD/sim/idl/output/whi_uv_unc.dat', usecols = [1])
fu = np.loadtxt('/mnt/SSD/sim/idl/output/whi_uv_unc.dat', usecols = [3])

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

for i in range(len(ax)):

    ax[i].set_xlim(70, 170)
    ax[i].set_xlim(70, 170)

    ax[i].xaxis.set_minor_locator(AutoMinorLocator(5))
    ax[i].xaxis.set_major_locator(MultipleLocator(10))

ax[0].plot(wo, fo, color = 'g', linewidth = 0.5, label = 'Old')
ax[0].plot(wr, fr, color = 'b',                  label = 'Rayleigh')
ax[0].plot(wn, fn, color = 'r', linewidth = 0.5, label = 'No_Rayleigh')

ax[0].fill_between(ww, fl, fu, color = 'grey', alpha = 0.5, label = 'WHI')

ax[1].plot(wr, fr / fn)

ax[1].set_xlabel('Wavelength, [nm]')
ax[0].set_ylabel('Flux, [W / m$^2$ / nm]')
ax[1].set_ylabel('Rayleigh / No_Rayleigh')

ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[0].set_yscale('log')

leg = ax[0].legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 20})

for handle in leg.legendHandles: handle.set_visible(False)

lo, lr, ln, lw = leg.get_texts()

plb.setp(lo, color = 'g')
plb.setp(lr, color = 'b')
plb.setp(ln, color = 'r')
plb.setp(lw, color = 'grey')

plt.show()
