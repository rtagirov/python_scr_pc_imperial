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

wvl, I1 = nessy.read_spec(paths.it1f + 'atlodf/old/fal/base/',                            )
wvl, I2 = nessy.read_spec(paths.it1f + 'atlodf/old/fal/base_no_opac_in_chromo/',          )
wvl, I3 = nessy.read_spec(paths.it1f + 'atlodf/new/fal/nopre_no_opac_in_chromo/',         )
wvl, I4 = nessy.read_spec(paths.it1f + 'atlodf/new/fal/nopre_noH_nopif_no_opac_in_chromo/')

wvl, idx1, h1 = nessy.read_tau(paths.it1f + 'atlodf/old/fal/base/',                            )
wvl, idx2, h2 = nessy.read_tau(paths.it1f + 'atlodf/old/fal/base_no_opac_in_chromo/',          )
wvl, idx3, h3 = nessy.read_tau(paths.it1f + 'atlodf/new/fal/nopre_no_opac_in_chromo/',         )
wvl, idx4, h4 = nessy.read_tau(paths.it1f + 'atlodf/new/fal/nopre_noH_nopif_no_opac_in_chromo/')

wvl = wvl / 10.0

temp = np.loadtxt('/mnt/SSD/sim/nessy/inp/atm/fal/FAL99_C', usecols = [1])

t1 = temp[idx1]
t2 = temp[idx2]
t3 = temp[idx3]
t4 = temp[idx4]

wvls, ts1 = spec.weigh_within_delta(wvl, t1, I1, 1.0)
wvls, ts2 = spec.weigh_within_delta(wvl, t2, I2, 1.0)
wvls, ts3 = spec.weigh_within_delta(wvl, t3, I3, 1.0)
wvls, ts4 = spec.weigh_within_delta(wvl, t4, I4, 1.0)

outname = 'atl_odf_ft'

np.savez(paths.npz + outname, w = wvls, t1 = ts1, t2 = ts2, t3 = ts3, t4 = ts4)

specs = np.load(paths.npz + outname + '.npz')

w =  specs['w']
t1 = specs['t1']
t2 = specs['t2']
t3 = specs['t3']
t4 = specs['t4']

plt.close('all')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (6.0, 3.0))

fig.tight_layout()

#ax.set_xlim(210, 500)
#ax.set_xlim(300, 320)
ax.set_xlim(100, 500)

ax.xaxis.set_minor_locator(AutoMinorLocator(10))
ax.xaxis.set_major_locator(MultipleLocator(100))

#ax.set_ylim(0.0, 800)

ax.plot(w, t1, label = 'FIOSS',                                 color = 'k', linewidth = 0.5)
ax.plot(w, t2, label = 'FIOSS (no opac in chromosphere)',       color = 'b', linewidth = 0.5)
ax.plot(w, t3, label = 'ATLAS (no preselection)',               color = 'r', linewidth = 0.5)
ax.plot(w, t4, label = 'ATLAS (no preselection, no H, no PIF)', color = 'g', linewidth = 0.5)

#ax.axhline(y = 1.0, color = 'k', linestyle = '--', linewidth = 0.5)

ax.set_xlabel('Wavelength, nm')
ax.set_ylabel('Temperature, K')

#ax.set_yscale('log')

leg = ax.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size':7.5})

for obj in leg.legendHandles:

    obj.set_linewidth(3.0)

auxplt.savepdf('spec_comp_no_opac_in_chromo_ft')
