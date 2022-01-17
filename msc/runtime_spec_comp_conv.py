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

#wvl, flu1 =     nessy.read_spec(paths.it1f + 'runtime/def',                      wvl1 = 1005, wvl2 = 4000)
#wvl, flu2 =     nessy.read_spec(paths.it1f + 'runtime/nlte11',                   wvl1 = 1005, wvl2 = 4000)
#wvl, flu3 =     nessy.read_spec(paths.it1f + 'runtime/nlte11_fgrid_red_5',       wvl1 = 1005, wvl2 = 4000)
#wvl, flu4 =     nessy.read_spec(paths.it1f + 'runtime/nlte11_fgrid_red_5_eps_3', wvl1 = 1005, wvl2 = 4000)
#wvl, flu5 =     nessy.read_spec(paths.it1f + 'runtime/nlte11_fgrid_red_5_eps_2', wvl1 = 1005, wvl2 = 4000)
#wvl, flu6 =     nessy.read_spec(paths.it1f + 'runtime/nlte11_fgrid_red_5_eps_1', wvl1 = 1005, wvl2 = 4000)

#wvl = wvl / 10.0

#wvls, flus1 =     spec.mean_within_delta(wvl, flu1, 1.0)
#wvls, flus2 =     spec.mean_within_delta(wvl, flu2, 1.0)
#wvls, flus3 =     spec.mean_within_delta(wvl, flu3, 1.0)
#wvls, flus4 =     spec.mean_within_delta(wvl, flu4, 1.0)
#wvls, flus5 =     spec.mean_within_delta(wvl, flu5, 1.0)
#wvls, flus6 =     spec.mean_within_delta(wvl, flu6, 1.0)

outname = 'runtime_specomp_conv'

#np.savez(paths.npz + outname, w = wvls, f1 = flus1,
#                                        f2 = flus2,
#                                        f3 = flus3,
#                                        f4 = flus4,
#                                        f5 = flus5,
#                                        f6 = flus6)

specs = np.load(paths.npz + outname + '.npz')

w =  specs['w']
f1 = specs['f1']
f2 = specs['f2']
f3 = specs['f3']
f4 = specs['f4']
f5 = specs['f5']
f6 = specs['f6']

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6.0, 6.0))

fig.tight_layout()

fig.suptitle('Convergence timing', y = 1.01)

ax[0].plot(w, f1, label = '(31m 00s) default',                 color = 'k')
ax[0].plot(w, f2, label = '(24m 15s) nlte 11',                 color = 'b',      linewidth = 0.5)
ax[0].plot(w, f3, label = '(12m 24s) nlte 11, cont 5, conv 4', color = 'r',      linewidth = 0.5)
ax[0].plot(w, f4, label = '(10m 32s) nlte 11, cont 5, conv 3', color = 'g',      linewidth = 0.25)
ax[0].plot(w, f5, label = '(07m 19s) nlte 11, cont 5, conv 2', color = 'purple', linewidth = 0.25)
ax[0].plot(w, f6, label = '(05m 22s) nlte 11, cont 5, conv 1', color = 'orange', linewidth = 0.25)

ax[1].plot(w, f1 / f2, color = 'b',      linewidth = 1)
ax[1].plot(w, f1 / f3, color = 'r',      linewidth = 1)
ax[1].plot(w, f1 / f4, color = 'g',      linewidth = 0.5)
ax[1].plot(w, f1 / f5, color = 'purple', linewidth = 0.5)
ax[1].plot(w, f1 / f6, color = 'orange', linewidth = 0.5)

for j in [0, 1]:

    ax[j].set_xlim(100, 400)

    ax[j].xaxis.set_minor_locator(AutoMinorLocator(10))
    ax[j].xaxis.set_major_locator(MultipleLocator(100))

#ax.set_ylim(0.95, 1.08)

ax[1].axhline(y = 1.0, color = 'k', linestyle = '--', linewidth = 0.5)

ax[1].set_xlabel('Wavelength, [nm]')
ax[0].set_ylabel('Flux, [W / m$^2$ / nm]')
ax[1].set_ylabel('Ratio')

#ax[0].set_yscale('log')

leg = ax[0].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size':7.5})

for obj in leg.legendHandles:

    obj.set_linewidth(3.0)

auxplt.savepdf('runtime_spec_comp_conv')
