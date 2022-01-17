import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

import importlib
import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)

if len(sys.argv) == 1: smoothing = '1'
if len(sys.argv) == 2: smoothing = sys.argv[1]

prefix0 = paths.it0f
prefix1 = paths.it1f

#wvl, C11h = nessy.read_spec(prefix1 + 'contrasts/11/C', wvl1 = 1005, wvl2 = 10000)
#wvl, P11h = nessy.read_spec(prefix1 + 'contrasts/11/P', wvl1 = 1005, wvl2 = 10000)
#wvl, S11h = nessy.read_spec(prefix1 + 'contrasts/11/S', wvl1 = 1005, wvl2 = 10000)
#wvl, C30h = nessy.read_spec(prefix1 + 'contrasts/30/C', wvl1 = 1005, wvl2 = 10000)
#wvl, P30h = nessy.read_spec(prefix1 + 'contrasts/30/P', wvl1 = 1005, wvl2 = 10000)
#wvl, S30h = nessy.read_spec(prefix1 + 'contrasts/30/S', wvl1 = 1005, wvl2 = 10000)
#wvl, Cleh = nessy.read_spec(prefix0 + 'contrasts/30/C', wvl1 = 1005, wvl2 = 10000)
#wvl, Pleh = nessy.read_spec(prefix0 + 'contrasts/30/P', wvl1 = 1005, wvl2 = 10000)
#wvl, Sleh = nessy.read_spec(prefix0 + 'contrasts/30/S', wvl1 = 1005, wvl2 = 10000)

#wvl = wvl / 10.0

#wvls, C11 = spec.running_mean(wvl, C11h, float(smoothing))
#wvls, P11 = spec.running_mean(wvl, P11h, float(smoothing))
#wvls, S11 = spec.running_mean(wvl, S11h, float(smoothing))
#wvls, C30 = spec.running_mean(wvl, C30h, float(smoothing))
#wvls, P30 = spec.running_mean(wvl, P30h, float(smoothing))
#wvls, S30 = spec.running_mean(wvl, S30h, float(smoothing))
#wvls, Cle = spec.running_mean(wvl, Cleh, float(smoothing))
#wvls, Ple = spec.running_mean(wvl, Pleh, float(smoothing))
#wvls, Sle = spec.running_mean(wvl, Sleh, float(smoothing))

#np.savez(paths.npz + 'contr_' + smoothing, w = wvls, fc11 = C11,\
#                                                     fp11 = P11,\
#                                                     fs11 = S11,\
#                                                     fc30 = C30,\
#                                                     fp30 = P30,\
#                                                     fs30 = S30,\
#                                                     fcle = Cle,\
#                                                     fple = Ple,\
#                                                     fsle = Sle,)

contr = np.load(paths.npz + 'contr_' + smoothing + '.npz')

w =   contr['w']

C11 = contr['fc11']
P11 = contr['fp11']
S11 = contr['fs11']
C30 = contr['fc30']
P30 = contr['fp30']
S30 = contr['fs30']
Cle = contr['fcle']
Ple = contr['fple']
Sle = contr['fsle']

PC11 = P11 - C11
PC30 = P30 - C30
PCle = Ple - Cle

SC11 = S11 - C11
SC30 = S30 - C30
SCle = Sle - Cle

RDP1130 = (PC11 - PC30) * 100.0 / PC30
RDS1130 = (SC11 - SC30) * 100.0 / SC30

RDPle30 = (PCle - PC30) * 100.0 / PC30
RDSle30 = (SCle - SC30) * 100.0 / SC30

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

#pltaux.figpar()

fig.tight_layout()

ax[0].plot(w, np.zeros(len(PC11)), 'k--')

ax[0].set_xlim(90, 1010)

if smoothing == '1':

    ls = ':'; lw = 1

if smoothing == '10':

    ls = '--'; lw = 2

ax[0].plot(w, PC30, color = 'r', linewidth = lw, label = 'Facula')
ax[0].plot(w, SC30, color = 'b', linewidth = lw, label = 'Spot')

ax[1].plot(w, abs(RDS1130), color = 'b', linewidth = lw, label = r'NLTE$_{11}$ vs. NLTE$_{30}$')
ax[1].plot(w, abs(RDP1130), color = 'r', linewidth = lw)

ax[1].plot(w, abs(RDSle30), color = 'b', linewidth = lw, linestyle = ls, label = r'LTE vs. NLTE$_{30}$')
ax[1].plot(w, abs(RDPle30), color = 'r', linewidth = lw, linestyle = ls)

ax[1].plot(w, np.ones(len(RDS1130)), 'k')

ax[1].set_yscale('log')

ax[1].set_xlim(90, 1010)
ax[1].set_ylim(1e-6, 1e+6)

for i in range(0, 2):

    ax[i].xaxis.set_major_locator(MultipleLocator(100))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))

ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[1].yaxis.set_major_locator(LogLocator(10))
ax[1].yaxis.set_minor_locator(LogLocator(base = 10.0, subs = (2, 3, 4, 5, 6, 7, 8, 9)))

ax[1].set_xlabel('Wavelength, [nm]',                              fontsize = 12.5)
ax[0].set_ylabel('Contrast (30 NLTE elements), [W / m$^2$ / nm]', fontsize = 12.5)
ax[1].set_ylabel('Relative difference, [%]',                      fontsize = 12.5)

leg0 = ax[0].legend(framealpha = 1, loc = 3, handletextpad = 1, prop = {'size': 20.5})
leg1 = ax[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('contr_' + smoothing)

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

#pltaux.figpar()

fig.tight_layout()

ax[0].plot(w, PC11 - PC30, color = 'm', label = r'$(P - C)_{11} - (P - C)_{30}$')
ax[0].plot(w, PC30,        color = 'y', label = r'$(P - C)_{30}$')

ax[0].set_ylabel('Difference, [W / m$^2$ / nm]', fontsize = 12.5)

ax[0].set_xlim(90, 1010)

ax[1].plot(w, PC11 - PC30, color = 'm', label = r'$(P - C)_{11} - (P - C)_{30}$')
ax[1].plot(w, PC30,        color = 'y', label = r'$(P - C)_{30}$')

ax[1].set_xlim(90, 1010)

ax[1].set_yscale('log')

ax[1].set_xlabel('Wavelength, [nm]',             fontsize = 12.5)
ax[1].set_ylabel('Difference, [W / m$^2$ / nm]', fontsize = 12.5)

leg = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('numbers_' + smoothing)
