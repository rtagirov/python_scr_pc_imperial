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

#1 - NESSY LTE  FAL
#2 - NESSY NLTE FAL

#3 - NESSY LTE  KUR
#4 - NESSY NLTE KUR

#wn, C1h = nessy.read_spec(prefix0 + 'var/Q/fal/', wvl1 = 1005, wvl2 = 10000)
#wn, P1h = nessy.read_spec(prefix0 + 'var/F/fal/', wvl1 = 1005, wvl2 = 10000)
#wn, C2h = nessy.read_spec(prefix1 + 'var/Q/fal/', wvl1 = 1005, wvl2 = 10000)
#wn, P2h = nessy.read_spec(prefix1 + 'var/F/fal/', wvl1 = 1005, wvl2 = 10000)

#wn, C3h = nessy.read_spec(prefix0 + 'var/Q/kur/', wvl1 = 1005, wvl2 = 10000)
#wn, P3h = nessy.read_spec(prefix0 + 'var/F/kur/', wvl1 = 1005, wvl2 = 10000)
#wn, C4h = nessy.read_spec(prefix1 + 'var/Q/kur/', wvl1 = 1005, wvl2 = 10000)
#wn, P4h = nessy.read_spec(prefix1 + 'var/F/kur/', wvl1 = 1005, wvl2 = 10000)

#wn = wn / 10.0

#wns, C1 = spec.mean(wn, C1h, float(smoothing))
#wns, P1 = spec.mean(wn, P1h, float(smoothing))
#wns, C2 = spec.mean(wn, C2h, float(smoothing))
#wns, P2 = spec.mean(wn, P2h, float(smoothing))

#wns, C3 = spec.mean(wn, C3h, float(smoothing))
#wns, P3 = spec.mean(wn, P3h, float(smoothing))
#wns, C4 = spec.mean(wn, C4h, float(smoothing))
#wns, P4 = spec.mean(wn, P4h, float(smoothing))

#np.savez(paths.npz + 'contr_var_1_' + smoothing, w = wns, c1 = C1,\
#                                                          p1 = P1,\
#                                                          c2 = C2,\
#                                                          p2 = P2,\
#                                                          c3 = C3,\
#                                                          p3 = P3,\
#                                                          c4 = C4,\
#                                                          p4 = P4)

contr = np.load(paths.npz + 'contr_var_1_' + smoothing + '.npz')

w =  contr['w']

C1 = contr['c1']
P1 = contr['p1']
C2 = contr['c2']
P2 = contr['p2']
C3 = contr['c3']
P3 = contr['p3']
C4 = contr['c4']
P4 = contr['p4']

PC1 = P1 - C1
PC2 = P2 - C2
PC3 = P3 - C3
PC4 = P4 - C4

RDPC12 = (PC2 - PC1) * 100.0 / PC1
RDPC34 = (PC4 - PC3) * 100.0 / PC3

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

#pltaux.figpar()

fig.tight_layout()

fig.suptitle('NESSY: FAL vs. Kurucz facular contrasts', y = 1.01)

ax[0].plot(w, np.zeros(len(PC1)), 'k--')

ax[0].set_xlim(100, 1000)

if smoothing == '1':

    ls = ':'; lw = 1

if smoothing == '5':

    ls = ':'; lw = 1.5

if smoothing == '10':

    ls = '--'; lw = 2

ax[0].plot(w, PC1, color = 'b', linewidth = lw, label = 'LTE   (FAL)')
ax[0].plot(w, PC2, color = 'r', linewidth = lw, label = 'NLTE (FAL)')
ax[0].plot(w, PC3, color = 'b', linewidth = lw, label = 'LTE   (KUR)', linestyle = ls)
ax[0].plot(w, PC4, color = 'r', linewidth = lw, label = 'NLTE (KUR)',  linestyle = ls)

ax[1].plot(w, abs(RDPC12), color = 'k', linewidth = lw, label = 'FAL')                #, label = r'')
ax[1].plot(w, abs(RDPC34), color = 'k', linewidth = lw, label = 'KUR', linestyle = ls)#, label = r'')

#ax[1].plot(w, np.ones(len(RDPC12)), 'k')

ax[1].set_yscale('log')

ax[1].set_xlim(100, 1000)
ax[1].set_ylim(1e-0, 1e+2)

for i in range(0, 2):

    ax[i].xaxis.set_major_locator(MultipleLocator(100))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))

ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[1].yaxis.set_major_locator(LogLocator(10))
ax[1].yaxis.set_minor_locator(LogLocator(base = 10.0, subs = (2, 3, 4, 5, 6, 7, 8, 9)))

ax[0].set_ylabel('Facular Contrast, [W / m$^2$ / nm]', fontsize = 12.5)
ax[1].set_ylabel('(NLTE - LTE) / LTE, [%]',            fontsize = 12.5)
ax[1].set_xlabel('Wavelength, [nm]',                   fontsize = 12.5)

leg0 = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})
leg1 = ax[1].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 20.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('fcontr_nes_falkur')
