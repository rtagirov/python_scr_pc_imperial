import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

import importlib
import math
import sys
import scipy.io

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import phys;   importlib.reload(phys)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)

prefix0 = paths.it0f
prefix1 = paths.it1f

#1 - NESSY LTE
#2 - NESSY NLTE
#3 - ATLAS

#4 - NESSY LTE  FAL
#5 - NESSY NLTE FAL

#wn, Q1h = nessy.read_spec(prefix0 + 'var/Q/kur_cardsdef', wvl1 = 1005, wvl2 = 16000)
#wn, F1h = nessy.read_spec(prefix0 + 'var/F/kur_cardsdef', wvl1 = 1005, wvl2 = 16000)

#wn, Q2h = nessy.read_spec(prefix1 + 'var/Q/kur_cardsdef', wvl1 = 1005, wvl2 = 16000)
#wn, F2h = nessy.read_spec(prefix1 + 'var/F/kur_cardsdef', wvl1 = 1005, wvl2 = 16000)

#wn, Q4h = nessy.read_spec(prefix0 + 'var/Q/fal_cardsdef', wvl1 = 1005, wvl2 = 16000)
wn, Q4h = nessy.read_spec(prefix0 + 'VAR/Q_kur_old', wvl1 = 1005, wvl2 = 16000)
#wn, F4h = nessy.read_spec(prefix0 + 'var/F/fal_cardsdef', wvl1 = 1005, wvl2 = 16000)

#wn, Q5h = nessy.read_spec(prefix1 + 'var/Q/fal_cardsdef', wvl1 = 1005, wvl2 = 16000)
wn, Q5h = nessy.read_spec(prefix1 + 'VAR/Q_kur_old', wvl1 = 1005, wvl2 = 16000)
#wn, F5h = nessy.read_spec(prefix1 + 'var/F/fal_cardsdef', wvl1 = 1005, wvl2 = 16000)

wn = wn / 10.0

#wns, Q1 = spec.mean_within_delta(wn, Q1h, 10)
#wns, F1 = spec.mean_within_delta(wn, F1h, 10)

#wns, Q2 = spec.mean_within_delta(wn, Q2h, 10)
#wns, F2 = spec.mean_within_delta(wn, F2h, 10)

wns, Q4 = spec.mean_within_delta(wn, Q4h, 10)
#wns, F4 = spec.mean_within_delta(wn, F4h, 10)

wns, Q5 = spec.mean_within_delta(wn, Q5h, 10)
#wns, F5 = spec.mean_within_delta(wn, F5h, 10)

np.savez(paths.npz + 'spec_var_whi_fallte_C_new', w = wns,
#                                            q1 = Q1,\
#                                            f1 = F1,\
#                                            q2 = Q2,\
#                                            f2 = F2,\
                                            q4 = Q4,\
#                                            f4 = F4,\
                                            q5 = Q5)\
#                                            f5 = F5,)

contr = np.load(paths.npz + 'spec_var_whi_fallte_C_new.npz')

w = contr['w']

#Q1 = contr['q1']
#F1 = contr['f1']

#Q2 = contr['q2']
#F2 = contr['f2']

Q4 = contr['q4']
#F4 = contr['f4']

Q5 = contr['q5']
#F5 = contr['f5']

plt.close('all')

#fig, ax = plt.subplots(nrows = 3, ncols = 1, figsize = (6.0, 12.06))
fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6.0, 8.00))

bbox = dict(boxstyle = 'round', ec = (1.0, 0.5, 0.5), fc = (1.0, 0.8, 0.8),)

auxplt.figpar(3, 3, 20)

fig.tight_layout()

plt.subplots_adjust(hspace = 0.15)

ls = ':'
lw = 1.0

#ax[0].plot(w, Q1, color = 'm', linewidth = lw, label = 'NESSY (LTE, U99)')
#ax[0].plot(w, Q2, color = 'g', linewidth = lw, label = 'NESSY (NLTE, U99)')
ax[0].plot(w, Q4, color = 'k', linewidth = lw, label = 'NESSY (LTE, FAL99)', linestyle = '--')
ax[0].plot(w, Q5, color = 'k', linewidth = lw, label = 'NESSY (NLTE, FAL99)')

#ax[0].text(720, 1.62, 'Quiet Sun', bbox = bbox)

#ax[1].plot(w, F1, color = 'm', linewidth = lw, label = 'NESSY (LTE, U99)')
#ax[1].plot(w, F2, color = 'g', linewidth = lw, label = 'NESSY (NLTE, U99)')

#ax[0].plot(w, F4, color = 'r', linewidth = lw, label = 'NESSY (LTE, FAL99)', linestyle = '--')
#ax[0].plot(w, F5, color = 'r', linewidth = lw, label = 'NESSY (NLTE, FAL99)')

#ax[0].text(700, 1.80, 'Facula',   bbox = bbox)

#ax[1].plot(w, Q2 / Q1, label = 'Quiet Sun, U99', color = 'orange')
ax[1].plot(w, Q5 / Q4, label = 'Quiet Sun, FAL99', color = 'k')
#ax[1].plot(w, F2 / F1, label = 'Facula, U99', linestyle = '--', color = 'orange')

#ax[1].plot(w, F5 / F4, label = 'Facula, FAL99',  color = 'r')

ax[1].axhline(y = 1.0, color = 'k')

ax[1].set_xlim(170, 1600)
ax[1].set_ylim(0.84, 1.22)
ax[0].set_ylim(bottom = 0.0)

ax[1].xaxis.set_major_locator(MultipleLocator(300))
#ax[1].xaxis.set_minor_locator(AutoMinorLocator(5))

ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[1].set_ylabel('NLTE / LTE')

for i in [0, 1]:

    ax[i].set_xlim(170, 1600)

    ax[i].xaxis.set_major_locator(MultipleLocator(300))
#    ax[i].xaxis.set_minor_locator(AutoMinorLocator(5))

    ax[i].yaxis.set_minor_locator(AutoMinorLocator(5))

ax[0].set_ylabel(r'Flux, [W / m$^2$ / nm]', fontsize = 20)

ax[1].set_xlabel('Wavelength, [nm]', fontsize = 20)

#leg0 = ax[0].legend(framealpha = 1, loc = 3, bbox_to_anchor = (0.539, 0.03), handletextpad = 1, prop = {'size': 22.0})
#leg2 = ax[2].legend(framealpha = 1, loc = 1,                                 handletextpad = 1, prop = {'size': 22.0})

#for obj in leg0.legendHandles: obj.set_linewidth(5.0)
#for obj in leg2.legendHandles: obj.set_linewidth(5.0)

auxplt.savepdf('var/Q_kur_old')
