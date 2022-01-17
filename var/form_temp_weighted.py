import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import more_itertools as mit

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import importlib
import auxplt
import auxsys
import paths
import nessy

importlib.reload(auxplt)
importlib.reload(auxsys)
importlib.reload(paths)
importlib.reload(nessy)

def plot_gr(wvl, contr, indices, col, label = ''):

    groups = [list(gr) for gr in mit.consecutive_groups(indices[0].tolist())]
    
    for i, g in enumerate(groups):

        idx = (np.array(g), )

        if len(label) != 0:

            if i == 0: plt.plot(wvl[idx], contr[idx], color = col, label = label)
            if i != 0: plt.plot(wvl[idx], contr[idx], color = col)

        else:

            plt.plot(wvl[idx], -contr[idx], color = col, linestyle = '--')

wvl, fhq, fTq = nessy.weighted_form_temp(paths.it0f + 'var/Q/kur', 1005, 11000)
wvl, fhf, fTf = nessy.weighted_form_temp(paths.it0f + 'var/F/kur', 1005, 11000)

contr0 = (fTf - fTq)

wvl, fhq, fTq = nessy.weighted_form_temp(paths.it1f + 'var/Q/kur', 1005, 11000)
wvl, fhf, fTf = nessy.weighted_form_temp(paths.it1f + 'var/F/kur', 1005, 11000)

contr1 = (fTf - fTq)

wvl, fhq, fTq = nessy.weighted_form_temp(paths.it1f + 'var/Q/fal', 1005, 11000)
wvl, fhf, fTf = nessy.weighted_form_temp(paths.it1f + 'var/F/fal', 1005, 11000)

contr2 = (fTf - fTq)

wvl, fhqr, fTqr = nessy.weighted_form_temp(paths.it1f + 'var/Q/fal_ring', 1005, 11000)
wvl, fhfr, fTfr = nessy.weighted_form_temp(paths.it1f + 'var/F/fal_ring', 1005, 11000)

contr3 = (fTfr - fTqr)

np.savez(paths.npz + 'ftemp', wvl = wvl,
                              contr0 = contr0,
                              contr1 = contr1,
                              contr2 = contr2,
                              contr3 = contr3)

ftemp = np.load(paths.npz + 'ftemp.npz')

wvl = ftemp['wvl'] / 10.0

contr0 = ftemp['contr0']# * 100
contr1 = ftemp['contr1']# * 100
contr2 = ftemp['contr2']# * 100
contr3 = ftemp['contr3']# * 100

plt.close('all')

#fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (const.golden * 5, 5))
fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (18, 6))

auxplt.figpar(3, 3, 15)

fig.tight_layout()

idx0_pos = np.where(contr0 >= 0.0)
idx1_pos = np.where(contr1 >= 0.0)
idx2_pos = np.where(contr2 >= 0.0)
idx3_pos = np.where(contr3 >= 0.0)

idx0_neg = np.where(contr0 < 0.0)
idx1_neg = np.where(contr1 < 0.0)
idx2_neg = np.where(contr2 < 0.0)
idx3_neg = np.where(contr3 < 0.0)

plot_gr(wvl, contr0, idx0_pos, 'm', '(NESSY, LTE, U99)')
plot_gr(wvl, contr1, idx1_pos, 'g', '(NESSY, NLTE, U99)')
plot_gr(wvl, contr2, idx2_pos, 'r', '(NESSY, NLTE, FAL99)')
                             
plot_gr(wvl, contr0, idx0_neg, 'm')
plot_gr(wvl, contr1, idx1_neg, 'g')
plot_gr(wvl, contr2, idx2_neg, 'r')

plt.xlim(100, 1100)

plt.yscale('log')

plt.ylim(top = 1.1e+3)

ax.xaxis.set_major_locator(MultipleLocator(100))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))

plt.xlabel('Wavelength, [nm]')
#plt.ylabel(r'$(T^\mathrm{form}_f - T^\mathrm{form}_q) / T^\mathrm{form}_q, [\%]$')
plt.ylabel(r'$T^\mathrm{form}_f - T^\mathrm{form}_q, [K]$')

leg = plt.legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 17.0})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/ftemp_weighted')

sys.exit()

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (const.golden * 5, 5))

auxplt.figpar(3, 3, 15)

fig.tight_layout()

plt.plot(wvl / 10.0, fhq, label = 'Quiet Sun (NESSY-FAL99)')
plt.plot(wvl / 10.0, fhf, label = 'Facula (NESSY-FAL99)')

plt.xlim(100, 1000)

plt.xlabel('Wavelength, [nm]')
plt.ylabel('Formation height, [km]')

leg = plt.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 17.0})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/form_height')
