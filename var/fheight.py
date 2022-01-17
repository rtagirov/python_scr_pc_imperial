import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import more_itertools as mit
import math as m

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

from tqdm import tqdm

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import importlib
import auxplt
import auxsys
import auxfunc
import paths
import nessy
import spec

importlib.reload(auxplt)
importlib.reload(auxsys)
importlib.reload(auxfunc)
importlib.reload(paths)
importlib.reload(nessy)
importlib.reload(spec)

cases = ['Q kur 0', 'F kur 0', 'Q kur 1', 'F kur 1', 'Q fal 1', 'F fal 1']

wvl = np.linspace(1005, 11005, 1001)

delta_h = 5

nws_h = int(m.ceil((wvl[len(wvl) - 1] - wvl[0]) / delta_h / 10))

fhs = np.zeros((len(cases), nws_h))

for i, case in enumerate(cases):

    ar = case.split()[0]
    mo = case.split()[1]
    it = case.split()[2]

    path = ''

    if it == '0': path = paths.it0f
    if it == '1': path = paths.it1f

    path += 'var/' + ar + '/rings/' + mo

    wvl, fh = nessy.weighted_form_height(path + '/1.0', wvl[0], wvl[len(wvl) - 1] - 5)

    np.savez(paths.npz + ar + '_' + mo + '_' + it + '.fhei_center.npz', fh = fh)

    fh = np.load(paths.npz + ar + '_' + mo + '_' + it + '.fhei_center.npz')['fh']

    wvls_h, heis = spec.mean_within_delta(wvl / 10, fh, delta_h)

    fhs[i, :] = heis

plt.close('all')

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (6.0, 3.375))

auxplt.figpar(3, 3, 15)

ax.set_yscale('log')
ax.set_xlim(100, 1100)
ax.set_ylim(1e+0, 2e+3)

ax.xaxis.set_major_locator(MultipleLocator(200))
ax.xaxis.set_minor_locator(AutoMinorLocator(4))

ax.fill_between(wvls_h, 140 * np.ones(len(wvls_h)), 400 * np.ones(len(wvls_h)), color = 'gray', alpha = 0.5)

ax.plot(wvls_h, fhs[0, :] - 134.237, color = 'm')
ax.plot(wvls_h, fhs[1, :] - 104.736, color = 'm', linestyle = '--')
ax.plot(wvls_h, 1e+6 * fhs[2, :] - 134.237, color = 'k', label = 'Quiet Sun')
ax.plot(wvls_h, 1e+6 * fhs[3, :] - 104.736, color = 'k', linestyle = '--', label = 'Facula')
#ax.plot(wvls_h, fhs[2, :] - 134.237, color = 'k', label = 'Quiet Sun')
#ax.plot(wvls_h, fhs[3, :] - 104.736, color = 'k', linestyle = '--', label = 'Facula')
ax.plot(wvls_h, fhs[4, :] - 114.227, 'r-', alpha = 0.5)
ax.plot(wvls_h, fhs[5, :] - 106.725, 'r--', alpha = 0.5)

ax.set_xlabel('Wavelength, [nm]')
ax.set_ylabel(r'Formation height ($\mu = 1$), [km]')

ax.axvline(x = 210, linestyle = ':', color = 'k')
ax.axvline(x = 400, linestyle = ':', color = 'k')

leg1 = ax.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 12.0})

for obj in leg1.legendHandles: obj.set_linewidth(3.0)

ax.text(450, 9.0e+2, 'NESSY, LTE, U99',    color = 'm', fontsize = 11)
ax.text(450, 4.0e+2, 'NESSY, NLTE, FAL99', color = 'r', fontsize = 11)

auxplt.savepdf('var/fheight')
