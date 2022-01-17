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

def plot_gr(wvl, contr, indices, col, label = ''):

    groups = [list(gr) for gr in mit.consecutive_groups(indices[0].tolist())]
    
    for i, g in enumerate(groups):

        idx = (np.array(g), )

        if len(label) != 0:

            if i == 0: ax[0].plot(wvl[idx], contr[idx], color = col, label = label)
            if i != 0: ax[0].plot(wvl[idx], contr[idx], color = col)

        else:

            ax[0].plot(wvl[idx], -contr[idx], color = col, linestyle = '--')

def formation_height(path, wvl1, wvl2, mode, num):

    nf =  2000 # number of lines in each .tau and .mdisp file

    ivl = 10   # length of the spectral interval of each .tau and .mdisp file (angstroems)
    mid = 5    # distance to the middle of each .tau and .mdisp file (angstroems)

    nint = int(ivl)

    wmin = int(wvl1)
    wmax = int(wvl2)

    imin = wmin - ((wmin - mid) % nint)
    imax = wmax - ((wmax - mid) % nint) + nint

    na = (imax - imin) / nint + 1 # number of arrays

    nw = na * nf

    wvl = np.array([])
    frm = np.array([])

    for i in tqdm(range(int(na)), ncols = auxfunc.term_width(), desc = 'mode ' + str(mode) + '; run ' + str(num) + '; ' + path):

        idx = str(imin + i * nint)

        f1 = path + '/tau/' + idx + '.tau'

        w, h, x, x, x = np.loadtxt(f1, unpack = True)

        wvl = np.concatenate((wvl, w))
        frm = np.concatenate((frm, h))

    frm[np.where(np.isnan(frm))] = 0.0

    return wvl, frm

def ring_weights(mu, mode):

    p = np.sqrt(1 - mu**2)

    w = np.zeros(len(p))

    if mode >= 0 and mode <= 10:

        w[mode] = 1.0

    if mode == 11:

        for i in range(1, len(p) - 1):

            w[i] = p[i + 1]**2 - p[i - 1]**2

        w[0] = p[1]**2

        w[len(p) - 1] = 1 - p[len(p) - 2]**2

    if mode == 12:

        for i in range(len(p) - 1):

            w[i] = p[i + 1]**2 - p[i]**2

        w[len(p) - 1] = 1 - p[len(p) - 1]**2

    if mode == 13:

        w = np.ones(len(p)) * 0.1

    return w / sum(w)

cases = ['Q fal 1', 'F fal 1']

p = np.sqrt(np.arange(0, 11)) / np.sqrt(10)

p_mid = np.zeros(len(p) - 1)

for i in range(len(p_mid)):

    p_mid[i] = (p[i + 1] + p[i]) / 2

mu = np.sqrt(1 - p_mid**2)

mode = int(sys.argv[1])

w = ring_weights(mu, mode)

wvl = np.linspace(1005, 11005, 1001)

fh = np.zeros((len(cases),  2000 * len(wvl)))

delta = 5

num = 1

for i, case in enumerate(cases):

    ar = case.split()[0]
    mo = case.split()[1]
    it = case.split()[2]

    path = ''

    if it == '0': path = paths.it0f
    if it == '1': path = paths.it1f

    path += 'var_od/' + ar + '/rings/' + mo

#    for j in range(len(mu)):

##        if mode >= 0 and mode <= 10 and j != mode: continue

#        wvl_h, fhr = formation_height(path + '/' + str(mu[j])[:4], wvl[0], wvl[len(wvl) - 1] - 5, mode, num)

#        fh[i, :] += fhr * mu[j] * w[j]

#        num += 1
    
#    np.savez(paths.npz + ar + '_' + mo + '_' + it + '.fhei_od_high_res.npz', fh = fh[i, :])

    fh[i, :] = np.load(paths.npz + ar + '_' + mo + '_' + it + '.fhei_od_high_res.npz')['fh']

#np.savez(paths.npz + ar + '_' + mo + '_' + it + '.fhei_od_high_res_wvl.npz', wh = wvl_h)

wh = np.load(paths.npz + ar + '_' + mo + '_' + it + '.fhei_od_high_res_wvl.npz')['wh']

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6, 7))

auxplt.figpar(3, 3, 15)

ax[0].set_xlim(5248, 5253)
ax[0].set_ylim(20, 2e+3)
ax[1].set_xlim(7697, 7703)
ax[1].set_ylim(45, 60)

ax[0].set_yscale('log')

ax[0].xaxis.set_major_locator(MultipleLocator(1))
ax[0].xaxis.set_minor_locator(AutoMinorLocator(4))
ax[1].xaxis.set_major_locator(MultipleLocator(1))
ax[1].xaxis.set_minor_locator(AutoMinorLocator(4))

ax[0].axvline(x = 5251.461390185789, color = 'g', linestyle = '-')
ax[1].axvline(x = 7701.099088671894, color = 'g', linestyle = '-')

ax[0].axvline(x = 5250.00, color = 'g', linestyle = '--')
ax[1].axvline(x = 7698.98, color = 'g', linestyle = '--')

ax[0].plot(wh, fh[0, :] - 114.388, 'k',                  label = 'Quiet Sun (FAL99-C)')
ax[0].plot(wh, fh[1, :] - 106.802, 'r', linewidth = 0.5, label = 'Facula (FAL99-P)')
ax[1].plot(wh, fh[0, :] - 114.388, 'k')
ax[1].plot(wh, fh[1, :] - 106.802, 'r', linewidth = 0.5)

ax[0].set_ylabel('Formation height, [km]')
ax[1].set_ylabel('Formation height, [km]')
ax[1].set_xlabel('Wavelength, [nm]')

leg1 = ax[0].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 12.0})

for obj in leg1.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/fheight_weighted_rings_2rows_od_high_res')
