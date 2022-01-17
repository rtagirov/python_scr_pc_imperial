import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import os

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)
import phys;   importlib.reload(phys)

from tqdm import tqdm

#wvln, opan = nessy.read_lopa(paths.it0f + 'runtime/def', wvl1 = 1005., wvl2 = 2100)

#wvln /= 10.0

#np.savez(paths.npz + 'unsorted_opac_100_200', w = wvln, o = opan)

opac = np.load(paths.npz + 'unsorted_opac_100_200.npz')

wvln = opac['w']
opan = opac['o']

wvla, opaa = np.loadtxt('/mnt/SSD/sim/nessy/inp/odf/high_res_100_210/fort.92', unpack = True)

n = np.loadtxt(paths.it0f + '/runtime/def/atm.inp', usecols = [3])

apm = 2.137995438028139e-024

opaa *= n[54] * apm

wvln_s, opan_s = spec.mean_within_delta(wvln, opan[54, :], 0.05)
wvla_s, opaa_s = spec.mean_within_delta(wvla, opaa,        0.05)

#wvl_n = wvln
#wvl_a = wvla
#opa_n = opan[54, :]
#opa_a = opaa
wvl_n = wvln_s
wvl_a = wvla_s
opa_n = opan_s
opa_a = opaa_s

plt.close('all')

xlim = [[100, 110],
        [110, 120],
        [120, 130],
        [130, 140],
        [140, 150],
        [150, 160],
        [160, 170],
        [170, 180],
        [180, 190],
        [190, 200],
        [200, 210]]

fig, ax = plt.subplots(nrows = len(xlim), ncols = 2, figsize = (12.0, 22.0))

for i in range(len(ax)):

    ax[i, 0].plot(wvl_n, opa_n, color = 'k')
    ax[i, 0].plot(wvl_a, opa_a, color = 'r')

    wvln_s = np.sort(wvl_n[np.where((wvl_n >= xlim[i][0]) & (wvl_n <= xlim[i][1]))])
    wvla_s = np.sort(wvl_a[np.where((wvl_a >= xlim[i][0]) & (wvl_a <= xlim[i][1]))])

    opan_s = np.sort(opa_n[np.where((wvl_n >= xlim[i][0]) & (wvl_n <= xlim[i][1]))])
    opaa_s = np.sort(opa_a[np.where((wvl_a >= xlim[i][0]) & (wvl_a <= xlim[i][1]))])

    ax[i, 1].plot(wvln_s, opan_s, color = 'k', label = 'NESSY')
    ax[i, 1].plot(wvla_s, opaa_s, color = 'r', label = 'ATLAS')

    ax[i, 0].set_yscale('log')
    ax[i, 1].set_yscale('log')

    ax[i, 0].set_xlim(xlim[i][0], xlim[i][1])
    ax[i, 1].set_xlim(xlim[i][0], xlim[i][1])

    ax[i, 0].set_ylabel('Opacity')

ax[len(ax) - 1, 0].set_xlabel('Wavelength, nm')
ax[len(ax) - 1, 1].set_xlabel('Wavelength, nm')

#leg = ax[0, 1].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 7.5}, bbox_to_anchor=(0, 1.08))
leg = ax[0, 1].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 7.5})

for obj in leg.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('unsorted_100_200', paths.figdir + 'plt_opac/')
