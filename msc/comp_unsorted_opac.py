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

wvl3 = np.arange(600, 621, 1) * 10.0

opa3 = np.array([-4134, -3763, -3575, -3433, -3305, -3170, -3040, -2881, -2645, -956,
                 -4054, -3708, -3506, -3357, -3224, -3081, -2912, -2689, -2239, -179, 0.0])

opa3 = 10.0**(opa3 / 1000.0)

#wvl1, opa1 = nessy.read_lopa(paths.it0f + '/atlodf/old/fal/base', wvl1 = 6000., wvl2 = 6200)
wvl2, opa2 = np.loadtxt('600_620v2.dat', unpack = True)

#np.savez(paths.npz + 'unsorted_opac', w = wvl1, o = opa1)

opac = np.load(paths.npz + 'unsorted_opac.npz')

wvl1 = opac['w']
opa1 = opac['o']

n = np.loadtxt(paths.it0f + '/atlodf/old/fal/base/atm.inp', usecols = [3])

apm = 2.137995438028139e-024

opa2 *= n[54] * apm
opa3 *= n[54] * apm

wvl2 *= 10.0

#wvl2 = phys.vac_to_air(wvl2)

step = 100

xbin =  [[i, i + step]  for i in np.arange(6000, 6200, step)]

pdfs = ''

os.system('rm ' + paths.figdir + 'plt_opac/*.pdf')

for i in tqdm(range(len(xbin))):

    name = str(xbin[i][0]) + '_' + str(xbin[i][1])

    pdfs += name + '.pdf '

    plt.close('all')

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 10))

#    ax.plot(wvl1, opa1[54, :], label = 'FIOSS', color = 'k',      alpha = 0.5)
#    ax.plot(wvl2, opa2,        label = 'ATLAS', color = 'orange', alpha = 0.5)

    idx1 = np.where((wvl1 >= xbin[i][0]) & (wvl1 <= xbin[i][1]))
    idx2 = np.where((wvl2 >= xbin[i][0]) & (wvl2 <= xbin[i][1]))

    wvl1s = wvl1[idx1]
    wvl2s = wvl2[idx2]

    opa1s = np.sort(opa1[54, idx1])
    opa2s = np.sort(opa2[idx2])

    xbina = [[i, i + 10] for i in np.arange(xbin[i][0], xbin[i][1], 10)]

    wvlm = np.zeros(len(xbina) + 1)

    opa1m = np.zeros(len(xbina) + 1)
    opa2m = np.zeros(len(xbina) + 1)

    for j in range(len(xbina)):

        idxa1 = np.where((wvl1s >= xbina[j][0]) & (wvl1s <= xbina[j][1]))
        idxa2 = np.where((wvl2s >= xbina[j][0]) & (wvl2s <= xbina[j][1]))

#        wvlm[j] = (xbina[j][0] + xbina[j][1]) / 2.0
        wvlm[j] = xbina[j][0]

        opa1m[j] = np.mean(opa1s[0, idxa1])
        opa2m[j] = np.mean(opa2s[idxa2])

        if j == len(xbina) - 1:

            wvlm[len(xbina)] = xbina[j][1]

#    plt.step(wvlm, opa1m, where = 'post', color = 'k')
    plt.step(wvlm, opa2m, where = 'post', color = 'r')

    plt.step(wvl3, opa3, where = 'post', color = 'purple')

#    plt.plot(wvlm, opa1m, color = 'k')
#    plt.plot(wvlm, opa2m, color = 'r')

#    ax.plot(wvl1s, opa1s[0, :], color = 'k')
#    ax.plot(wvl2s, opa2s,       color = 'r')

#    opa1m = np.mean(opa1[54, idx1])
#    opa2m = np.mean(opa2[idx2])

#    ax.axhline(y = opa1m, color = 'k', linestyle = '--')
#    ax.axhline(y = opa2m, color = 'r', linestyle = '--')

    ax.set_xlim(xbin[i][0], xbin[i][1])

    ax.set_xlabel('Wavelength, A')
    ax.set_ylabel(r'Opacity, cm$^{-1}$')

    ax.set_yscale('log')

    ax.xaxis.set_major_locator(MultipleLocator(10))

#    leg = ax.legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 7.5}, bbox_to_anchor=(0, 1.08))

#    for obj in leg.legendHandles: obj.set_linewidth(3.0)

    auxplt.savepdf(name, paths.figdir + 'plt_opac/')

os.chdir(paths.figdir + 'plt_opac/')

os.system('pdftk ' + pdfs + ' output overall.pdf')

os.chdir(paths.mscdir)
