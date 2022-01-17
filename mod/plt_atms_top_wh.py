import numpy as np

import random

import sys
import os

if not '../aux/' in sys.path: sys.path.append('../aux/')

from tqdm import tqdm

import matplotlib.pyplot as plt

import auxplt
import paths

f = open('inputfile', 'r')

nr = int(f.readline())
nd = int(f.readline())

h, T, p, d = np.loadtxt('inputfile', skiprows = 2, unpack = True)

ho = h.reshape(nr, nd)
To = T.reshape(nr, nd)
po = p.reshape(nr, nd)
do = d.reshape(nr, nd)

f.close()

for i in range(nr):

    ho[i, :] += abs(min(ho[i, :]))

    ho[i, :] = np.flip(ho[i, :])

tau = np.arange(5, 10)
#tau = np.array([9])

nf = len(tau)

hc = np.zeros((nf, nr, nd))
Tc = np.zeros((nf, nr, nd))
pc = np.zeros((nf, nr, nd))
dc = np.zeros((nf, nr, nd))

for k in range(len(tau)):

    f = open('outputfile_tau_top_' + str(tau[k]), 'r')

    f.readline()
    f.readline()

    for i in tqdm(range(nr)):

        for j in range(nd):

            values = f.readline().split()

            hc[k, i, j] = float(values[1]) / 1.0e+5
            Tc[k, i, j] = float(values[2])
            pc[k, i, j] = float(values[3])
            dc[k, i, j] = float(values[4])

        for j in range(3): f.readline()

        hc[k, i, :] = -hc[k, i, :]

        hc[k, i, :] = np.flip(hc[k, i, :])

        hc[k, i, :] -= min(hc[k, i, :])

        idxmin = abs(To[i, :] - Tc[k, i, nd - 1]).argmin()

        hc[k, i, :] += ho[i, idxmin]

    f.close()

plt.close('all')

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

pdfs = ''

#for i in tqdm(range(0, nr, 256)):
#for i in tqdm(range(0, nr, 1)):
for i in tqdm(range(0, nr, 4)):
#for i in tqdm(range(0, nr, 32)):

    fig, ax = plt.subplots(nrows = 3, ncols = 1, figsize = (7.5, 10))

    pdfs += str(i + 1) + '.pdf '
#    j = random.choice(np.arange(0, 512))

    ax[0].plot(po[i, :], To[i, :], color = 'k', label = 'extracted', linewidth = 1.0)
    ax[1].plot(do[i, :], To[i, :], color = 'k', linewidth = 1.0)
    ax[2].plot(po[i, :], do[i, :], color = 'k', linewidth = 1.0)

    for k in range(len(tau)):

        ax[0].plot(pc[k, i, :], Tc[k, i, :], label = 'interpolated (' + r'$\tau_{500}^\mathrm{top}$ = 1.0e-' + str(tau[k]) + ')', linewidth = 1.1 + (len(tau) - 2 * k) * 0.5)
        ax[1].plot(dc[k, i, :], Tc[k, i, :], linewidth = 1.1 + (len(tau) - 2 * k) * 0.5)
        ax[2].plot(pc[k, i, :], dc[k, i, :], linewidth = 1.1 + (len(tau) - 2 * k) * 0.5)

#    ax[1].set_yscale('log')
#    ax[2].set_yscale('log')
    ax[2].set_yscale('log')

#    for j in [0, 1, 2]: ax[j].set_xlim(1600, 3600)

    ax[0].set_xlim(0.0, 300000)
    ax[1].set_xlim(0.0, 0.0000004)
    ax[2].set_xlim(0.0, 300000)

#    ax[0].set_ylim(top = 14000)

    ax[0].set_xlabel('Pressure')
    ax[0].set_ylabel('Temperature')
    ax[1].set_xlabel('Density')
    ax[1].set_ylabel('Temperature')
    ax[2].set_xlabel('Pressure')
    ax[2].set_ylabel('Density')
#    ax[2].set_xlabel('Height')

    leg = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 4.5})

    for obj in leg.legendHandles: obj.set_linewidth(3.0)

    ax[0].text(0.12, 0.12, str(i + 1), transform=ax[0].transAxes, fontsize=10,
               verticalalignment='top', bbox=props)

    ax[0].text(0.12, 0.95, r'$\tau_{500}^\mathrm{bot}$ = 2.5', transform=ax[0].transAxes, fontsize=10,
               verticalalignment='top', bbox=props)

    auxplt.savepdf(str(i + 1), './plt_atms/')

    plt.close('all')

os.chdir('./plt_atms/')

os.system('pdftk ' + pdfs + ' output atms_top_wh.pdf')

os.chdir('../')
