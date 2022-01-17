import numpy             as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths
import importlib
import glob
import auxplt
import auxfile
import spec
import os

from tqdm import tqdm

dir_n = '/mnt/SSD/sim/runs/hminus/IT1/atlodf/new/fal/nopre_noH_nopif/'
dir_o = '/mnt/SSD/sim/runs/hminus/IT1/atlodf/old/fal/base/lbkg/'

file_n = dir_n + 'linop.out'

T = np.loadtxt(dir_n + 'atm.inp', usecols = [1])
n = np.loadtxt(dir_n + 'atm.inp', usecols = [3])

tmdp = np.argmin(T) + 1

ndp = auxfile.num_lines(dir_n + 'atm.inp')

nwv = auxfile.num_lines(file_n)

wvl_o = np.arange(1000, 10000, 10)

opa_o = np.zeros((len(wvl_o), ndp))

for i, wvl in enumerate(wvl_o):

    opa_o[i, :] = np.loadtxt(dir_o + str(wvl + 5) + '.lbkg', skiprows = 1)

lam_n = np.loadtxt(file_n, usecols = [0])
kap_n = np.loadtxt(file_n, usecols = [4])

nwv = int(nwv / ndp)

wvl_n = np.zeros(nwv)
opa_n = kap_n.reshape((nwv, ndp))

j = 0

for i, elem in enumerate(lam_n):

    if i % ndp == 0:

        wvl_n[j] = elem

        j += 1

pdfs = ''

apm = 2.137995438028139e-024

wv1_v = np.arange(200, 260, 1)
wv2_v = np.arange(590, 670, 1)
wv2_u = np.arange(600, 620, 1)
wv3_v = np.arange(800, 860, 1)

op1_v = 10.0**(np.loadtxt(paths.inp + 'interpolated_200_250_vturb_2.dat', usecols = [3]) / 1000.0)
op2_v = 10.0**(np.loadtxt(paths.inp + 'interpolated_590_670_vturb_2.dat', usecols = [3]) / 1000.0)
op3_v = 10.0**(np.loadtxt(paths.inp + 'interpolated_800_850_vturb_2.dat', usecols = [3]) / 1000.0)

op1_v = op1_v.reshape((6, 10, 82))
op1_v = np.flip(op1_v, axis = 1)

op2_v = op2_v.reshape((8, 10, 82))
op2_v = np.flip(op2_v, axis = 1)

op3_v = op3_v.reshape((6, 10, 82))
op3_v = np.flip(op3_v, axis = 1)

#op2_v = np.loadtxt(paths.inp + 'interpolated_590_660_vturb_2.dat', usecols = [2])
#op2_v = op2_v.reshape(7, 10)
#op2_v = np.fliplr(op2_v)
#op2_v = op2_v.reshape(70)

op2_u = 10.0**(np.array([-4134, -3763, -3575, -3433, -3305, -3170, -3040, -2881, -2645, -956,
                         -4054, -3708, -3506, -3357, -3224, -3081, -2912, -2689, -2239, -179]) / 1000.0)

os.system('rm ' + paths.figdir + 'plt_opac/*.pdf')

#label1 = 'NESSY (sorting and averaging)'
#label2 = 'ATLAS (sorting and averaging)'
#label3 = 'NESSY (interpolation)'
#label4 = 'ATLAS (interpolation)'

for i in tqdm(range(ndp)):
#for i in [54]:

    pdfs += str(i + 1) + '.pdf '

    plt.close('all')

    fig, ax = plt.subplots(nrows = 3, ncols = 1, figsize = (10, 10))

    for j in range(len(ax)):

        if j == 0:

            wv = wv1_v
            op = op1_v[:, :, i].reshape(60)

        if j == 1:

            wv = wv2_v
            op = op2_v[:, :, i].reshape(80)

        if j == 2:

            wv = wv3_v
            op = op3_v[:, :, i].reshape(60)

        ax[j].step(wvl_o / 10.0, opa_o[:, i],     label = 'NESSY (sorting and averaging)', where = 'post', color = 'k')

        if i == 54 and j == 1:

#            ax[1].step(wv2_v, apm * n[i] * op2_v, label = 'ATLAS (interpolation)', where = 'post', color = 'g',      linewidth = 0.8)
#            ax[1].step(wv2_u, apm * n[i] * op2_u, label = 'ATLAS (sorting and averaging)', where = 'post', color = 'purple', linewidth = 0.8)
#            ax[1].step(wv2_v, apm * n[i] * op2_v, where = 'post', color = 'g',      linewidth = 0.8)
            ax[1].step(wv2_u, apm * n[i] * op2_u, label = 'ATLAS (sorting and averaging)', where = 'post', color = 'purple', linewidth = 0.8)

        ax[j].step(wv,           apm * n[i] * op, label = 'ATLAS (interpolation)', where = 'post', color = 'g', linewidth = 0.8)
        ax[j].step(wvl_n / 10.0, opa_n[:, i],     label = 'NESSY (interpolation)', where = 'post', color = 'r', linewidth = 0.5)

#        ax[j].plot(wvl_o / 10.0, opa_o[:, i], label = 'NESSY (sorting and averaging)', color = 'k')
#        ax[j].plot(wvl_n / 10.0, opa_n[:, i], label = 'NESSY (interpolation)',         color = 'r', linewidth = 0.5)

        ax[j].set_ylabel(r'Opacity, cm$^{-1}$')

        ax[j].set_yscale('log')

        ax[j].xaxis.set_major_locator(MultipleLocator(10))
        ax[j].xaxis.set_minor_locator(AutoMinorLocator(10))

    ax[0].set_xlim(200, 260)
    ax[1].set_xlim(590, 670)
    ax[2].set_xlim(800, 860)

    ax[2].set_xlabel('Wavelength, nm')

    title = 'depth point ' + str(i + 1)

    if i + 1 == 1:    title +=  ', top'
    if i + 1 == tmdp: title += r', $T_\mathrm{min}$'
    if i + 1 == ndp:  title +=  ', bottom'

    ax[0].set_title(title)

#    if i != 54: leg = ax[0].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 7.5}, bbox_to_anchor=(0, 1.20))
    leg = ax[0].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 7.5}, bbox_to_anchor=(0, 1.20))
    if i == 54: leg = ax[1].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 7.5})

    for obj in leg.legendHandles: obj.set_linewidth(3.0)

    auxplt.savepdf(str(i + 1), paths.figdir + 'plt_opac/')

os.chdir(paths.figdir + 'plt_opac/')

os.system('pdftk ' + pdfs + ' output overall.pdf')

os.chdir(paths.mscdir)
