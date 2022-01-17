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

#new_dir1 = '/mnt/SSD/sim/runs/hminus/IT1/atlodf/new/fal/nopre/'
#new_dir2 = '/mnt/SSD/sim/runs/hminus/IT1/atlodf/new/fal/nopre_noH_nopif/'
#new_dir1 = '/mnt/SSD/sim/runs/hminus/IT1/atlodf/new/kur/nopre/'
#new_dir2 = '/mnt/SSD/sim/runs/hminus/IT1/atlodf/new/kur/nopre_noH_nopif/'
new_dir1 = '/mnt/SSD/sim/runs/hminus/IT1/atlodf/new/fal/vdop1/'
new_dir2 = '/mnt/SSD/sim/runs/hminus/IT1/atlodf/new/fal/vdop2/'

#fo = '/mnt/SSD/sim/runs/hminus/IT1/atlodf/old/kur/base/lbkg/'
fo = '/mnt/SSD/sim/runs/hminus/IT1/atlodf/old/fal/base/lbkg/'

fn1 = new_dir1 + 'linop.out'
fn2 = new_dir2 + 'linop.out'

T = np.loadtxt(new_dir1 + 'atm.inp', usecols = [1])

tmdp = np.argmin(T) + 1

ndp = auxfile.num_lines(fo + '1005.lbkg') - 1

#wvl_o = np.arange(1005, 10000, 10)
wvl_o = np.arange(1000, 10000, 10)
vertl = np.arange(100, 1000, 10)

opa_o = np.zeros((len(wvl_o), ndp))

for i, wvl in enumerate(wvl_o):

#    opa_o[i, :] = np.loadtxt('/mnt/SSD/sim/runs/hminus/IT1/test_atlodf_old_totop/lbkg/' + str(wvl + 5) + '.lbkg', skiprows = 1)
    opa_o[i, :] = np.loadtxt(fo + str(wvl + 5) + '.lbkg', skiprows = 1)
#    opa_o[i, :] = np.loadtxt('./lbkg/' + str(wvl) + '.lbkg', skiprows = 1)

lam_n = np.loadtxt(fn1, usecols = [0])

kap_n1 = np.loadtxt(fn1, usecols = [4])
kap_n2 = np.loadtxt(fn2, usecols = [4])

#wvl_n = np.zeros(924)
wvl_n = np.zeros(941)

#opa_n1 = kap_n1.reshape((924, ndp))
#opa_n2 = kap_n2.reshape((924, ndp))
opa_n1 = kap_n1.reshape((941, ndp))
opa_n2 = kap_n2.reshape((941, ndp))

j = 0

for i, elem in enumerate(lam_n):

    if i % ndp == 0:

        wvl_n[j] = elem

        j += 1

pdfs = ''

os.system('rm ' + paths.figdir + 'plt_opac/*.pdf')

wvl = np.arange(105, 1000, 10)

wvl_ai = np.arange(590, 660, 1)

opac_ai = np.loadtxt('../../inp/interpolated_590_660_vturb_2.dat', usecols = [2])

opac_ai = opac_ai.reshape(7, 10)

opac_ai = np.fliplr(opac_ai)

opac_ai = opac_ai.reshape(70)

n = np.loadtxt(paths.it0f + '/atlodf/old/fal/base/atm.inp', usecols = [3])

apm = 2.137995438028139e-024

opac_ai *= n[54] * apm

wvl3 = np.arange(600, 621, 1)

opa3 = np.array([-4134, -3763, -3575, -3433, -3305, -3170, -3040, -2881, -2645, -956,
                 -4054, -3708, -3506, -3357, -3224, -3081, -2912, -2689, -2239, -179, 0.0])

opa3 = 10.0**(opa3 / 1000.0)

opa3 *= n[54] * apm

for i in tqdm(range(ndp)):
#for i in [54]:

#    opacn = np.interp(wvl_o, wvl_n, opa_n[:, i])

    pdfs += str(i + 1) + '.pdf '

    plt.close('all')

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (10, 10))

#    ax.step(wvl_o / 10.0, opa_o[:, i], where = 'post', label = 'FIOSS', color = 'k')

#    ax.plot(wvl_n / 10.0, opa_n1[:, i], label = 'ATLAS (nopre)',               color = 'r')
#    ax.plot(wvl_n / 10.0, opa_n2[:, i], label = 'ATLAS (no H, no PIF)', color = 'g', linewidth = 0.5)
#    ax.step(wvl_n / 10.0, opa_n1[:, i], where = 'post', label = 'ATLAS (nopre)',               color = 'r')
#    ax.step(wvl_n / 10.0, opa_n2[:, i], where = 'post', label = 'ATLAS (nopre, no H, no PIF)', color = 'g', linewidth = 0.5)

#    ax.plot(wvl_o / 10.0, opa_o[:, i], label = 'FIOSS', color = 'k', linewidth = 0.5)
#    ax.plot(wvl_ai, opac_ai, label = 'Interpolated', color = 'g', linewidth = 2)
#    ax.plot(wvl_n / 10.0, opa_n2[:, i], label = 'Interpolated by Rinat', color = 'r', linewidth = 0.5)
#    ax.plot(wvl3, opa3, label = 'Sorted and averaged', color = 'r', linewidth = 0.5)

#    ax.plot(wvl_n / 10.0, opa_n1[:, i], label = 'vdop = 1 km / s', color = 'r')
#    ax.plot(wvl_n / 10.0, opa_n2[:, i], label = 'vdop = 2 km / s', color = 'g', linewidth = 0.5)

#    ax.bar(wvl_o / 10.0, opa_o[:, i], width = 1, label = 'FIOSS', color = 'k')
#    ax.bar(wvl_o / 10.0, opacn,       width = 1, label = 'ATLAS', alpha = 0.5, color = 'r')

#    opao = spec.mean_within_delta_over_grid(opa_o[:, i], wvl_o / 10.0, 10.0, wvl)

    opan1 = spec.mean_within_delta_over_grid(opa_n1[:, i], wvl_n / 10.0, 10.0, wvl)
    opan2 = spec.mean_within_delta_over_grid(opa_n2[:, i], wvl_n / 10.0, 10.0, wvl)

#    ax.step(wvl, opao,  label = 'FIOSS', color = 'k')

#    ax.step(wvl, opan1, label = 'ATLAS (no preselection)',               alpha = 1.0, color = 'r')
#    ax.step(wvl, opan2, label = 'ATLAS (no preselection, no H, no PIF)', alpha = 1.0, color = 'g')
    ax.step(wvl * 10 - 50, opan1, label = 'vdop = 1 km / s', alpha = 1.0, color = 'r', where = 'post')
    ax.step(wvl * 10 - 50, opan2, label = 'vdop = 2 km / s', alpha = 1.0, color = 'g', where = 'post')

    ax.set_xlabel('Wavelength, nm')
    ax.set_ylabel(r'Opacity, cm$^{-1}$')

    ax.set_yscale('log')

#    ax.set_xlim(100, 1000)
    ax.set_xlim(1000, 4500)
#    ax.set_xlim(600, 620)
#    ax.set_xlim(590, 660)
#    ax.set_ylim(1e-35, 1e+0)
#    ax.set_ylim(1e-17, 1e-1)
    ax.set_ylim(1e-18, 1e-2)
#    ax.set_ylim(1e-14, 1e-8)

#    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_locator(MultipleLocator(1000))
    ax.xaxis.set_minor_locator(AutoMinorLocator(10))

    title = 'depth point ' + str(i + 1)

    if i + 1 == 1:    title +=  ', top'
    if i + 1 == tmdp: title += r', $T_\mathrm{min}$'
    if i + 1 == ndp:  title +=  ', bottom'

#    for elem in vertl:

#        ax.axvline(x = elem, linewidth = 0.3, color = 'g')

    ax.set_title(title)

    leg = ax.legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 7.5}, bbox_to_anchor=(0, 1.08))

    for obj in leg.legendHandles: obj.set_linewidth(3.0)

    auxplt.savepdf(str(i + 1), paths.figdir + 'plt_opac/')

#    sys.exit()

os.chdir(paths.figdir + 'plt_opac/')

os.system('pdftk ' + pdfs + ' output overall.pdf')

os.chdir(paths.mscdir)
