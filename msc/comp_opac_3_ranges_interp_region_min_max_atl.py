import numpy             as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths
import phys
import importlib
import glob
import auxplt
import auxfile
import spec
import os

from tqdm import tqdm

os.system('rm ' + paths.figdir + 'plt_opac/*.pdf')

odf_grid = np.loadtxt(paths.inp + 'odf.table.grid')

numt = int(odf_grid[0] )
nump = int(odf_grid[1])

tabt = odf_grid[2 : 2 + numt]
tabp = odf_grid[2 + numt : ]

nbins = 90
nsubbins = 10

#wvlgrid = np.zeros(nbins + 1)
#odf = np.zeros((nbins, nsubbins, nump, numt))

#f = open(paths.inp + 'odf.table.nopre_noH_nopif', 'r')

#for inu in range(nbins):

#    line = f.readline().rstrip("\n")

#    wvlgrid[inu] =     int(float(line.split()[0]))
#    wvlgrid[inu + 1] = int(float(line.split()[1]))

#    for ip in range(nump):

#        for it in range(numt):

#            line = f.readline().rstrip("\n")

#            odf[inu, 0 : nsubbins, ip, it] = np.array(line.split()).astype(int)

#f.close()

#np.savez(paths.npz + 'odf.table', odf = odf)

odf_table = np.load(paths.npz + 'odf.table.npz')

odf = odf_table['odf']

otw = np.arange(100, 1000, 1)

odf = 10.0**(odf.reshape((900, nump, numt)) / 1000.0)

dir_n = '/mnt/SSD/sim/runs/hminus/IT1/atlodf/new/fal/nopre_noH_nopif/'
dir_o = '/mnt/SSD/sim/runs/hminus/IT1/atlodf/old/fal/base/lbkg/'

file_n = dir_n + 'linop.out'

T =  np.loadtxt(dir_n + 'atm.inp', usecols = [1])
n =  np.loadtxt(dir_n + 'atm.inp', usecols = [3])

ne = np.loadtxt(dir_n + 'inine.out') * n

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

apm = 2.137995438028139e-024

#dp = np.array([1, 2, 3, 4, 25, 37, 38, 39, 40, 41, 42, 53, 54, 55, 56, 57, 58, 70, 71, 72, 73, 74])
dp = np.array([1, 2, 3, 4, 25, 37, 38, 39, 40, 41, 42, 53, 54, 55, 56, 57, 58, 70, 71, 72, 73])

wv1_v = np.arange(240, 310, 1)
wv2_v = np.arange(600, 660, 1)
wv3_v = np.arange(790, 850, 1)

op1_v = 10.0**(np.loadtxt(paths.inp + 'send/interpolated_new/interpolated_240_310.dat', usecols = [3]) / 1000.0)
op2_v = 10.0**(np.loadtxt(paths.inp + 'send/interpolated_new/interpolated_600_660.dat', usecols = [3]) / 1000.0)
op3_v = 10.0**(np.loadtxt(paths.inp + 'send/interpolated_new/interpolated_790_850.dat', usecols = [3]) / 1000.0)

op1_v = op1_v.reshape((7, 10, 22))
op1_v = np.flip(op1_v, axis = 1)

op2_v = op2_v.reshape((6, 10, 22))
op2_v = np.flip(op2_v, axis = 1)

op3_v = op3_v.reshape((6, 10, 22))
op3_v = np.flip(op3_v, axis = 1)

wv_u = np.arange(240, 860, 1)

op_u = np.zeros((620, 22))

for k in range(len(dp)):

    i = dp[k]

    if i < 10:  number = '0' + str(i)
    if i >= 10: number =       str(i)

    f = open(paths.inp + 'send/calculated' + number + '/p00big3.bdf', 'r')

    j = 0

    opac = np.array([])

    for line in f:

        j += 1

        if j % 2 != 0:

            continue

        if j % 2 == 0:

            opac = np.concatenate((opac, np.array(line.split()).astype(int)))

    op_u[:, k] = 10.0**(opac / 1000.0)

pdfs = ''

odf_min = np.zeros((ndp, len(otw)))
odf_max = np.zeros((ndp, len(otw)))

for i in range(ndp):

    p = (n[i] + ne[i]) * phys.boltz * T[i]

    it = np.searchsorted(tabt, T[i])
    ip = np.searchsorted(tabp, p)

    for j in range(len(otw)):

        odf_min[i, j] = apm * n[i] * min(odf[j, ip, it], odf[j, ip - 1, it - 1], odf[j, ip - 1, it], odf[j, ip, it - 1])
        odf_max[i, j] = apm * n[i] * max(odf[j, ip, it], odf[j, ip - 1, it - 1], odf[j, ip - 1, it], odf[j, ip, it - 1])

xmin = np.array([240, 600, 790])
xmax = np.array([300, 660, 850])

for k in tqdm(range(len(dp))):

    if k < 17:  m = k
    if k >= 17: m = k + 1

    i = dp[k] - 1

    p = (n[i] + ne[i]) * phys.boltz * T[i]

    it = np.searchsorted(tabt, T[i])
    ip = np.searchsorted(tabp, p)

    pdfs += str(i + 1) + '.pdf '

    plt.close('all')

    fig, ax = plt.subplots(nrows = 3, ncols = 1, figsize = (10, 10))

    for j in range(len(ax)):

        ax[j].set_xlim(xmin[j], xmax[j])

        if j == 0:

            wv = wv1_v
            op = op1_v[:, :, m].reshape(70)

        if j == 1:

            wv = wv2_v
            op = op2_v[:, :, m].reshape(60)

        if j == 2:

            wv = wv3_v
            op = op3_v[:, :, m].reshape(60)

#        ax[j].fill_between(otw, odf_min[i, :], odf_max[i, :], label = 'odf table', step = 'post', facecolor = 'grey', edgecolor = 'grey')

#        ax[j].step(wvl_o / 10.0, opa_o[:, i],                 label = 'NESSY (sorted and averaged)', where = 'post', color = 'k')
        ax[j].step(wv_u,         apm * n[i] * op_u[:, k],     label = 'ATLAS (sorted and averaged)', where = 'post', color = 'purple', linewidth = 1.5)
#        ax[j].step(wvl_n / 10.0, opa_n[:, i],                 label = 'NESSY (interpolation)',       where = 'post', color = 'r',      linewidth = 0.5)
        ax[j].step(wv,           apm * n[i] * op,             label = 'ATLAS (interpolation)',       where = 'post', color = 'g',      linewidth = 0.8)
#        ax[j].step(wv_u,         op_u[:, k],     label = 'ATLAS (sorted and averaged)', where = 'post', color = 'purple', linewidth = 1.5)
#        ax[j].step(wv,           op,             label = 'ATLAS (interpolation)',       where = 'post', color = 'g',      linewidth = 0.8)

#        idx = np.where((wv_u >= wv[0]) & (wv_u <= wv[len(wv) - 1]))[0]

#        wvl = wv_u[idx]

#        ax[j].step(wvl, ((op_u[idx, k] / op) - 1) * 100, where = 'post', color = 'k', linewidth = 1.5)

        ax[j].set_ylabel(r'Opacity, cm$^{-1}$')
#        ax[j].set_ylabel('Relative difference, %')

        ax[j].set_yscale('log')

        ax[j].xaxis.set_major_locator(MultipleLocator(10))
        ax[j].xaxis.set_minor_locator(AutoMinorLocator(10))

    ax[2].set_xlabel('Wavelength, nm')

    title = r'$L = $' + str(i + 1)

    if i + 1 == 1:    title +=  ', top'
    if i + 1 == tmdp: title += r', $T_\mathrm{min}$'
    if i + 1 == ndp:  title +=  ', bottom'

    leg0 = ax[0].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 7.5}, bbox_to_anchor=(0, 1.35))

    for obj in leg0.legendHandles: obj.set_linewidth(3.0)

    if tabp[ip - 1] > 0.00: presstr = r'$j = %.2i$' % (ip + 1, ) + r': $p_{j-1} = %.4f$' % (tabp[ip - 1], ) + r', $p_L = %.4f$' % (p, )    + r', $p_j = %.4f$' % (tabp[ip], )
    if tabp[ip - 1] > 1e+0: presstr = r'$j = %.2i$' % (ip + 1, ) + r': $p_{j-1} = %.4f$' % (tabp[ip - 1], ) + r', $p_L = %.4f$' % (p, )    + r', $p_j = %.4f$' % (tabp[ip], )
    if tabp[ip - 1] > 1e+1: presstr = r'$j = %.2i$' % (ip + 1, ) + r': $p_{j-1} = %.3f$' % (tabp[ip - 1], ) + r', $p_L = %.3f$' % (p, )    + r', $p_j = %.3f$' % (tabp[ip], )
    if tabp[ip - 1] > 1e+2: presstr = r'$j = %.2i$' % (ip + 1, ) + r': $p_{j-1} = %.2f$' % (tabp[ip - 1], ) + r', $p_L = %.2f$' % (p, )    + r', $p_j = %.2f$' % (tabp[ip], )
    if tabp[ip - 1] > 1e+3: presstr = r'$j = %.2i$' % (ip + 1, ) + r': $p_{j-1} = %.1f$' % (tabp[ip - 1], ) + r', $p_L = %.1f$' % (p, )    + r', $p_j = %.1f$' % (tabp[ip], )
    if tabp[ip - 1] > 1e+4: presstr = r'$j = %.2i$' % (ip + 1, ) + r': $p_{j-1} = %.0f$' % (tabp[ip - 1], ) + r', $p_L = %.0f$' % (p, )    + r', $p_j = %.0f$' % (tabp[ip], )
    if tabt[it - 1] > 0.00: tempstr = r'$i = %.2i$' % (it + 1, ) + r': $T_{i-1} = %.4f$' % (tabt[it - 1], ) + r', $T_L = %.4f$' % (T[i], ) + r', $T_i = %.4f$' % (tabt[it], )
    if tabt[it - 1] > 1e+0: tempstr = r'$i = %.2i$' % (it + 1, ) + r': $T_{i-1} = %.4f$' % (tabt[it - 1], ) + r', $T_L = %.4f$' % (T[i], ) + r', $T_i = %.4f$' % (tabt[it], )
    if tabt[it - 1] > 1e+1: tempstr = r'$i = %.2i$' % (it + 1, ) + r': $T_{i-1} = %.3f$' % (tabt[it - 1], ) + r', $T_L = %.3f$' % (T[i], ) + r', $T_i = %.3f$' % (tabt[it], )
    if tabt[it - 1] > 1e+2: tempstr = r'$i = %.2i$' % (it + 1, ) + r': $T_{i-1} = %.2f$' % (tabt[it - 1], ) + r', $T_L = %.2f$' % (T[i], ) + r', $T_i = %.2f$' % (tabt[it], )
    if tabt[it - 1] > 1e+3: tempstr = r'$i = %.2i$' % (it + 1, ) + r': $T_{i-1} = %.1f$' % (tabt[it - 1], ) + r', $T_L = %.1f$' % (T[i], ) + r', $T_i = %.1f$' % (tabt[it], )
    if tabt[it - 1] > 1e+4: tempstr = r'$i = %.2i$' % (it + 1, ) + r': $T_{i-1} = %.0f$' % (tabt[it - 1], ) + r', $T_L = %.0f$' % (T[i], ) + r', $T_i = %.0f$' % (tabt[it], )

    textstr = r'$\quad\quad\quad\quad\quad\quad\quad\quad\quad\quad\quad$' + title + '\n' + tempstr + '\n' + presstr

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    ax[0].text(0.55, 1.28, textstr, transform=ax[0].transAxes, fontsize=10,
               verticalalignment='top', bbox=props)

    auxplt.savepdf(str(i + 1), paths.figdir + 'plt_opac/')

os.chdir(paths.figdir + 'plt_opac/')

os.system('pdftk ' + pdfs + ' output overall.pdf')

os.chdir(paths.mscdir)
