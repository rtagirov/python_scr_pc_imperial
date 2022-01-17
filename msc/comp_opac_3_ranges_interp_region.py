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

l, it, ip, Tl, T, Tu, pl, p, pu = np.loadtxt(paths.inp + 'interpolation.txt', unpack = True)

os.system('rm ' + paths.figdir + 'plt_opac/*.pdf')

plt.close('all')

fig, ax = plt.subplots(nrows = 3, ncols = 1, figsize = (10, 10))

plt.subplots_adjust(hspace = 0.1)

for i in range(len(ax)):

    ax[i].set_xlim(1, 82)

    ax[i].xaxis.set_minor_locator(AutoMinorLocator(10))

ax[0].plot(l, it, color = 'blue', label = 'temperature')
ax[0].plot(l, ip, color = 'brown', label = 'pressure')

ax[1].plot(l, T, color = 'blue', label = 'FAL-C', linewidth = 0.5)

ax[1].fill_between(l, Tl, Tu, color = 'grey', label = 'ODF grid')

ax[1].set_yscale('log')
ax[2].set_yscale('log')

ax[2].plot(l, p, color = 'brown', label = 'FAL-C', linewidth = 0.5)
ax[2].fill_between(l, pl, pu, color = 'grey', label = 'ODF grid')

ax[0].yaxis.set_major_locator(MultipleLocator(10))
ax[0].yaxis.set_minor_locator(AutoMinorLocator(10))

ax[0].set_ylim(20, 100)

ax[0].set_ylabel('ODF T- and p-grid indices')
ax[1].set_ylabel('Temperature, K')
ax[2].set_ylabel('Pressure, cgs')
ax[2].set_xlabel('FAL-C depth index')

leg = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 10.5})
leg = ax[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 10.5})
leg = ax[2].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 10.5})

auxplt.savepdf('interp_info', paths.figdir + 'plt_opac/')

#sys.exit()

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

op2_u = 10.0**(np.array([-4134, -3763, -3575, -3433, -3305, -3170, -3040, -2881, -2645, -956,
                         -4054, -3708, -3506, -3357, -3224, -3081, -2912, -2689, -2239, -179]) / 1000.0)

pdfs = ''

#for i in tqdm(range(ndp)):
for i in range(ndp):
#for i in [54]:

    p = (n[i] + ne[i]) * phys.boltz * T[i]

    it = np.searchsorted(tabt, T[i])
    ip = np.searchsorted(tabp, p)

    print(i + 1, ' ', it + 1, ' ', ip + 1, ' ', tabt[it - 1], ' ', T[i], ' ', tabt[it],
                                           ' ', tabp[ip - 1], ' ', p,    ' ', tabp[ip])

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

        ax[j].fill_between(otw, apm * n[i] * odf[:, ip - 1, it - 1], apm * n[i] * odf[:, ip, it], step = 'post', facecolor = 'grey', label = 'odf table', edgecolor = 'grey')

        ax[j].step(wvl_o / 10.0, opa_o[:, i],     label = 'NESSY (sorting and averaging)', where = 'post', color = 'k')

        if i == 54 and j == 1:

            ax[1].step(wv2_u, apm * n[i] * op2_u, label = 'ATLAS (sorting and averaging)', where = 'post', color = 'purple', linewidth = 0.8)

        ax[j].step(wv,           apm * n[i] * op, label = 'ATLAS (interpolation)', where = 'post', color = 'g', linewidth = 0.8)
        ax[j].step(wvl_n / 10.0, opa_n[:, i],     label = 'NESSY (interpolation)', where = 'post', color = 'r', linewidth = 0.5)

        ax[j].set_ylabel(r'Opacity, cm$^{-1}$')

        ax[j].set_yscale('log')

        ax[j].xaxis.set_major_locator(MultipleLocator(10))
        ax[j].xaxis.set_minor_locator(AutoMinorLocator(10))

    ax[0].set_xlim(200, 260)
    ax[1].set_xlim(590, 670)
    ax[2].set_xlim(800, 860)

    ax[2].set_xlabel('Wavelength, nm')

    title = r'$L = $' + str(i + 1)

    if i + 1 == 1:    title +=  ', top'
    if i + 1 == tmdp: title += r', $T_\mathrm{min}$'
    if i + 1 == ndp:  title +=  ', bottom'

#    ax[0].set_title(title)

    leg = ax[0].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 7.5}, bbox_to_anchor=(0, 1.35))
    if i == 54: leg = ax[1].legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 7.5})

    for obj in leg.legendHandles: obj.set_linewidth(3.0)

    if tabp[ip - 1] > 0.00: presstr = r'$j = %.2i$' % (ip, ) + r': $p_{j-1} = %.4f$' % (tabp[ip - 1], ) + r', $p_L = %.4f$' % (p, )    + r', $p_j = %.4f$' % (tabp[ip], )
    if tabp[ip - 1] > 1e+0: presstr = r'$j = %.2i$' % (ip, ) + r': $p_{j-1} = %.4f$' % (tabp[ip - 1], ) + r', $p_L = %.4f$' % (p, )    + r', $p_j = %.4f$' % (tabp[ip], )
    if tabp[ip - 1] > 1e+1: presstr = r'$j = %.2i$' % (ip, ) + r': $p_{j-1} = %.3f$' % (tabp[ip - 1], ) + r', $p_L = %.3f$' % (p, )    + r', $p_j = %.3f$' % (tabp[ip], )
    if tabp[ip - 1] > 1e+2: presstr = r'$j = %.2i$' % (ip, ) + r': $p_{j-1} = %.2f$' % (tabp[ip - 1], ) + r', $p_L = %.2f$' % (p, )    + r', $p_j = %.2f$' % (tabp[ip], )
    if tabp[ip - 1] > 1e+3: presstr = r'$j = %.2i$' % (ip, ) + r': $p_{j-1} = %.1f$' % (tabp[ip - 1], ) + r', $p_L = %.1f$' % (p, )    + r', $p_j = %.1f$' % (tabp[ip], )
    if tabp[ip - 1] > 1e+4: presstr = r'$j = %.2i$' % (ip, ) + r': $p_{j-1} = %.0f$' % (tabp[ip - 1], ) + r', $p_L = %.0f$' % (p, )    + r', $p_j = %.0f$' % (tabp[ip], )
    if tabt[it - 1] > 0.00: tempstr = r'$i = %.2i$' % (it, ) + r': $T_{i-1} = %.4f$' % (tabt[it - 1], ) + r', $T_L = %.4f$' % (T[i], ) + r', $T_i = %.4f$' % (tabt[it], )
    if tabt[it - 1] > 1e+0: tempstr = r'$i = %.2i$' % (it, ) + r': $T_{i-1} = %.4f$' % (tabt[it - 1], ) + r', $T_L = %.4f$' % (T[i], ) + r', $T_i = %.4f$' % (tabt[it], )
    if tabt[it - 1] > 1e+1: tempstr = r'$i = %.2i$' % (it, ) + r': $T_{i-1} = %.3f$' % (tabt[it - 1], ) + r', $T_L = %.3f$' % (T[i], ) + r', $T_i = %.3f$' % (tabt[it], )
    if tabt[it - 1] > 1e+2: tempstr = r'$i = %.2i$' % (it, ) + r': $T_{i-1} = %.2f$' % (tabt[it - 1], ) + r', $T_L = %.2f$' % (T[i], ) + r', $T_i = %.2f$' % (tabt[it], )
    if tabt[it - 1] > 1e+3: tempstr = r'$i = %.2i$' % (it, ) + r': $T_{i-1} = %.1f$' % (tabt[it - 1], ) + r', $T_L = %.1f$' % (T[i], ) + r', $T_i = %.1f$' % (tabt[it], )
    if tabt[it - 1] > 1e+4: tempstr = r'$i = %.2i$' % (it, ) + r': $T_{i-1} = %.0f$' % (tabt[it - 1], ) + r', $T_L = %.0f$' % (T[i], ) + r', $T_i = %.0f$' % (tabt[it], )

    textstr = r'$\quad\quad\quad\quad\quad\quad\quad\quad\quad\quad\quad$' + title + '\n' + tempstr + '\n' + presstr

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    ax[0].text(0.55, 1.28, textstr, transform=ax[0].transAxes, fontsize=10,
               verticalalignment='top', bbox=props)

    auxplt.savepdf(str(i + 1), paths.figdir + 'plt_opac/')

os.chdir(paths.figdir + 'plt_opac/')

os.system('pdftk ' + pdfs + ' output overall.pdf')

os.chdir(paths.mscdir)
