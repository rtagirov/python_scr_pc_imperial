import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import phys
import auxplt

import itertools

#wk, nk, gfk, elk, x, euk, x, grk, gsk, gvk, x = np.loadtxt('/mnt/SSD/sim/nessy/inp/lin/kurucz_comb_air_parts/Kurucz_comb_air_6900to11000.dat', unpack = True)
#wv, nv, gfv, elv, x, euv, x, grv, gsv, gvv, x = np.loadtxt('/home/rtagirov/VALD_Alexander/vald_in_kurucz_format.dat', unpack = True)

elem = 'Cr'; istage = 'II'; elemistage = elem + istage

#if istage =='I':   istagedec = 0.00
#if istage =='II':  istagedec = 0.01
#if istage =='III': istagedec = 0.02
#if istage == 'IV': istagedec = 0.03
#if istage == 'V':  istagedec = 0.04
#if istage == 'VI': istagedec = 0.05

#idxk = np.where(nk == phys.ptable[elem]['atnum'] + istagedec)[0]
#idxv = np.where(nv == phys.ptable[elem]['atnum'] + istagedec)[0]

#wk = wk[idxk]
#wv = wv[idxv]

#gfk = gfk[idxk]
#gfv = gfv[idxv]

#elk = elk[idxk]
#elv = elv[idxv]

#euk = euk[idxk]
#euv = euv[idxv]

#grk = grk[idxk]
#grv = grv[idxv]

#gsk = gsk[idxk]
#gsv = gsv[idxv]

#gvk = gvk[idxk]
#gvv = gvv[idxv]

#np.savez('vald_kur_' + elemistage + '_comp', wk = wk,
#                                             wv = wv,
#                                             gfk = gfk,
#                                             gfv = gfv,
#                                             elk = elk,
#                                             elv = elv,
#                                             euk = euk,
#                                             euv = euv,
#                                             grk = grk,
#                                             grv = grv,
#                                             gsk = gsk,
#                                             gsv = gsv,
#                                             gvk = gvk,
#                                             gvv = gvv)

wk = np.load('vald_kur_' + elemistage + '_comp.npz')['wk']
wv = np.load('vald_kur_' + elemistage + '_comp.npz')['wv']

gfk = np.load('vald_kur_' + elemistage + '_comp.npz')['gfk']
gfv = np.load('vald_kur_' + elemistage + '_comp.npz')['gfv']

elk = np.load('vald_kur_' + elemistage + '_comp.npz')['elk']
elv = np.load('vald_kur_' + elemistage + '_comp.npz')['elv']

euk = np.load('vald_kur_' + elemistage + '_comp.npz')['euk']
euv = np.load('vald_kur_' + elemistage + '_comp.npz')['euv']

grk = np.load('vald_kur_' + elemistage + '_comp.npz')['grk']
grv = np.load('vald_kur_' + elemistage + '_comp.npz')['grv']

gsk = np.load('vald_kur_' + elemistage + '_comp.npz')['gsk']
gsv = np.load('vald_kur_' + elemistage + '_comp.npz')['gsv']

gvk = np.load('vald_kur_' + elemistage + '_comp.npz')['gvk']
gvv = np.load('vald_kur_' + elemistage + '_comp.npz')['gvv']

plt.close('all')

#fig, ax = plt.subplots(nrows = 6, ncols = 1, figsize = (6.0, 40.5))
fig, ax = plt.subplots(nrows = 2, ncols = 3, figsize = (18, 13.5))

auxplt.figpar(3, 3, 15)

#fig.tight_layout()

#fig.suptitle(elemistage, y = 1.00002)

ax[0, 0].scatter(wk, gfk, s = 5.00, label = 'Kurucz', color = 'k', marker = 'x')
ax[0, 0].scatter(wv, gfv, s = 0.75, label = 'VALD',   color = 'r')
ax[0, 1].scatter(wk, elk, s = 5.00, label = 'Kurucz', color = 'k', marker = 'x')
ax[0, 1].scatter(wv, elv, s = 0.75, label = 'VALD',   color = 'r')
ax[0, 2].scatter(wk, euk, s = 5.00, label = 'Kurucz', color = 'k', marker = 'x')
ax[0, 2].scatter(wv, euv, s = 0.75, label = 'VALD',   color = 'r')
ax[1, 0].scatter(wk, grk, s = 5.00, label = 'Kurucz', color = 'k', marker = 'x')
ax[1, 0].scatter(wv, grv, s = 0.75, label = 'VALD',   color = 'r')
ax[1, 1].scatter(wk, gsk, s = 5.00, label = 'Kurucz', color = 'k', marker = 'x')
ax[1, 1].scatter(wv, gsv, s = 0.75, label = 'VALD',   color = 'r')
ax[1, 2].scatter(wk, gvk, s = 5.00, label = 'Kurucz', color = 'k', marker = 'x')
ax[1, 2].scatter(wv, gvv, s = 0.75, label = 'VALD',   color = 'r')

#ax[0, 0].set_ylim(-16, 1)
ax[0, 0].set_ylim(-20, 20)

ax[0, 1].ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))
ax[0, 2].ticklabel_format(style = 'sci', axis = 'y', scilimits = (0, 0))

for i, j in itertools.product(range(2), range(3)):

#    ax[i, j].set_xlim(780, 850)
    ax[i, j].set_xlim(740, 750)

#    ax[i, j].xaxis.set_major_locator(MultipleLocator(10))
    ax[i, j].xaxis.set_major_locator(MultipleLocator(2))

#    leg = ax[i, j].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size': 12.0})

#    for obj in leg.legendHandles: obj.set_linewidth(3.0)

    ax[i, j].set_xlabel('Wavelength, [nm]')  

#ax[1, 0].set_xlabel('Wavelength, [nm]')
#ax[1, 1].set_xlabel('Wavelength, [nm]')
#ax[1, 2].set_xlabel('Wavelength, [nm]')

ax[0, 0].set_ylabel(r'$\log(gf)$')
ax[0, 1].set_ylabel(r'Lower level energy, [$\mathrm{cm}^{-1}$]')
ax[0, 2].set_ylabel(r'Upper level energy, [$\mathrm{cm}^{-1}$]')
ax[1, 0].set_ylabel(r'$\log(\gamma_\mathrm{rad})$')
ax[1, 1].set_ylabel(r'$\log(\gamma_\mathrm{stark})$')
ax[1, 2].set_ylabel(r'$\log(\gamma_\mathrm{Van-der-Waals})$')

auxplt.savepdf('vald_kur_' + elemistage + '_comp')
