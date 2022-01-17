import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import pylab             as plb
import numpy             as np
import oper_file         as of

import os
import itertools

import nessy_spec; reload(nessy_spec)
import spec;       reload(spec)
import pltaux;     reload(pltaux)
import sysaux;     reload(sysaux)
import paths;      reload(paths)
import auxfunc;    reload(auxfunc)

from progressbar import *

wvl, irr = of.read2(paths.inp + 'atlas3.dat')

irr_l = []; irr_u = []

for ind in range(len(wvl)):

    if wvl[ind] <= 121.0:

       irr_l.append(irr[ind] * (1.0 - 30.0 / 100.0))
       irr_u.append(irr[ind] * (1.0 + 30.0 / 100.0))

    if wvl[ind] > 121.0 and wvl[ind] <= 200.0:

       irr_l.append(irr[ind] * (1.0 - 3.5 / 100.0))
       irr_u.append(irr[ind] * (1.0 + 3.5 / 100.0))

    if wvl[ind] > 200.0 and wvl[ind] <= 400.0:

       irr_l.append(irr[ind] * (1.0 - 4.0 / 100.0))
       irr_u.append(irr[ind] * (1.0 + 4.0 / 100.0))

    if wvl[ind] > 400.0 and wvl[ind] <= 870.0:

       irr_l.append(irr[ind] * (1.0 - 2.0 / 100.0))
       irr_u.append(irr[ind] * (1.0 + 2.0 / 100.0))

    if wvl[ind] > 870.0:

       irr_l.append(irr[ind] * (1.0 - 3.0 / 100.0))
       irr_u.append(irr[ind] * (1.0 + 3.0 / 100.0))

of.write4(paths.out + 'atl_unc.dat', wvl, irr_l, irr, irr_u)

atlas = np.loadtxt(paths.out + 'atl_unc.dat')

wvl_atl = atlas[:, 0]
irr_atl = atlas[:, 2] * 1e-3

irr_atl_l = atlas[:, 1] * 1e-3
irr_atl_u = atlas[:, 3] * 1e-3

#atl_unc_u = (irr_atl_u / irr_atl - 1.0) * 100.
#atl_unc_l = (irr_atl_l / irr_atl - 1.0) * 100.

#idx1 = np.where((wvl_atl >= 170.0) & (wvl_atl <= 400.0))
#idx2 = np.where((wvl_atl >= 400.0) & (wvl_atl <= 1000.0))

#wvl_nes, irr_nes = nessy_spec.read(paths.it3f + 'NEW_SPEC_ALI_COL_HMI_STI_GAU', swav = 1700.0, ewav = 24000.0)

#irr_nes = nessy_spec.atl_con(paths.it3f + 'NEW_SPEC_ALI_COL_HMI_STI_GAU', wvl_nes, irr_nes)

#np.savez(paths.npz + 'new_spec', irr_nes = irr_nes)

irr_nes = np.load(paths.npz + 'new_spec.npz')['irr_nes']

delta = 1.5

wvls, irrs_atl = spec.running_mean(wvl_atl, irr_atl, delta)
wvls, irrs_nes = spec.running_mean(wvl_atl, irr_nes, delta)

wvls, irrs_atl_u = spec.running_mean(wvl_atl, irr_atl_u, delta)
wvls, irrs_atl_l = spec.running_mean(wvl_atl, irr_atl_l, delta)

atl_unc_u = (irrs_atl_u / irrs_atl - 1.0) * 100.
atl_unc_l = (irrs_atl_l / irrs_atl - 1.0) * 100.

idx1 = np.where((wvls >= 170.0) & (wvls <= 400.0))
idx2 = np.where((wvls >= 400.0) & (wvls <= 2400.0))

rel_dev = (irrs_nes / irrs_atl - 1.0) * 100.0

sysaux.clean_dir(paths.figdir, mode = 'verbose')

fontsize = 30

pltaux.figpar(xtick_maj_pad = 20, ytick_maj_pad = 15, fontsize = fontsize)

fig, ax = plt.subplots(nrows = 2, ncols = 2, figsize = (20.0, 15.0))

fig.tight_layout()

plt.subplots_adjust(wspace = 0.20, hspace = 0.1)

for i, j in itertools.product(range(2), range(2)):

    if (i == 0 and j == 0) or (i == 1 and j == 0):

        ax[i, j].set_xlim(170, 400)
        ax[i, j].xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax[i, j].xaxis.set_major_locator(ticker.MultipleLocator(50))
	   
    if (i == 0 and j == 1) or (i == 1 and j == 1):

        ax[i, j].set_xlim(400, 2400)
        ax[i, j].xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
        ax[i, j].xaxis.set_major_locator(ticker.MultipleLocator(400))

    if (i == 1 and j == 0) or (i == 1 and j == 1):

        ax[i, j].yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
        ax[i, j].set_xlabel(r'$\mathrm{Wavelength}, [\mathrm{nm}]$', labelpad = 20, fontsize = fontsize)

    if (i == 0 and j == 0) or (i == 0 and j == 1):

        ax[i, j].tick_params(labelbottom='off')
    
    ax[i, j].tick_params(axis='x', length = 10,  which = 'major')
    ax[i, j].tick_params(axis='x', length = 5,   which = 'minor')
    ax[i, j].tick_params(axis='y', length = 9,   which = 'major')
    ax[i, j].tick_params(axis='y', length = 4.5, which = 'minor')

ax[0, 0].fill_between(wvls[idx1], irrs_atl_l[idx1], irrs_atl_u[idx1], color = 'grey', label = 'ATLAS3', alpha = 0.5)

ax[0, 0].plot(wvls, irrs_nes, 'r', linewidth = 1.0, label = 'NESSY')

ax[0, 0].set_ylim(0.0, 1.5)

ax[0, 0].set_ylabel(r'$\mathrm{Irradiance}, [\mathrm{W} / \mathrm{m}^2 / \mathrm{nm}]$', labelpad = 20, fontsize = fontsize)

ax[0, 0].yaxis.set_minor_locator(ticker.AutoMinorLocator(4))

ax[0, 1].fill_between(wvls[idx2], irrs_atl_l[idx2], irrs_atl_u[idx2], color = 'grey', alpha = 0.5)

ax[0, 1].plot(wvls, irrs_nes, 'r', linewidth = 1.0)

ax[0, 1].set_ylim(0.0, 2.3)

ax[0, 1].yaxis.set_minor_locator(ticker.AutoMinorLocator(5))

#------------ inset plot start ---------------

idx3 = np.where((wvls >= 400.0) & (wvls <= 600.0))

left, bottom, width, height = [0.704, 0.695, 0.24, 0.25] # These are in unitless percentages of the figure size. (0,0 is bottom left)

ax_inset = fig.add_axes([left, bottom, width, height])

ax_inset.fill_between(wvls[idx3], irrs_atl_l[idx3], irrs_atl_u[idx3], color = 'grey', alpha = 0.5)

ax_inset.plot(wvls, irrs_nes, 'r', linewidth = 1.0)

ax_inset.set_xlim(400, 600)
ax_inset.set_ylim(1.2, 2.3)

ax_inset.tick_params(axis='x', length = 10,  which = 'major')
ax_inset.tick_params(axis='x', length = 5,   which = 'minor')
ax_inset.tick_params(axis='y', length = 9,   which = 'major')
ax_inset.tick_params(axis='y', length = 4.5, which = 'minor')

ax_inset.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
ax_inset.xaxis.set_major_locator(ticker.MultipleLocator(50))

ax_inset.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
ax_inset.yaxis.set_major_locator(ticker.MultipleLocator(0.4))

plt.xticks(fontsize = 19)
plt.yticks(fontsize = 19)

#------------ inset plot end ---------------

ax[1, 0].fill_between(wvls[idx1], atl_unc_l[idx1], atl_unc_u[idx1], color = 'grey', alpha = 0.5)

ax[1, 0].plot(wvls, np.zeros(len(wvls)), 'k--', linewidth = 1.0)
ax[1, 0].plot(wvls, rel_dev,             'r',   linewidth = 1.0)

ax[1, 0].set_ylim(-25, 25)

ax[1, 0].set_ylabel(r'$\mathrm{Deviation}, [\%]$', labelpad = 20, fontsize = fontsize)

ax[1, 1].fill_between(wvls[idx2], atl_unc_l[idx2], atl_unc_u[idx2], color = 'grey', alpha = 0.5)

ax[1, 1].plot(wvls, np.zeros(len(wvls)), 'k--', linewidth = 1.0)
ax[1, 1].plot(wvls, rel_dev,             'r',   linewidth = 1.0)

ax[1, 1].set_ylim(-10, 10)
ax[1, 1].yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
ax[1, 1].yaxis.set_major_locator(ticker.MultipleLocator(2))

leg = ax[0, 0].legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 50})

for handle in leg.legendHandles: handle.set_visible(False)

text_nes, text_atl = leg.get_texts()

plb.setp(text_nes, color = 'r')
plb.setp(text_atl, color = 'grey')

pltaux.savepdf(paths.figdir, 'nessy_vs_atl')
