#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker

import pylab as plb

import nessy_spec
import pltaux
import sysaux
import path
import os

icl_keg_spec =              np.load(path.npz + 'icl_keg_spec.npz')
icl_keg_spec_atl_con =      np.load(path.npz + 'icl_keg_spec_atl_con.npz')
icl_keg_spec_lodf =         np.load(path.npz + 'icl_keg_spec_lodf.npz')
icl_keg_spec_lodf_atl_con = np.load(path.npz + 'icl_keg_spec_lodf_atl_con.npz')

wvl =  icl_keg_spec_lodf['wvl_icl']
wvls = icl_keg_spec_lodf_atl_con['wvls_icl']

flu_it3 =  icl_keg_spec['flu_icl']
flus_it3 = icl_keg_spec_atl_con['flus_icl']

flu_lod =  icl_keg_spec_lodf['flu_icl']
flus_lod = icl_keg_spec_lodf_atl_con['flus_icl']

sysaux.clean_dir(path.figdir)

pltaux.figpar()

fig, ax = plt.subplots(nrows = 5, ncols = 1, figsize = (10.0, 14.5))

fig.tight_layout()

rat =  (flu_it3  - flu_lod)  * 100 / flu_it3
rats = (flus_it3 - flus_lod) * 100 / flus_it3

xmin = [125, 160, 240, 320, 400]
xmax = [160, 240, 320, 400, 480]

for i in range(len(ax)):

    ax[i].plot(wvl,  rat,  'k', linewidth = 0.5, label = 'high resolution')
    ax[i].plot(wvls, rats, 'r', linewidth = 2.5, label = 'atlas convolution')

    ax[i].set_ylabel('(IT3 - LODF) / IT3, [%]', labelpad = 20, fontsize = 15)

    ax[i].set_xlim(xmin[i], xmax[i])

    idx = np.where((wvl > xmin[i]) & (wvl < xmax[i]))

    ax[i].set_ylim(min(rat[idx]), max(rat[idx]))

    if i == 0: minloc_x = ticker.AutoMinorLocator(5)
    if i != 0: minloc_x = ticker.AutoMinorLocator(4)

#    if i <= 2: minloc_y = ticker.AutoMinorLocator(4)
    minloc_y = ticker.AutoMinorLocator(5)

    ax[i].xaxis.set_minor_locator(minloc_x)
    ax[i].yaxis.set_minor_locator(minloc_y)

    if i != 0: ax[i].xaxis.set_major_locator(ticker.MultipleLocator(20))
#    if i == 1: ax[i].yaxis.set_major_locator(ticker.MultipleLocator(0.02))

    ax[i].tick_params(axis='x', length = 10,  which = 'major')
    ax[i].tick_params(axis='x', length = 5,   which = 'minor')
    ax[i].tick_params(axis='y', length = 9,   which = 'major')
    ax[i].tick_params(axis='y', length = 4.5, which = 'minor')

    if i == len(ax) - 1: ax[i].set_xlabel('Wavelength, [nm]', labelpad = 20)

leg = ax[0].legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 20})

for handle in leg.legendHandles:
    handle.set_visible(False)

text_hir, text_atc = leg.get_texts()

plb.setp(text_hir, color = 'k')
plb.setp(text_atc, color = 'r')

pltaux.savepdf(path.figdir, 'it3_lod_1')

fig, ax = plt.subplots(nrows = 5, ncols = 1, figsize = (10.0, 14.5))

fig.tight_layout()

xmin = [480, 560, 640, 720, 800]
xmax = [560, 640, 720, 800, 880]

for i in range(len(ax)):

    ax[i].plot(wvl,  rat,  'k', linewidth = 0.5, label = 'high resolution')
    ax[i].plot(wvls, rats, 'r', linewidth = 2.5, label = 'atlas convolution')

    ax[i].set_ylabel('(IT3 - LODF) / IT3, [%]', labelpad = 20, fontsize = 15)

    ax[i].set_xlim(xmin[i], xmax[i])

    idx = np.where((wvl > xmin[i]) & (wvl < xmax[i]))

    ax[i].set_ylim(min(rat[idx]), max(rat[idx]))

    if i == 0: minloc_x = ticker.AutoMinorLocator(5)
    if i != 0: minloc_x = ticker.AutoMinorLocator(4)

#    if i <= 2: minloc_y = ticker.AutoMinorLocator(4)
    minloc_y = ticker.AutoMinorLocator(5)

    ax[i].xaxis.set_minor_locator(minloc_x)
    ax[i].yaxis.set_minor_locator(minloc_y)

    if i != 0: ax[i].xaxis.set_major_locator(ticker.MultipleLocator(20))
#    if i == 1: ax[i].yaxis.set_major_locator(ticker.MultipleLocator(0.02))

    ax[i].tick_params(axis='x', length = 10,  which = 'major')
    ax[i].tick_params(axis='x', length = 5,   which = 'minor')
    ax[i].tick_params(axis='y', length = 9,   which = 'major')
    ax[i].tick_params(axis='y', length = 4.5, which = 'minor')

    if i == len(ax) - 1: ax[i].set_xlabel('Wavelength, [nm]', labelpad = 20)

leg = ax[0].legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 20})

for handle in leg.legendHandles:
    handle.set_visible(False)

text_hir, text_atc = leg.get_texts()

plb.setp(text_hir, color = 'k')
plb.setp(text_atc, color = 'r')

pltaux.savepdf(path.figdir, 'it3_lod_2')

os.system('pdftk ' + path.figdir + '*' + ' output ' + path.figdir + 'it3_lod.pdf')
