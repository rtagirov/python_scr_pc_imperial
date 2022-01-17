#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import paths as p
import matplotlib.ticker as ticker

import pylab as plb

import nessy_spec
import pltaux

#wvl_icl, flu_icl = nessy_spec.read(p.it3f + 'TST_ICL', 1005., 9995.)
#wvl_keg, flu_keg = nessy_spec.read(p.it3f + 'TST_KEG', 1005., 9995.)

#wvl_icl = wvl_icl / 10.0
#wvl_keg = wvl_keg / 10.0

#np.savez(p.npz + 'icl_keg_spec', wvl_icl = wvl_icl, flu_icl = flu_icl, \
#                                 wvl_keg = wvl_keg, flu_keg = flu_keg)

#icl_keg_spec = np.load(p.npz + 'icl_keg_spec.npz')

#wvl_icl = icl_keg_spec['wvl_icl']
#flu_icl = icl_keg_spec['flu_icl']
#wvl_keg = icl_keg_spec['wvl_keg']
#flu_keg = icl_keg_spec['flu_keg']

#wvls_icl, flus_icl = nessy_spec.atl_con(wvl_icl * 10.0, flu_icl, 0.6 * 3.5)
#wvls_keg, flus_keg = nessy_spec.atl_con(wvl_keg * 10.0, flu_keg, 0.6 * 3.5)

#wvls_icl = wvls_icl / 10.0
#wvls_keg = wvls_keg / 10.0

#np.savez(p.npz + 'icl_keg_spec_atl_con', wvls_icl = wvls_icl, flus_icl = flus_icl, \
#                                         wvls_keg = wvls_keg, flus_keg = flus_keg)

icl_keg_spec_atl_con = np.load(p.npz + 'icl_keg_spec_atl_con.npz')

wvls_icl = icl_keg_spec_atl_con['wvls_icl']
flus_icl = icl_keg_spec_atl_con['flus_icl']
wvls_keg = icl_keg_spec_atl_con['wvls_keg']
flus_keg = icl_keg_spec_atl_con['flus_keg']

plb.rcParams['xtick.major.pad'] = 15
plb.rcParams['ytick.major.pad'] = 10

#plt.rcParams['font.size'] = 20
#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = 'Times'
#plt.rcParams['font.monospace'] = 'Ubuntu Mono'

plt.rc('font', size = 20)
plt.rc('font', family = 'serif')
#plt.rc('font', serif = 'Times New Roman')

fig, ax = plt.subplots(nrows = 5, ncols = 1, figsize = (10.0, 14.5))

fig.tight_layout()

x1 = wvls_keg
x2 = wvls_icl

xmin = [125, 160, 240, 320, 400]
xmax = [160, 240, 320, 400, 480]

ymax = [5.00, 0.08, 1.00, 2.00, 2.50]

for i in range(len(ax)):

    if i == 0: 

       y1 = 1e+4 * flus_keg
       y2 = 1e+4 * flus_icl

    if i != 0: 

       y1 = flus_keg
       y2 = flus_icl

    ax[i].plot(x2, y2, 'r', linewidth = 2.5, label = 'imperial')
    ax[i].plot(x1, y1, 'k', linewidth = 0.5, label = 'pmod/wrc')

    if i == 0: ax[i].set_ylabel(r'$\mathrm{Flux}\times 10^{4}, [\mathrm{W} / \mathrm{m}^2 / \mathrm{nm}]$', \
                                labelpad = 20, \
                                fontsize = 15)

    if i != 0: ax[i].set_ylabel(r'$\mathrm{Flux}, [\mathrm{W} / \mathrm{m}^2 / \mathrm{nm}]$', \
                                labelpad = 20, \
                                fontsize = 20)

    ax[i].set_xlim(xmin[i], xmax[i])
    ax[i].set_ylim(0.0,     ymax[i])

    if i == 0: minloc_x = ticker.AutoMinorLocator(5)
    if i != 0: minloc_x = ticker.AutoMinorLocator(4)

    if i <= 2: minloc_y = ticker.AutoMinorLocator(4)
    if i >  2: minloc_y = ticker.AutoMinorLocator(5)

    ax[i].xaxis.set_minor_locator(minloc_x)
    ax[i].yaxis.set_minor_locator(minloc_y)

    if i != 0: ax[i].xaxis.set_major_locator(ticker.MultipleLocator(20))
    if i == 1: ax[i].yaxis.set_major_locator(ticker.MultipleLocator(0.02))

    ax[i].tick_params(axis='x', length = 10,  which = 'major')
    ax[i].tick_params(axis='x', length = 5,   which = 'minor')
    ax[i].tick_params(axis='y', length = 9,   which = 'major')
    ax[i].tick_params(axis='y', length = 4.5, which = 'minor')

    if i == 4: ax[i].set_xlabel('Wavelength, [nm]', labelpad = 20)

leg = ax[0].legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 20})

for handle in leg.legendHandles:
    handle.set_visible(False)

text_icl, text_pmo = leg.get_texts()

plb.setp(text_icl, color = 'r')
plb.setp(text_pmo, color = 'k')

pltaux.savepdf(p.figdir, 'icl_keg_spec_IT3_1')

fig, ax = plt.subplots(nrows = 5, ncols = 1, figsize = (10.0, 14.5))

fig.tight_layout()

xmin = [480, 560, 640, 720, 800]
xmax = [560, 640, 720, 800, 880]

ymax = [2.5, 2.0, 2.0, 1.5, 1.2]

y1 = flus_keg
y2 = flus_icl

for i in range(len(ax)):

    ax[i].plot(x2, y2, 'r', linewidth = 2.5, label = 'imperial')
    ax[i].plot(x1, y1, 'k', linewidth = 0.5, label = 'pmod/wrc')

    ax[i].set_ylabel(r'$\mathrm{Flux}, [\mathrm{W} / \mathrm{m}^2 / \mathrm{nm}]$', \
                     labelpad = 20, \
                     fontsize = 20)

    ax[i].set_xlim(xmin[i], xmax[i])
    ax[i].set_ylim(0.0,     ymax[i])

    minloc_x = ticker.AutoMinorLocator(4)

    if i != 4: minloc_y = ticker.AutoMinorLocator(5)
    if i == 4: minloc_y = ticker.AutoMinorLocator(4)

    ax[i].xaxis.set_minor_locator(minloc_x)
    ax[i].yaxis.set_minor_locator(minloc_y)

    ax[i].xaxis.set_major_locator(ticker.MultipleLocator(20))
    if i == 3: ax[i].yaxis.set_major_locator(ticker.MultipleLocator(0.5))

    ax[i].tick_params(axis='x', length = 10,  which = 'major')
    ax[i].tick_params(axis='x', length = 5,   which = 'minor')
    ax[i].tick_params(axis='y', length = 9,   which = 'major')
    ax[i].tick_params(axis='y', length = 4.5, which = 'minor')

    if i == 4: ax[i].set_xlabel('Wavelength, [nm]', labelpad = 20)

leg = ax[0].legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 20})

for handle in leg.legendHandles:
    handle.set_visible(False)

text_icl, text_pmo = leg.get_texts()

plb.setp(text_icl, color = 'r')
plb.setp(text_pmo, color = 'k')

pltaux.savepdf(p.figdir, 'icl_keg_spec_IT3_2')
