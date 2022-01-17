import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pylab             as plb
import numpy             as np

import os
import importlib

#import nessy_spec; importlib.reload(nessy_spec)
import pltaux;     importlib.reload(pltaux)
import sysaux;     importlib.reload(sysaux)
import paths;      importlib.reload(paths)
import auxfunc;    importlib.reload(auxfunc)

from   scipy.io.idl import readsav

data_o  = readsav(paths.idlout + 'INTS_O.sav')

data_n1 = readsav(paths.idlout + 'INTS_N1.sav')
data_n2 = readsav(paths.idlout + 'INTS_N2.sav')
data_n3 = readsav(paths.idlout + 'INTS_N3.sav')
data_n4 = readsav(paths.idlout + 'INTS_N4.sav')

wavs_o = data_o['wavs_o']
ints_o = data_o['ints_o']

wavs_n1 = data_n1['wavs_n1']
ints_n1 = data_n1['ints_n1']

wavs_n2 = data_n2['wavs_n2']
ints_n2 = data_n2['ints_n2']

wavs_n3 = data_n3['wavs_n3']
ints_n3 = data_n3['ints_n3']

wavs_n4 = data_n4['wavs_n4']
ints_n4 = data_n4['ints_n4']

wavs_n1 = wavs_n1 / 10.0
wavs_n2 = wavs_n2 / 10.0
wavs_n3 = wavs_n3 / 10.0
wavs_n4 = wavs_n4 / 10.0

r1 = (ints_n1 / ints_o - 1.0) * 100.0
r2 = (ints_n2 / ints_o - 1.0) * 100.0
r3 = (ints_n3 / ints_o - 1.0) * 100.0
r4 = (ints_n4 / ints_o - 1.0) * 100.0

z = np.zeros(len(wavs_n1))

idx = np.where((wavs_n1 >= 90.0) & (wavs_n1 <= 125.0))

f = np.ones(len(wavs_n1[idx])) * 5.0
h = np.ones(len(wavs_n1[idx])) * (-100.0)

sysaux.clean_dir(paths.figdir, mode = 'verbose')

pltaux.figpar()

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

plt.subplots_adjust(hspace = 0.30)

#ax[0].plot(wavs_n1, r1, 'k',   linewidth = 1.5, label = 'ALI')

ax[0].fill_between(wavs_n1[idx], h, f, color = 'grey', alpha = 0.5)

ax[0].plot(wavs_n2, r2, 'r', linewidth = 1.5, label = 'ALI')
ax[0].plot(wavs_n3, r3, 'g', linewidth = 1.5, label = 'ALI + COL')
ax[0].plot(wavs_n4, r4, 'b', linewidth = 1.5, label = 'ALI + COL + C.1')

ax[0].plot(wavs_n1, z,  'k--', linewidth = 0.5)

ax[0].set_xlim( 90,  160)
ax[0].set_ylim(-100, 5)

ax[0].xaxis.set_minor_locator(ticker.AutoMinorLocator(4))
ax[0].yaxis.set_minor_locator(ticker.AutoMinorLocator(4))

ax[0].set_xlabel(r'$\mathrm{Wavelength}, [\mathrm{nm}]$', labelpad = 10)
ax[0].set_ylabel(r'$\mathrm{Deviation},  [\mathrm{\%}]$', labelpad = 10)

#ax[1].plot(wavs_n1, r1, 'k', linewidth = 1.5, label = 'ALI')
ax[1].plot(wavs_n4, r4, 'b', linewidth = 0.5, label = 'ALI + STI + COL + HMI')
ax[1].plot(wavs_n3, r3, 'g', linewidth = 1.0, label = 'ALI + STI + COL')
ax[1].plot(wavs_n2, r2, 'r', linewidth = 0.5, label = 'ALI + STI')

ax[1].plot(wavs_n1, z,  'k--', linewidth = 0.5)

ax[1].set_xlim(160, 4000)
ax[1].set_ylim(-3,  6)

ax[1].set_xscale('log')

ax[1].set_xticks([160, 300, 500, 700, 1000, 2000, 3000, 4000])
ax[1].get_xaxis().set_major_formatter(ticker.ScalarFormatter())

ax[1].yaxis.set_minor_locator(ticker.AutoMinorLocator(4))

ax[1].set_xlabel(r'$\mathrm{Wavelength}, [\mathrm{nm}]$', labelpad = 10)
ax[1].set_ylabel(r'$\mathrm{Deviation},  [\mathrm{\%}]$', labelpad = 10)

#plt.xscale('log')

for i in range(2):

    ax[i].tick_params(axis='x', length = 8, direction = 'in', top = 'on', right = 'on', which = 'major')
    ax[i].tick_params(axis='x', length = 5, direction = 'in', top = 'on', right = 'on', which = 'minor')
    ax[i].tick_params(axis='y', length = 8, direction = 'in', top = 'on', right = 'on', which = 'major')
    ax[i].tick_params(axis='y', length = 5, direction = 'in', top = 'on', right = 'on', which = 'minor')

leg = ax[0].legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 25})

for handle in leg.legendHandles: handle.set_visible(False)

#text_a, text_as, text_asc, text_asch = leg.get_texts()
text_as, text_asc, text_asch = leg.get_texts()

#plb.setp(text_a,    color = 'k')
plb.setp(text_as,   color = 'r')
plb.setp(text_asc,  color = 'g')
plb.setp(text_asch, color = 'b')

pltaux.savepdf(paths.figdir, 'new_vs_old')
