import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import pylab             as plb
import numpy             as np

import importlib

from sys import exit

import pltaux;     importlib.reload(pltaux)
import sysaux;     importlib.reload(sysaux)
import paths;      importlib.reload(paths)
import oper_file;  importlib.reload(oper_file)

def cat(col1, col2, col3):

    col = np.zeros(len(col1) + len(col2) + len(col3))

    for i in range(len(col2)):

        j = 3 * i

        col[j] =     col1[i]
        col[j + 1] = col2[i]
        col[j + 2] = col3[i]

    col[len(col) - 1] = col1[len(col1) - 1]

    return col

col1_bol, col2_bol, col3_bol = oper_file.read3(paths.inp + 'rne_bol')
col1_sbe, col2_sbe, col3_sbe = oper_file.read3(paths.inp + 'rne_sbe')

rne_bol = cat(col1_bol, col2_bol, col3_bol)
rne_sbe = cat(col1_sbe, col2_sbe, col3_sbe)

rat = np.zeros(len(rne_sbe))

rat = rne_bol / rne_sbe

z = np.loadtxt(paths.inp + 'FCHHT_B', skiprows = 0, usecols = [0]) / 1000.0
T = np.loadtxt(paths.inp + 'FCHHT_B', skiprows = 0, usecols = [1])

fontsize = 25

pltaux.figpar(xtick_maj_pad = 20, ytick_maj_pad = 15, fontsize = fontsize)

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (20.0, 15.0))

fig.tight_layout()

plt.subplots_adjust(hspace = 0.02)

ax[0].plot(z, rne_bol, color = 'r', linewidth = 3, label = 'BOL')
ax[0].plot(z, rne_sbe, color = 'b', linewidth = 3, label = 'SBE')

ax[0].set_ylabel(r'$n_\mathrm{e} / n_\mathrm{tot}$', labelpad = 20)

axt = ax[0].twinx()

axt.plot(z, T, color = 'k', linewidth = 3)

axt.set_yscale('log')

axt.set_ylabel(r'$\mathrm{Temperature}, [\mathrm{K}]$', labelpad = 20)

axt.set_ylim(3e+3, 3e+5)

axt.tick_params(axis='y', length = 9,   width = 2, which = 'major')
axt.tick_params(axis='y', length = 4.5, width = 2, which = 'minor')

ax[1].plot(z, rat, color = 'k', linewidth = 3)

ax[1].set_ylabel('ratio', labelpad = 20)

for i in range(len(ax)):

    ax[i].set_xlim(min(z), 2.15)

    ax[i].tick_params(axis='x', length = 12.5, width = 2, which = 'major')
    ax[i].tick_params(axis='x', length = 7.5,  width = 2, which = 'minor')
    ax[i].tick_params(axis='y', length = 9,    width = 2, which = 'major')
    ax[i].tick_params(axis='y', length = 4.5,  width = 2, which = 'minor')

    ax[i].xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

    if i == 0:

        ax[i].set_yscale('log')
        ax[i].set_ylim(3e-5, 2e+0)
        ax[i].tick_params(labelbottom = 'off')

    if i == 1:

        ax[i].set_ylim(0.0, 1.7)

        ax[i].yaxis.set_minor_locator(ticker.AutoMinorLocator(4))

ax[1].set_xlabel('Height, [Mm]', labelpad = 20)

leg = ax[0].legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 40})

for handle in leg.legendHandles: handle.set_visible(False)

text_bol, text_sbe = leg.get_texts()

plb.setp(text_bol, color = 'r')
plb.setp(text_sbe, color = 'b')

pltaux.savepdf(paths.figdir, 'rne')
