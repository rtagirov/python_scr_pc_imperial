import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import pylab             as plb
import numpy             as np

import importlib
import itertools

from sys import exit

import pltaux;     importlib.reload(pltaux)
import paths;      importlib.reload(paths)

h1 = np.loadtxt(paths.it0h + 'redmur/mur/1/ATM_STR',  skiprows = 2, usecols = [1])
T1 = np.loadtxt(paths.it0h + 'redmur/mur/1/ATM_STR',  skiprows = 2, usecols = [2])
n1 = np.loadtxt(paths.it0h + 'redmur/mur/1/ATM_STR',  skiprows = 2, usecols = [3])
i1 = np.loadtxt(paths.it0h + 'redmur/mur/1/CONV/ALL', skiprows = 2, usecols = [0])
c1 = np.loadtxt(paths.it0h + 'redmur/mur/1/CONV/ALL', skiprows = 2, usecols = [4])

h2 = np.loadtxt(paths.it0h + 'redmur/mur/2/ATM_STR',  skiprows = 2, usecols = [1])
T2 = np.loadtxt(paths.it0h + 'redmur/mur/2/ATM_STR',  skiprows = 2, usecols = [2])
n2 = np.loadtxt(paths.it0h + 'redmur/mur/2/ATM_STR',  skiprows = 2, usecols = [3])
i2 = np.loadtxt(paths.it0h + 'redmur/mur/2/CONV/ALL', skiprows = 2, usecols = [0])
c2 = np.loadtxt(paths.it0h + 'redmur/mur/2/CONV/ALL', skiprows = 2, usecols = [4])

h3 = np.loadtxt(paths.it0h + 'redmur/mur/3/ATM_STR',  skiprows = 2, usecols = [1])
T3 = np.loadtxt(paths.it0h + 'redmur/mur/3/ATM_STR',  skiprows = 2, usecols = [2])
n3 = np.loadtxt(paths.it0h + 'redmur/mur/3/ATM_STR',  skiprows = 2, usecols = [3])
i3 = np.loadtxt(paths.it0h + 'redmur/mur/3/CONV/ALL', skiprows = 2, usecols = [0])
c3 = np.loadtxt(paths.it0h + 'redmur/mur/3/CONV/ALL', skiprows = 2, usecols = [4])

h4 = np.loadtxt(paths.it0h + 'redmur/mur/4/ATM_STR',  skiprows = 2, usecols = [1])
T4 = np.loadtxt(paths.it0h + 'redmur/mur/4/ATM_STR',  skiprows = 2, usecols = [2])
n4 = np.loadtxt(paths.it0h + 'redmur/mur/4/ATM_STR',  skiprows = 2, usecols = [3])
i4 = np.loadtxt(paths.it0h + 'redmur/mur/4/CONV/ALL', skiprows = 2, usecols = [0])
c4 = np.loadtxt(paths.it0h + 'redmur/mur/4/CONV/ALL', skiprows = 2, usecols = [4])

h5 = np.loadtxt(paths.it0h + 'redmur/mur/5/ATM_STR',  skiprows = 2, usecols = [1])
T5 = np.loadtxt(paths.it0h + 'redmur/mur/5/ATM_STR',  skiprows = 2, usecols = [2])
n5 = np.loadtxt(paths.it0h + 'redmur/mur/5/ATM_STR',  skiprows = 2, usecols = [3])
i5 = np.loadtxt(paths.it0h + 'redmur/mur/5/CONV/ALL', skiprows = 2, usecols = [0])
c5 = np.loadtxt(paths.it0h + 'redmur/mur/5/CONV/ALL', skiprows = 2, usecols = [4])

cl = np.ones(len(i1)) * 1e-4

fontsize = 10

pltaux.figpar(xtick_maj_pad = 5, ytick_maj_pad = 5, fontsize = fontsize)

fig = plt.figure(figsize = (7.0, 6.3))

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0), colspan = 2)

plt.subplots_adjust(wspace = 0.45, hspace = 0.35)

#fig.tight_layout()

ax1.plot(h1, T1, linewidth = 0.5)
ax1.plot(h2, T2, linewidth = 0.5)
ax1.plot(h3, T3, linewidth = 0.5)
ax1.plot(h4, T4, linewidth = 0.5)
ax1.plot(h5, T5, linewidth = 0.5)

#ax1.set_xlim()
#ax1.set_ylim()

ax1.set_xlabel('Height, [km]',     labelpad = 5)
ax1.set_ylabel('Temperature, [K]', labelpad = 5)

ax2.plot(h1, n1, linewidth = 0.5)
ax2.plot(h2, n2, linewidth = 0.5)
ax2.plot(h3, n3, linewidth = 0.5)
ax2.plot(h4, n4, linewidth = 0.5)
ax2.plot(h5, n5, linewidth = 0.5)

#ax2.set_xlim()
#ax2.set_ylim()

ax2.set_yscale('log')

ax2.minorticks_off()

ax2.set_xlabel('Height, [km]',                labelpad = 5)
ax2.set_ylabel('Number density, [cm$^{-3}$]', labelpad = 5)

ax3.plot(i1, c1, linewidth = 0.5)
ax3.plot(i2, c2, linewidth = 0.5)
ax3.plot(i3, c3, linewidth = 0.5)
ax3.plot(i4, c4, linewidth = 0.5)
ax3.plot(i5, c5, linewidth = 0.5)

ax3.plot(i1, cl, 'k--', linewidth = 0.5)

ax3.set_xlim(0, 320)
ax3.set_ylim(1e-10, 1e+1)

ax3.set_yscale('log')

ax3.minorticks_off()

ax3.set_xlabel('Iteration number', labelpad = 5)
ax3.set_ylabel('CORMAX',  labelpad = 5)

ax1.tick_params(axis='x', length = 6.25, width = 1, which = 'major')
ax1.tick_params(axis='x', length = 3.12, width = 1, which = 'minor')
ax1.tick_params(axis='y', length = 4.5,  width = 1, which = 'major')
ax1.tick_params(axis='y', length = 4.5,  width = 1, which = 'minor')

ax2.tick_params(axis='x', length = 6.25, width = 1, which = 'major')
ax2.tick_params(axis='x', length = 3.12, width = 1, which = 'minor')
ax2.tick_params(axis='y', length = 4.5,  width = 1, which = 'major')
ax2.tick_params(axis='y', length = 4.5,  width = 1, which = 'minor')

ax3.tick_params(axis='x', length = 6.25, width = 1, which = 'major')
ax3.tick_params(axis='x', length = 3.12, width = 1, which = 'minor')
ax3.tick_params(axis='y', length = 4.5,  width = 1, which = 'major')
ax3.tick_params(axis='y', length = 4.5,  width = 1, which = 'minor')

ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
ax3.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))

#leg1 = ax1.legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 10})
#leg2 = ax3.legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 10})

#for handle in leg1.legendHandles: handle.set_visible(False)
#for handle in leg2.legendHandles: handle.set_visible(False)

#text_c, text_kn = leg1.get_texts()

#plb.setp(text_c,  color = 'g')
#plb.setp(text_kn, color = 'r')

#text_kn, text_ctm, text_c = leg2.get_texts()

#plb.setp(text_kn,  color = 'r')
#plb.setp(text_c,   color = 'g')
#plb.setp(text_ctm, color = 'k')

pltaux.savepdf(paths.figdir, 'mur_conv_mur')
