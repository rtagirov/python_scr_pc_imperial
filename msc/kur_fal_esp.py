import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import pylab             as plb
import numpy             as np

import importlib
import itertools

from sys import exit

import pltaux;     importlib.reload(pltaux)
import paths;      importlib.reload(paths)
import oper_file;  importlib.reload(oper_file)

# FAL99-C
z_c_full = np.loadtxt(paths.it0h + 'esc_prob/c_full/FAL_VD', usecols = [0]) / 1000.
T_c_full = np.loadtxt(paths.it0h + 'esc_prob/c_full/FAL_VD', usecols = [1])
d_c_full = np.loadtxt(paths.it0h + 'esc_prob/c_full/FAL_VD', usecols = [3])

z_c_utmi = np.loadtxt(paths.it0h + 'esc_prob/c_utmi/FAL_VD', usecols = [0]) / 1000.
T_c_utmi = np.loadtxt(paths.it0h + 'esc_prob/c_utmi/FAL_VD', usecols = [1])
d_c_utmi = np.loadtxt(paths.it0h + 'esc_prob/c_utmi/FAL_VD', usecols = [3])

# FAL99-C convergence
i_c_full = np.loadtxt(paths.it0h + 'esc_prob/c_full/CONV/ALL', skiprows = 2, usecols = [0])
c_c_full = np.loadtxt(paths.it0h + 'esc_prob/c_full/CONV/ALL', skiprows = 2, usecols = [3])

i_c_utmi = np.loadtxt(paths.it0h + 'esc_prob/c_utmi/CONV/ALL', skiprows = 2, usecols = [0])
c_c_utmi = np.loadtxt(paths.it0h + 'esc_prob/c_utmi/CONV/ALL', skiprows = 2, usecols = [3])

# Kurucz
z_k = np.loadtxt(paths.it0h + 'esc_prob/k/kur_atm', skiprows = 7, usecols = [2]) / 1000.
T_k = np.loadtxt(paths.it0h + 'esc_prob/k/kur_atm', skiprows = 7, usecols = [7])

d_k = 10.0**np.loadtxt(paths.it0h + 'esc_prob/k/kur_atm', skiprows = 7, usecols = [6])

#Kurucz convergence
i_k = np.loadtxt(paths.it0h + 'esc_prob/k/CONV/ALL', skiprows = 2, usecols = [0])
c_k = np.loadtxt(paths.it0h + 'esc_prob/k/CONV/ALL', skiprows = 2, usecols = [3])

cl = np.ones(len(i_c_utmi)) * 1e-4

fontsize = 10

pltaux.figpar(xtick_maj_pad = 5, ytick_maj_pad = 5, fontsize = fontsize)

fig = plt.figure(figsize = (7.0, 6.3))

#ax[0, 0].plot(z_c_utmi, T_c_utmi, color = 'm', linewidth = 2.5)

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax2 = plt.subplot2grid((2, 2), (0, 1))
ax3 = plt.subplot2grid((2, 2), (1, 0), colspan = 2)

plt.subplots_adjust(wspace = 0.45, hspace = 0.35)

#fig.tight_layout()

ax1.plot(z_c_utmi, T_c_utmi, color = 'k', linewidth = 2.5)
ax1.plot(z_c_full, T_c_full, color = 'g', linewidth = 0.5, label = 'FAL99-C')
ax1.plot(z_k,      T_k,  color = 'r', linewidth = 0.5, label = 'Kurucz (ICL)')

ax1.set_xlim(0, 2.3)
ax1.set_ylim(3200, 10000)

ax1.set_xlabel('Height, [km]',     labelpad = 5)
ax1.set_ylabel('Temperature, [K]', labelpad = 5)

ax2.plot(z_c_utmi, d_c_utmi, color = 'k', linewidth = 2.5)
ax2.plot(z_c_full, d_c_full, color = 'g', linewidth = 0.5)
ax2.plot(z_k,      d_k,      color = 'r', linewidth = 0.5)

ax2.set_xlim(0, 2.3)
ax2.set_ylim(1e+10, 1e+18)

ax2.set_yscale('log')

ax2.minorticks_off()

ax2.set_xlabel('Height, [km]',                labelpad = 5)
ax2.set_ylabel('Number density, [cm$^{-3}$]', labelpad = 5)

ax3.plot(i_k, c_k, color = 'r', linewidth = 0.5, label = 'Kurucz (ICL)')

ax3.plot(i_c_utmi, c_c_utmi, color = 'k', linewidth = 0.5, label = 'FAL99-C (T$_\mathrm{min}$)')
ax3.plot(i_c_full, c_c_full, color = 'g', linewidth = 0.5, label = 'FAL99-C')

ax3.plot(i_c_utmi, cl, 'k--', linewidth = 0.5)

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

leg1 = ax1.legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 10})
leg2 = ax3.legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 10})

for handle in leg1.legendHandles: handle.set_visible(False)
for handle in leg2.legendHandles: handle.set_visible(False)

text_c, text_kn = leg1.get_texts()

plb.setp(text_c,  color = 'g')
plb.setp(text_kn, color = 'r')

text_kn, text_ctm, text_c = leg2.get_texts()

plb.setp(text_kn,  color = 'r')
plb.setp(text_c,   color = 'g')
plb.setp(text_ctm, color = 'k')

pltaux.savepdf(paths.figdir, 'fal99c_vs_kurucz_esp_nhn')
