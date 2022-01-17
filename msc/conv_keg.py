import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import pylab             as plb
import numpy             as np

import importlib

from sys import exit

import pltaux;     importlib.reload(pltaux)
#import sysaux;     importlib.reload(sysaux)
import paths;      importlib.reload(paths)
import oper_file;  importlib.reload(oper_file)

# temperature C
z_c_full = np.loadtxt(paths.it0h + 'conv_keg/c_full/FAL_VD',   usecols = [0])
T_c_full = np.loadtxt(paths.it0h + 'conv_keg/c_full/FAL_VD',   usecols = [1])

z_c_utmi = np.loadtxt(paths.it0h + 'conv_keg/c_utmi/FAL_VD',   usecols = [0])
T_c_utmi = np.loadtxt(paths.it0h + 'conv_keg/c_utmi/FAL_VD',   usecols = [1])

# convergence C
i_c_full = np.loadtxt(paths.it0h + 'conv_keg/c_full/CONV/ALL', skiprows = 2, usecols = [0])
c_c_full = np.loadtxt(paths.it0h + 'conv_keg/c_full/CONV/ALL', skiprows = 2, usecols = [4])

i_c_utmi = np.loadtxt(paths.it0h + 'conv_keg/c_utmi/CONV/ALL', skiprows = 2, usecols = [0])
c_c_utmi = np.loadtxt(paths.it0h + 'conv_keg/c_utmi/CONV/ALL', skiprows = 2, usecols = [4])

# temperature S
z_s_full = np.loadtxt(paths.it0h + 'conv_keg/s_full/FAL_VD',   usecols = [0])
T_s_full = np.loadtxt(paths.it0h + 'conv_keg/s_full/FAL_VD',   usecols = [1])

z_s_utmi = np.loadtxt(paths.it0h + 'conv_keg/s_utmi/FAL_VD',   usecols = [0])
T_s_utmi = np.loadtxt(paths.it0h + 'conv_keg/s_utmi/FAL_VD',   usecols = [1])

# convergence S
i_s_full = np.loadtxt(paths.it0h + 'conv_keg/s_full/CONV/ALL', skiprows = 2, usecols = [0])
c_s_full = np.loadtxt(paths.it0h + 'conv_keg/s_full/CONV/ALL', skiprows = 2, usecols = [4])

i_s_utmi = np.loadtxt(paths.it0h + 'conv_keg/s_utmi/CONV/ALL', skiprows = 2, usecols = [0])
c_s_utmi = np.loadtxt(paths.it0h + 'conv_keg/s_utmi/CONV/ALL', skiprows = 2, usecols = [4])

# temperature Kurucz
z_k = np.loadtxt(paths.it0h + 'conv_keg/k/kur_atm', skiprows = 7, usecols = [2])
T_k = np.loadtxt(paths.it0h + 'conv_keg/k/kur_atm', skiprows = 7, usecols = [7])

# convergence Kurucz
i_k = np.loadtxt(paths.it0h + 'conv_keg/k/CONV/ALL', skiprows = 2, usecols = [0])
c_k = np.loadtxt(paths.it0h + 'conv_keg/k/CONV/ALL', skiprows = 2, usecols = [4])

fontsize = 10

pltaux.figpar(xtick_maj_pad = 10, ytick_maj_pad = 7.5, fontsize = fontsize)

fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (10.0, 5.0))

fig.tight_layout()

plt.subplots_adjust(wspace = 0.2)

ax[0].plot(z_c_utmi, T_c_utmi, color = 'm', linewidth = 2.5)
ax[0].plot(z_s_utmi, T_s_utmi, color = 'k', linewidth = 2.5)

ax[0].plot(z_c_full, T_c_full, color = 'g', linewidth = 0.5, label = 'FAL99-C')
ax[0].plot(z_s_full, T_s_full, color = 'b', linewidth = 0.5, label = 'FAL99-S')
ax[0].plot(z_k,      T_k,      color = 'r', linewidth = 0.5, label = 'Kurucz')

ax[0].set_xlim(0, 2300)
ax[0].set_ylim(3200, 10000)

ax[0].set_xlabel('Height, [km]',     labelpad = 10)
ax[0].set_ylabel('Temperature, [K]', labelpad = 10)

ax[0].set_title('Temperature', fontsize = 15)

ax[1].plot(i_c_utmi, c_c_utmi, color = 'm', linewidth = 1.5)
ax[1].plot(i_s_utmi, c_s_utmi, color = 'k', linewidth = 1.5)

ax[1].plot(i_c_full, c_c_full, color = 'g', linewidth = 0.5, label = 'FAL99-C')
ax[1].plot(i_s_full, c_s_full, color = 'b', linewidth = 0.5, label = 'FAL99-S')
ax[1].plot(i_k,      c_k,      color = 'r', linewidth = 0.5, label = 'Kurucz')

ax[1].set_xlim(0, 100)
ax[1].set_ylim(1e-4, 1e+9)

ax[1].set_yscale('log')

ax[1].set_xlabel('Iteration number', labelpad = 10)
ax[1].set_ylabel('Hydrogen CORMAX',  labelpad = 10)

ax[1].set_title('Hydrogen convergence', fontsize = 15)

for i in range(len(ax)):

    ax[i].tick_params(axis='x', length = 6.25, width = 1, which = 'major')
    ax[i].tick_params(axis='x', length = 3.12, width = 1, which = 'minor')
    ax[i].tick_params(axis='y', length = 4.5,  width = 1, which = 'major')
    ax[i].tick_params(axis='y', length = 4.5,  width = 1, which = 'minor')

    if i == 0: ax[i].xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
    if i == 1: ax[i].xaxis.set_major_locator(ticker.MultipleLocator(10))

leg = ax[0].legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 20})

for handle in leg.legendHandles: handle.set_visible(False)

text_c, text_s, text_k = leg.get_texts()

plb.setp(text_c, color = 'g')
plb.setp(text_s, color = 'b')
plb.setp(text_k, color = 'r')

pltaux.savepdf(paths.figdir, 'conv_keg')
