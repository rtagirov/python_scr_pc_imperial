import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import FixedLocator

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import auxplt
import paths

ei_old = np.loadtxt(paths.inp + 'eion_old.out', usecols = [6])
ei_new = np.loadtxt(paths.inp + 'eion_new.out', usecols = [6])

el_old = np.loadtxt(paths.inp + 'eion_old.out', usecols = [5])
el_new = np.loadtxt(paths.inp + 'eion_new.out', usecols = [5])

ch = np.loadtxt(paths.inp + 'eion_new.out', usecols = [3])
ln = np.loadtxt(paths.inp + 'eion_new.out', usecols = [0])

lnam = []

f = open(paths.inp + 'eion_new.out', 'r')

for line in f:

    lnam.append(line.split()[1])

for i in range(len(lnam)):

    lnam[i] = lnam[i].replace('.', '').replace('-','').replace('I', '')

el_old_1 = el_old[np.where(ch == 1)]
el_new_1 = el_new[np.where(ch == 1)]

j = 0

for i in range(len(el_old)):

    if ch[i] == 1:

        continue

    el_old[i] += el_old_1[j]
    el_new[i] += el_new_1[j]

    if ch[i] != ch[i + 1]:

        j += 1

T = np.loadtxt(paths.nessy + '/inp/atm/fal/FAL99_C', usecols = [1])
h = np.loadtxt(paths.nessy + '/inp/atm/fal/FAL99_C', usecols = [0])

delta_E = el_new[np.where(ch == 0)] - el_old[np.where(ch == 0)]

delta_E[np.where(delta_E == 0.0)] = np.nan

const = 1.4388

lnn = ln[np.where(ch == 0)]

lnamn = []

for i in range(len(lnam)):

    if ch[i] == 1:

        continue

    lnamn.append(lnam[i])

plt.close('all')

fig, ax1 = plt.subplots(nrows = 2, ncols = 1, figsize = (28.0, 20))

fig.tight_layout()

ax1[0].scatter(lnn, delta_E, s = 200)

for i in range(len(lnn)):

    if np.isnan(delta_E[i]):

        continue

    ax1[0].axvline(x = lnn[i], linewidth = 0.5)

    if delta_E[i] > 2:

        ax1[0].text(lnn[i] - 0.3, -20, lnamn[i], fontsize = 12)

        if lnamn[i][0 : 2] != 'Fe' and lnamn[i][0 : 2] != 'Si':

            ax1[1].plot(h, np.exp(const * delta_E[i] / T), color = 'k')

        if lnamn[i][0 : 2] == 'Si':

            ax1[1].plot(h, np.exp(const * delta_E[i] / T), color = 'r', linewidth = 4)

        if lnamn[i][0 : 2] == 'Fe':

            ax1[1].plot(h, np.exp(const * delta_E[i] / T), color = 'g', linewidth = 4)

ax1[0].text(41, 750, r'$\chi_i = \chi_0 - \chi_{0, i} + \chi_{1, 1}$', bbox = dict(facecolor = 'red', alpha = 0.5), fontsize = 30)
ax1[0].text(41, 650, r'$\chi_i^\dagger = \chi_0^\dagger - \chi_{0, i}^\dagger + \chi_{1, 1}^\dagger$', bbox = dict(facecolor = 'red', alpha = 0.5), fontsize = 30)
ax1[0].text(79, 740, r'$T_\mathrm{min} = ' + str(min(T)) + '\ K$', bbox = dict(facecolor = 'red', alpha = 0.5), fontsize = 30)

ax1[1].text(2000, 1.20, 'Iron',    fontsize = 50, color = 'g')
ax1[1].text(2000, 1.23, 'Silicon', fontsize = 50, color = 'r')

ax1[0].xaxis.set_major_locator(MultipleLocator(5))
ax1[0].xaxis.set_minor_locator(AutoMinorLocator(5))

ax1[1].xaxis.set_major_locator(MultipleLocator(250))
ax1[1].xaxis.set_minor_locator(AutoMinorLocator(5))
ax1[1].yaxis.set_minor_locator(AutoMinorLocator(5))

ax1[0].tick_params(labelbottom = 'off')

ax1[0].set_xlim(29, 110)
ax1[0].set_ylim(0, 800)

ax1[1].set_ylim(1, 1.25)
ax1[1].set_xlim(0, 2282.66)

ax1[0].yaxis.grid(True)

ax2 = ax1[0].twinx()

ticks = ax1[0].get_yticks()

ax2.set_yticks(ticks)

ticklabels = []

for i in range(len(ticks)):

    ticklabels.append(str(np.exp(const * ticks[i] / min(T)))[0 : 4])

ax2.set_yticklabels(ticklabels)

ax1[0].tick_params(axis = 'both', labelsize = 19)
ax2.tick_params(axis = 'both', labelsize = 19)

ax1[0].set_ylabel(r'$\chi_i - \chi_i^\dagger, [\mathrm{cm}^{-1}]$', fontsize = 30)
ax2.set_ylabel(r'$\exp{\left(\frac{hc(\chi_i - \chi_i^\dagger)}{kT_\mathrm{min}}\right)}$', fontsize = 30)

ax1[1].set_ylabel(r'$\exp{\left(\frac{hc(\chi_i - \chi_i^\dagger)}{kT}\right)}$', fontsize = 30)
ax1[1].set_xlabel('Height, [km]', fontsize = 30)

auxplt.savepdf('delta_E')
