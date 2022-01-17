import numpy             as np
import matplotlib.pyplot as plt

import os
import importlib

import paths;  importlib.reload(paths)
import pltaux; importlib.reload(pltaux)

plt.close('all')

l, dm, eftd = np.loadtxt(paths.it0h + 'tst_eftd/eftd.out', unpack = True)

h = np.zeros(len(l))

height = np.loadtxt(paths.it0h + 'tst_eftd/ATM_MOD', usecols = [0])

lidx = range(0, len(l))

for ind in lidx: h[ind] = height[int(l[ind]) - 1]

fig, ax1 = plt.subplots(ncols = 1, nrows = 1, figsize = (20.0, 15.0))

pltaux.figpar(fontsize = 25)

fig.tight_layout()

ax1.set_ylim(0, 1)

ax1.scatter(h, -eftd, color = 'k', label = 'Correction')

ax1.set_xlabel('Height, [km]')

ax1.set_ylabel('Correction')

ax2 = ax1.twinx()

ax2.set_ylim(1e-10, 1e+1)

ax2.set_yscale('log')

ax2.set_ylabel('Electron concentration')

pop_en = np.loadtxt(paths.it0h + 'tst_eftd/NLTE/LEV/ELECTR', usecols = [3], skiprows = 2)
pop_el = np.loadtxt(paths.it0h + 'tst_eftd/NLTE/LEV/ELECTR', usecols = [2], skiprows = 2)
pop_hn = np.loadtxt(paths.it0h + 'tst_eftd/NLTE/LEV/HMINUS', usecols = [3], skiprows = 2)
pop_hl = np.loadtxt(paths.it0h + 'tst_eftd/NLTE/LEV/HMINUS', usecols = [2], skiprows = 2)

ax2.plot(height, pop_en, color = 'r',   label = 'Electrons (NLTE)')
ax2.plot(height, pop_el,         'r--', label = 'Electrons (LTE)')
ax2.plot(height, pop_hn, color = 'b',   label = r'$H^-$ (NLTE)')
ax2.plot(height, pop_hl,         'b--', label = r'$H^-$ (NLTE)')

#fig.suptitle(f, y = 1.02)

leg1 = ax1.legend(loc = 2, handlelength = 1, handletextpad=1, prop={'size': 30})
leg2 = ax2.legend(loc = 4, handlelength = 1, handletextpad=1, prop={'size': 30})

pltaux.savepdf('eftd')

#files = os.listdir(paths.it0h + 'tst_eftd/NLTE/LEV')

#for f in files:

#    plt.close('all')

#    print('Plotting ' + f)

#    fig, ax1 = plt.subplots(ncols = 1, nrows = 1, figsize = (20.0, 15.0))

#    pltaux.figpar(fontsize = 25)

#    fig.tight_layout()

#    ax1.set_ylim(0, 1)

#    ax1.scatter(h, -eftd)

#    ax1.set_xlabel('Height, [km]')

#    ax1.set_ylabel('Correction')

#    ax2 = ax1.twinx(); ax2.set_yscale('log'); ax2.set_ylim(1e-30, 1e+1)

#    pop = np.loadtxt(paths.it0h + 'tst_eftd/NLTE/LEV/' + f, usecols = [3], skiprows = 2)

#    ax2.plot(height, pop)

#    fig.suptitle(f, y = 1.03)

#    pltaux.savepdf('etfd/' + f)
