import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pylab             as plb
import numpy             as np
import paths             as p
import os

from matplotlib.ticker import AutoMinorLocator

import pltaux

z = np.loadtxt('/mnt/SSD/sim/nessy/inp/atm/FAL99_C', skiprows = 0, usecols = [0])
T = np.loadtxt('/mnt/SSD/sim/nessy/inp/atm/FAL99_C', skiprows = 0, usecols = [1])
N = np.loadtxt('/mnt/SSD/sim/nessy/inp/atm/FAL99_C', skiprows = 0, usecols = [3])

z = z / 1000.0

os.system('rm -fv ' + p.figdir + '*.pdf')

plt.rc('font', size = 20)
plt.rc('font', family = 'serif')

pltaux.figpar(xtick_maj_pad = 10, ytick_maj_pad = 5)

fig, ax1 = plt.subplots()

fig.tight_layout()

ax1.plot(z, T, 'g', linewidth = 2.0, label = 'Temperature')

plt.xlabel('Height, [Mm]')
plt.ylabel(r'$\mathrm{Temperature}, [\mathrm{K}]$')

plt.ylim(4000, 10000)

ax2 = ax1.twinx()
ax2.plot(z, N, 'b', linewidth = 2.0, label = 'Number density')

plt.ylabel(r'$\mathrm{Number density}, \left[\mathrm{cm}^{-3}\right]$')

plt.yscale('log')

plt.minorticks_off()

plt.xlim(min(z), max(z))

fig.suptitle(r'$T(z)$ and $n(z)$ are known', size = 35, y = 1.02)

leg1 = ax1.legend(framealpha = 0, loc = 1, handlelength = 0, handletextpad=0, prop={'size': 20})
leg2 = ax2.legend(framealpha = 0, loc = 2, handlelength = 0, handletextpad=0, prop={'size': 20})

for handle in leg1.legendHandles: handle.set_visible(False)
for handle in leg2.legendHandles: handle.set_visible(False)

text_t = leg1.get_texts()
text_n = leg2.get_texts()

plb.setp(text_t, color = 'g')
plb.setp(text_n, color = 'b')

ax1.tick_params(axis='x', length = 9,   which = 'major')
ax1.tick_params(axis='y', length = 9,   which = 'major')
ax2.tick_params(axis='y', length = 9,   which = 'major')

ax1.tick_params(axis='x', length = 4.5, which = 'minor')
ax1.tick_params(axis='y', length = 4.5, which = 'minor')
ax2.tick_params(axis='y', length = 4.5, which = 'minor')

ax1.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))

pltaux.savepdf(p.figdir, 'falc')
