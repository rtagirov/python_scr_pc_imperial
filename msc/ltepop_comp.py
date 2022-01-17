import matplotlib.pyplot as plt
import pylab             as pl
import numpy             as np
import paths             as p
import os

from matplotlib.ticker import AutoMinorLocator

import pltaux

dir_icl = 'TST_ICL_LTEPOP/'

dir_keg = 'TST_KEG_LTEPOP/'

z = np.loadtxt(p.it0h + dir_icl + 'FAL_VD', skiprows = 0, usecols = [0])

z = z / 1000.0

#os.system('pkill evince')
#os.system('rm -fv ' + p.figdir + 'overall.pdf')
os.system('rm -fv ' + p.figdir + '*.pdf')

for fname in os.listdir(p.it0h + dir_icl + p.lev):

    n_icl = np.loadtxt(p.it0h + dir_icl + p.lev + fname, skiprows = 2, usecols = [2])
    n_keg = np.loadtxt(p.it0h + dir_keg + p.lev + fname, skiprows = 2, usecols = [2])

    plt.rc('font', size = 20)
    plt.rc('font', family = 'serif')

    fig, ax = plt.subplots()

    ax.plot(z, n_icl, 'r', linewidth = 3.0, label = 'imperial')
    ax.plot(z, n_keg, 'k', linewidth = 0.5, label = 'pmod/wrc')

#    leg = ax.legend(framealpha = 0, loc = 4, handlelength = 0, handletextpad=0, prop={'size':20})
    leg = ax.legend(framealpha = 0, loc = 'best', handlelength = 0, handletextpad=0, prop={'size': 20})

    for handle in leg.legendHandles:
        handle.set_visible(False)

    text_icl, text_pmo = leg.get_texts()

    pl.setp(text_icl, color = 'r')
    pl.setp(text_pmo, color = 'k')

    plt.yscale('log')

    plt.xlabel('Height, [Mm]')
#    plt.ylabel('$n / n_\mathrm{tot}$')
    plt.ylabel(r'LTE population, [$\mathrm{cm}^{-3}$]')

    fig.suptitle(fname, fontsize=20)

    plt.grid(True)

    plt.xlim(min(z), max(z))
#    plt.ylim(1.0E-1, 2.0E+2)
    plt.ylim(min(n_icl), max(n_icl))

    minorloc = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minorloc)

    pltaux.savepdf(p.figdir, fname)

os.system('pdftk ' + p.figdir + '*' + ' output ' + p.figdir + 'overall_ltepop.pdf')
#os.system('evince' + ' ' + p.figdir + 'overall.pdf' + ' &')
