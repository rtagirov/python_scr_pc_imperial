#!/usr/bin/env python

import matplotlib.pyplot as plt
import pylab             as pl
import numpy             as np
import paths             as p
import os

from matplotlib.ticker import AutoMinorLocator

import pltaux

dir_icl = 'TST_ICL_FCHHT_B/'

dir_keg = 'TST_KEG_FCHHT_B/'

z = np.loadtxt(p.it0h + dir_icl + 'FAL_VD', skiprows = 0, usecols = [0])

z = z / 1000.0

#os.system('pkill evince')
#os.system('rm -fv ' + p.figdir + 'overall.pdf')
os.system('rm -fv ' + p.figdir + '*.pdf')

for fname in os.listdir(p.it0h + dir_icl + p.lev):

    b_icl = np.loadtxt(p.it0h + dir_icl + p.lev + fname, skiprows = 2, usecols = [4])
    b_keg = np.loadtxt(p.it0h + dir_keg + p.lev + fname, skiprows = 2, usecols = [4])

    plt.rc('font', size = 20)
    plt.rc('font', family = 'serif')

    fig, ax = plt.subplots()

    ax.plot(z, b_icl, 'r', linewidth = 3.0, label = 'imperial')
    ax.plot(z, b_keg, 'k', linewidth = 0.5, label = 'pmod/wrc')

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
    plt.ylabel('NLTE departure coefficient')

    fig.suptitle(fname, fontsize=20)

    plt.grid(True)

    plt.xlim(min(z), max(z))
#    plt.ylim(1.0E-1, 2.0E+2)
    plt.ylim(min(b_icl), max(b_icl))

    minorloc = AutoMinorLocator(5)
    ax.xaxis.set_minor_locator(minorloc)

    pltaux.savepdf(p.figdir, fname)

os.system('pdftk ' + p.figdir + '*' + ' output ' + p.figdir + 'overall_dep.pdf')
#os.system('evince' + ' ' + p.figdir + 'overall.pdf' + ' &')
