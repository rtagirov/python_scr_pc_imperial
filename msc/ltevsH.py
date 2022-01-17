import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import LogLocator

import importlib
import sys
import os

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)
import phys;   importlib.reload(phys)

import math as m

rnames = 'H\
          H_C\
          H_C_Mg\
          H_C_Mg_Al\
          H_C_Mg_Al_Si\
          H_C_Mg_Al_Si_Fe\
          H_Mg\
          H_Al\
          H_Si\
          H_Fe'

wn,  nes_lte_hr = nessy.read_spec(paths.it0f + 'ltevsH/lte', wvl1 = 1005, wvl2 = 10000)
wns, nes_lte    = spec.mean(wn / 10.0, nes_lte_hr, 1.0)
np.savez(paths.npz + 'ltevsH_lte', w = wns, lte = nes_lte)

lte = np.load(paths.npz + 'ltevsH_lte.npz'); nes_lte = lte['lte']; wn = lte['w']; 

wa =      np.loadtxt(paths.atlruns + 'ltevsH/lte/spec.out', skiprows = 2, usecols = [0])
atl_lte = np.loadtxt(paths.atlruns + 'ltevsH/lte/spec.out', skiprows = 2, usecols = [1])

#atl_lte = atl_lte * phys.c / (wa * 1.0e-8)**2.0 * 1.0e-7 * (phys.r_sun / phys.au)**2.0 * 1.0e-3 * m.pi
atl_lte = atl_lte * phys.c / (wa * 1.0e-7)**2.0 * 1.0e-7 * (phys.r_sun / phys.au)**2.0 * 1.0e-3 * m.pi

for name in rnames.split():

    w, nes_hnl_hr = nessy.read_spec(paths.it1f + 'ltevsH/' + name, wvl1 = 1005, wvl2 = 10000)
    ws, nes_hnl = spec.mean(w / 10.0, nes_hnl_hr, 1.0)
    np.savez(paths.npz + 'ltevsH_' + name, hnl = nes_hnl)

    ltevsH = np.load(paths.npz + 'ltevsH_' + name + '.npz')

    nes_hnl = ltevsH['hnl']

    atl_hnl = np.loadtxt(paths.atlruns + 'ltevsH/' + name + '/spec.out', skiprows = 2, usecols = [1])
    atl_hnl = atl_hnl * phys.c / (wa * 1.0e-8)**2.0 * 1.0e-7 * (phys.r_sun / phys.au)**2.0 * 1.0e-3 * m.pi

    plt.close('all')

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig.tight_layout()

    fig.suptitle(name.replace('_', ', '), y = 1.02)

    ax.plot(wn, nes_hnl / nes_lte, color = 'k', label = 'NESSY')
    ax.plot(wa, atl_hnl / atl_lte, color = 'r', label = 'ATLAS')

    ax.plot(wn, np.ones(len(wn)), 'k--')

    ax.set_xlim(110.5, 400)

    if name == 'H_Si':            ax.set_yscale('log')
    if name == 'H_C_Mg_Al_Si':    ax.set_yscale('log')
    if name == 'H_C_Mg_Al_Si_Fe': ax.set_yscale('log')

    if name == 'H'            or name == 'H_C':                    ax.set_ylim(0.98, 1.08)
    if name == 'H_C_Mg'       or name == 'H_C_Mg_Al':              ax.set_ylim(0.95, 4.50)
    if name == 'H_C_Mg_Al_Si' or name == 'H_C_Mg_Al_Si_Fe':        ax.set_ylim(0.10, 6.00)

#    if name == 'H_Mg'         or name == 'H_Al' or name == 'H_Fe': ax.set_ylim(0.95, 3.50)
    if name == 'H_Mg'         or name == 'H_Al': ax.set_ylim(0.95, 3.50)

#    if name == 'H_Fe':

#        for i in range(len(nes_hnl) - 1, 0, -1):

#            print(wn[i], nes_hnl[i])

#        sys.exit()

    ax.set_xlabel('Wavelength, [nm]', fontsize = 12.5)
    ax.set_ylabel('NLTE / LTE',       fontsize = 12.5)

    leg = ax.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 20.5})

    for obj in leg.legendHandles: obj.set_linewidth(3.0)

    auxplt.savepdf('ltevsH/' + name)

os.chdir(paths.figdir + 'ltevsH')

names = rnames.split()

pdfnames = ''

for name in rnames.split():

    pdfname = name + '.pdf'

    pdfnames = pdfnames + pdfname + ' '

os.system('pdftk ' + pdfnames + ' output ' + 'spec.pdf')

os.chdir(paths.pydir + 'msc')
