import numpy as np

import matplotlib.pyplot as plt

import importlib

import paths;     importlib.reload(paths)
import auxfunc;   importlib.reload(auxfunc)
import pltaux;    importlib.reload(pltaux)
import sysaux;    importlib.reload(sysaux)
import oper_file; importlib.reload(oper_file)

import sys

import os

from tqdm import tqdm

#number = sys.argv[1]

#fgr = np.loadtxt(paths.it0h + 'test_grid/FGRID') 

#xjc_1d = np.loadtxt(paths.it0h + 'test_grid/xjc.out', usecols = [3])

wvl = np.loadtxt(paths.it1h + 'fgrid_new/sigma.out', skiprows = 560 * 113, usecols = [2]) / 10.0

sig_1d = np.loadtxt(paths.it1h + 'fgrid_new/sigma.out', usecols = [3])

sig_1d[sig_1d == 0.0] = np.nan

#xjc = np.reshape(xjc_1d, (82, 1948))
sig = np.reshape(sig_1d, (114, 1948))
#sig = np.reshape(sig_1d, (114, 560))

plt.close('all')

#f1 = plt.figure(1)

#for i in range(0, len(fgr)): plt.axvline(fgr[i], color = 'r',                   linewidth = 0.5)
#for i in range(0, len(wvl)): plt.axvline(wvl[i], color = 'k', linestyle = '--', linewidth = 0.1)

#plt.plot(wvl, xjc[81, :])
#plt.plot(wvl, xjc[77, :])
#plt.plot(wvl, xjc[75, :])
#plt.plot(wvl, xjc[70, :])
#plt.plot(wvl, xjc[65, :])
#plt.plot(wvl, xjc[0, :])

#plt.xscale('log')

#f1.show()

H =  []
He = []
Li = []
Be = []
B =  []
C =  []
N =  []
O =  []
F =  []
Ne = []
Na = []
Mg = []
Al = []
Si = []
P =  []
S =  []
Cl = []
Ar = []
K =  []
Ca = []
Sc = []
Ti = []
V =  []
Cr = []
Mn = []
Fe = []
Co = []
Ni = []
Cu = []
Zn = []

#lev = [] * oper_file.num_lines(paths.it0h + 'test_grid/levels.out')
lev = [] * oper_file.num_lines(paths.it1h + 'fgrid/new/levels.out')

i = 0

#for line in open(paths.it0h + 'test_grid/levels.out'):
for line in open(paths.it1h + 'fgrid/new/levels.out'):

    parts = line.split()

    lev.append(parts[1])

    if lev[i][0 : 2] == 'H' or lev[i][0 : 2] == 'HE': lev[i] = parts[1] + parts[2]

    if lev[i][0 : 2] == 'HM': H.append(i)
    if lev[i][0 : 2] == 'HI': H.append(i)
    if lev[i][0 : 2] == 'HE': He.append(i)
    if lev[i][0 : 2] == 'Li': Li.append(i)
    if lev[i][0 : 2] == 'Be': Be.append(i)
    if lev[i][0 : 2] == 'BI': B.append(i)
    if lev[i][0 : 2] == 'CI': C.append(i)
    if lev[i][0 : 2] == 'NI': N.append(i)
    if lev[i][0 : 2] == 'OI': O.append(i)
    if lev[i][0 : 2] == 'FI': F.append(i)
    if lev[i][0 : 2] == 'Ne': Ne.append(i)
    if lev[i][0 : 2] == 'Na': Na.append(i)
    if lev[i][0 : 2] == 'Mg': Mg.append(i)
    if lev[i][0 : 2] == 'Al': Al.append(i)
    if lev[i][0 : 2] == 'Si': Si.append(i)
    if lev[i][0 : 2] == 'PI': P.append(i)
    if lev[i][0 : 2] == 'SI': S.append(i)
    if lev[i][0 : 2] == 'Cl': Cl.append(i)
    if lev[i][0 : 2] == 'Ar': Ar.append(i)
    if lev[i][0 : 2] == 'KI': K.append(i)
    if lev[i][0 : 2] == 'Ca': Ca.append(i)
    if lev[i][0 : 2] == 'Sc': Sc.append(i)
    if lev[i][0 : 2] == 'Ti': Ti.append(i)
    if lev[i][0 : 2] == 'VI': V.append(i)
    if lev[i][0 : 2] == 'Cr': Cr.append(i)
    if lev[i][0 : 2] == 'Mn': Mn.append(i)
    if lev[i][0 : 2] == 'Fe': Fe.append(i)
    if lev[i][0 : 2] == 'Co': Co.append(i)
    if lev[i][0 : 2] == 'Ni': Ni.append(i)
    if lev[i][0 : 2] == 'Cu': Cu.append(i)
    if lev[i][0 : 2] == 'Zn': Zn.append(i)

    i = i + 1

fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

fig.tight_layout()

ax.set_xlim(20, 1e+4)

for m in range(0, len(wvl)): ax.axvline(wvl[m], color = 'k', linewidth = 0.05)

for i in tqdm(range(1, 31), ncols = auxfunc.term_width(), desc = 'Plotting'):

    if i == 1:  elem = H;  title = 'H'
    if i == 2:  elem = He; title = 'He'
    if i == 3:  elem = Li; title = 'Li'
    if i == 4:  elem = Be; title = 'Be'
    if i == 5:  elem = B;  title = 'B'
    if i == 6:  elem = C;  title = 'C'
    if i == 7:  elem = N;  title = 'N'
    if i == 8:  elem = O;  title = 'O'
    if i == 9:  elem = F;  title = 'F'
    if i == 10: elem = Ne; title = 'Ne'
    if i == 11: elem = Na; title = 'Na'
    if i == 12: elem = Mg; title = 'Mg'
    if i == 13: elem = Al; title = 'Al'
    if i == 14: elem = Si; title = 'Si'
    if i == 15: elem = P;  title = 'P'
    if i == 16: elem = S;  title = 'S'
    if i == 17: elem = Cl; title = 'Cl'
    if i == 18: elem = Ar; title = 'Ar'
    if i == 19: elem = K;  title = 'K'
    if i == 20: elem = Ca; title = 'Ca'
    if i == 21: elem = Sc; title = 'Sc'
    if i == 22: elem = Ti; title = 'Ti'
    if i == 23: elem = V;  title = 'V'
    if i == 24: elem = Cr; title = 'Cr'
    if i == 25: elem = Mn; title = 'Mn'
    if i == 26: elem = Fe; title = 'Fe'
    if i == 27: elem = Co; title = 'Co'
    if i == 28: elem = Ni; title = 'Ni'
    if i == 29: elem = Cu; title = 'Cu'
    if i == 30: elem = Zn; title = 'Zn'

    for l in elem:

        if l != elem[len(elem) - 1]: ax.plot(wvl, sig[l, :], linewidth = 3, label = lev[l])

#    if i != 2:

#        leg = ax.legend(framealpha = 1, loc = 'best',        handletextpad = 1, prop = {'size': 20})

#    else:

#        leg = ax.legend(framealpha = 1, loc = 'lower right', handletextpad = 1, prop = {'size': 20})

#    for obj in leg.legendHandles: obj.set_linewidth(4.0)

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('Wavelength, [nm]')
    ax.set_ylabel(r'Cross section, [cm$^{-2}$]')

pltaux.savepdf('cs_fgrid_new_2020')
