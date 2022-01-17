import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import os
import glob

import paths;  importlib.reload(paths)
import pltaux; importlib.reload(pltaux)
import sysaux; importlib.reload(sysaux)

def get_elems(filename):

    f = open(filename, 'r')

    elems = []; anums = []

    for line in f:

        if line[0 : 7] == 'ELEMENT':

            parts = line.split()

            elem = parts[2].strip(' () ')

            if parts[3] == ')': anum = int(float(parts[4]))

            if parts[3] != ')': anum = int(float(parts[3]))

            if elem == 'HE':
 
                elem = 'He'

                anum = 2

            elems.append(elem)

            anums.append(anum)

    f.close()

    return elems, anums

def get_time(directory, lin, lis):

    alltstamps = np.loadtxt(directory + '/CONV/ALL', dtype = 'str', skiprows = 2, usecols = [1])

    time = []

    i = 0

    for stamp in alltstamps:

        i = i + 1

        if i % lis != 0: continue

        time_str = stamp.strip("b'")

        time_parts = time_str.split(':')

        hours =   time_parts[0]
                         
        minutes = time_parts[1]

        if hours[0] == '0': hours = hours[1]

        time.append(hours + ':' + minutes)

    return time

if len(sys.argv) == 2: sysaux.abort('List of directory names has to be provided. Abort.')

if len(sys.argv) == 1: sysaux.abort('Mode has to be provided. Abort.')

mode = sys.argv[1]

if mode != 's' and mode != 'p' and mode != 'pj': sysaux.abort('Mode is not understood. Abort.')

#nltedir = 'nltemark_rd/'
nltedir = 'nltemark_noder/'

os.system('mkdir -p ' + paths.figdir + nltedir)

sysaux.clean_dir(paths.figdir + nltedir, mode = 'noverbose')

inpstr = sys.argv[2]

flag = 0

for char in inpstr:

    if char == '*':

        flag = 1

        break

if flag == 0:

    path = paths.it0h + nltedir

    dirs = ''

    for name in inpstr.split(): dirs = dirs + path + name + ' '

    dirs = dirs.split()

if flag == 1: dirs = glob.glob(paths.it0h + nltedir + inpstr)

names = ''

for d in dirs:

    pe = d.split('/')

    name = pe[len(pe) - 1]

    names = names + name + ' '

    conv = np.loadtxt(d + '/CONV/ELEM', skiprows = 2)

    lin = conv[:, 0]

    time = get_time(d, lin, 25)

    elems,      anums =      get_elems(d + '/DATOM_FULL')

    nlte_elems, nlte_anums = get_elems(d + '/DATOM_NLTE')

    plt.close('all')

    fig1, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize = (12.0, 10.0))

    fig1.suptitle(name, y = 1.02)

    if mode == 'p' or mode == 'pj': fig1.tight_layout()

    k = 0

    for i in anums:

        restart = False

        for j in nlte_anums:

            if (i == j):

                restart = True

                break

        if restart: continue

        if k == 0: ax1.plot(lin, conv[:, i], color = 'k', label = 'LTE')

        if k != 0: ax1.plot(lin, conv[:, i], color = 'k')

        k = k + 1

    for j in nlte_anums: ax1.plot(lin, conv[:, j], label = str(j) + ' (' + elems[j - 1] + ')')

    ax1.plot(lin, conv[:, 31], color = 'r', label = 'Electrons')

    ax1.plot(np.arange(1, 321), 1e-4 * np.ones(320), '--', color = 'k', linewidth = 0.5)

    ax1.set_xlabel(r'$\Lambda$-iteration',  labelpad = 10.5)
    ax1.set_ylabel('Convergence parameter', labelpad = 10.5)

    ax1.set_ylim(1e-8, 1e+5)
    ax1.set_xlim(0,    400)

    #ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))

    ax1.set_yscale('log')

    ax1.xaxis.grid(True)

    ax1.xaxis.set_major_locator(MultipleLocator(25))

#    ax2 = ax1.twiny()

#    ax2.set_xlim(ax1.get_xlim())

#    ax2.set_xticks(lin[np.where(lin % 25 == 0)])

#    ax2.set_xticklabels(time)

#    ax2.set_xlabel('time, [h]', labelpad = 10.5)

    leg = ax1.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 10.5})

    for obj in leg.legendHandles: obj.set_linewidth(3.0)

    if mode == 's': fig1.show()

#    if mode == 'p' or mode == 'pj': pltaux.savepdf('nltered_rd/' + name)
    if mode == 'p' or mode == 'pj': pltaux.savepdf(nltedir + name)

if mode == 'pj':

    pdf_names = ''

    for name in names.split(): pdf_names = pdf_names + name + '.pdf '

#    os.chdir(paths.figdir + 'nltered_rd/')
    os.chdir(paths.figdir + nltedir)

    os.system('pdftk ' + pdf_names + 'output ' + 'joined.pdf')

    os.chdir(paths.pydir)
