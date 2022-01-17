import numpy             as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from scipy             import interpolate

import importlib
import sys
import glob
import auxplt
import os

from tqdm import tqdm

#wvl_o = np.arange(1005, 10000, 10)
wvl_o = np.arange(1000, 10000, 10)

opa_o = np.zeros((len(wvl_o), 82))

for i, wvl in enumerate(wvl_o):

    opa_o[i, :] = np.loadtxt('./lbkg/' + str(wvl + 5) + '.lbkg', skiprows = 1)

lam_n = np.loadtxt('linop.out', usecols = [0])
kap_n = np.loadtxt('linop.out', usecols = [4])

wvl_n = np.zeros(924)
opa_n = kap_n.reshape((924, 82))

j = 0

for i, elem in enumerate(lam_n):

    if i % 82 == 0: 

        wvl_n[j] = elem

        j += 1

pdfs = ''

os.system('rm ./plt/*.pdf')

ranges = [[x, x + 100] for x in range(100, 1000, 100)]

for pair in ranges:

    for i in tqdm(range(54, 82)):

#        opacn = np.interp(wvl_o, wvl_n, opa_n[:, i])

        pdf_name = 'range_' + str(pair[0]) + '_' + str(pair[1]) + '_depth_' + str(i + 1)

        pdfs += pdf_name + '.pdf '

        plt.close('all')

        fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (6.0, 6.0))

        ax.plot(wvl_o / 10.0, opa_o[:, i], label = 'FIOSS', color = 'k')
#        ax.plot(wvl_o / 10.0, opacn,       label = 'ATLAS', color = 'r')
        ax.plot(wvl_n / 10.0, opa_n[:, i], label = 'ATLAS', color = 'r')

#        ax.bar(wvl_o / 10.0, opa_o[:, i], width = 1, label = 'FIOSS', color = 'k')
#        ax.bar(wvl_o / 10.0, opacn,       width = 1, label = 'ATLAS', alpha = 0.5, color = 'r')

        ax.set_xlabel('Wavelength, nm')
        ax.set_ylabel(r'Opacity, cm$^{-1}$')

        ax.set_yscale('log')

        ax.xaxis.grid(True)

        ax.set_xlim(pair[0], pair[1])
        ax.set_ylim(1e-14, 1e-1)

        title = '       range_' + str(pair[0]) + '_' + str(pair[1]) + '_depth_' + str(i + 1)

        if i + 1 == 1:  title +=  ', top'
        if i + 1 == 55: title += r', $T_\mathrm{min}$'
        if i + 1 == 82: title +=  ', bottom'

        ax.set_title(title)

        ax.xaxis.set_major_locator(MultipleLocator(10))

        leg = ax.legend(framealpha = 1, loc = 2, handletextpad = 1, prop = {'size': 10.5}, bbox_to_anchor=(0, 1.15))

        for obj in leg.legendHandles: obj.set_linewidth(3.0)

        auxplt.savepdf(pdf_name, './plt/')

os.chdir('./plt')

os.system('pdftk ' + pdfs + ' output overall.pdf')

os.chdir('../')
