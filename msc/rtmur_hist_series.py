import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import os

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)

from tqdm import tqdm

#x, fhi_fal, y = nessy.read_tau(paths.it1f + 'runtime/def/',     wvl1 = 1005, wvl2 = 3000)
#x, fhi_kur, y = nessy.read_tau(paths.it1f + 'runtime/def_kur/', wvl1 = 1005, wvl2 = 3000)

outname = 'rtmur_def_frmh_hr_fal_kur'

#np.savez(paths.npz + outname, f1 = fhi_fal, f2 = fhi_kur)

frmh_hs = np.load(paths.npz + outname + '.npz')

fhi_fal = frmh_hs['f1']
fhi_kur = frmh_hs['f2']

N =  41
NL = 80

#fhi = np.zeros((N, 1802000))
#fhi = np.zeros((N, 244000))
fhi = np.zeros((N, 402000))

#for i in range(N):

#    wvlh_hr, fhi[i, :], x = nessy.read_tau(paths.it1f + 'rtmur/def/' + str(i) + '/', wvl1 = 1005, wvl2 = 3000)

#    wvlh_hr /= 10.0

#    outname = 'rtmur_def_frmh_hr_' + str(i)

#    np.savez(paths.npz + outname, wv = wvlh_hr, fh = fhi[i, :])

b = np.arange(100, 300, 1)
c = np.arange(100, 300, 1) + 0.5

n =     np.zeros((N, len(c), NL))
n_kur = np.zeros((len(c), NL))

bbox = dict(boxstyle = 'round', ec = (1.0, 0.5, 0.5), fc = (1.0, 0.8, 0.8),)

for i in range(N):

    outname = 'rtmur_def_frmh_hr_' + str(i)

    frmh_hs = np.load(paths.npz + outname + '.npz')

    wvlh_hr =   frmh_hs['wv']
    fhi[i, :] = frmh_hs['fh']

os.system('rm ' + paths.figdir + 'plt_opac/*.pdf')

pdfs = ''

pph = np.zeros(NL)

for k in tqdm(range(NL)):

    pdfs += str(k + 2) + '.pdf '

    plt.close('all')

    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (6.0, 3.0))

    fig.tight_layout()

    for i in range(N):

        wvl = wvlh_hr[np.where(fhi[i, :] == k + 1)]

        for j in range(len(n[i, :, k])):

            n[i, j, k] = len(np.where((wvl >= c[j] - 0.5) & (wvl <= c[j] + 0.5))[0])

        if i == 0: ax.step(b, n[i, :, k], where = 'post', color = 'gray', label = 'MURaM')
        if i != 0: ax.step(b, n[i, :, k], where = 'post', color = 'gray')

        wvl_kur = wvlh_hr[np.where(fhi_kur == k + 1)]

        for j in range(len(n_kur[:, k])):

            n_kur[j, k] = len(np.where((wvl_kur >= c[j] - 0.5) & (wvl_kur <= c[j] + 0.5))[0])

    ax.step(b, np.mean(n[:, :, k], axis = 0), where = 'post', color = 'k', label = 'MURaM average')
    ax.step(b, n_kur[:, k],                   where = 'post', color = 'r', label = 'Kurucz')

    pph[k] = np.sum(np.mean(n[:, :, k], axis = 0)) * 100 / (len(c) * 2000)

#    ax.text(240, max(50, np.amax(n[:, :, k])) - 50,
    ax.text(240, 300, str(pph[k])[0 : 4] + '%', color = 'k', fontsize = 10, bbox = bbox)

#    ax.text(240, 100, str(np.sum(np.mean(n[:, :, k], axis = 0)) * 100 / (len(c) * 2000))[0 : 4] + '%',
#            color = 'k', fontsize = 10, bbox = bbox)

#    ax.set_xlim(180, 300)
    ax.set_xlim(100, 300)
#    ax.set_ylim(0,   max(50, np.amax(n[:, :, k])) + 50)
#    ax.set_ylim(0,   1100)
#    ax.set_ylim(0,   400)
    ax.set_ylim(0,   2100)

    ax.yaxis.set_major_locator(MultipleLocator(200))
#    ax.yaxis.set_major_locator(MultipleLocator(50))
#    ax.yaxis.set_major_locator(MultipleLocator(100))

    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
#    ax.xaxis.set_major_locator(MultipleLocator(10))
    ax.xaxis.set_major_locator(MultipleLocator(25))
#    ax.xaxis.set_minor_locator(AutoMinorLocator(4))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    ax.set_ylabel('Layer ' + str(k + 2))

    ax.set_xlabel('Wavelength, nm')

    leg = ax.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size':7.5})

    for obj in leg.legendHandles: obj.set_linewidth(3.0)

    auxplt.savepdf(str(k + 2), paths.figdir + 'plt_opac/')

os.chdir(paths.figdir + 'plt_opac/')

os.system('pdftk ' + pdfs + ' output rtmur_hist_series.pdf')

os.chdir(paths.mscdir)

print('sum 2 - 2 is: ', np.sum(pph[0 : 1]))
print('sum 2 - 3 is: ', np.sum(pph[0 : 2]))
print('sum 2 - 4 is: ', np.sum(pph[0 : 3]))
print('sum 2 - 5 is: ', np.sum(pph[0 : 4]))
print('sum 2 - 6 is: ', np.sum(pph[0 : 5]))
print('sum 2 - 7 is: ', np.sum(pph[0 : 6]))
print('sum 2 - 8 is: ', np.sum(pph[0 : 7]))
print('sum 2 - 9 is: ', np.sum(pph[0 : 8]))
print('sum 2 - 10 is: ', np.sum(pph[0 : 9]))
print('sum 2 - 11 is: ', np.sum(pph[0 : 10]))
print('sum total is: ',  np.sum(pph))
