import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)

#x, fhi_fal, y = nessy.read_tau(paths.it1f + 'runtime/def/',     wvl1 = 1800, wvl2 = 3000)
#x, fhi_kur, y = nessy.read_tau(paths.it1f + 'runtime/def_kur/', wvl1 = 1800, wvl2 = 3000)

outname = 'rtmur_def_frmh_hr_fal_kur'

#np.savez(paths.npz + outname, f1 = fhi_fal, f2 = fhi_kur)

frmh_hs = np.load(paths.npz + outname + '.npz')

fhi_fal = frmh_hs['f1']
fhi_kur = frmh_hs['f2']

N = 41

#fhi = np.zeros((N, 1802000))
fhi = np.zeros((N, 244000))

#for i in range(N):

#    wvlh_hr, fhi[i, :], x = nessy.read_tau(paths.it1f + 'rtmur/def/' + str(i) + '/', wvl1 = 1800, wvl2 = 3000)

#    wvlh_hr /= 10.0

#    outname = 'rtmur_def_frmh_hr_' + str(i)

#    np.savez(paths.npz + outname, wv = wvlh_hr, fh = fhi[i, :])

b = np.arange(180, 300, 1)
c = np.arange(180, 300, 1) + 0.5

n = np.zeros((N, len(c), 9))

n_kur = np.zeros((len(c), 9))

plt.close('all')

fig, ax = plt.subplots(nrows = 9, ncols = 1, figsize = (6.0, 13.0))

fig.tight_layout()

bbox = dict(boxstyle = 'round', ec = (1.0, 0.5, 0.5), fc = (1.0, 0.8, 0.8),)

for i in range(N):

    outname = 'rtmur_def_frmh_hr_' + str(i)

    frmh_hs = np.load(paths.npz + outname + '.npz')

    wvlh_hr =   frmh_hs['wv']
    fhi[i, :] = frmh_hs['fh']

for k in range(9):

    for i in range(N):

        wvl = wvlh_hr[np.where(fhi[i, :] == k + 1)]

        for j in range(len(n[i, :, k])):

            n[i, j, k] = len(np.where((wvl >= c[j] - 0.5) & (wvl <= c[j] + 0.5))[0])

        if i == 0 and k == 0:

            ax[k].step(b, n[i, :, k], where = 'post', color = 'gray', label = 'MURaM')

        else:

            ax[k].step(b, n[i, :, k], where = 'post', color = 'gray')

        wvl_kur = wvlh_hr[np.where(fhi_kur == k + 1)]

        for j in range(len(n_kur[:, k])):

            n_kur[j, k] = len(np.where((wvl_kur >= c[j] - 0.5) & (wvl_kur <= c[j] + 0.5))[0])

    if k == 0:

        ax[k].step(b, np.mean(n[:, :, k], axis = 0), where = 'post', color = 'k', label = 'MURaM average')
        ax[k].step(b, n_kur[:, k],                   where = 'post', color = 'r', label = 'Kurucz')

    if k != 0:

        ax[k].step(b, np.mean(n[:, :, k], axis = 0), where = 'post', color = 'k')
        ax[k].step(b, n_kur[:, k],                   where = 'post', color = 'r')

    if k == 0: ax[k].text(240, 250, str(np.sum(np.mean(n[:, :, k], axis = 0)) * 100 / (len(c) * 2000))[0 : 4] + '%',
                                    color = 'k', fontsize = 10, bbox = bbox)

    if k != 0: ax[k].text(240, 100, str(np.sum(np.mean(n[:, :, k], axis = 0)) * 100 / (len(c) * 2000))[0 : 4] + '%',
                                    color = 'k', fontsize = 10, bbox = bbox)

    ax[k].set_xlim(180, 300)

    if k == 0: ax[k].set_ylim(0, 400)
    if k != 0: ax[k].set_ylim(0, 150)
#    if k != 0: ax[k].set_ylim(0, 400)
    if k == 7: ax[k].set_ylim(0, 175)
#    if k == 7: ax[k].set_ylim(0, 400)

    if k == 0: ax[k].yaxis.set_major_locator(MultipleLocator(100))
    if k != 0: ax[k].yaxis.set_major_locator(MultipleLocator(50))

    ax[k].yaxis.set_minor_locator(AutoMinorLocator(4))
    ax[k].xaxis.set_minor_locator(AutoMinorLocator(4))

    ax[k].set_ylabel('Layer ' + str(k + 2))

    if k == 0:

        leg = ax[k].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size':7.5})

        for obj in leg.legendHandles: obj.set_linewidth(3.0)

ax[8].set_xlabel('Wavelength, nm')

auxplt.savepdf('rtmur_hist')
