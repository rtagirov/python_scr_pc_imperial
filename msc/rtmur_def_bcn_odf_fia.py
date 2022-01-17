import numpy as np
import matplotlib.pyplot as plt

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import importlib
import sys
import math

if not '../aux/' in sys.path: sys.path.append('../aux/')

import paths;  importlib.reload(paths)
import spec;   importlib.reload(spec)
import nessy;  importlib.reload(nessy)
import auxsys; importlib.reload(auxsys)
import auxplt; importlib.reload(auxplt)

l1 = 1005
l2 = 3000

Nw = int(l2 / 10.0) - math.floor(l1 / 10.0) + 1

#wvl, faldef = nessy.read_spec(paths.it1f + 'runtime/def',                             wvl1 = l1, wvl2 = l2)
#wvl, kurdef = nessy.read_spec(paths.it1f + 'runtime/def_kur',                         wvl1 = l1, wvl2 = l2)
#wvl, falbcn = nessy.read_spec(paths.it1f + 'runtime/nlte11_fgrid_red_5_20_eps_2',     wvl1 = l1, wvl2 = l2)
#wvl, kurbcn = nessy.read_spec(paths.it1f + 'runtime/nlte11_fgrid_red_5_20_eps_2_kur', wvl1 = l1, wvl2 = l2)
#wvl, falodf = nessy.read_spec(paths.it1f + 'runtime/def_atlodf_nesab',                wvl1 = l1, wvl2 = l2)
#wvl, kurodf = nessy.read_spec(paths.it1f + 'runtime/def_kur_atlodf_nesab',            wvl1 = l1, wvl2 = l2)
#wvl, falfia = nessy.read_spec(paths.it1f + 'runtime/def_fioss_acc',                   wvl1 = l1, wvl2 = l2)
#wvl, kurfia = nessy.read_spec(paths.it1f + 'runtime/def_kur_fioss_acc',               wvl1 = l1, wvl2 = l2)
#wvl, falall = nessy.read_spec(paths.it1f + 'runtime/fal_all_mod',                     wvl1 = l1, wvl2 = l2)
#wvl, kurall = nessy.read_spec(paths.it1f + 'runtime/kur_all_mod',                     wvl1 = l1, wvl2 = l2)

#wvl = wvl / 10.0

#wvls, faldef_s = spec.mean_within_delta(wvl, faldef, 1.0)
#wvls, kurdef_s = spec.mean_within_delta(wvl, kurdef, 1.0)
#wvls, falbcn_s = spec.mean_within_delta(wvl, falbcn, 1.0)
#wvls, kurbcn_s = spec.mean_within_delta(wvl, kurbcn, 1.0)
#wvls, falodf_s = spec.mean_within_delta(wvl, falodf, 1.0)
#wvls, kurodf_s = spec.mean_within_delta(wvl, kurodf, 1.0)
#wvls, falfia_s = spec.mean_within_delta(wvl, falfia, 1.0)
#wvls, kurfia_s = spec.mean_within_delta(wvl, kurfia, 1.0)
#wvls, falall_s = spec.mean_within_delta(wvl, falall, 1.0)
#wvls, kurall_s = spec.mean_within_delta(wvl, kurall, 1.0)

outname = 'rtmur_def_bcn_odf_fia'

#np.savez(paths.npz + outname, w = wvls, f1 =  faldef_s,
#                                        f2 =  kurdef_s,
#                                        f3 =  falbcn_s,
#                                        f4 =  kurbcn_s,
#                                        f5 =  falodf_s,
#                                        f6 =  kurodf_s,
#                                        f7 =  falfia_s,
#                                        f8 =  kurfia_s,
#                                        f9 =  falall_s,
#                                        f10 = kurall_s)

specs = np.load(paths.npz + outname + '.npz')

faldef = specs['f1']
kurdef = specs['f2']
falbcn = specs['f3']
kurbcn = specs['f4']
falodf = specs['f5']
kurodf = specs['f6']
falfia = specs['f7']
kurfia = specs['f8']
falall = specs['f9']
kurall = specs['f10']

N = 41
#N = 15
#N = 30

fd = np.zeros((N, Nw))
fb = np.zeros((N, Nw))
fo = np.zeros((N, Nw))
fa = np.zeros((N, Nw))
ft = np.zeros((N, Nw))

dev_b = np.zeros((N, Nw))
dev_o = np.zeros((N, Nw))
dev_a = np.zeros((N, Nw))
dev_t = np.zeros((N, Nw))

for i in range(N):

#    wvl, flu_def = nessy.read_spec(paths.it1f + 'rtmur/def/' + str(i), wvl1 = l1, wvl2 = l2)
#    wvl, flu_bcn = nessy.read_spec(paths.it1f + 'rtmur/bcn/' + str(i), wvl1 = l1, wvl2 = l2)
#    wvl, flu_odf = nessy.read_spec(paths.it1f + 'rtmur/odf/' + str(i), wvl1 = l1, wvl2 = l2)
#    wvl, flu_fia = nessy.read_spec(paths.it1f + 'rtmur/fia/' + str(i), wvl1 = l1, wvl2 = l2)
#    wvl, flu_all = nessy.read_spec(paths.it1f + 'rtmur/all/' + str(i), wvl1 = l1, wvl2 = l2)

#    wvl = wvl / 10.0

#    wvls, fd[i, :] = spec.mean_within_delta(wvl, flu_def, 1.0)
#    wvls, fb[i, :] = spec.mean_within_delta(wvl, flu_bcn, 1.0)
#    wvls, fo[i, :] = spec.mean_within_delta(wvl, flu_odf, 1.0)
#    wvls, fa[i, :] = spec.mean_within_delta(wvl, flu_fia, 1.0)
#    wvls, ft[i, :] = spec.mean_within_delta(wvl, flu_all, 1.0)

    outname = 'rtmur_def_bcn_odf_fia_' + str(i)

#    np.savez(paths.npz + outname, wv = wvls, fd = fd[i, :],
#                                             fb = fb[i, :],
#                                             fo = fo[i, :],
#                                             fa = fa[i, :],
#                                             ft = ft[i, :])

    specs = np.load(paths.npz + outname + '.npz')

    wv =       specs['wv']
    fd[i, :] = specs['fd']
    fb[i, :] = specs['fb']
    fo[i, :] = specs['fo']
    fa[i, :] = specs['fa']
    ft[i, :] = specs['ft']

plt.close('all')

fig, ax = plt.subplots(nrows = 5, ncols = 1, figsize = (6.0, 10.0))

fig.tight_layout()

for i in range(N):

    dev_b[i, :] = fd[i, :] / fb[i, :]
    dev_o[i, :] = fd[i, :] / fo[i, :]
    dev_a[i, :] = fd[i, :] / fa[i, :]
    dev_t[i, :] = fd[i, :] / ft[i, :]

    if i == 0: ax[0].plot(wv, dev_b[i, :], color = 'gray', label = 'MURaM')
    if i != 0: ax[0].plot(wv, dev_b[i, :], color = 'gray')

    ax[1].plot(wv, dev_o[i, :], color = 'gray')
    ax[2].plot(wv, dev_a[i, :], color = 'gray')
    ax[3].plot(wv, dev_t[i, :], color = 'gray')

mean_fd = np.mean(fd[:, :], axis = 0)
mean_fb = np.mean(fb[:, :], axis = 0)
mean_fo = np.mean(fo[:, :], axis = 0)
mean_fa = np.mean(fa[:, :], axis = 0)
mean_ft = np.mean(ft[:, :], axis = 0)

mean_dev_b = mean_fd / mean_fb
mean_dev_o = mean_fd / mean_fo
mean_dev_a = mean_fd / mean_fa
mean_dev_t = mean_fd / mean_ft

ax[0].plot(wv, mean_dev_b,      color = 'k', label = 'MURaM average')
ax[0].plot(wv, faldef / falbcn, color = 'b', label = 'FAL99-C')
ax[0].plot(wv, kurdef / kurbcn, color = 'r', label = 'Kurucz')

ax[1].plot(wv, faldef / falodf, color = 'b')
ax[1].plot(wv, kurdef / kurodf, color = 'r')
ax[1].plot(wv, mean_dev_o,      color = 'k')

ax[2].plot(wv, faldef / falfia, color = 'b')
ax[2].plot(wv, kurdef / kurfia, color = 'r')
ax[2].plot(wv, mean_dev_a,      color = 'k')

ax[3].plot(wv, faldef / falall, color = 'b')
ax[3].plot(wv, kurdef / kurall, color = 'r')
ax[3].plot(wv, mean_dev_t,      color = 'k')

ax[4].plot(wv, mean_dev_b,      color = 'k', label = 'mod')
ax[4].plot(wv, mean_dev_o,      color = 'b', label = 'odf')
ax[4].plot(wv, mean_dev_a,      color = 'g', label = 'ssa')
ax[4].plot(wv, mean_dev_t,      color = 'r', label = 'tot')

#for i in [0, 1]:

#    ax[i].axvline(x = 170,  color = 'k', linestyle = '--', linewidth = 1.0)
#    ax[i].axhline(y = 1.00, color = 'k', linestyle = '--', linewidth = 1.0)
#    ax[i].axhline(y = 1.07, color = 'k', linestyle = '--', linewidth = 1.0)
#    ax[i].axhline(y = 0.93, color = 'k', linestyle = '--', linewidth = 1.0)

for i in [0, 1, 2, 3, 4]:

    ax[i].set_xlim(100, 300)
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(5))

ax[0].set_ylabel('Ratio (mod)')
ax[1].set_ylabel('Ratio (odf)')
ax[2].set_ylabel('Ratio (ssa)')
ax[3].set_ylabel('Ratio (tot)')
ax[4].set_ylabel('Ratio')

ax[4].set_xlabel('Wavelength, nm')

ax[0].yaxis.set_major_locator(MultipleLocator(0.025))
ax[1].yaxis.set_major_locator(MultipleLocator(0.1))
ax[3].yaxis.set_major_locator(MultipleLocator(0.1))

leg0 = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size':7.5})
leg4 = ax[4].legend(framealpha = 1, loc = 4, handletextpad = 1, prop = {'size':7.5})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)
for obj in leg4.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('rtmur_def_bcn_odf_fia')
