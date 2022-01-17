import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import importlib
import auxplt
import paths
import spec

importlib.reload(auxplt)
importlib.reload(paths)
importlib.reload(spec)

def P(A, mu):

    P = 0.0

    for i, a in enumerate(A):

        P += a * mu**i

    return P

def find_bracket(wvl, w):

    i = np.abs(wvl - w).argmin()

    if w >= wvl[i] and w <= wvl[i + 1]:

        return [i, i + 1], np.array([wvl[i], wvl[i + 1]])

    if w >= wvl[i - 1] and w <= wvl[i]:

        return [i - 1, i], np.array([wvl[i - 1], wvl[i]])

def interp_clv(idx, wvl, I, w):

    II = np.zeros(len(I[0, :]))

    for j in range(len(I[0, :])):

        w0 = wvl[0]
        w1 = wvl[1]

        I0 = I[idx[0], j]
        I1 = I[idx[1], j]

        II[j] = I0 + (w - w0) * (I1 - I0) / (w1 - w0)

    return II / II[0]

mu = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05])

w = [329.897, 401.970, 445.125, 519.930, 669.400, 1046.600]

A = np.array([[0.09188, 0.92459,  0.19604, -0.39546,  0.23599, -0.05303],  # 329.897
              [0.12323, 1.08648, -0.43974,  0.45912, -0.32759,  0.09850],  # 401.970
              [0.15248, 1.38517, -1.49615,  1.99886, -1.48155,  0.44119],  # 445.125
              [0.23695, 1.29927, -1.28034,  1.37760, -0.85054,  0.21706],  # 519.930
              [0.34685, 1.37539, -2.04425,  2.70493, -1.94290,  0.55999],  # 669.400
              [0.49870, 1.21429, -2.06976,  2.80703, -2.05247,  0.60221]]) # 1046.600

clv_at = np.loadtxt(paths.atlruns + 'var_m/Q/spec.out')

#clv_n0 = np.loadtxt(paths.it0f    + 'var_od/Q/kur/CLV_UVI')
#clv_n1 = np.loadtxt(paths.it1f    + 'var_od/Q/kur/CLV_UVI')
#clv_n2 = np.loadtxt(paths.it1f    + 'var_od/Q/fal/CLV_UVI')

clv_n0 = np.loadtxt(paths.it0f    + 'var_od/Q/kur_cont/CLV_UVI')
clv_n1 = np.loadtxt(paths.it1f    + 'var_od/Q/kur_cont/CLV_UVI')
clv_n2 = np.loadtxt(paths.it1f    + 'var_od/Q/fal_cont/CLV_UVI')
clv_n3 = np.loadtxt(paths.it1f    + 'var_od/Q/fal/CLV_UVI')

np.savez(paths.npz + 'clv_var', clv_at = clv_at,
                                clv_n0 = clv_n0,
                                clv_n1 = clv_n1,
                                clv_n2 = clv_n2,
                                clv_n3 = clv_n3)

clv_var = np.load(paths.npz + 'clv_var.npz')

clv_at = clv_var['clv_at']
clv_n0 = clv_var['clv_n0']
clv_n1 = clv_var['clv_n1']
clv_n2 = clv_var['clv_n2']
clv_n3 = clv_var['clv_n3']

wvl_at = clv_at[:, 1]
int_at = clv_at[:, 4 : ]

wvl_n0 = clv_n0[:, 0] / 10.0
wvl_n1 = clv_n1[:, 0] / 10.0
wvl_n2 = clv_n2[:, 0] / 10.0
wvl_n3 = clv_n3[:, 0] / 10.0

int_n0 = clv_n0[:, 1 : ]
int_n1 = clv_n1[:, 1 : ]
int_n2 = clv_n2[:, 1 : ]
int_n3 = clv_n3[:, 1 : ]

plt.close('all')

fig, ax = plt.subplots(nrows = len(w), ncols = 1, figsize = (6.0, 6.75))

auxplt.figpar(3, 3, 15)

fig.tight_layout()

plt.subplots_adjust(hspace = 0.15)

bbox = dict(boxstyle = 'round', ec = (1.0, 0.5, 0.5), fc = (1.0, 0.8, 0.8),)

atl_pairs = []
nes_pairs = []

for i in range(len(w)):

    neckel = P(A[i, :], mu)

    idx_atl, wvl_atl = find_bracket(wvl_at, w[i])

    idx_ne0, wvl_ne0 = find_bracket(wvl_n0, w[i])
    idx_ne1, wvl_ne1 = find_bracket(wvl_n1, w[i])
    idx_ne2, wvl_ne2 = find_bracket(wvl_n2, w[i])

    print('atl = ', wvl_atl, 'ne0 = ', wvl_ne0)

    atl_pairs.append(wvl_atl)
    nes_pairs.append(wvl_ne0)

    clv_atl = interp_clv(idx_atl, wvl_atl, int_at, w[i])

    clv_ne0 = interp_clv(idx_ne0, wvl_ne0, int_n0, w[i])
    clv_ne1 = interp_clv(idx_ne1, wvl_ne1, int_n1, w[i])
    clv_ne2 = interp_clv(idx_ne2, wvl_ne2, int_n2, w[i])

    ax[i].plot(mu, (clv_atl - neckel) * 100 / neckel, color = 'k', linewidth = 1.8, label = 'ATLAS9, LTE, U99')
    ax[i].plot(mu, (clv_ne0 - neckel) * 100 / neckel, color = 'm', linewidth = 0.9, label = 'NESSY, LTE, U99')
    ax[i].plot(mu, (clv_ne1 - neckel) * 100 / neckel, color = 'g', linewidth = 0.9, label = 'NESSY, NLTE, U99')
    ax[i].plot(mu, (clv_ne2 - neckel) * 100 / neckel, color = 'r', linewidth = 0.9, label = 'NESSY, NLTE, FAL99')

    ax[i].set_xlim(1.0, 0.05)

    if i == 3: ax[i].axes.set_ylabel(r'$(\mathrm{CLV} - \mathrm{CLV}_\mathrm{Neckel}) / \mathrm{CLV}_\mathrm{Neckel}$, [\%]')

    if i != 5: ax[i].axhline(y = 0.0, color = 'k', linestyle = '--', linewidth = 0.7)

    ax[i].xaxis.set_major_locator(MultipleLocator(0.1))

    ax[i].yaxis.set_minor_locator(AutoMinorLocator(5))
    ax[i].xaxis.set_minor_locator(AutoMinorLocator(4))

    xt = 0.12

    if i in [0, 1, 2]: yt = 0.73
    if i in [3, 4, 5]: yt = 0.27

    ax[i].text(xt, yt, str(round(w[i])) + ' nm', ha = 'center', va = 'center', transform = ax[i].transAxes, bbox = bbox)

    if i != len(w) - 1:

        ax[i].xaxis.set_tick_params(direction = 'in', which = 'both')
        ax[i].tick_params(labelbottom = 'off')

ax[len(w) - 1].set_ylim(top = 0.0)

ax[len(w) - 1].set_xlabel(r'$\mu$')

ax[3].yaxis.set_label_coords(-0.10, 1.1)

shift = 0.030

ax[0].text(1.0 - shift,   15.5, 'ATLAS9, LTE, U99', fontsize = 10)
ax[0].text(0.775 - shift, 15.5, 'NESSY, LTE, U99', color = 'm', fontsize = 10)
ax[0].text(0.565 - shift, 15.5, 'NESSY, NLTE, U99', color = 'g', fontsize = 10)
ax[0].text(0.335 - shift, 15.5, 'NESSY, NLTE, FAL99', color = 'r', fontsize = 10)

auxplt.savepdf('var/neckel_clv_comp_7_interp_od_cont')

int_fd_at = clv_at[:, 3]

int_fd_ne_line = np.zeros(len(int_n3[:, 0]))
int_fd_ne_cont = np.zeros(len(int_n2[:, 0]))

ip = np.sqrt(1 - mu**2)

area = np.zeros(len(ip) - 1)
iprc = np.zeros(len(ip) - 1)

for j in range(len(area)):

    area[j] = ip[j + 1]**2 - ip[j]**2

    iprc[j] = (ip[j + 1] + ip[j]) / 2

for i, elem in enumerate(wvl_n2):

    int_ne_rc = np.interp(iprc, ip, int_n2[i, :])

    int_fd_ne_cont[i] = sum(area * int_ne_rc)

for i, elem in enumerate(wvl_n3):

    int_ne_rc = np.interp(iprc, ip, int_n3[i, :])

    int_fd_ne_line[i] = sum(area * int_ne_rc)

int_fd_ne_line_red = spec.mean_over_grid(int_fd_ne_line, wvl_n3, wvl_at)

plt.close('all')

#fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6.0, 6.75))
fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize = (6.0, 6.75))

ax.plot(wvl_n3, int_fd_ne_line * 1e+5,     color = 'gray')
ax.plot(wvl_n2, int_fd_ne_cont * 1e+5,     color = 'black')
ax.plot(wvl_at, int_fd_at * 1e+5,          color = 'yellow')
ax.plot(wvl_at, int_fd_ne_line_red * 1e+5, color = 'purple')

for i, elem in enumerate(atl_pairs):

    ax.axvspan(atl_pairs[i][0], atl_pairs[i][1], alpha = 0.5, color = 'red')
    ax.axvspan(nes_pairs[i][0], nes_pairs[i][1], alpha = 1.0, color = 'red')

#for i in range(len(ax)):

#    ax[i].set_xlim(300, 1100)
ax.set_xlim(400, 1100)
ax.set_ylim(0.0, 3.7)

#    if i != 0: ax[i].set_xlabel('Wavelength, [nm]')
#    ax[i].set_ylabel(r'Flux $\times 10^5$, [W / m$^2$ / nm]')

ax.set_xlabel('Wavelength, [nm]')
ax.set_ylabel(r'Flux $\times 10^5$, [W / m$^2$ / nm]')

plt.show()
