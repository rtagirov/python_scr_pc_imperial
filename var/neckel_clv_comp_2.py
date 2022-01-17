import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import importlib
import auxplt
import paths

importlib.reload(auxplt)
importlib.reload(paths)

def P(A, mu):

    P = 0.0

    for i, a in enumerate(A):

        P += a * mu**i

    return P

def wvl_selection(wvl, wvl_neckel):

    idxs = []
    wvls = []

    for i, w in enumerate(wvl_neckel):

        idx = np.abs(wvl - w).argmin()

        idxs.append(idx)
        wvls.append(wvl[idx])

    return idxs, wvls

mu = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05])

w = [303.327, 329.897, 401.970, 445.125, 519.930, 669.400, 1046.600]

A = np.array([[0.08011, 0.70695,  0.49910, -0.31080, -0.02177,  0.04642],  # 303.327
              [0.09188, 0.92459,  0.19604, -0.39546,  0.23599, -0.05303],  # 329.897
              [0.12323, 1.08648, -0.43974,  0.45912, -0.32759,  0.09850],  # 401.970
              [0.15248, 1.38517, -1.49615,  1.99886, -1.48155,  0.44119],  # 445.125
              [0.23695, 1.29927, -1.28034,  1.37760, -0.85054,  0.21706],  # 519.930
              [0.34685, 1.37539, -2.04425,  2.70493, -1.94290,  0.55999],  # 669.400
              [0.49870, 1.21429, -2.06976,  2.80703, -2.05247,  0.60221]]) # 1046.600

clv_at = np.loadtxt(paths.atlruns + 'var_m/Q/spec.out')
clv_n0 = np.loadtxt(paths.it0f + 'var/Q/kur/CLV_UVI')
clv_n1 = np.loadtxt(paths.it1f + 'var/Q/kur/CLV_UVI')
clv_n2 = np.loadtxt(paths.it1f + 'var/Q/fal/CLV_UVI')

wvl_at = clv_at[:, 1]
int_at = clv_at[:, 4 : ]

wvl_n0 = clv_n0[:, 0] / 10.0
wvl_n1 = clv_n1[:, 0] / 10.0
wvl_n2 = clv_n2[:, 0] / 10.0

int_n0 = clv_n0[:, 1 : ]
int_n1 = clv_n1[:, 1 : ]
int_n2 = clv_n2[:, 1 : ]

idx_at, wvl_atlas = wvl_selection(wvl_at, np.array(w))

idx_n0, wvl_nessy = wvl_selection(wvl_n0, np.array(w))
idx_n1, wvl_nessy = wvl_selection(wvl_n1, np.array(w))
idx_n2, wvl_nessy = wvl_selection(wvl_n2, np.array(w))

plt.close('all')

fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (6.0, 6.75))

#col = ['g', 'r', 'c', 'm', 'b', 'y', 'k']

col = ['k', 'y', 'b', 'm', 'c', 'r', 'g']

auxplt.figpar(3, 3, 15)

fig.tight_layout()

plt.subplots_adjust(hspace = 0.00)

ss = 10.1

for i in range(len(w) - 1, -1, -1):

    neckel = P(A[i, :], mu)

    ax[0].plot(mu, neckel, label = str(w[i]) + ' nm', color = col[i])

    clv_atl = int_at[idx_at[i], :] / int_at[idx_at[i], 0]

    clv_ne0 = int_n0[idx_n0[i], :] / int_n0[idx_n0[i], 0]
    clv_ne1 = int_n1[idx_n1[i], :] / int_n1[idx_n1[i], 0]
    clv_ne2 = int_n2[idx_n2[i], :] / int_n2[idx_n2[i], 0]

    if i == 0: break

    ax[1].scatter(mu, (neckel - clv_atl) * 100 / neckel, s = ss, marker = '*', color = col[i])

    ax[1].scatter(mu, (neckel - clv_ne0) * 100 / neckel, s = ss * 2, marker = 'o', facecolors = 'none', edgecolors = col[i])
    ax[1].scatter(mu, (neckel - clv_ne1) * 100 / neckel, s = ss * 2, marker = 's', facecolors = 'none', edgecolors = col[i])
    ax[1].scatter(mu, (neckel - clv_ne2) * 100 / neckel, s = ss * 2, marker = '+', color = col[i])

ax[1].scatter(mu, (neckel - clv_atl) * 100 / neckel, s = ss, marker = '*', color = col[i], label = 'ATLAS9, LTE, U99')

ax[1].scatter(mu, (neckel - clv_ne0) * 100 / neckel, s = ss * 2, marker = 'o', facecolors = 'none', edgecolors = col[i], label = 'NESSY, LTE, U99')
ax[1].scatter(mu, (neckel - clv_ne1) * 100 / neckel, s = ss * 2, marker = 's', facecolors = 'none', edgecolors = col[i], label = 'NESSY, NLTE, U99')
ax[1].scatter(mu, (neckel - clv_ne2) * 100 / neckel, s = ss * 2, marker = '+', color = col[i], label = 'NESSY, NLTE, FAL99')

for i in range(len(ax)):

    ax[i].set_xlim(1.0, 0.0)
    ax[i].set_ylabel(r'$I(\mu) / I_c$')

ax[1].set_xlabel(r'$\mu$')

leg0 = ax[0].legend(framealpha = 1, loc = 3, handletextpad = 1, prop = {'size': 17.0})
leg1 = ax[1].legend(framealpha = 1, loc = 3, handletextpad = 1, prop = {'size': 17.0})

for obj in leg0.legendHandles: obj.set_linewidth(3.0)

auxplt.savepdf('var/neckel_clv_comp_2')
