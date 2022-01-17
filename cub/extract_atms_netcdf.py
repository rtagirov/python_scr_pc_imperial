import netCDF4

import matplotlib.pyplot as plt
import numpy as np

import random
import itertools

import sys
import os
import pathlib

from tqdm import tqdm

snapshot = str(int(np.loadtxt('snapshot.inp')))

r1 = random.sample(range(512), 30)
r2 = random.sample(range(512), 30)
#r1 = [439, 241, 395, 237, 205]
#r2 = [197, 418, 495, 327, 46]

dpns = [30, 40, 50, 60, 70, 80, 82, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210]
#dpns = [40, 50, 60, 70, 80, 82, 90, 120]
#dpns = [40, 50, 60, 70, 80, 90, 120]
#dpns = np.arange(0, 11)

for n in tqdm(dpns):

#    z = netCDF4.Dataset('./nc/Z_onTau.'   + snapshot + '.nc' + '.' + str(n))['Z']
#    T = netCDF4.Dataset('./nc/T_onTau.'   + snapshot + '.nc' + '.' + str(n))['T']
#    p = netCDF4.Dataset('./nc/P_onTau.'   + snapshot + '.nc' + '.' + str(n))['P']
#    d = netCDF4.Dataset('./nc/rho_onTau.' + snapshot + '.nc' + '.' + str(n))['R']
    z = netCDF4.Dataset('./nc/Z_onTau.'   + snapshot + '.' + str(n) + '.nc')['Z']
    T = netCDF4.Dataset('./nc/T_onTau.'   + snapshot + '.' + str(n) + '.nc')['T']
    p = netCDF4.Dataset('./nc/P_onTau.'   + snapshot + '.' + str(n) + '.nc')['P']
    d = netCDF4.Dataset('./nc/rho_onTau.' + snapshot + '.' + str(n) + '.nc')['R']

    z = np.array(z) / 1e+5
    T = np.array(T)
    p = np.array(p)
    d = np.array(d)

    if not os.path.isdir('./atms/' + str(n)):

        pathlib.Path('./atms/' + str(n)).mkdir(parents=True, exist_ok=True)

    k = 1

    for i, j in itertools.product(r1, r2):

        zk = z[:, j, i]
        Tk = T[:, j, i]
        pk = p[:, j, i]
        dk = d[:, j, i]

        np.savetxt('./atms/' + str(n) + '/atm.' + str(k), \
                   np.column_stack([zk, Tk, pk, dk]), \
                   fmt = ('%6.1f', '%7.1f', '%7.5e', '%7.5e'), delimiter = '  ')

        k += 1

T = np.fromfile('eosT.'     + snapshot + '.bin', dtype = np.float32)
p = np.fromfile('eosP.'     + snapshot + '.bin', dtype = np.float32)
d = np.fromfile('result_0.' + snapshot + '.bin', dtype = np.float32)

dims = np.loadtxt('dims.inp', comments = '!')

Nz = int(dims[0, 0])

T = T.reshape(512, 512, Nz)
p = p.reshape(512, 512, Nz)
d = d.reshape(512, 512, Nz)

dz = dims[1, 0] / 1e+5

z = np.arange(Nz - 1, -1, -1) * dz

if not os.path.isdir('./atms/' + str(Nz)):

    os.mkdir('./atms/' + str(Nz))

k = 1

for i, j in itertools.product(r1, r2):

    Tk = np.flip(T[i, j, :])
    pk = np.flip(p[i, j, :])
    dk = np.flip(d[i, j, :])

    np.savetxt('./atms/' + str(Nz) + '/atm.' + str(k), \
               np.column_stack([z, Tk, pk, dk]), \
               fmt = ('%6.1f', '%7.1f', '%7.5e', '%7.5e'), delimiter = '  ')

    k += 1
