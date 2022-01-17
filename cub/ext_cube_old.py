# program to extract data from FITS files into files to
# repressent a box of the stellar atmosphere
# including turbulent velocity

#import pyfits
import numpy as np
import os.path
import os
import sys
import time
import glob

from astropy.io import fits
from tqdm import tqdm

import importlib
import phys

importlib.reload(phys)

def getdata(f):

    hdulist = fits.open(f)

    hdu = hdulist[0]

    a = hdu.data

    return a

def swap_and_flip(a):

    a = np.swapaxes(a, 0, 2)
    a = np.swapaxes(a, 1, 2)

    a = a[:, :, ::-1]

    return a

"""
dz = height resolution: 
1.25 km for F3V
10 km for G2V
6 km for K0V
4 km for M0V
3.2 km for M2V

for the width and length, dx=dy, the resolution is
58.5938 km for F3V,
17.5781 km for G2V,
11.7188 km for K0V, 
4.8828 km for M0V,
3.0469 km for M2V.

apn - added points number (number of points added when extending an extracted atmospheric structure)
"""

num = '123000'
dz = 10.0
#apn = 24
apn = 0

abund_wh = np.loadtxt('abund.inp')

abund = np.zeros(30)

abund[0] = abund_wh[0]

abund[1] = 1.0 - sum(abund_wh)

abund[2 : len(abund)] = abund_wh[1 : len(abund_wh)]

apm = phys.average_particle_mass(abund)

T = getdata('eosT.'     + num + '.fits')
p = getdata('eosP.'     + num + '.fits')
d = getdata('result_0.' + num + '.fits')

T = swap_and_flip(T)
p = swap_and_flip(p)
d = swap_and_flip(d)

T_top = T[:, :, 0]
p_top = p[:, :, 0]

h1d = -np.arange(1, apn + 1) * dz * 1.0e+5

h = np.broadcast_to(h1d, T_top.shape + h1d.shape)

T_add = np.broadcast_to(T_top[...,None], T_top.shape + (apn,))
p_top = np.broadcast_to(p_top[...,None], p_top.shape + (apn,))

H = phys.boltz * T_add / apm / phys.grav_sun

p_add = p_top * np.exp(h / H)

p_add = p_add[:, :, ::-1]

d_add = apm * p_add / phys.boltz / T_add

T = np.concatenate((T_add, T), axis = 2)
p = np.concatenate((p_add, p), axis = 2)
d = np.concatenate((d_add, d), axis = 2)

T = T[:, :, ::-1]
p = p[:, :, ::-1]
d = d[:, :, ::-1]

T = T.flatten().astype(np.float32)
p = p.flatten().astype(np.float32)
d = d.flatten().astype(np.float32)

#T.tofile('eosT.'     + num + '.ext.' + str(apn) + '.bin')
#p.tofile('eosP.'     + num + '.ext.' + str(apn) + '.bin')
#d.tofile('result_0.' + num + '.ext.' + str(apn) + '.bin')
T.tofile('eosT.'     + num + '.bin')
p.tofile('eosP.'     + num + '.bin')
d.tofile('result_0.' + num + '.bin')
