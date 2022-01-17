from astropy.io    import fits
from tqdm          import tqdm
from scipy.ndimage import filters

import numpy             as np
import matplotlib.pyplot as plt

import os
import sys
import importlib
import itertools

import paths;   importlib.reload(paths)
import auxfunc; importlib.reload(auxfunc)
import mask;    importlib.reload(mask)

int_file_names = os.listdir(paths.inp + 'int/')

alpha_s = alpha_f = np.array([])

for i in tqdm(range(len(int_file_names)), ncols = auxfunc.term_width()):

    int_file_name = int_file_names[i]

    date = int_file_name[0 : 10]

    filename_int = 'int/' + date + '.fits'
    filename_mag = 'mag/' + date + '.fits'

    hdulist_int = fits.open(paths.inp + filename_int)
    hdulist_mag = fits.open(paths.inp + filename_mag)

    hdu_int = hdulist_int[0]
    hdu_mag = hdulist_mag[0]

    x0 =   hdu_int.header['X0']
    y0 =   hdu_int.header['Y0']
    rsun = hdu_int.header['R_SUN']

    r = 0.9 * rsun

    gs = 10

    I = filters.gaussian_filter(hdu_int.data, gs)
    B = filters.gaussian_filter(hdu_mag.data, gs)

    hdulist_int.close()
    hdulist_mag.close()

    rad = np.zeros((4096, 4096))

    for y, x in itertools.product(range(4096), range(4096)):

        rad[y, x] = (x - x0)**2.0 + (y - y0)**2.0

    rad = np.sqrt(rad)

    spt_mask, a_s = mask.spt(I, x0, y0, r, rsun, rad)
    mag_mask, a_m = mask.mag(B,         r, rsun, rad)

    fac_mask = mag_mask - spt_mask

    fac_mask[np.where(fac_mask == -1)] = 0

    a_f = len(fac_mask[np.where(fac_mask == 1.0)]) / len(fac_mask[np.where(rad <= r)])

    alpha_s = np.concatenate((alpha_s, np.array([a_s])))
    alpha_f = np.concatenate((alpha_f, np.array([a_f])))

np.savez(paths.npz + 'ff', alpha_s = alpha_s, alpha_f = alpha_f)

ff = np.load(paths.npz + 'ff.npz')

alpha_s = ff['alpha_s']
alpha_f = ff['alpha_f']

plt.figure(1)

plt.scatter(alpha_f, alpha_s)

plt.show()
