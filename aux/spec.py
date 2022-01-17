import numpy as np
import math as m

import importlib

import auxfunc; importlib.reload(auxfunc)

from tqdm import tqdm

def weigh_within_delta(wvl, val, I, delta, message = 'Averaging'):

    wvl_min = wvl[0]
    wvl_max = wvl[len(wvl) - 1]

    nws = int(m.ceil((wvl_max - wvl_min) / delta))

    wvls = np.zeros(nws)
    vals = np.zeros(nws)

    for i in tqdm(range(nws), ncols = auxfunc.term_width(), desc = message):

        wvls[i] = wvl_min + (wvl_max - wvl_min) * i / (nws - 1)

        idx = np.where((wvl > wvls[i] - delta / 2.0) & (wvl < wvls[i] + delta / 2.0))

        vals[i] = np.sum(val[idx] * I[idx]) / np.sum(I[idx])

    return wvls, vals

def mean_within_delta(wvl, flu, delta, message = 'Averaging'):

    wvl_min = wvl[0]
    wvl_max = wvl[len(wvl) - 1]

    nws = int(m.ceil((wvl_max - wvl_min) / delta))

    wvls = np.zeros(nws)
    flus = np.zeros(nws)

#    for i in tqdm(range(nws), ncols = auxfunc.term_width(), desc = message):
    for i in tqdm(range(nws), desc = message):

        wvls[i] = wvl_min + (wvl_max - wvl_min) * i / (nws - 1)

        idx = np.where((wvl > wvls[i] - delta / 2.0) & (wvl < wvls[i] + delta / 2.0))

        flus[i] = np.mean(flu[idx])

    return wvls, flus

def mean_over_grid(fh, wh, wl, message = 'Averaging'):

    nws = len(wl)

    fl = np.zeros(nws)

    for i in tqdm(range(nws), ncols = auxfunc.term_width(), desc = message):

        if i == 0:      

            idx = np.where((wh >= wl[i]) & (wh <= wl[i + 1]))

        elif i == nws - 1:

            idx = np.where((wh >= wl[i - 1]) & (wh <= wl[i]))

        else:

            idx = np.where((wh >= wl[i - 1]) & (wh <= wl[i + 1]))

        fl[i] = np.mean(fh[idx])

    return fl

def mean_within_delta_over_grid(fh, wh, delta, wl, message = 'Averaging'):

    nws = len(wl)

    fl = np.zeros(nws)

    for i in tqdm(range(nws), ncols = auxfunc.term_width(), desc = message):

        idx = np.where((wh >= wl[i] - delta / 2.0) & (wh <= wl[i] + delta / 2.0))

        fl[i] = np.mean(fh[idx])

    return fl

def sort_within_delta_over_grid(f, w, delta, wb, message = 'Averaging'):

    nws = len(wb)

    fs = np.array([])

    for i in tqdm(range(nws), ncols = auxfunc.term_width(), desc = message):

        idx = np.where((w >= wb[i] - delta / 2.0) & (w <= wb[i] + delta / 2.0))

        fs = np.concatenate((fs, np.sort(f[idx])))

    return fs
