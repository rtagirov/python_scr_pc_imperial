import numpy as np

import importlib
import julian; importlib.reload(julian)
import sys
import math

from tqdm import tqdm

if not '../aux/' in sys.path: sys.path.append('../aux/')

import auxfunc; importlib.reload(auxfunc)
import auxfile; importlib.reload(auxfile)
import paths;   importlib.reload(paths)
import phys;    importlib.reload(phys)


def time(jdate):

    leap_years = [2012, 2016]

    t = np.zeros(len(jdate))

    for i in range(len(jdate)):

        cdate = julian.from_jd(jdate[i])

        dd = 31.0

        if cdate.month % 2 == 0:

            dd = 30.0

            if cdate.month == 2:

                dd = 28.0

                for year in leap_years:

                    if year == cdate.year:

                        dd = 29.0

        t[i] = cdate.year + (cdate.month + cdate.day / dd) / 12.0

    return t

def bvar(ssi, w, w1, w2):

    nd = len(ssi[0, :])

    bv = np.zeros(nd)

    idx = np.where((w >= w1) & (w <= w2))[0]

    for i in tqdm(range(nd), ncols = auxfunc.term_width(), desc = 'Integrating over the band'):

        for j in idx:

            bv[i] += (ssi[j, i] + ssi[j + 1, i]) * (w[j + 1] - w[j]) / 2.0

    return bv

def svar(ssi, t):

    nw = len(ssi[:, 0])

    nd = len(ssi[0, :])

    sv = np.zeros(nw)

    for i in tqdm(range(nw), ncols = auxfunc.term_width(), desc = 'Integrating over time'):

        for j in range(nd - 1):

            sv[i] += (ssi[i, j] + ssi[i, j + 1]) * (t[j + 1] - t[j]) / 2.0

    sv = sv / (t[len(t) - 1] - t[0]) / 365.0

    return sv

def mmmvar(ssi, t, t1, t2, t3, t4):

    nw = len(ssi[:, 0])

    mv = np.zeros(nw)

    idx12 = np.where((t >= t1) & (t <= t2))[0]
    idx34 = np.where((t >= t3) & (t <= t4))[0]

    for i in tqdm(range(nw), ncols = auxfunc.term_width(), desc = 'Averaging over time'):

        ma = np.mean(ssi[i, idx34])

        mi = np.mean(ssi[i, idx12])

        av = (ma + mi) / 2.0

#        mv[i] = (ma - mi) * 100.0 / av
        mv[i] = (ma - mi) * 100.0 / mi

#    mv[np.where(np.isnan(mv))] = 1.0e-99

    return mv

def spec(mode, model):

    if mode != 'lte' and mode != 'nlt': auxsys.abort('varfunc.py: spec: mode is not recognized. Abort.')

    if mode == 'nlt': path = paths.satnlt
    if mode == 'lte': path = paths.satlte

    spec_dir = 'intensity_models/'

    path = path + spec_dir + model

    if mode == 'nlt': nw = auxfile.num_lines(path)
    if mode == 'lte': nw = auxfile.num_lines(path) // 3

    w = np.zeros(nw)

    its = np.zeros((nw, 11))

    f = open(path, 'r')

    if mode == 'lte':

        lc = 1

        wc = 0

        for l in f:

            a = l.split()

            if lc % 3 == 1:

                w[wc] = float(a[1])

            if lc % 3 == 2:

               for j in range(len(a)):

                   its[wc, j] = float(a[j])

            if lc % 3 == 0:

               for j in range(len(a)):

                   its[wc, j + 8] = float(a[j])

               wc += 1

            lc += 1

    if mode == 'nlt':

        wc = 0

        for l in f:

            a = l.split()

            w[wc] = float(a[1])

            for j in range(3, len(a)):

                its[wc, j - 3] = float(a[j])

            wc += 1

    f.close()

    mu = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05])

    p = np.sqrt(1.0 - mu**2.0)

    dii = np.zeros(nw)

    for i in range(nw):

        for j in range(1, len(p)):

            dii[i] += 0.5 * (its[i, j - 1] + its[i, j]) * (p[j]**2.0 - p[j - 1]**2.0)

    flux = dii * phys.c / (w * 1.0e-7)**2.0 * 1.0e-7 * (phys.r_sun / phys.au)**2.0 * 1.0e-3 * math.pi

    return w, flux

#def readssi(path):

#    nl = 16727310

#    nw = 1070

#    nd = 15633

#    wvl = np.zeros(nw)

#    t = np.zeros(nd)

#    jd, w1, w2, f, s = np.loadtxt(path, skiprows = 28, unpack = True)

#    for i in range(nw):

#        wvl[i] = (w1[i] + w2[i]) / 2.0

#    t[0] = jd[0]

#    tc = 1

#    for i in range(1, nl):

#        if jd[i] != jd[i - 1]:

#            t[tc] = jd[i]

#            tc += 1

#    tau = time(t)

#    ssi = np.reshape(f, (nw, nd))
#    ssi = np.reshape(f, (nd, nw))

#    return wvl, tau, ssi

def readssi(path):

    nw = 1070

#    nd = 15633
    nd = 2287

    wvl = np.zeros(nw)

    t = np.zeros(nd)

    ssi = np.zeros((nw, nd))

    f = open(path, 'r')

    j = -1

    wc = 0

    lc = 1

    for line in f:

        if line[0] != ' ':

            continue

        va = line.split()

        if int(va[0]) < 2455317 or int(va[0]) > 2457603:

            continue

        if lc <= nw:

            wvl[wc] = (float(va[1]) + float(va[2])) / 2.0

            wc += 1

        if lc % nw == 1:

            j += 1

            i = 0

            t[j] = float(va[0])

        ssi[i, j] = float(va[3])

        i += 1

        lc += 1

    tau = time(t)

    f.close()

    return wvl, tau, ssi
