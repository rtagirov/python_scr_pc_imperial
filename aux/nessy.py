import numpy as np
import math as m

import os
import importlib

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import auxfunc;   importlib.reload(auxfunc)
import paths;     importlib.reload(paths)
import spec;      importlib.reload(spec)
import auxfile;   importlib.reload(auxfile)
import phys;      importlib.reload(phys)

from tqdm import tqdm

def read_spec(frpath, wvl1 = 1005., wvl2 = 9995., mode = 'dir'):

    if mode == 'dir':

       spath = frpath + '/mdisp/'

       nf =  2000 # number of lines in each .mdisp file

       ivl = 10   # length of the spectral interval of each .mdisp file (angstroems)
       mid = 5    # distance to the middle of each .mdisp file (angstroems)

       nint = int(ivl)

       wmin = int(wvl1)
       wmax = int(wvl2)

       imin = wmin - ((wmin - mid) % nint)
       imax = wmax - ((wmax - mid) % nint) + nint

       na = (imax - imin) / nint + 1 # number of arrays

       nw = na * nf

       wvl = np.array([])
       its = np.array([]) 

#       for i in tqdm(range(int(na)), ncols = auxfunc.term_width(), desc = 'Reading ' + frpath):
       for i in tqdm(range(int(na)), desc = 'Reading ' + frpath):

           idx = str(imin + i * nint)

           f = spath + idx + '.mdisp'

           spec = np.loadtxt(f)

           wvl_i = spec[:, 0]; wvl = np.concatenate((wvl, wvl_i))
           its_i = spec[:, 1]; its = np.concatenate((its, its_i))

    if mode == 'file':

        spath = frpath + '/SPE_UVI'

        spec = np.loadtxt(spath)

        wvl = spec[:, 0]
        its = spec[:, 1]

#   conversion to flux at 1 AU in W / m^2 / nm
    flu = its * phys.c / (wvl * 1.0e-8)**2.0 * 1.0e-7 * (phys.r_sun / phys.au)**2.0 * 1.0e-3 * m.pi
#    flu = its

    return wvl, flu

def read_lopa(frpath, wvl1 = 1005., wvl2 = 9995.):

    ndp = auxfile.num_lines(frpath + '/atm.inp')

    spath = frpath + '/lopa/'

    ivl = 10 # length of the spectral interval of each .tau file (angstroems)
    mid = 5  # distance to the middle of each .tau file (angstroems)

    nint = int(ivl)

    wmin = int(wvl1)
    wmax = int(wvl2)

    imin = wmin - ((wmin - mid) % nint)
    imax = wmax - ((wmax - mid) % nint) + nint

    idx = np.arange(imin, imax + ivl, ivl)

    for i in tqdm(range(len(idx)), ncols = auxfunc.term_width(), desc = 'Reading ' + frpath):

        wvl, opa = np.loadtxt(spath + str(idx[i]) + '.lopa', unpack = True)

        nw = int(wvl[0])

        wvl = wvl[1 : nw + 1]

        opa = opa.reshape(ndp, nw + 1) # 1D -> 2D conversion

        opa = opa[:, 1 : ]             # removing the depth point numbers

        if i == 0:

            wvls = wvl
            opas = opa

        if i != 0:

            wvls = np.concatenate((wvls, wvl))

            opas = np.concatenate((opas, opa), axis = 1)

    return wvls, opas

def read_tau(frpath, wvl1 = 1005., wvl2 = 9995., mode = 'dir'):

    ndp = auxfile.num_lines(frpath + 'atm.inp')

    if mode == 'dir':

       spath = frpath + '/tau/'

       nf =  2000 # number of lines in each .tau file

       ivl = 10   # length of the spectral interval of each .tau file (angstroems)
       mid = 5    # distance to the middle of each .tau file (angstroems)

       nint = int(ivl)

       wmin = int(wvl1)
       wmax = int(wvl2)

       imin = wmin - ((wmin - mid) % nint)
       imax = wmax - ((wmax - mid) % nint) + nint

       na = (imax - imin) / nint + 1 # number of arrays

       nw = na * nf

       wvl = np.array([])
       fid = np.array([]) 
       hei = np.array([]) 

       for i in tqdm(range(int(na)), ncols = auxfunc.term_width(), desc = 'Reading ' + frpath):

           idx = str(imin + i * nint)

           f = spath + idx + '.tau'

           tau = np.loadtxt(f)

           wvl_i = tau[:, 0]; wvl = np.concatenate((wvl, wvl_i))
           fid_i = tau[:, 4]; fid = np.concatenate((fid, fid_i))
           hei_i = tau[:, 1]; hei = np.concatenate((hei, hei_i))

    if mode == 'file':

        spath = frpath + '/TAU_UVI'

        tau = np.loadtxt(spath)

        wvl = tau[:, 0]
        fid = tau[:, 4]
        hei = tau[:, 1]

    hei[np.where(fid == ndp + 1)] = 0.0
    fid[np.where(fid == ndp + 1)] = ndp

    return wvl, fid.astype(int) - 1, hei

def weighted_form_temp(path, wvl1 = 1005., wvl2 = 9995.):

    h = np.loadtxt(path + '/atm.inp', usecols = [0])
    T = np.loadtxt(path + '/atm.inp', usecols = [1])

    ndp = len(T)

    nf =  2000 # number of lines in each .tau and .mdisp file

    ivl = 10   # length of the spectral interval of each .tau and .mdisp file (angstroems)
    mid = 5    # distance to the middle of each .tau and .mdisp file (angstroems)

    nint = int(ivl)

    wmin = int(wvl1)
    wmax = int(wvl2)

    imin = wmin - ((wmin - mid) % nint)
    imax = wmax - ((wmax - mid) % nint) + nint

    na = (imax - imin) / nint + 1 # number of arrays

    nw = na * nf

    wvl = []
    fhw = []
    fTw = []

    for i in tqdm(range(int(na)), ncols = auxfunc.term_width(), desc = 'Reading ' + path):

        idx = str(imin + i * nint)

        f1 = path + '/tau/' + idx + '.tau'
        f2 = path + '/mdisp/' + idx + '.mdisp'

        fid = np.loadtxt(f1, usecols = [4])
        I = np.loadtxt(f2, usecols = [1])

        fid[np.where(fid == ndp + 1)] = ndp

        fid -= 1

        fh = h[fid.astype(int)]
        fT = T[fid.astype(int)]

        wvl.append(float(idx))

        fhw.append(sum(I * fh) / sum(I))
        fTw.append(sum(I * fT) / sum(I))

    return np.array(wvl), np.array(fhw), np.array(fTw)

def weighted_form_height(path, wvl1 = 1005., wvl2 = 9995.):

    nf =  2000 # number of lines in each .tau and .mdisp file

    ivl = 10   # length of the spectral interval of each .tau and .mdisp file (angstroems)
    mid = 5    # distance to the middle of each .tau and .mdisp file (angstroems)

    nint = int(ivl)

    wmin = int(wvl1)
    wmax = int(wvl2)

    imin = wmin - ((wmin - mid) % nint)
    imax = wmax - ((wmax - mid) % nint) + nint

    na = (imax - imin) / nint + 1 # number of arrays

    nw = na * nf

    wvl = []
    fhw = []

#    for i in tqdm(range(int(na)), ncols = auxfunc.term_width(), desc = 'Reading ' + path):
    for i in tqdm(range(int(na)), desc = 'Reading ' + path):

        idx = str(imin + i * nint)

        f1 = path + '/tau/' + idx + '.tau'
        f2 = path + '/mdisp/' + idx + '.mdisp'

        h = np.loadtxt(f1, usecols = [1])
        I = np.loadtxt(f2, usecols = [1])

        wvl.append(float(idx))

        fhw.append(sum(I * h) / sum(I))

    return np.array(wvl), np.array(fhw)

def atl_con(path, wvl, flu, delta = 0.6 * 3.5):

#    wvl_min = wvl[0]
#    wvl_max = wvl[len(wvl) - 1]

#    nws = 2 * int(m.ceil(wvl_max - wvl_min))

#	nws = len(wvl_atl)

#    wvls = np.zeros(nws)

	wvl_atl = np.loadtxt(paths.inp + 'atlas3.dat', usecols = [0]) * 10.0

	flus = np.zeros(len(wvl_atl))

	pb = ProgressBar(widgets = auxfunc.pbar_widgets('Performing the ATLAS convolution of ' + path + '/'), term_width = auxfunc.terminal_width())

	for i in pb(range(len(wvl_atl))):

#        wvls[i] = wvl_min + (wvl_max - wvl_min) * i / (nws - 1)

		idx = np.where((wvl > wvl_atl[i] - 3.0 * delta) & (wvl < wvl_atl[i] + 3.0 * delta))

		for j in range(len(idx[0]) - 1):

			flus[i] = flus[i] + flu[idx[0][j]] * phi(wvl_atl[i], wvl[idx[0][j]], delta) * \
                                                 (wvl[idx[0][j + 1]] - wvl[idx[0][j]])

	return flus


def phi(wvl0, wvl, delta):

    return m.exp(-((wvl0 - wvl) / delta)**2.0) / m.sqrt(m.pi) / delta


def get_elems(filename):

    f = open(filename, 'r')

    elems = []; anums = []

    for line in f:

        if line[0 : 7] == 'ELEMENT':

            parts = line.split()

            elem = parts[2].strip(' () ')

            if parts[3] == ')': anum = int(float(parts[4]))

            if parts[3] != ')': anum = int(float(parts[3]))

            if elem == 'HE':
 
                elem = 'He'

                anum = 2

            elems.append(elem)

            anums.append(anum)

    f.close()

    return elems, anums


def read_popnum(path):

    srcdir = os.getcwd()

    os.chdir(path)

    os.system('grep -w LEVEL datom.inp > levels.out')

    dpn = auxfile.num_lines('atm.inp')

    ln = auxfile.num_lines('levels.out')

    rne = np.zeros(dpn)

    popnum = np.zeros((ln, dpn))

    dpn3 = m.ceil(dpn / 3.0)

    dpnln3 = int(dpn * ln / 3)

    head = int(10 + 2 * dpn3)

    os.system('head -n ' + str(head) + ' POPNUM | tail -n ' + str(dpn3) + ' > rne.out')

    head = int(head + 1 + dpnln3)

    os.system('head -n ' + str(head) + ' POPNUM | tail -n ' + str(dpnln3) + ' > popnum.out')

    lev = []

    f = open('levels.out', 'r')

    for line in f:

        parts = line.split()

        lev_name = parts[1]

        if (lev_name[0] == 'H'): lev_name = parts[1] + parts[2]

        lev_name = lev_name.replace('-', '')

        if (lev_name[0 : 2] == 'HI' or lev_name[0 : 2] == 'HM' or lev_name[0 : 4] == 'HEII'):

            lev_name = lev_name.replace('.', '')

        if (lev_name == 'HMINUS1'): lev_name = 'HMINUS'

        lev.append(lev_name)

    f.close()

    f = open('rne.out', 'r')

    i = 0

    for line in f:

        parts = line.split()

        for part in parts:

            rne[i] = float(part)

            i = i + 1

    f.close()

    f = open('popnum.out', 'r')

    i = 0; j = 0

    for line in f:

        parts = line.split()

        for part in parts:

            popnum[j, i] = float(part)

            i = i + 1

            if i != 0 and i % dpn == 0:

                i = 0

                j = j + 1

    f.close()

    os.system('rm levels.out rne.out popnum.out')

    os.chdir(srcdir)

    return lev, rne, popnum
