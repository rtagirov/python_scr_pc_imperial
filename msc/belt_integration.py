import numpy as np
import math

import sys

if not '../aux/' in sys.path: sys.path.append('../aux/')

import importlib

import auxplt; importlib.reload(auxplt)
import paths;  importlib.reload(paths)
import phys;   importlib.reload(phys)

from tqdm import tqdm

def find_bracket(mu_grid, mu):

    i = np.abs(mu_grid - mu).argmin()

    if mu <= mu_grid[i - 1] and mu >= mu_grid[i]:

        return [i - 1, i], np.array([mu_grid[i - 1], mu_grid[i]])

    if mu <= mu_grid[i] and mu >= mu_grid[i + 1]:

        return [i, i + 1], np.array([mu_grid[i], mu_grid[i + 1]])

def interp_clv(idx, mu_bracket, I, mu):

    mu0 = mu_bracket[0]
    mu1 = mu_bracket[1]

    I0 = I[idx[0]]
    I1 = I[idx[1]]

    Imu = I0 + (mu - mu0) * (I1 - I0) / (mu1 - mu0)

    return Imu

mu_grid = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05])

p = np.sqrt(1.0 - mu_grid**2.0)

p_mid = np.zeros(len(p) - 1)
p_rar = np.zeros(len(p) - 1)

for i in range(len(p_mid)):

    p_mid[i] = (p[i] + p[i + 1]) / 2.0
    p_rar[i] = p[i + 1]**2.0 - p[i]**2.0

N = 100

theta_b = 5.0

#models = ['A', 'B', 'C', 'F', 'P', 'S']
models = ['B', 'C', 'F', 'P', 'S']
#models = ['A']

for model in models:

    clv = np.loadtxt(paths.it1f + 'tatiana/' + model + '/CLV_UVI')

    wvl = clv[:, 0] / 10.0

    I = clv[:, 1 : ]

    for theta_u in [30.0, 40.0]:

        yu = np.sqrt(1.0 - np.cos(theta_u * math.pi / 180)**2.0)
        yb = np.sqrt(1.0 - np.cos(theta_b * math.pi / 180)**2.0)

        xu = np.sqrt(1.0 - yu**2.0)
        xb = np.sqrt(1.0 - yb**2.0)

        Ar = (yu - yb) * xu
        At = (yu - yb) * (xb - xu) / 2.0

        phi_u = np.arctan(yu / xu)
        phi_b = np.arctan(yb / xb)

        theta = phi_u - phi_b

        As = (theta - np.sin(theta)) / 2.0

        Ab = 4.0 * (Ar + At + As)

#       print('Area fraction ', 4.0 * (At + As) / Ab, Ab / np.pi)

#       continue

        x_grid = np.linspace(0,  1, N)

        delta_x = x_grid[1] - x_grid[0]

        x_grid = np.delete(x_grid + delta_x / 2.0, len(x_grid) - 1)

        I_belt = np.zeros(len(wvl))

        for k in tqdm(range(len(wvl)), desc = 'model ' + model + ', ' + 'theta_u = ' + str(theta_u)):
#       for k in range(len(wvl)):

#            disc_counter = 0

            for i, x in enumerate(x_grid):

                yt = min(yu, np.sqrt(1.0 - x**2.0))

                y_grid = np.linspace(yb, yt, N)

                delta_y = y_grid[1] - y_grid[0]

                y_grid = np.delete(y_grid + delta_y / 2.0, len(y_grid) - 1)

                for j, y in enumerate(y_grid):

                    if x > np.sqrt(1.0 - yu**2.0) and (y + delta_y / 2.0 > yt or x + delta_x / 2.0 > np.sqrt(1.0 - y**2.0)):

#                        disc_counter += 1

                        continue

                    if 1.0 - x**2.0 - y**2.0 < 0.0:

                        continue

                    else:

                        mu = np.sqrt(1.0 - x**2.0 - y**2.0)

                    if mu <= 0.05: continue

                    idx, mu_bracket = find_bracket(mu_grid, mu)

                    Ixy = interp_clv(idx, mu_bracket, I[k, :], mu)

                    I_belt[k] += Ixy * delta_x * delta_y / math.pi

#            print('lalala ', disc_counter)

#            sys.exit()

        I_belt *= math.pi * (phys.r_sun / phys.au)**2.0
        I_belt *= phys.c / (wvl * 1.0e-7)**2.0
        I_belt *= 1.0e-10
        I_belt *= 4.0

        if theta_u == 30.0: np.savetxt('sbelt_' + model, np.transpose((wvl, I_belt)), fmt = ('%10.4f', '%10.4E'), delimiter = '  ')
        if theta_u == 40.0: np.savetxt('fbelt_' + model, np.transpose((wvl, I_belt)), fmt = ('%10.4f', '%10.4E'), delimiter = '  ')

        if theta_u == 40.0:

            I_disk = np.zeros(len(wvl))

            for k in tqdm(range(len(wvl)), desc = 'full disk'):

                for i, p in enumerate(p_mid):

                    mu_mid = np.sqrt(1.0 - p**2.0)

                    idx, mu_bracket = find_bracket(mu_grid, mu_mid)

                    I_mid = interp_clv(idx, mu_bracket, I[k, :], mu_mid)

                    I_disk[k] += I_mid * p_rar[i]

            I_disk *= math.pi * (phys.r_sun / phys.au)**2.0
            I_disk *= phys.c / (wvl * 1.0e-7)**2.0
            I_disk *= 1.0e-10

            np.savetxt('disk_' + model, np.transpose((wvl, I_disk)), fmt = ('%10.4f', '%10.4E'), delimiter = '  ')
