import numpy as np

import subprocess
import os
import sys

if not '/mnt/SSD/sim/python/src/aux/' in sys.path: sys.path.append('/mnt/SSD/sim/python/src/aux/')

import paths

#p = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95])
#mu = np.array([1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05])
#mu = np.array([0.87])

p = np.sqrt(np.arange(0, 11)) / np.sqrt(10)

p_mid = np.zeros(len(p) - 1)

for i in range(len(p_mid)):

    p_mid[i] = (p[i + 1] + p[i]) / 2

mu = np.sqrt(1 - p_mid**2)

ar = sys.argv[1]

model = sys.argv[2]

it = sys.argv[3]

if it == '0': firuns = paths.it0f
if it == '1': firuns = paths.it1f

#gen_dir = firuns + 'var/' + ar + '/' + model + '_for_rings'
gen_dir = firuns + 'var_od/' + ar + '/' + model + '_for_rings'
#hmi_dir = paths.it0h + 'var/' + ar + '/' + model
hmi_dir = paths.it0h + 'var_od/' + ar + '/' + model

os.system('cp ' + hmi_dir + '/ATM_STR ' + gen_dir)

prev_dir = os.getcwd()

os.chdir(paths.nessy)

os.system('cp -r ./obj/fioss.exe ./src -t ' + gen_dir)
os.system('cp    ./inp/car/fioss.cards.template.central.ray ' + gen_dir + '/CARDS.TEMPLATE')

os.chdir(prev_dir)

remove = ['JOB-*', 'clv_uvi', 'CLV', 'CLV_RAD', 'CLV_UVI', 'contr.txt', \
          'form.height', 'JOB-*', 'log.rad', 'lopa', 'radio', 'output', 'title', 'abemlin']

for i in range(len(mu)):

    run_name = str(mu[i])[:4]

#    run_dir = firuns + 'var/' + ar + '/rings/' + model + '/' + run_name
    run_dir = firuns + 'var_od/' + ar + '/rings/' + model + '/' + run_name

    if os.path.exists(run_dir): os.system('rm -rf ' + run_dir)

    os.system('cp -r ' + gen_dir + ' ' + run_dir)

    f = open(run_dir + '/mu', 'w')

    f.write(str(mu[i]))

    f.close()

#    command = '/mnt/SSD/sim/nessy/scr/run_nessy var/' + ar + '/rings/' + model + '/' + run_name + ' -i ' + it + ' -P 16 -F uvi &'
    command = '/mnt/SSD/sim/nessy/scr/run_nessy var_od/' + ar + '/rings/' + model + '/' + run_name + ' -i ' + it + ' -P 12 -F uvi &'

    proc = subprocess.Popen(command.split())

    proc.wait()

    for j, elem in enumerate(remove):

        os.system('rm -rf ' + run_dir + '/' + elem)
