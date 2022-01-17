import subprocess
import os
import sys

if not '/mnt/SSD/sim/python/src/aux/' in sys.path: sys.path.append('/mnt/SSD/sim/python/src/aux/')

import paths

#runs = ['Q kur 0', 'F kur 0', 'Q kur 1', 'F kur 1', 'Q fal 1', 'F fal 1']
runs = ['Q kur 0', 'F kur 0', 'Q fal 1', 'F fal 1']
#runs = ['Q fal 1', 'F fal 1']

os.system('rm -f report')

for run in runs:

    args = run.split()

    command = 'python ' + paths.pydir + 'var/rings.py ' + args[0] + ' ' + args[1] + ' ' + args[2] + ' &'

    proc = subprocess.Popen(command.split())

    proc.wait()

    f = open('report', 'a+')

    f.write(run + '\n')

    f.close()
