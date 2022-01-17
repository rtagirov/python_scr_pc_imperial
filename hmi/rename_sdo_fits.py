import glob
import os
import sys

names = glob.glob('*_fit')

for name in names:

    year =  name[0 : 4]

    month = name[5 : 7]

    day =   name[8 : 10]

#    kind =  name[24 : 27]

#    os.system('mv ' + name + ' ' + year + '_' + month + '_' + day + '_' + kind + '.fits')
    os.system('mv ' + name + ' ' + year + '_' + month + '_' + day + '.fits')
