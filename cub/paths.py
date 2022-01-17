import socket

hostname = socket.gethostname()

if hostname == 'ph-rtagirov':  sim = '/mnt/SSD/sim'
if hostname == 'Tagirov-SLM8': sim = '/mnt/SSD/sim'

if hostname == 'schiphol':    sim = '/Users/rinattagirov/sim'

npz =    sim + '/python/npz/'

it0h =   sim + '/runs/hminus/IT0/'
it1h =   sim + '/runs/hminus/IT1/'
it2h =   sim + '/runs/hminus/IT2/'
it3h =   sim + '/runs/hminus/IT3/'

it0f =   sim + '/runs/fioss/IT0/'
it1f =   sim + '/runs/fioss/IT1/'
it2f =   sim + '/runs/fioss/IT2/'
it3f =   sim + '/runs/fioss/IT3/'

figdir = sim + '/python/fig/'
mscdir = sim + '/python/src/msc/'

inp =    sim + '/python/inp/'
out =    sim + '/python/out/'

idlinp = sim + '/idl/input/'
idlout = sim + '/idl/output/'

pydir =  sim + '/python/src/'

atlruns = sim + '/runs/atlas9/'

lev = 'NLTE/LEV/'

satnlt = sim + '/satire_nlte/'
satlte = sim + '/satire_lte/'
satrnm = sim + '/satire/'
satnss = sim + '/satire_nlte_ss/'

nessy = sim + '/nessy'
