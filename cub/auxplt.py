from pylab                           import rcParams
from matplotlib                      import rc
from matplotlib.backends.backend_pdf import PdfPages

#import seaborn as sns
import numpy as np

#import importlib

import paths

reload(paths)

#importlib.reload(paths)

def figpar(xtick_maj_pad = 15, ytick_maj_pad = 10, fontsize = 10):

#    sns.set(rc={"figure.figsize": figsize})
#    np.random.seed(sum(map(ord, "palettes")))

#    seaborn.palplot(seaborn.color_palette('cubehelix', 11))

#    seaborn.palplot(seaborn.diverging_palette(255, 133, l=60, n=11, center='dark'))

#    sns.palplot(sns.color_palette('dark'))

#    sns.set(font_scale = 1.5)

#    sns.set_style(style = 'white')

    rcParams['xtick.major.pad'] = xtick_maj_pad
    rcParams['ytick.major.pad'] = ytick_maj_pad

    rc('font', size = fontsize)
    rc('font', family = 'serif')
    rc('font', serif = 'Times')
    rc('text', usetex = True)

def savepdf(filename, figdir = paths.figdir):

    pp = PdfPages(figdir + filename + '.pdf')

    pp.savefig(bbox_inches = 'tight')
    pp.close()
