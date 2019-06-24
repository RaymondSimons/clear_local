import astropy
from astropy.io import fits
import glob
from glob import glob
import numpy
from numpy import *
import matplotlib.pyplot as plt
plt.ioff()
plt.close('all')
figdir = '/Users/rsimons/Desktop/clear/figures/continuum'


fields = ['GS4']

for field in fields:
    fls = glob('/Volumes/pegasus/clear/grizli_extractions/%s/*/Prep/*GrismFLT.fits'%field)

    for f, fl in enumerate(fls):

        print (f, fl)
        a = fits.open(fl)
        fig, axes = plt.subplots(2,2, figsize = (10,10))

        xmn = 500
        xmx = 1100
        vmn = 0.0001
        vmx = 0.08

        vmn_im = 0.
        vmx_im = 1.e-20
        im = a['DREF'].data
        #im-=im.ravel().min()
        axes[0,0].imshow(im, cmap = 'Greys_r', vmin = vmn_im, vmax = vmx_im)
        axes[0,1].imshow(a['GSCI'].data, cmap = 'viridis', vmin = vmn, vmax = vmx)
        axes[1,0].imshow(a['MODEL'].data, cmap = 'viridis',  vmin = vmn, vmax = vmx)
        axes[1,1].imshow(a['GSCI'].data - a['MODEL'].data, cmap = 'viridis',  vmin = vmn, vmax = vmx)

        for ax in axes.ravel(): 
            ax.axis('off')
            ax.set_xlim(xmn, xmx)
            ax.set_ylim(xmn, xmx)

        outfile = fl.split('/')[-1].replace('GrismFLT.fits', '%s.png'%field)
        fig.subplots_adjust(hspace = 0.0, wspace = 0.0, top = 1.0, bottom = 0.0, right = 1.0, left = 0.0)
        fig.savefig(figdir + '/' + outfile)




