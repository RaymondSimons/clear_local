import astropy
from astropy.io import fits
import glob
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy import *

plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['xtick.labelsize'] = 12


field, di = 'ERSPRIME', 40192

metal_dir = '/Volumes/pegasus/clear/metal_maps/local_testing'

fls = glob(metal_dir + '/%s_%.5i_metals.fits'%(field, di))
for f, fl in enumerate(fls):
    a = fits.open(fl)
    cal_array = array(['r3', 'r2', 'r23', 'o32', 'all'])
    cal_array = array(['r23', 'all'])

    xs = arange(35, 45)
    ys = arange(35, 45)

    for i in xs:
        for j in ys:
            fig, axes = plt.subplots(1,len(cal_array), figsize = (3*len(cal_array), 3))
            for n, calib in enumerate(cal_array):
                print (calib)
                cube = a['Z_%s_full'%calib.upper()].data
                ratio_calib = calib.upper()
                if calib == 'all': ratio_calib = 'R23'
                R = a['%s'%ratio_calib].data
                eR = a['e%s'%ratio_calib].data/R/log(10)

                R = log10(R)
                post = cube[i,j,:]
                R_ij  = R[i,j] 
                eR_ij = eR[i,j] 
                if ~isnan(post[0]):
                    axes[n].hist(post, normed = True, histtype = 'step', color = 'black', alpha = 0.3, linewidth = 0.2, bins = linspace(6.8, 9.2, 100))
                    axes[n].annotate(r'%.2f $\pm$ %.2f'%(R_ij, eR_ij), (0.7, 0.9), xycoords = 'axes fraction')
                axes[n].set_xlabel(r'$\log$ (O/H)$_{\text{%s}}$ + 12'%calib.upper(), fontsize = 15)



            axes[n].set_xlim(6.8, 9.5)


            fig.tight_layout()
            fig.savefig('/Users/rsimons/Desktop/clear/figures/metal_pixels_posteriors/%s_%.5i_%i_%i.png'%(field, di, i, j), dpi = 300)

            plt.close('all')