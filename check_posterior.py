import astropy
from astropy.io import fits
import glob
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.convolution import Box1DKernel, convolve_fft
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages
from numpy import *

plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['xtick.labelsize'] = 12


field, di = 'ERSPRIME', 40192
field, di = 'GN5', 32616
field, di = 'GN7', 13197
field, di = 'GN7', 17477
field, di = 'GS3', 38912
field, di = 'GN3', 35204
#field, di = 'GN3', 30204
#field, di = 'GN1', 37395


metal_dir = '/Volumes/pegasus/clear/metal_maps'

fls = glob(metal_dir + '/%s_%.5i_metals.fits'%(field, di))

kern = Box1DKernel(2)
for f, fl in enumerate(fls):
    a = fits.open(fl)
    #cal_array = array([ 'r2', 'r3', 'r23', 'o32', 'all'])
    #cal_array = array(['r23', 'all'])
    cal_array = array(['r3', 'all'])

    xs = arange(28, 52)
    ys = arange(39, 40)

    with PdfPages('/Users/rsimons/Desktop/clear/figures/metal_pixels_posteriors/%s_%.5i.pdf'%(field, di)) as pdf:
        for i in xs:
            print (i)
            for j in ys:
                fig, axes = plt.subplots(1,len(cal_array), figsize = (3*len(cal_array), 3))
                for n, calib in enumerate(cal_array):
                    Z = a['Z_%s'%calib.upper()].data
                    cube = a['Z_%s_full'%calib.upper()].data
                    ratio_calib = calib.upper()
                    if calib == 'all': ratio_calib = 'R3'
                    R = a['%s'%ratio_calib].data
                    eR = a['e%s'%ratio_calib].data/R/log(10)

                    R = log10(R)
                    post = cube[i,j,:]
                    R_ij  = R[i,j] 
                    eR_ij = eR[i,j]
                    if eR_ij < 1./2./log(10): 
                        if ~isnan(post[0]):
                            print R_ij
                            axes[n].hist(post, normed = True, histtype = 'step', color = 'black', alpha = 1.0, linewidth = 0.2, bins = linspace(6.8, 9.2, 100))
                            #axes[n].annotate(r'%.2f $\pm$ %.2f'%(R_ij, eR_ij), (0.7, 0.9), xycoords = 'axes fraction')

                            hst, xedges = histogram(post, normed = True, bins = arange(6.8, 9.2, 0.1))
                            xedges = mean((xedges[:-1], xedges[1:]), axis = 0)

                            hst = convolve_fft(hst, kern)
                            axes[n].plot(xedges, hst, 'r-', linewidth = 0.3)

                            gd_mx = xedges[argmax(hst)]
                            axes[n].axvline(gd_mx, color = 'red', linestyle = '-', alpha = 0.4)


                        z50 = Z[i,j][0]
                        z84 = z50 + Z[i,j][1]
                        z16 = z50 - Z[i,j][2]

                        axes[n].axvline(z50, color = 'blue', linestyle = '-', alpha = 0.4)
                        axes[n].axvline(z84, color = 'blue', linestyle = '--', alpha = 0.4)
                        axes[n].axvline(z16, color = 'blue', linestyle = '--', alpha = 0.4)



                    axes[n].set_xlabel(r'$\log$ (O/H)$_{\text{%s}}$ + 12'%calib.upper(), fontsize = 15)



                axes[n].set_xlim(6.8, 9.5)


                fig.tight_layout()
                #fig.savefig('/Users/rsimons/Desktop/clear/figures/metal_pixels_posteriors/%s_%.5i_%i_%i.png'%(field, di, i, j), dpi = 300)
                pdf.savefig()
                plt.close('all')