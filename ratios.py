import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import click
import astropy
from astropy.io import fits
import glob
from glob import glob
import numpy as np
import emcee
import joblib
from joblib import Parallel, delayed
from numpy import *
from matplotlib.backends.backend_pdf import PdfPages
import photutils
plt.ioff()
plt.close('all')
from photutils import detect_sources
#mpl.rcParams['text.usetex'] = True


metal_maps_dir = '/Volumes/pegasus/clear/metal_maps'


fls = glob(metal_maps_dir + '/*metals.fits')



def lnprob(theta, x, y, yerr):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, yerr)


def lnprior(theta):
    m, b = theta
    if -5.0 < m < 5. and -5 < b < 5:
        return 0.0
    return -np.inf

def lnlike(theta, x, y, yerr):
    m, b = theta
    model = m * x + b
    inv_sigma2 = 1.0/(yerr**2)
    return -0.5*(np.sum((y-model)**2*inv_sigma2 - np.log(inv_sigma2)))






def make_fits(rat, line1, line2, line3):
    x = np.arange(-40, 40, 1)
    y = np.arange(-40, 40, 1)
    xv, yv = np.meshgrid(x, y)

    rmax_arc = 1.4
    xmn_pix = 26
    xmx_pix = 54
    with PdfPages('/Volumes/pegasus/clear/line_ratio_profiles/%s.pdf'%rat) as pdf:
        for f, fl in enumerate(fls):
            rv = 0.1*np.sqrt(xv**2. + yv**2.)
            rv = rv[xmn_pix:xmx_pix, xmn_pix:xmx_pix]

            save_name =  fl.split('/')[-1].strip('_metals.fits')
            if True:#save_name == 'GS2_45378':
                print (f)
                a = fits.open(fl)
                ext_info = a.info(output = False)
                ext_names = [ext[1] for ext in ext_info]
                if rat in ext_names:
                    cmap = mpl.cm.viridis
                    cmap.set_bad('k', 1.)

                    if line3 is not None:
                        fig = plt.figure(figsize = (10,10*3./4.))
                        ax = plt.subplot2grid((3,4), (1,0), colspan = 4, rowspan = 2)
                        axl1 = plt.subplot2grid((3,4), (0,0), colspan = 1)
                        axl2 = plt.subplot2grid((3,4), (0,1), colspan = 1)
                        axl3 = plt.subplot2grid((3,4), (0,2), colspan = 1)
                        axR = plt.subplot2grid((3,4), (0,3), colspan = 1)




                        for axx in [axl1, axl2, axl3, axR]:
                            axx.axis('off')

                    if line3 is None:
                        fig = plt.figure(figsize = (10,10*3./3.))
                        ax = plt.subplot2grid((3,3), (1,0), colspan = 4, rowspan = 2)
                        axl1 = plt.subplot2grid((3,3), (0,0), colspan = 1)
                        axl2 = plt.subplot2grid((3,3), (0,1), colspan = 1)
                        axR = plt.subplot2grid((3,3), (0,2), colspan = 1)

                        for axx in [axl1, axl2, axR]:
                            axx.axis('off')




                    imR = a[rat].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix]
                    eimR = a['e'+rat].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix]
                    imR[(imR<0) | (eimR/imR/log(10) > 0.4)] = nan


                    segm = detect_sources(imR, -99, npixels=5)

                    lbl_interest = array(segm.data)[14, 14]
                    if lbl_interest == 0:
                        small_box = array(segm.data)[12:16, 12:16].ravel() 
                        if len(small_box[small_box > 0]) > 0:                    
                            lbl_interest = min(small_box[small_box > 0])
                        else:
                            imR[:,:] = nan




                    imR[segm.data != lbl_interest] = nan


                    #return imR
                    #imR[] = nan
                    axR.imshow(log10(imR), cmap = 'viridis', interpolation = 'nearest')
                    axR.plot([(xmx_pix - xmn_pix)/2.], [(xmx_pix - xmn_pix)/2.], 'x', color = 'grey', markersize = 10)





                    im_line1 = a[line1].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix]
                    im_line2 = a[line2].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix]


                    def get_vmin(im):
                        im_rvl = im.ravel()
                        im_rvl = im_rvl[~isnan(im_rvl)]
                        im_rvl = sort(im_rvl)
                        vmin = im_rvl[int(0.10*len(im_rvl))]
                        return vmin

                    def get_vmax(im):
                        im_rvl = im.ravel()
                        im_rvl = im_rvl[~isnan(im_rvl)]
                        im_rvl = sort(im_rvl)
                        vmax = im_rvl[int(0.98*len(im_rvl))]
                        return vmax


                    axl1.imshow(im_line1, vmin = get_vmin(im_line1), vmax = get_vmax(im_line1), cmap = 'Greys_r')
                    axl2.imshow(im_line2, vmin = get_vmin(im_line2), vmax = get_vmax(im_line2), cmap = 'Greys_r')

                    if line3 is not None:
                        im_line3 = a[line3].data[xmn_pix:xmx_pix, xmn_pix:xmx_pix]
                        axl3.imshow(im_line3, vmin = get_vmin(im_line3),  vmax = get_vmax(im_line3), cmap = 'Greys_r')





                    et = eimR.ravel()/imR.ravel()/log(10)
                    dt = log10(imR.ravel())
                    rvt = rv.ravel()

                    gd = where(~isnan(dt))[0]
                    

                    ax.errorbar(rvt[gd], dt[gd], yerr = et[gd], fmt = 'o', markersize = 0.4, linewidth = 0.2, color = 'black')

                    x = rvt[gd]
                    y = dt[gd]
                    yerr = et[gd]

                    ndim, nwalkers = 2, 100
                    
                    import scipy.optimize as op
                    nll = lambda *args: -lnlike(*args)
                    result = op.minimize(nll, [0., 0.], args=(x, y, yerr))
                    m_ml, b_ml = result["x"]
                    pos = [result["x"] + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]

                    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, yerr))
                    sampler.run_mcmc(pos, 500)
                    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

                    xl = np.array([0, rmax_arc])
                    for m, b in samples[np.random.randint(len(samples), size=100)]:
                        ax.plot(xl, m*xl+b, color="k", alpha=0.1)
                    ax.errorbar(x, y, yerr=yerr, fmt=".k")
                    #ax = axes.ravel()[f]
                    #ax.imshow(a[rat].data)                
                    m_mcmc, b_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                                 zip(*np.percentile(samples, [16, 50, 84],
                                                                    axis=0)))

                    #print ('%.4f  %.4f  %.4f'%(m_mcmc[0], m_mcmc[1], m_mcmc[2]))
                    #ax.set_ylim(10**-1.5,10**1.5)
                    ax.set_ylim(-1.5,1.5)
                    ax.set_xlim(0, rmax_arc)
                    ax.set_ylabel(r'$\log$ %s'%rat, fontsize = 15)
                    ax.annotate('m = %.3f (-%.3f, +%.3f)'%(m_mcmc[0], m_mcmc[1], m_mcmc[2]), (0.05, 0.12), xycoords = 'axes fraction', fontsize = 12)
                    ax.annotate('b = %.3f (-%.3f, +%.3f)'%(b_mcmc[0], b_mcmc[1], b_mcmc[2]), (0.05, 0.05), xycoords = 'axes fraction', fontsize = 12)
                    ann_name = fl.split('/')[-1].strip('_metals.fits').replace('_', ' ')  
                    ax.annotate('%s'%ann_name, (0.7, 0.92), xycoords = 'axes fraction', fontsize = 15)

                    fit_results = {}

                    fit_results['r'] = x
                    fit_results['l'+rat] = y
                    fit_results['e_l' + rat] = yerr
                    fit_results['m_mcmc'] = [m_mcmc[0], m_mcmc[1], m_mcmc[2]]
                    fit_results['b_mcmc'] = [b_mcmc[0], b_mcmc[1], b_mcmc[2]]

                    save_name =  fl.split('/')[-1].strip('_metals.fits')
                    np.save('/Volumes/pegasus/clear/line_ratio_profiles/%s/%s_%s_profile.npy'%(rat, save_name, rat), fit_results)

                    #ax.set_xlim(0, 10)
                    ax.set_xlabel('distance from center (arcsec)', fontsize = 15)
                    #ax.set_yscale('log')
                    fig.tight_layout()
                    #fig.subplots_adjust(left = 0.15, top = 0.95, right = 0.95, wspace = 0.0)
                    pdf.savefig(fig)
                    plt.close()

#ratios = ['O32', 'R2', 'R3', 'R23']
#ratios = ['O32', 'R2', 'R3', 'R23']
#ratios = ['O32', 'R2', 'R3', 'R23']
#ratios = ['O32', 'R2', 'R3', 'R23']

#Parallel(n_jobs = 1)(delayed(make_fits)(rat = rat) for rat in ratios)

import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--ratio', default = 'R23')
args = vars(parser.parse_args())

rat = args['ratio']

if __name__ == '__main__':

    if rat == 'R23': line1, line2, line3 = 'OII', 'OIII', 'HB'
    if rat == 'R2': line1, line2, line3 = 'OII', 'HB', None
    if rat == 'R3': line1, line2, line3 = 'OIII', 'HB', None
    if rat == 'O32': line1, line2, line3 = 'OII', 'OIII', None



    make_fits(rat = rat, line1 = line1, line2 = line2, line3 = line3)





















