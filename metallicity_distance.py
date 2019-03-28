import time
import os
import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
import glob
import importlib
import os
from astropy.table import Table
from matplotlib.colors import LogNorm
from IPython.display import Image
from numpy import *
import photutils
from astropy.cosmology import Planck15 as cosmo
from matplotlib.backends.backend_pdf import PdfPages
import glob
from glob import glob
from astropy.convolution import Gaussian2DKernel, convolve_fft, Box2DKernel
from scipy.interpolate import interp1d
import joblib
from joblib import Parallel, delayed
from astropy.stats import sigma_clip
import emcee
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 14
plt.ioff()




def O32_OH(O3, O2, eO3, eO2):
    O32_arr = zeros(O3.shape) * nan
    O32_e_arr = zeros(O3.shape) * nan

    OH_arr = zeros(O3.shape) * nan
    OH_e_arr = zeros(O3.shape) * nan

    for i in arange(O3.shape[0]):
        for j in arange(O3.shape[1]):
            O2_arr_temp = np.random.normal(O2[i,j], eO2[i,j], 200)
            O3_arr_temp = np.random.normal(O3[i,j], eO3[i,j], 200)
            O32_arr[i,j] = O3[i,j]/O2[i,j]
            O32_e_arr[i,j] = (O3[i,j]/O2[i,j]) * sqrt((eO3[i,j]/O3[i,j])**2. + (eO2[i,j]/O2[i,j])**2.)
            OH_arr[i,j]   = np.nanmean(8.54 - 0.59 * O3_arr_temp/O2_arr_temp)
            OH_e_arr[i,j] = np.nanstd(8.54 - 0.59 * O3_arr_temp/O2_arr_temp)
    
    return O32_arr, O32_e_arr, OH_arr, OH_e_arr#np.std(OH_z_arr)


def R23_OH(O3, O2, Hb, eO3, eO2, eHb):
    O32_arr = zeros(O3.shape) * nan
    O32_e_arr = zeros(O3.shape) * nan

    OH_arr = zeros(O3.shape) * nan
    OH_e_arr = zeros(O3.shape) * nan

    for i in arange(O3.shape[0]):
        for j in arange(O3.shape[1]):
            O2_arr_temp = np.random.normal(O2[i,j], eO2[i,j], 200)
            O3_arr_temp = np.random.normal(O3[i,j], eO3[i,j], 200)
            O32_arr[i,j] = O3[i,j]/O2[i,j]
            O32_e_arr[i,j] = (O3[i,j]/O2[i,j]) * sqrt((eO3[i,j]/O3[i,j])**2. + (eO2[i,j]/O2[i,j])**2.)
            OH_arr[i,j]   = np.nanmean(8.54 - 0.59 * O3_arr_temp/O2_arr_temp)
            OH_e_arr[i,j] = np.nanstd(8.54 - 0.59 * O3_arr_temp/O2_arr_temp)
    
    return O32_arr, O32_e_arr, OH_arr, OH_e_arr#np.std(OH_z_arr)






def O32_OH_profile(r, O3, O2, eO3, eO2, r_min, r_max, dr):
    rrs = arange(r_min + dr/2., r_max + dr/2., dr)
    O32    = zeros((len(rrs), 2)) * nan
    O3_bin = zeros((len(rrs), 2)) * nan
    O2_bin = zeros((len(rrs), 2)) * nan
    OH     = zeros((len(rrs), 2)) * nan  
    for i, rr in enumerate(rrs):
        gd = where((r > (rr - dr/2.)) & (r < (rr + dr/2.)))[0]
        O3_rr = sum(O3[gd]/eO3[gd]**2.)/sum(1./eO3[gd]**2.)
        O2_rr = sum(O2[gd]/eO2[gd]**2.)/sum(1./eO2[gd]**2.)
        eO3_rr = sqrt(sum(eO3[gd]**2.))/len(gd)
        eO2_rr = sqrt(sum(eO2[gd]**2.))/len(gd)

        O3_bin[i,0] = O3_rr
        O2_bin[i,0] = O2_rr

        O3_bin[i,1] = eO3_rr
        O2_bin[i,1] = eO2_rr


        O2s = array([])
        O3s = array([])
        for g in gd:
            O3s = concatenate((O3s, np.random.normal(O3[g], eO3[g], 1000)))
            O2s = concatenate((O2s, np.random.normal(O2[g], eO2[g], 1000)))




        #O32.append(nanmean(O3s/O2s))
        #eO32.append(nanstd(O3s/O2s))
        O3s = np.random.normal(O3_rr, eO3_rr, 1000)
        O2s = np.random.normal(O2_rr, eO2_rr, 1000)


        O32_rr = sigma_clip(O3s/O2s, sigma = 2)
        O32_rr = O32_rr.data[O32_rr.mask == False]


        O32_rr_mean = nanmean(O32_rr)
        O32_rr_std = nanstd(O32_rr)

        O32[i,0] = O32_rr_mean
        O32[i,1] = O32_rr_std


        O32_bs = np.random.normal(O32_rr_mean, O32_rr_std, 1000)
        OH[i,0] = np.nanmean(8.54 - 0.59 * O32_bs)
        OH[i,1] = np.nanstd(8.54 - 0.59 * O32_bs)




    return rrs, O3_bin, O2_bin, O32, OH



def load_galfit(field, id_fit,  ra, dec, gfit_cat_gdn, gfit_cat_gds):
    if 'GN' in field: gfit_cat = gfit_cat_gdn
    elif 'GS' in field: 
        gfit_cat = gfit_cat_gds
    else: return nan, nan

    diff_arc = np.sqrt((ra - gfit_cat[:,1])**2. + (dec - gfit_cat[:,2])**2.)*3600.

    if min(diff_arc) > 0.3:
        return nan, nan
    else:
        re    = gfit_cat[np.argmin(diff_arc), 6]
        tht   = gfit_cat[np.argmin(diff_arc), 12]
        ethr  = gfit_cat[np.argmin(diff_arc), 13]
        ab    = gfit_cat[np.argmin(diff_arc), 10]
        eab   = gfit_cat[np.argmin(diff_arc), 11]

    return tht, ab






def metallicity_distance(field, id_fit, gfit_cat_gdn, gfit_cat_gds, rmx = 1.0):
    fits_file = glob(PATH_TO_GE + '/%s/j*/Prep/{0}_{1:05d}.full.fits'.format(field, field, id_fit))[0]
    #fits_file = PATH_TO_PREP + '/{0}_{1:05d}.full.fits'.format(field, id_fit)
    if os.path.isfile(fits_file):

        fit_hdu = fits.open(fits_file)
        pix_scale = abs(fit_hdu['DSCI'].header['CD1_1'] * 60. * 60.)
        ra, dec = fit_hdu[0].header['ra'], fit_hdu[0].header['dec']
        tht, ab = load_galfit(field, id_fit, ra, dec, gfit_cat_gdn, gfit_cat_gds)
        tht_rad = tht*pi/180.
        try:
            if (tht != -999) & (~isnan(tht)):
                with PdfPages('/Users/rsimons/Dropbox/rcs_clear/z_radius_plots/%s_%i.pdf'%(field, id_fit)) as pdf:
                    if prt: print('/Users/rsimons/Dropbox/rcs_clear/z_radius_plots/%s_%i.pdf'%(field, id_fit))
                    fig, axes = plt.subplots(len(lines)+2,2, figsize = (14, 5 * (len(lines)+2)))

                    for ax in axes[:,0]:
                        ax.set_xticklabels([])
                        ax.set_yticklabels([])
                    for ax in axes[:,1]:
                        ax.axhline(y = 0.0, color = 'grey', alpha = 0.3)




                    direct_im = fit_hdu['DSCI'].data
                    seg_im = fit_hdu['SEG'].data
                    direct_im_temp = direct_im.copy()
                    direct_im_temp[seg_im != seg_im[int(shape(seg_im)[0]/2.), int(shape(seg_im)[0]/2.)]] = nan
                    x1, y1 = photutils.centroid_2dg(direct_im_temp)

                    x = (np.arange(0, shape(direct_im)[0]) - x1 + 0.5) * pix_scale
                    y = (np.arange(0, shape(direct_im)[1]) - y1 + 0.5) * pix_scale

                    xv, yv = np.meshgrid(x, y)
                    r = sqrt(xv**2. + yv**2.)

                    axes[0,0].plot(x1, y1, marker = 'x', color = 'Grey', markersize = 30)



                    srt_rvl = np.sort(direct_im.ravel())
                    vmn = srt_rvl[int(0.1*len(srt_rvl))]
                    vmx = srt_rvl[int(0.9999*len(srt_rvl))]



                    derr = 1./np.sqrt(fit_hdu['DWHT'].data)
                    ymn, ymx =  y1 - 1.25/pix_scale, y1 + 1.25/pix_scale
                    xmn, xmx =  x1 - 1.25/pix_scale, x1 + 1.25/pix_scale


                    ymn_i, ymx_i = int(ymn), int(ymx)
                    xmn_i, xmx_i = int(xmn), int(xmx)



                    srt_rvl = np.sort(direct_im[xmn_i:xmx_i, ymn_i:ymx_i].ravel())
                    vmn = srt_rvl[int(0.05*len(srt_rvl))]
                    vmx = srt_rvl[int(0.975*len(srt_rvl))]
                    
                    axes[0,0].imshow(direct_im,vmin = vmn, vmax = vmx, cmap = 'Greys_r')
                    axes[0,0].set_title('direct', fontsize = 40)
                    axes[0,1].axis('off')
                    kern = Box2DKernel(2)

                    for l, line in enumerate(lines):
                        line_im = fit_hdu['LINE', line].data
                        line_err = 1/np.sqrt(fit_hdu['LINEWHT', line].data)
                        line_im = convolve_fft(line_im, kern)
                        line_im[~isfinite(line_err)] = 0.
                        line_err[~isfinite(line_err)] = 0.

                        srt_rvl = np.sort(line_im.ravel())
                        vmn = srt_rvl[int(0.1*len(srt_rvl))]
                        vmx = srt_rvl[int(0.99*len(srt_rvl))]
                        

                        axes[l+1, 0].imshow(line_im, vmin = vmn, vmax = vmx)


                        axes[l+1,1].set_xlabel('distance from center (arcsec)', fontsize = 20)
                        axes[l+1,1].set_ylabel('surface brightness', fontsize = 20)

                        gd = where((line_im !=0.) & (r < rmx))

                        axes[l+1, 1].errorbar(r[gd].ravel(), line_im[gd].ravel(), yerr = line_err[gd].ravel(),  color = 'black', fmt = 'o', ms = 4, linewidth = 0.4, alpha = 0.15)
                        axes[l+1, 1].annotate(line, (0.75, 0.85), xycoords = 'axes fraction', color = 'black', fontweight = 'bold', fontsize = 40)

                    O2_im = convolve_fft(fit_hdu['LINE', 'OII'].data, kern)
                    O3_im = convolve_fft(fit_hdu['LINE', 'OIII'].data, kern)

                    O2_err = 1/np.sqrt(fit_hdu['LINEWHT', 'OII'].data)
                    O3_err = 1/np.sqrt(fit_hdu['LINEWHT', 'OIII'].data)                
                    #line_im[~isfinite(line_err)] = 0.
                    #line_err[~isfinite(line_err)] = 0.



                    O32_im, O32_eim, OH_z, eOH_z = O32_OH(O3 = O3_im, O2 = O2_im, eO3 = O3_err, eO2 = O3_err)
                    dr = 0.1
                    rrs, O3_bin, O2_bin,O32, OH = O32_OH_profile(r = r.ravel(), O3 = O3_im.ravel(), O2 = O2_im.ravel(), eO3 = O3_err.ravel(), eO2 = O3_err.ravel(), r_min = min(r.ravel()), r_max = rmx, dr = dr)


                    axes[1, 1].errorbar(rrs, O2_bin[:,0], yerr = O2_bin[:,1],  color = 'blue', fmt = 'o', ms = 4, linewidth = 0.4)
                    axes[2, 1].errorbar(rrs, O3_bin[:,0], yerr = O3_bin[:,1],  color = 'blue', fmt = 'o', ms = 4, linewidth = 0.4)


                    gd = where((O2_im!=0) & (O3_im!=0))# & (O32_im/O32_eim > 1/3.))

                    
                    #O32_im = O3_im/O2_im

                    #O32_im[O3_im < 0] = nan
                    #O32_im[O2_im < 0] = nan
                    
                    srt_rvl = np.sort(O32_im.ravel())

                    vmn = srt_rvl[int(0.1*len(srt_rvl))]
                    vmx = srt_rvl[int(0.99*len(srt_rvl))]
            
                    vmn = 0.
                    vmx = 3.
                    axes[len(lines)+1, 0].imshow(O32_im, vmin = vmn, vmax = vmx)


                    #axes[len(lines)+1, 1].errorbar(r[gd].ravel(), O32_arr[gd].ravel(), yerr = O32_e_arr[gd].ravel(),  color = 'black', fmt = 'o', ms = 4, linewidth = 0.4, alpha = 0.3)
                    #axes[len(lines)+1, 1].errorbar(r[gd].ravel(), O32_im[gd].ravel(), yerr = O32_eim[gd].ravel(),  color = 'black', fmt = 'o', ms = 4, linewidth = 0.4, alpha = 0.3)
                    
                    axes[len(lines)+1, 1].errorbar(rrs, O32[:,0], yerr = O32[:,1], color = 'blue', fmt = 'o', ms = 4, linewidth = 2., alpha = 1.0)
                    



                    axes[len(lines)+1, 1].annotate('OIII/OII', (0.58, 0.85), xycoords = 'axes fraction', color = 'black', fontweight = 'bold', fontsize = 40)
                    axes[len(lines)+1, 1].set_ylim(0.,5.)
                    axes[len(lines)+1,1].set_ylabel('OIII/OII', fontsize = 20)

                    for ax in axes[:,0]:
                        ax.plot(x1, y1, marker = 'x', color = 'Grey', markersize = 30)
                        ax.set_xlim(xmn, xmx)
                        ax.set_ylim(ymn, ymx)
                    for ax in axes[:,1]:
                        ax.set_xlim(0,rmx)
                        ax.set_xlabel('distance from center (arcsec)', fontsize = 20)

                    #axes[len(lines)+1, 0].set_xlim(x1 - 1.25/pix_scale, x1 + 1.25/pix_scale, )                
                    #axes[len(lines)+1, 0].set_ylim(y1 - 1.25/pix_scale, y1 + 1.25/pix_scale, )                





                    pdf.savefig()



                    fig3, ax3 = plt.subplots(1,1, figsize = (12, 7))

                    #gd = where((O3_im !=0.) & (O2_im !=0.) & (eOH_z > 1.) & (r < rmx))
                    #ax3.errorbar(r[gd].ravel(), OH_z[gd].ravel(), yerr = eOH_z[gd].ravel(),color = 'grey', fmt = 'o', alpha = 0.2, markersize = 4, linewidth = 0.05,zorder = 2)


                    #gd = where((O3_im != 0.) & (O2_im != 0.) & (O2_im/O2_err > 1.) & (O3_im/O3_err > 1.) & (eOH_z < 1.))
                    #ax3.errorbar(r[gd].ravel(), OH_z[gd].ravel(), yerr = eOH_z[gd].ravel(),color = 'blue', fmt = 'o', alpha = 1.0, markersize = 8, linewidth = 1.,zorder = 2)
                    gd = where(O32[:,0] > O32[:,1])[0]
                    ax3.errorbar(rrs, OH[:,0], yerr = OH[:,1],color = 'grey', fmt = 'x', alpha = 0.3, markersize = 8, linewidth = 0.3,zorder = 2)
                    ax3.errorbar(rrs[gd], OH[gd,0], yerr = OH[gd,1],color = 'blue', fmt = 'o', alpha = 1.0, markersize = 8, linewidth = 1.,zorder = 2)
                    z = fit_hdu[1].header['Z50']
                    if len(gd) > 4:
                        p, V = np.polyfit(rrs[gd], OH[gd,0], deg = 1., w = 1./OH[gd,1], cov = True)
                        plot_x = max(rrs[gd] + dr/2.)
                        x = np.linspace(0, plot_x, 200)
                        x_rest = np.linspace(plot_x, rmx, 200)
                        draws = np.random.multivariate_normal(p, V, size = 100)

                        for d in draws:
                            ax3.plot(x, x*d[0] + d[1], color = 'blue', alpha = 0.1)
                            ax3.plot(x_rest, x_rest*d[0] + d[1], color = 'blue',alpha = 0.03, linestyle = '--')

                        cat.write('%s   %.5i   %.3f   %.3f   %.3f   %.3f\n'%(field, id_fit, p[0], np.sqrt(V[0,0]), p[0]*cosmo.arcsec_per_kpc_proper(z).value, np.sqrt(V[0,0])*cosmo.arcsec_per_kpc_proper(z).value))

                        ax3.annotate(r'$\Delta$(O/H)/$\Delta$r = %.3f $\pm$ %.3f dex kpc$^{-1}$'%(p[0]*cosmo.arcsec_per_kpc_proper(z).value, np.sqrt(V[0,0])*cosmo.arcsec_per_kpc_proper(z).value), xy = (0.03, 0.9), color = 'blue', fontsize = 20, xycoords = 'axes fraction')

                    ax3.set_ylim(5, 10)

                    ax3.set_ylabel(r'12 + log(O/H)', fontsize = 20)
                    ax3.set_xlabel(r'distance from center (arcsec)', fontsize = 20)
                    kpc_ticks = np.arange(0, 25, 1)

                    arc_ticks = np.array([cosmo.arcsec_per_kpc_proper(z).value * k for k in kpc_ticks])
                    ax3_t = ax3.twiny()

                    ax3_t.set_xticks(arc_ticks)
                    for ax in [ax3, ax3_t]:
                        ax.set_xlim(0,rmx)



                    ax3_t.set_xticklabels(np.array(['%i'%k for k in kpc_ticks]))
                    ax3_t.set_xlabel('Semi-major axis radius [kpc]', fontsize = 20, labelpad = 12)

                    pdf.savefig()

                    return
            else:
                if prt: print 'bad'
                return
        except:
            if prt: print 'exception'
            return
    else:
        if prt: print '%s does not exist'%fits_file
        return
    plt.ioff()
    plt.close('all')




if __name__ == '__main__':
    plt.ioff()
    global PATH_TO_PREP, PATH_TO_GE
    #PATH_TO_PREP = '/Users/rsimons/Desktop/clear/grizli_v2.1/all_full'    
    PATH_TO_GE   = '/user/rsimons/grizli_extractions'

    lines = ['OII', 'OIII']

    #gfit_cat_dir = '/Users/rsimons/Desktop/clear/Catalogs/galfit'
    gfit_cat_dir = '/user/rsimons/grizli_extractions/Catalogs/galfit'
    gfit_cat_gdn = np.loadtxt(gfit_cat_dir + '/gn_all_candels_wfc3_f105w_060mas_v0.8_galfit.cat')
    gfit_cat_gds = np.loadtxt(gfit_cat_dir + '/gs_all_candels_ers_udf_f105w_v0.5_galfit.cat')


    diagnostic = 'O32'

    #objects = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/z_r_sample.cat', dtype = 'str')
    objects = np.loadtxt('/user/rsimons/grizli_extractions/Catalogs/sample_cats/%s_sample.cat'%diagnostic, dtype = 'str')
    objects = objects
    #objects = [('GN7', 17293)]


    #cat = open('/Users/rsimons/Desktop/clear/Catalogs/z_r.cat', 'w+')
    cat = open('/user/rsimons/grizli_extractions/Catalogs/z_r_%s.cat'%diagnostic, 'w+')

    rmx = 1.0

    for o, obj in enumerate(objects):    
        metallicity_distance(field = obj[0], id_fit = int(obj[1]), gfit_cat_gdn = gfit_cat_gdn, gfit_cat_gds = gfit_cat_gds,  rmx = rmx)
    #Parallel(n_jobs = 1, backend = 'threading')(delayed(metallicity_distance)(field = obj[0], id_fit = int(obj[1]), gfit_cat_gdn = gfit_cat_gdn, gfit_cat_gds = gfit_cat_gds, rmx = rmx, cat = cat_f) for o, obj in enumerate(objects))

    cat.close()

























