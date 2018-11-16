import time
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
import drizzlepac
import grizli
import glob
from grizli import utils
import importlib
from grizli.prep import process_direct_grism_visit
from hsaquery import query, overlaps
from grizli.pipeline import auto_script
from grizli.multifit import GroupFLT, MultiBeam, get_redshift_fit_defaults
import os
from grizli.pipeline import photoz
from astropy.table import Table
from matplotlib.colors import LogNorm
from IPython.display import Image
from numpy import *
import photutils
from astropy.cosmology import Planck15 as cosmo
from matplotlib.backends.backend_pdf import PdfPages
import glob
from glob import glob




def OH(O3, O2, eO3, eO2):
    
    O3_arr = np.random.normal(O3, eO3, 10000)
    O2_arr = np.random.normal(O2, eO2, 10000)
    
    OH_z_arr = 8.54 - 0.59 * O3_arr/O2_arr
    
    return np.mean(OH_z_arr), np.std(OH_z_arr)
    





def load_galfit(field, id_fit,  ra, dec, gfit_cat_gdn, gfit_cat_gds):
    if 'GN' in field: gfit_cat = gfit_cat_gdn
    elif 'GS' in field: gfit_cat = gfit_cat_gds
    else: return nan, nan

    diff_arc = np.sqrt((ra - gfit_cat[:,1])**2. + (dec - gfit_cat[:,2])**2.)*3600.

    if min(diff_arc) > 0.1:
        return nan, nan
    else:
        re    = gfit_cat[np.argmin(diff_arc), 6]
        tht   = gfit_cat[np.argmin(diff_arc), 12]
        ethr  = gfit_cat[np.argmin(diff_arc), 13]
        ab    = gfit_cat[np.argmin(diff_arc), 10]
        eab   = gfit_cat[np.argmin(diff_arc), 11]

    return tht, ab


mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 14

plt.ioff()
plt.close('all')


#objects = [('GN2', 18574), ('GN2', 21022), ('GN2', 14850)]

objects = [('GN1', 36108),
           ('GN1', 36854),
           ('GN1', 37031),
           ('GN1', 37208),
           ('GN1', 37636),
           ('GN1', 37525),
           ('GN1', 38091),
           ('GN2', 10512),
           ('GN2', 11339),
           ('GN2', 13453),
           ('GN2', 14319),
           ('GN2', 16758),
           ('GN2', 16867),
           ('GN2', 17898),
           ('GN2', 18150),
           ('GN2', 18197), 
           ('GN2', 18315),
           ('GN2', 18763),
           ('GN2', 18908),
           ('GN2', 21022),
           ('GN3', 26293)
           ]


objects = [('GS1', 43403)]



lines = ['OII', 'OIII']


mnx = 20
mxx = 60



gfit_cat_dir = '/Users/rsimons/Desktop/clear/Catalogs/galfit'
gfit_cat_gdn = np.loadtxt(gfit_cat_dir + '/gn_all_candels_wfc3_f105w_060mas_v0.8_galfit.cat')
gfit_cat_gds = np.loadtxt(gfit_cat_dir + '/gs_all_candels_ers_udf_f105w_v0.5_galfit.cat')


for o, obj in enumerate(objects):    
    try:
        field = obj[0]
        id_fit = obj[1]
        PATH_TO_PREP = '/Users/rsimons/Desktop/clear/for_hackday/Prep'
        #PATH_TO_PREP = glob('/Volumes/wd/clear/%s/*/Prep'%field)[0]
        fits_file = PATH_TO_PREP + '/{0}_{1:05d}.full.fits'.format(field, id_fit)
        fit_hdu = fits.open(fits_file)
        pix_scale = abs(fit_hdu['DSCI'].header['CD1_1'] * 60. * 60.)



        ra, dec = fit_hdu[0].header['ra'], fit_hdu[0].header['dec']

        tht, ab = load_galfit(field, id_fit, ra, dec, gfit_cat_gdn, gfit_cat_gds, )





        with PdfPages('/Users/rsimons/Dropbox/rcs_clear/z_radius_plots/%s_%i.pdf'%(field, id_fit)) as pdf:



            fig, axes = plt.subplots(len(lines)+1,2, figsize = (14, 5 * (len(lines)+1)))

            for ax in axes[:,0]:
                ax.set_xticklabels([])
                ax.set_yticklabels([])




            direct_im = fit_hdu['DSCI'].data[mnx:mxx, mnx:mxx]
            x1, y1 = photutils.centroid_2dg(direct_im)

            axes[0,0].plot(x1, y1, marker = 'x', color = 'Grey', markersize = 30)



            srt_rvl = np.sort(direct_im.ravel())
            vmn = srt_rvl[int(0.1*len(srt_rvl))]
            vmx = srt_rvl[int(0.9999*len(srt_rvl))]



            derr = 1./np.sqrt(fit_hdu['DWHT'].data)



            axes[0,0].imshow(direct_im)
            axes[0,0].set_title('direct')


            min_a = 0.01
            max_a = (mxx - mnx)/2. - 3.

            da = 1.

            a_arr = np.arange(min_a, max_a, da)
            try:
                for a_in in a_arr:
                    a_out = a_in + da
                    ea = photutils.EllipticalAnnulus(positions = (x1, y1), a_in = a_in, a_out = a_out, b_out = a_out/ab, theta = tht*pi/180.)
                    ap_sums = ea.do_photometry(direct_im, derr)
                    
                    
                    axes[0, 1].errorbar((a_in+da/2.) * pix_scale, ap_sums[0]/ea.area(), xerr = da/2. * pix_scale, yerr = ap_sums[1]/ea.area(), color = 'black', marker = 'o', ms = 10)
                    #for ax in axes[:, 0]:
                    ea.plot(ax = axes[0, 0])
            except:
                pass
            #ax2.set_yscale('log')
                   
                    
            aps = np.zeros((3, len(a_arr), len(lines)))
            clrs = ['blue', 'purple', 'red', 'green']
            for l, line in enumerate(lines):

                line_im = fit_hdu['LINE', line].data[mnx:mxx, mnx:mxx]
                line_err = 1/np.sqrt(fit_hdu['LINEWHT', line].data)[mnx:mxx, mnx:mxx]
                line_im[~isfinite(line_im)] = 0.
                line_err[~isfinite(line_err)] = 1000.



                srt_rvl = np.sort(line_im.ravel())
                vmn = srt_rvl[int(0.1*len(srt_rvl))]
                vmx = srt_rvl[int(0.9*len(srt_rvl))]
                
                axes[l+1, 0].plot(x1, y1, marker = 'x', color = 'Grey', markersize = 30)

                axes[l+1, 0].imshow(line_im, vmin = vmn, vmax = vmx)
                try:
                    for a, a_in in enumerate(a_arr):
                        a_out = a_in + da
                        ea = photutils.EllipticalAnnulus(positions = (x1, y1), a_in = a_in, a_out = a_out, b_out = a_out/ab, theta = tht*pi/180.)
                        ap_sums = ea.do_photometry(line_im, line_err)
                        
                        aps[0, a, l] = (a_in + da/2.) * pix_scale
                        aps[1, a, l] = ap_sums[0]/ea.area()
                        aps[2, a, l] = ap_sums[1]/ea.area()
                except:
                    pass

                axes[l+1, 1].annotate(line, (0.75, 0.85), xycoords = 'axes fraction', color = 'black', fontweight = 'bold', fontsize = 40)
                axes[l+1, 1].axhline(y = 0.0, color = 'black')

                axes[l+1, 1].plot(aps[0, :, l], aps[1, :, l], marker = 'o', color = 'black', linewidth = 0, markersize = 0.0)
                ylm = axes[l+1, 1].get_ylim()
                axes[l+1, 1].errorbar(aps[0, :, l], aps[1, :, l], xerr = da/2. * pix_scale, yerr = aps[2, :, l], color = 'black', fmt = 'o', ms = 10)
                ylm2 = axes[l+1, 1].get_ylim()
                axes[l+1, 1].set_ylim(max(ylm[0]*2.0, ylm2[0]), min(ylm[1]*2.0, ylm2[1]))

            fig.subplots_adjust(wspace = 0.0, hspace = 0.2)
            pdf.savefig()



                #if ap_sums[1]/ea.area() < 3.0:
                #axes[l+1, 1].errorbar((a_in+da/2.) * pix_scale, ap_sums[0]/ea.area(), xerr = da/2. * pix_scale, yerr = ap_sums[1]/ea.area(), color = 'black', marker = 'o', ms = 10)
                #else:
                #    axes[l+1, 1].errorbar((a_in+da/2.) * pix_scale, ap_sums[0]/ea.area(), xerr = da/2. * pix_scale, yerr = 0, color = 'grey', marker = 'x', ms = 10, alpha = 0.3)





            fig3, ax3 = plt.subplots(1,1, figsize = (15, 6))
            z = fit_hdu[1].header['Z50']
                


            to_fit = []
            for a in np.arange(len(a_arr)):
                OH_z, eOH_z = OH(O3 = aps[1, a, 1], O2 = aps[1, a, 0], eO3 = aps[2, a, 1], eO2 = aps[2, a, 0])


                if eOH_z/OH_z > 1/3.:
                    alp = 0.1
                    clr = 'grey'

                elif eOH_z/OH_z < 0.1:
                    alp = 1.0
                    clr = 'darkblue'
                else:
                    alp = 1.0
                    clr = 'black'            

                to_fit.append([aps[0,a,0], OH_z, eOH_z])
                ax3.errorbar(aps[0,a,0], OH_z, yerr = eOH_z,color = clr, fmt = 'o', alpha = alp, markersize = 10, zorder = 2)

                ax3.set_ylim(7., 9.5)
                ax3.set_ylabel(r'12 + log(O/H)', fontsize = 20)
                ax3.set_xlabel(r'Semi-major axis radius [arcsec]', fontsize = 20)
                ax3_t = ax3.twiny()


                kpc_ticks = np.arange(1, 10)

                arc_ticks = np.array([cosmo.arcsec_per_kpc_proper(z).value * k for k in kpc_ticks])






                ax3_t.set_xticks(arc_ticks)
                ax3_t.set_xticklabels(np.array(['%i'%k for k in kpc_ticks]))
                ax3_t.set_xlabel('Semi-major axis radius [kpc]', fontsize = 20)

            to_fit = np.array(to_fit)


            clrs = ['blue', 'darkblue', 'black']

            for b, ar in enumerate(np.array([10., 5., 3.])):  
                g = where(to_fit[:,1]/to_fit[:,2] > ar)[0]
                try:
                    p, V = np.polyfit(to_fit[g,0], to_fit[g,1], deg = 1., w = 1./to_fit[g,2], cov = True)

                    x = np.linspace(0, 2, 1000)

                    draws = np.random.multivariate_normal(p, V, size = 100)
                    ax3.annotate(r'$\Delta$(O/H)/$\Delta$r = %.3f $\pm$ %.3f dex kpc$^{-1}$'%(p[0]*cosmo.arcsec_per_kpc_proper(z).value, np.sqrt(V[0,0])*cosmo.arcsec_per_kpc_proper(z).value), xy = (0.03, 0.18 - b*0.07),color = clrs[b], fontsize = 20, xycoords = 'axes fraction')


                    for d in draws:
                        ax3.plot(x, x*d[0] + d[1], color = clrs[b], alpha = 0.1)

                except:
                    pass



            pdf.savefig()
    except:
        pass
























