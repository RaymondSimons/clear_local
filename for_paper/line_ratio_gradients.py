import astropy
from astropy.io import fits
import glob
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import joblib
from joblib import Parallel, delayed
import photutils
from photutils import detect_sources
from numpy import *
from astropy.visualization import AsinhStretch

plt.ioff()
plt.close('all')


def make_figure(fl):

    master_hdulist = []
    prihdr = fits.Header()
    prihdr['COMMENT'] = "Storing the line ratio maps in this FITS file."
    prihdu = fits.PrimaryHDU(header=prihdr)    
    master_hdulist.append(prihdu)


    di =  int(fl.split('/')[-1].split('_')[1])

    metal_fits = fits.open(fl)
    master_hdulist.append(metal_fits['DSCI'])
    master_hdulist.append(metal_fits['DWHT'])
    master_hdulist.append(metal_fits['SEG'])

    seg = metal_fits['SEG'].data[20:60, 20:60]

    hdu_names = [hdu.name for hdu in metal_fits]
    fig, axes = plt.subplots(13,4, figsize = (12,40))



    line_ratio_fig_dir = '/Users/rsimons/Desktop/clear/figures/line_ratio_maps'
    for ax in axes.ravel(): 
        ax.set_xticks([])
        ax.set_yticks([])

        ax.annotate('x', (0.5, 0.5), va = 'center', ha = 'center', xycoords = 'axes fraction', fontsize =20, color = 'grey', fontweight = 'bold')
        ax.imshow(np.zeros((40,40)), cmap = 'Greys_r')
        ax.imshow(np.zeros((40,40)), cmap = 'Greys_r')

    seg_masked = seg.copy()
    seg_masked = seg_masked.astype('float')
    seg_masked[seg !=di] = np.nan


    seg_maskzeros = seg.copy()
    seg_maskzeros = seg_maskzeros.astype('float')
    seg_maskzeros[seg ==0.] = np.nan


    dir_im = metal_fits['DSCI'].data[20:60, 20:60]

    dir_im -= nanmin(dir_im)
    dir_im/= nanmax(dir_im)
    print (dir_im.min(), dir_im.max())

    dir_im_ravel = sort(dir_im.ravel())



    vmin = dir_im_ravel[int(0.1*len(dir_im_ravel))]
    vmax = dir_im_ravel[int(0.999*len(dir_im_ravel))]



    dir_im_asinh = np.log10(dir_im)#AsinhStretch(values = dir_im, a = vmin)

    axes[0, 0].imshow(dir_im_asinh, vmin = log10(vmin), vmax = 0., cmap = 'Greys_r')
    #axes[0, 0].imshow(dir_im, vmin = vmin, vmax = vmax, cmap = 'Greys_r')

    axes[0, 1].imshow(seg_masked, vmin = di - 5, vmax = di+ 5)
    axes[0, 2].imshow(seg_maskzeros, vmin = di - 5, vmax = di+ 5)

    axes[0, 0].annotate('F105W', (0.05, 0.95), va = 'top', ha = 'left', xycoords = 'axes fraction', fontsize =20, color = 'white', fontweight = 'bold')
    axes[0, 1].annotate(r'Seg$_{F105W}$', (0.05, 0.95), va = 'top', ha = 'left', xycoords = 'axes fraction', fontsize =20, color = 'white', fontweight = 'bold')
    axes[0, 2].annotate(r'Seg$_{F105W, all}$', (0.05, 0.95), va = 'top', ha = 'left', xycoords = 'axes fraction', fontsize =20, color = 'white', fontweight = 'bold')


    lines = ['OII', 'OIII', 'HB']
    line_ratios = ['O32', 'R2', 'R3', 'R23']
    typs = ['', '_S', '_EC', '_S_EC']
    for t, typ in enumerate(typs):
        for l, line in enumerate(lines):
            if line == 'OII': line_str = '[OII]'


            if line == 'OIII':  line_str = '[OIII]'

            if line == 'HB': line_str = r'H$\beta$'

            axes[t*3 + 1,l].annotate(line_str, (0.05, 0.95), va = 'top', ha = 'left', xycoords = 'axes fraction', fontsize =20, color = 'white', fontweight = 'bold')
            if '%s%s'%(line, typ) not in hdu_names: 
                axes[t*3 + 1,l].imshow(np.zeros((40,40)), cmap = 'Greys_r')
                continue

            if line == 'OII': 
                line_str = '[OII]'
                O2, eO2 = metal_fits['OII%s'%typ].data,  metal_fits['eOII%s'%typ].data

            if line == 'OIII': 
                line_str = '[OIII]'
                O3, eO3 = metal_fits['OIII%s'%typ].data, metal_fits['eOIII%s'%typ].data


            if line == 'HB': 
                line_str = r'H$\beta$'
                HB, eHB = metal_fits['HB%s'%typ].data, metal_fits['eHB%s'%typ].data

            line_map = metal_fits['%s%s'%(line, typ)].data
            eline_map = metal_fits['E%s%s'%(line, typ)].data

            axes[t*3 + 1,l].imshow(line_map)



        
        cmap = plt.cm.viridis
        cmap.set_bad('k')

        def mask_LR_maps(imR, eimR, seg, di, SN_thresh = 1.):
            gd = where((imR/eimR < SN_thresh) | (imR < 0.) | (seg != di))
            imR[gd] = np.nan
            eimR[gd] = np.nan

            segm = detect_sources(imR, -99, npixels=5)
            if segm is not None:
                lbl_interest = array(segm.data)[20, 20]
                if lbl_interest == 0:
                    small_box = array(segm.data)[16:24, 16:24].ravel() 
                    if len(small_box[small_box > 0]) > 0:                    
                        lbl_interest = min(small_box[small_box > 0])
                    else:
                        imR[:,:] = nan
                        eimR[:,:] = nan


                imR[segm.data != lbl_interest] = nan
                eimR[segm.data != lbl_interest] = nan

            return imR, eimR

        for lr, line_ratio in enumerate(line_ratios):
            vmn = -1
            vmx = 1.
            axes[t*3 + 2,lr].annotate(line_ratio, (0.05, 0.95), va = 'top', ha = 'left', xycoords = 'axes fraction', fontsize =20, color = 'white', fontweight = 'bold')
            axes[t*3 + 3,lr].annotate(line_ratio, (0.05, 0.95), va = 'top', ha = 'left', xycoords = 'axes fraction', fontsize =20, color = 'white', fontweight = 'bold')

            if (line_ratio == 'O32') & ('OII%s'%typ in hdu_names) & ('OIII%s'%typ in hdu_names):
                R  = O3/O2
                eR = R * np.sqrt((eO2/O2)**2. + (eO3/O3)**2.)
            elif (line_ratio == 'R2') & ('OII%s'%typ in hdu_names) & ('HB%s'%typ in hdu_names):
                R = O2/HB
                eR = R * np.sqrt((eO2/O2)**2. + (eHB/HB)**2.)

            elif (line_ratio == 'R3') & ('OIII%s'%typ in hdu_names) & ('HB%s'%typ in hdu_names):
                R = O3/HB
                eR = R * np.sqrt((eO3/O3)**2. + (eHB/HB)**2.)
     
            elif (line_ratio == 'R23') & ('OII%s'%typ in hdu_names) & ('OIII%s'%typ in hdu_names) & ('HB%s'%typ in hdu_names):
                R = (O2 + O3)/HB
                etop = np.sqrt(eO3**2. + eO2**2.)
                eR = R * np.sqrt((etop/(O3 + O2))**2. + (eHB/HB)**2.)
            else:
                R = np.nan*zeros((40,40))
                eR = np.nan*zeros((40,40))

            axes[t*3 + 2,lr].imshow(np.log10(R),  vmin = vmn, vmax = vmx, cmap = cmap)
            if len(where(~np.isnan(R.ravel()))[0]) > 0:  R_mask, eR_mask = mask_LR_maps(R, eR, seg, di)
            else:                                        R_mask, eR_mask = R.copy(), R.copy()

            axes[t*3 + 3,lr].imshow(np.log10(R),  vmin = vmn, vmax = vmx, cmap = cmap)


            master_hdulist.append(fits.ImageHDU(data = stack((R_mask, eR_mask)), header = prihdr, name = '%s%s'%(line_ratio, typ)))


    axes[0,3].axis('off')
    axes[3,3].axis('off')
    axes[6,3].axis('off')
    axes[9,3].axis('off')

    fig.tight_layout()
    fig.savefig(line_ratio_fig_dir + '/' + fl.split('/')[-1].replace('_maps.fits', '_line_ratios.pdf'))
    plt.close(fig)



    out_dir = '/Users/rsimons/Desktop/clear/maps/line_ratio_maps'
    fits_name = out_dir + '/' +  fl.split('/')[-1].replace('_maps.fits', '_line_ratios.fits')
    print ('\tSaving to ' + fits_name)
    thdulist = fits.HDUList(master_hdulist)
    thdulist.writeto(fits_name, overwrite = True)






orig_maps_dir = '/Users/rsimons/Desktop/clear/maps/orig_maps'
fls = glob(orig_maps_dir + '/*_maps.fits')




Parallel(n_jobs = -1)(delayed(make_figure)(fl) for fl in fls)

#for fl in fls:     make_figure(fl)

plt.close('all')


























