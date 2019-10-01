import astropy
from astropy.io import fits, ascii
import glob
from glob import glob
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy.ma as ma
import photutils
from photutils import detect_sources
plt.ioff()
plt.close('all')


figdir = '/Users/rsimons/Desktop/clear/figures/izi_metal_maps/good_new'
indir = '/Users/rsimons/Dropbox/rcs_clear/data/izi_metal_maps'

izi_cat = np.loadtxt('/Users/rsimons/Dropbox/clear/Catalogs/good_izi.cat', dtype = 'str')

flds = izi_cat[:,0]
dis = izi_cat[:,1]


fls = glob(indir + '/*metals.fits')


fit_types = array(['', '_S', '_EC', '_S_EC'])
zmax = 0.5

for f, (fld, di) in enumerate(zip(flds, dis)):
    np.random.seed(1)
    di = int(di)
    #if di != 19913: continue
    #if f > 10: break
    lines = [('OII', 'oii3726;oii3729', 3727.),
             ('OIII', 'oiii4959;oiii5007', 5007.),
             ('Hb', 'hbeta', 4863.),
             ('Ha', 'nii6548;halpha;nii6584', 6563.),
             ('SII', 'sii6717;sii6731', 6725.)
            ]

    fl = indir + '/{}_{}_metals.fits'.format(fld, di)
    #di = int(fl.split('/')[-1].split('_')[1])

    mm_fits = fits.open(fl)
    
    line_cnt = 0
    for l, (line, izi_line, wave) in enumerate(lines):
        try: im = mm_fits['%s'%line].data
        except: continue
        line_cnt+=1



    nzcols = 7
    nrows = line_cnt + nzcols
    ncols = 4

    fig = plt.figure(figsize = (8., 2.*nrows))
    res = {}

    for ft, fit_type in enumerate(fit_types):
        ct = 0
        dir_im = mm_fits['DSCI'].data[20:60, 20:60]
        dir_seg = mm_fits['SEG'].data[20:60, 20:60]


        Zax1 = plt.subplot2grid((nrows, ncols), (0, ft))
        Zax2 = plt.subplot2grid((nrows, ncols), (1, ft))
        Zax3 = plt.subplot2grid((nrows, ncols), (2, ft))
        Zax4 = plt.subplot2grid((nrows, ncols), (3, ft))
        Zax5 = plt.subplot2grid((nrows, ncols), (4, ft))
        Dax1 = plt.subplot2grid((nrows, ncols), (5, ft))
        Sax1 = plt.subplot2grid((nrows, ncols), (6, ft))

        Zim = mm_fits['Z%s'%(fit_type)].data

        Z_mode = Zim[:,:,0].copy()
        Z_l68 = Zim[:,:,1]
        Z_u68 = Zim[:,:,2]
        Z_np = Zim[:,:,3]
        Z_lerr = abs(Z_l68 - Z_mode)
        Z_uerr = abs(Z_u68 - Z_mode)

        Z_modeforseg = Z_mode.copy()
        Z_modeforseg[(Z_lerr > zmax) | (Z_uerr > zmax)] = np.nan

        if True:
            #Make segmentation map off of the Zmap
            segm = detect_sources(Z_modeforseg, -99, npixels=5)
            midx = 20
            midy = 20

            lbl_interest = array(segm.data)[int(midx), int(midy)]
            if lbl_interest == 0:
                small_box = array(segm.data)[int(midx - 1):int(midx +1), int(midy - 1):int(midy +1)].ravel() 
                if len(small_box[small_box > 0]) > 0:  lbl_interest = min(small_box[small_box > 0])
                else: lbl_interest = 99

        di_seg = dir_seg[20,20]

        Z_show = ma.masked_where((segm.data != lbl_interest) | (Z_lerr > zmax) | (Z_uerr > zmax) | (np.isnan(Z_mode)), Z_mode)
        Z_seg = ma.masked_where((dir_seg != di_seg ) | (Z_lerr > zmax) | (Z_uerr > zmax) | (np.isnan(Z_mode)), Z_mode)
        Z_both = ma.masked_where((dir_seg != di_seg ) | (segm.data != lbl_interest) | (Z_lerr > zmax) | (Z_uerr > zmax) | (np.isnan(Z_mode)), Z_mode)
        eZl_both = ma.masked_where((dir_seg != di_seg ) | (segm.data != lbl_interest) | (Z_lerr > zmax) | (Z_uerr > zmax) | (np.isnan(Z_mode)), Z_lerr)
        eZu_both = ma.masked_where((dir_seg != di_seg ) | (segm.data != lbl_interest) | (Z_lerr > zmax) | (Z_uerr > zmax) | (np.isnan(Z_mode)), Z_uerr)
        Z_np_both = ma.masked_where((dir_seg != di_seg ) | (segm.data != lbl_interest) | (Z_lerr > zmax) | (Z_uerr > zmax) | (np.isnan(Z_mode)), Z_np)

        cmap = plt.cm.viridis
        cmap.set_bad('k')

        vmin = 7.8
        vmax = 9.2

        Zax1.imshow(Z_mode, cmap = cmap, vmin = vmin,  vmax = vmax)
        Zax2.imshow(Z_lerr, cmap = cmap, vmin = 0.,    vmax = 0.5)
        Zax3.imshow(Z_show, cmap = cmap, vmin = vmin,  vmax = vmax)
        Zax4.imshow(Z_seg,  cmap = cmap, vmin = vmin,  vmax = vmax)
        print (Z_both[20,20])
        Zax5.imshow(Z_both,  cmap = cmap, vmin = vmin, vmax = vmax)

        Dax1.imshow(dir_im,  cmap = 'Greys_r')
        dir_seg_clip = dir_seg.copy()
        dir_seg_clip[dir_seg_clip != di_seg ] = 0
        Sax1.imshow(dir_seg_clip, cmap = 'Greys_r')

        dir_im_seg = dir_im.copy()
        dir_im_seg[dir_seg != di_seg] = 0.
        midx, midy = photutils.centroid_2dg(dir_im_seg)

        fs = 15
        if ft == 0:
            Zax1.annotate('Z',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Zax2.annotate('Z$_{err}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Zax3.annotate('Z$_{Z, seg}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Zax4.annotate('Z$_{dir, seg}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Zax5.annotate('Z$_{both, seg}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Dax1.annotate('direct$_{F105W}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
            Sax1.annotate('seg$_{F105W}$',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)

        fs = 8

        if ft == 1: Zax1.annotate('smoothed',(0.1, 0.95),va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
        if ft == 2: Zax1.annotate('exctinction-corrected',(0.1, 0.95),va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)
        if ft == 3: Zax1.annotate('smoothed\nexctinction-corrected',(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)




        for ax in [Zax1, Zax2, Zax3, Zax4, Zax5, Dax1, Sax1]:
            ax.axis('off')
            ax.set_xticklabels([''])
            ax.set_yticklabels([''])
            ax.plot(midx, midy, 'x', color = 'black')


        x = (np.arange(0, shape(Z_both)[0]) - midx)
        y = (np.arange(0, shape(Z_both)[1]) - midy)

        xv, yv = np.meshgrid(x, y)
        r = sqrt(xv**2. + yv**2.)


        res['r{}'.format(fit_type)] = r
        res['Z{}'.format(fit_type)] = Z_both
        res['Npeaks{}'.format(fit_type)] = Z_np_both

        res['elZ{}'.format(fit_type)] = eZl_both
        res['euZ{}'.format(fit_type)] = eZu_both

        res['midx{}'.format(fit_type)] = midx
        res['midy{}'.format(fit_type)] = midy

        for l, (line, izi_line, wave) in enumerate(lines):
            try: im = mm_fits['%s%s'%(line, fit_type)].data
            except: continue
            ax = plt.subplot2grid((nrows, ncols), (ct+nzcols, ft))
            ct+=1
            ax.set_xticklabels([''])
            ax.set_yticklabels([''])
            ax.axis('off')
            #if ft == 0: 
            pl = ax.imshow(im, cmap = 'viridis')
            #else: ax.imshow(im, cmap = 'viridis', vmin = pl.get_clim()[0], vmax = pl.get_clim()[1])

            fs = 15
            if ft == 0:
                ax.annotate(  '%s'%(line),(0.1, 0.95), va = 'top', xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = fs)


    np.save('/Users/rsimons/Desktop/clear/izi_metal_profiles/{}_{}.npy'.format(fld, di), res)



    figname = figdir + '/' + fl.split('/')[-1].replace('fits', 'png')
    print ('saving %s'%figname)
    fig.tight_layout()
    fig.savefig(figname)

    plt.close('all')






