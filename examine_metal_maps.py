import astropy
from astropy.io import fits
from glob import glob
from astropy.convolution import Gaussian2DKernel, convolve_fft, Box2DKernel
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import numpy as np
from numpy import *
plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 

metal_dir = '/Volumes/pegasus/clear/metal_maps'

xmn = 28
xmx = 52

fls = glob(metal_dir + '/*.fits')#[0:10]

SN = 2.



def determine_zm(zmap_masked, xmn, xmx):

    for_zm_small = sort(zmap_masked.data[35:45,35:45][zmap_masked.mask[35:45,35:45] == False])
    for_zm_all = sort(zmap_masked.data[xmn:xmx,xmn:xmx][zmap_masked.mask[xmn:xmx,xmn:xmx] == False])

    if len(for_zm_small) > 1:
        zmn = round(math.floor(for_zm_small[int(0.03*len(for_zm_small))]*100)/100., 1)
        zmx = round(math.ceil(for_zm_small[int(0.95*len(for_zm_small))]*100)/100., 1)
        return zmn, zmx
    elif len(for_zm_all) > 1:
        zmn = round(math.floor(for_zm_all[int(0.03*len(for_zm_all))]*100)/100., 1)
        zmx = round(math.ceil(for_zm_all[int(0.95*len(for_zm_all))]*100)/100., 1)
        return zmn, zmx
    else:
        zmn = 6
        zmx = 9
        return zmn, zmx


for f, fl in enumerate(fls):
    if True:#'ERSPRIME_40192' in fl:
        #a = fits.open(metal_dir + '/' + fl)
        a = fits.open(fl)
        fld = fl.split('/')[-1].split('_')[0]
        di = fl.split('/')[-1].split('_')[1]


        full_fl = glob('/Volumes/pegasus/clear/grizli_extractions/%s/j*/Prep/%s_%s.full.fits'%(fld, fld, di))[0]

        full = fits.open(full_fl)
        fig = plt.figure(figsize = (7/4.*5, 5))

        axes_m1   = plt.subplot2grid((4, 7), (0,0), rowspan = 1, colspan = 1, fig = fig)
        axes_m2   = plt.subplot2grid((4, 7), (0,1), rowspan = 1, colspan = 1, fig = fig)
        axes_m3   = plt.subplot2grid((4, 7), (0,2), rowspan = 1, colspan = 1, fig = fig)


        axes_O32   = plt.subplot2grid((4, 7), (0,3), rowspan = 1, colspan = 1, fig = fig)
        axes_R2    = plt.subplot2grid((4, 7), (0,4), rowspan = 1, colspan = 1, fig = fig)
        axes_R3    = plt.subplot2grid((4, 7), (0,5), rowspan = 1, colspan = 1, fig = fig)
        axes_R23   = plt.subplot2grid((4, 7), (0,6), rowspan = 1, colspan = 1, fig = fig)


        axes_Z    = plt.subplot2grid((4, 7), (1,3), rowspan = 3, colspan = 3, fig = fig)
        axes_H    = plt.subplot2grid((4, 7), (1,0), rowspan = 3, colspan = 3, fig = fig)



        all_axes = [axes_m1, axes_m2, axes_m3, axes_O32, axes_R2, axes_R3, axes_R23, axes_Z, axes_H]


        for l, ln in enumerate(array(['[OII]', '[OIII]', 'Hb'])):
            try: ln_im = a[ln.replace('[', '').replace(']','')].data
            except: continue

            ax = all_axes[l]
            eln_im = a['E' + ln.replace('[', '').replace(']','')].data
            for_vm = sort(ln_im[xmn:xmx, xmn:xmx].ravel())
            vmn = for_vm[int(0.15*len(for_vm))]
            vmx = for_vm[int(0.995*len(for_vm))]
            ax.imshow(ln_im, vmin = vmn, vmax = vmx, cmap = 'viridis')
            ax.annotate(ln, (0.97, 0.03), va = 'bottom', ha = 'right', fontsize = 20, xycoords = 'axes fraction', fontweight = 'bold', color = 'white')            


        if True:
            eR_map_all = []
            for d, dg in enumerate(array(['O32', 'R2', 'R3', 'R23'])):
                try: zmap = a['Z_%s'%dg].data[:,:,0]
                except: continue
                cm = mpl.cm.cool
                ax = all_axes[3 + d]               
                ezmap_l = a['Z_%s'%dg].data[:,:,1]
                ezmap_u = a['Z_%s'%dg].data[:,:,2]

                R = a['%s'%dg].data
                eR = a['e%s'%dg].data
                eR = eR/R/log(10)
                R = log10(10)
                zmap_masked =np.ma.masked_where((eR > 1./SN/log(10)) |
                                                (isnan(zmap)), zmap)

                eR_map_all.append(eR)
                zmn, zmx = determine_zm(zmap_masked, xmn, xmx)
                mp = ax.imshow(zmap_masked, vmin = zmn, vmax = zmx, cmap = cm)
                cbar_coords = [0.45 + 0.143*d, 0.95, 0.06, 0.02]
                cbaxes = fig.add_axes(cbar_coords)
                cbr = plt.colorbar(mp, cax = cbaxes, orientation = 'horizontal')
                cbr.ax.axes.tick_params(color = 'white', labelcolor = 'white', labelsize = 8)
                if zmx - zmn > 1.51: tcks = 0.5  
                elif zmx - zmn > 0.61: tcks = 0.2    
                elif zmx - zmn > 0.21: tcks = 0.1    
                else: tcks = 0.05
                cbr.set_ticks([zmn, zmx])

                bbox_props = dict(boxstyle="square", fc="k", ec=None, alpha=0.4)
                ax.annotate(dg, (0.93, 0.05), va = 'bottom', ha = 'right', fontsize = 15, xycoords = 'axes fraction', fontweight = 'bold', color = 'white', bbox = bbox_props)            


        eR_map_all = array(eR_map_all)
        eR_map_min = amin(eR_map_all, axis = 0)


        if True:
            zmap = a['Z_ALL'].data[:,:,0]
            cm.set_bad('black', 1.)

            #Mask where none of the ratios make it above this threshold
            zmap_masked = np.ma.masked_where((eR_map_min > 1./SN/log(10))|
                                            (isnan(zmap)), zmap)

            cm = mpl.cm.cool
            cm.set_bad('black', 1.)

            zmn, zmx = determine_zm(zmap_masked, xmn, xmx)            

            cbar_coords = [0.45, 0.06, 0.18, 0.05]

            mp = axes_Z.imshow(zmap_masked, vmin = zmn, vmax = zmx, cmap = cm)
            cbaxes = fig.add_axes(cbar_coords)

            cbr = plt.colorbar(mp, cax = cbaxes, orientation = 'horizontal')
            cbr.ax.axes.tick_params(color = 'white', labelcolor = 'white')
            if zmx - zmn > 1.51: tcks = 0.5  
            elif zmx - zmn > 0.61: tcks = 0.2    
            elif zmx - zmn > 0.21: tcks = 0.1    
            else: tcks = 0.05



            cbr.set_ticks([zmn, zmx])
            axes_Z.annotate(r'12 + $\log$(O/H)', (cbar_coords[0], cbar_coords[1] + cbar_coords[3]), xycoords = 'figure fraction', ha = 'left', va = 'bottom', fontsize = 18, fontweight = 'bold', color = 'white')
            bbox_props = dict(boxstyle="square", fc="k", ec=None, alpha=0.4)
            axes_Z.annotate('All', (0.93, 0.04), va = 'bottom', ha = 'right', fontsize = 25, xycoords = 'axes fraction', fontweight = 'bold', color = 'white', bbox = bbox_props)            




        im = full['DSCI', 'F105W'].data
        im-=min(im[xmn:xmx, xmn:xmx].ravel()) + 0.01


        for_vm = sort(log10(im[xmn:xmx, xmn:xmx].ravel()))
        vmn = for_vm[int(0.30*len(for_vm))]
        vmx = max(for_vm)#for_vm[int(0.995*len(for_vm))]

        axes_H.imshow(log10(im), cmap = 'Greys_r', vmin = vmn, vmax = vmx)


        axes_H.annotate('F105W', (0.93, 0.05), va = 'bottom', ha = 'right', fontsize = 20, xycoords = 'axes fraction', fontweight = 'bold', color = 'white', bbox = bbox_props)            


        for ax in all_axes:
            ax.set_xlim(xmn, xmx)
            ax.set_ylim(xmn, xmx)
            ax.axis('off')


        figdir = '/Users/rsimons/Desktop/clear/figures/metal_maps'

        #figdir = '/Users/rsimons/Desktop'

        fig.subplots_adjust(left = 0.0, right = 1.0, top = 1.0, bottom = 0.00, wspace = 0.0, hspace = 0.01)

        if not os.path.isdir('%s/SN_cut_%.1f'%(figdir, SN)):
            os.system('mkdir %s/SN_cut_%.1f/'%(figdir, SN))
        fig.savefig('%s/SN_cut_%.1f/'%(figdir, SN) + fl.split('/')[-1].replace('.fits', '.png'), dpi = 300)

        plt.close('all')














