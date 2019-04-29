import astropy
from astropy.io import fits
from glob import glob
from astropy.convolution import Gaussian2DKernel, convolve_fft, Box2DKernel

plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 


metal_dir = '/Users/rsimons/Desktop/clear/metal_maps'



fls = glob(metal_dir + '/*.fits')



xmn = 28
xmx = 52


sn_limit_l = np.inf
sn_limit_u = np.inf



fls = array(['ERSPRIME_40078_metals.fits',
'ERSPRIME_40192_metals.fits',
'ERSPRIME_40776_metals.fits',
'ERSPRIME_42109_metals.fits',
'GS1_45069_metals.fits',
'GS4_24615_metals.fits',
'GS4_24851_metals.fits',
'GS4_25745_metals.fits',
'GS4_26021_metals.fits',
'GS4_26087_metals.fits',
'GS4_26550_metals.fits',
'GS4_28794_metals.fits',
'GS4_29349_metals.fits',
'GS4_29460_metals.fits',
'GS4_30427_metals.fits',
'GS4_31311_metals.fits',
'GS5_37637_metals.fits',
'GS5_38616_metals.fits',
'GS5_38513_metals.fits',
'GS5_41761_metals.fits'])

#fls = array(['ERSPRIME_40776_metals.fits'])
#fls = array(['GS4_26087_metals.fits'])
#fls = array(['GS4_25745_metals.fits'])


for f, fl in enumerate(fls):
    if True:#'ERSPRIME_40192' in fl:
        a = fits.open(metal_dir + '/' + fl)
        fig = figure(figsize = (5, 5))

        axes_m1   = plt.subplot2grid((4, 4), (0,0), rowspan = 1, colspan = 1, fig = fig)
        axes_m2   = plt.subplot2grid((4, 4), (0,1), rowspan = 1, colspan = 1, fig = fig)
        axes_m3   = plt.subplot2grid((4, 4), (0,2), rowspan = 1, colspan = 1, fig = fig)


        axes_O32   = plt.subplot2grid((4, 4), (0,3), rowspan = 1, colspan = 1, fig = fig)
        axes_R2    = plt.subplot2grid((4, 4), (1,3), rowspan = 1, colspan = 1, fig = fig)
        axes_R3    = plt.subplot2grid((4, 4), (2,3), rowspan = 1, colspan = 1, fig = fig)
        axes_R23   = plt.subplot2grid((4, 4), (3,3), rowspan = 1, colspan = 1, fig = fig)


        axes_Z    = plt.subplot2grid((4, 4), (1,0), rowspan = 3, colspan = 3, fig = fig)

        all_axes = [axes_m1, axes_m2, axes_m3, axes_O32, axes_R2, axes_R3, axes_R23, axes_Z]


        for l, ln in enumerate(array(['[OII]', '[OIII]', 'Hb'])):
            ax = all_axes[l]
            try:
                ln_im = a[ln.replace('[', '').replace(']','')].data
                eln_im = a['E' + ln.replace('[', '').replace(']','')].data

                for_vm = sort(ln_im[xmn:xmx, xmn:xmx].ravel())
                vmn = for_vm[int(0.15*len(for_vm))]
                vmx = for_vm[int(0.995*len(for_vm))]
                ax.imshow(ln_im, vmin = vmn, vmax = vmx, cmap = 'viridis')
                ax.annotate(ln, (0.97, 0.03), va = 'bottom', ha = 'right', fontsize = 20, xycoords = 'axes fraction', fontweight = 'bold', color = 'white')            
            except:
                print 'broke'
                pass





        zmap = a['Z_ALL'].data[:,:,0]
        ezmap_l = a['Z_ALL'].data[:,:,1]
        ezmap_u = a['Z_ALL'].data[:,:,2]

        zmap_masked =np.ma.masked_where((ezmap_l > sn_limit_l) |
                                        (ezmap_u > sn_limit_u) |
                                        (isnan(zmap)), zmap)
        cm = mpl.cm.cool
        cm.set_bad('black', 1.)


        for_zm = sort(zmap_masked.data[35:45,35:45][zmap_masked.mask[35:45,35:45] == False])

        zmn = round(math.floor(for_zm[int(0.03*len(for_zm))]*100)/100., 1)
        zmx = round(math.ceil(for_zm[int(0.95*len(for_zm))]*100)/100., 1)
        #zmn = 6.8



        mp = axes_Z.imshow(zmap_masked, vmin = zmn, vmax = zmx, cmap = cm)


        cbar_coords = [0.05, 0.06, 0.3, 0.05]

        cbaxes = fig.add_axes(cbar_coords)
        cbr = plt.colorbar(mp, cax = cbaxes, orientation = 'horizontal')
        cbr.ax.axes.tick_params(color = 'white', labelcolor = 'white')
        if zmx - zmn > 1.51: tcks = 0.5  
        elif zmx - zmn > 0.61: tcks = 0.2    
        elif zmx - zmn > 0.21: tcks = 0.1    
        else: tcks = 0.05

        for ax in all_axes:
            ax.set_xlim(xmn, xmx)
            ax.set_ylim(xmn, xmx)
            ax.axis('off')



        #cbr.set_ticks(arange(zmn, zmx)]+0.2, 0.2))
        cbr.set_ticks([zmn, zmx])

        #cbr.set_ticks(arange(6., 10, tcks))
        axes_Z.annotate(r'12 + $\log$(O/H)', (cbar_coords[0], cbar_coords[1] + cbar_coords[3]), xycoords = 'figure fraction', ha = 'left', va = 'bottom', fontsize = 18, fontweight = 'bold', color = 'white')
        bbox_props = dict(boxstyle="square", fc="k", ec=None, alpha=0.4)
        axes_Z.annotate('All', (0.93, 0.04), va = 'bottom', ha = 'right', fontsize = 25, xycoords = 'axes fraction', fontweight = 'bold', color = 'white', bbox = bbox_props)            




        for d, dg in enumerate(array(['O32', 'R2', 'R3', 'R23'])):
            ax = all_axes[3 + d]
            try:
                ln_im = a[ln.replace('[', '').replace(']','')].data
                zmap = a['Z_%s'%dg].data[:,:,0]
                ezmap_l = a['Z_%s'%dg].data[:,:,1]
                ezmap_u = a['Z_%s'%dg].data[:,:,2]
                zmap_masked =np.ma.masked_where((ezmap_l > sn_limit_l) |
                                                (ezmap_u > sn_limit_u) |
                                                (isnan(zmap)), zmap)

                for_zm = sort(zmap_masked.data[35:45,35:45][zmap_masked.mask[35:45,35:45] == False])    

                zmn = round(math.floor(for_zm[int(0.05*len(for_zm))]*100)/100., 1)
                zmx = round(math.ceil(for_zm[int(0.98*len(for_zm))]*100)/100., 1)
                '''
                if dg == 'O32':  
                    zmx = 8.0

                if dg == 'R2':  
                    #zmn = 7.4
                    zmx = 8.5
                    #print zmx
                '''
                mp = ax.imshow(zmap_masked, vmin = zmn, vmax = zmx, cmap = cm)
                cbar_coords = [0.78, 0.95 - 0.25*d, 0.15, 0.02]
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



            except:
                pass




        fig.subplots_adjust(left = 0.0, right = 1.0, top = 1.0, bottom = 0.00, wspace = 0.0, hspace = 0.01)
        fig.savefig('/Users/rsimons/Desktop/clear/figures/metal_maps_alma/' + fl.split('/')[-1].replace('.fits', '.png'), dpi = 300)
        plt.close('all')














