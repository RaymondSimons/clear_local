import astropy
from astropy.io import fits, ascii
from astropy.convolution import Gaussian2DKernel, convolve_fft, Box2DKernel
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 





plt.ioff()
plt.close('all')

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
            OH_arr[i,j]   = np.nanmean(8.54 - 0.59 * O3_arr_temp/O2_arr_temp)
    
    return O32_arr, O32_e_arr, OH_arr, OH_e_arr#np.std(OH_z_arr)




fig = figure(figsize = (11.0,4))

axes_O2   = plt.subplot2grid((2,13), (0,0), rowspan = 1, colspan = 2, fig = fig)
axes_O3   = plt.subplot2grid((2,13), (1,0), rowspan = 1, colspan = 2, fig = fig)
axes_Z    = plt.subplot2grid((2,13), (0,2), rowspan = 2, colspan = 4, fig = fig)
axes_p    = plt.subplot2grid((2,13), (0,7), rowspan = 2, colspan = 6, fig = fig)




field, di = 'GN3', 35204
field, di = 'GS2', 46938
#field, di = 'GS4', 20651

full = fits.open('/Users/rsimons/Desktop/random/%s_%i.full.fits'%(field, di))



lineo2 = full['LINE', 'OII'].data
lineo3 = full['LINE', 'OIII'].data

elineo2 = 1./np.sqrt(full['LINEWHT', 'OII'].data)
elineo3 = 1./np.sqrt(full['LINEWHT', 'OIII'].data)


kern = Box2DKernel(2.)
#kern = Gaussian2DKernel(0.8)

O2 = convolve_fft(lineo2, kern)
O3 = convolve_fft(lineo3, kern)





Z = 8.54 - 0.59 * O3/O2


#elineo2[:,:] = np.inf



Z_masked =np.ma.masked_where((O2/elineo2 < 0.8) | (O3/elineo3 < 0.8), Z)
cm = mpl.cm.cool
cm.set_bad('black', 1.)

axes_O2.imshow(flipud(O2), vmin = -0.1, vmax = 0.4)
axes_O3.imshow(flipud(O3), vmin = -0.1, vmax = 0.5)
mp = axes_Z.imshow(flipud(Z_masked), vmin = 6.8, vmax = 8.1, cmap = cm)





for ax in [axes_O2, axes_O3, axes_Z]: 
    ax.set_xticks([])
    ax.set_yticks([])
    dx = -1
    dy = -1
    ax.set_xlim(32 + dx, 48 + dx)
    ax.set_ylim(32 + dy, 48 + dy)

axes_O2.annotate('[OII]', (0.64, 0.06), fontsize = 20, xycoords = 'axes fraction', fontweight = 'bold', color = 'white')
axes_O3.annotate('[OIII]', (0.64, 0.06),fontsize = 20,  xycoords = 'axes fraction', fontweight = 'bold', color = 'white')
axes_Z.annotate(r'12 + $\log$(O/H)', (0.48, 0.13),fontsize = 20,  xycoords = 'axes fraction', fontweight = 'bold', color = 'white')

axes_p.set_ylabel(r'12 + log(O/H)', fontsize = 18)
axes_p.set_xlabel('distance from center (arcsec)', fontsize = 18)


cbaxes = fig.add_axes([0.27, 0.22, 0.18, 0.02]) 
cbr = plt.colorbar(mp, cax = cbaxes, orientation = 'horizontal')
cbr.ax.axes.tick_params(color = 'white', labelcolor = 'white')
cbr.set_ticks([7, 7.5, 8, 8.5])

cat = ascii.read('/Users/rsimons/Dropbox/rcs_clear/catalogs/z_r_O32.cat')


#fig.tight_layout()
fig.subplots_adjust(hspace = 0.0, left = 0.01, right = 0.99, top = 0.97, bottom = 0.15)
fig.savefig('/Users/rsimons/Desktop/random/R23_z.png', dpi = 300)


plt.close('all')