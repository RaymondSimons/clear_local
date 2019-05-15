import astropy
from astropy.io import fits
import glob
from glob import glob
plt.ioff()
plt.close('all')

metal_dir = '/Users/rsimons/Desktop/clear/metal_maps'






if False:
    fls = glob(metal_dir + '/*metals.fits')
    fls = glob(metal_dir + '/ERSPRIME_40192_metals.fits')
    z_r2  = {}
    z_r3  = {}
    z_r23 = {}
    z_o32 = {}
    z_all = {}


    for z_ in [z_r3, z_r2, z_r23, z_o32, z_all]:
        z_['z'] = []
        z_['uez'] = []
        z_['lez'] = []
        z_['r'] = []
        z_['er'] = []


    for f, fl in enumerate(fls):
        a = fits.open(fl)
        b = a.info(False)
        headernames = array([bb[1] for bb in b])
        for (cdict, calib) in array([(z_r3, 'r3'), (z_r2, 'r2'), (z_r23, 'r23'), (z_o32, 'o32'), (z_all, 'all')]):
            if len(where(headernames == calib.upper())[0] > 0):
                cdict['z'].extend(a['Z_' + calib.upper()].data[:,:,0].ravel())
                cdict['uez'].extend(a['Z_' + calib.upper()].data[:,:,1].ravel())
                cdict['lez'].extend(a['Z_' + calib.upper()].data[:,:,2].ravel())
                if calib != 'all':
                    cdict['r'].extend(a[calib.upper()].data[:,:].ravel())
                    cdict['er'].extend(a['E' + calib.upper()].data[:,:].ravel())

                else:
                    cdict['r'].extend(nan*a['Z_' + calib.upper()].data[:,:,0].ravel())
                    cdict['er'].extend(nan*a['Z_' + calib.upper()].data[:,:,0].ravel())



if False:
    master_hdulist = []
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)    
    master_hdulist.append(prihdu)
    for (cdict, calib) in array([(z_r3, 'r3'), (z_r2, 'r2'), (z_r23, 'r23'), (z_o32, 'o32'), (z_all, 'all')]):
        cols = []
        gd = where(~isnan(cdict['z']))[0]
        cols.append(fits.Column(name = 'z', array = np.array(cdict['z'])[gd], format = 'D'))
        cols.append(fits.Column(name = 'uez', array = np.array(cdict['uez'])[gd], format = 'D'))
        cols.append(fits.Column(name = 'lez', array = np.array(cdict['lez'])[gd], format = 'D'))    
        cols.append(fits.Column(name = 'r', array = np.array(cdict['r'])[gd], format = 'D'))
        cols.append(fits.Column(name = 'er', array = np.array(cdict['er'])[gd], format = 'D'))    

        cols = fits.ColDefs(cols)
        
        master_hdulist.append(fits.BinTableHDU.from_columns(cols, name = calib))
    thdulist = fits.HDUList(master_hdulist)
    thdulist.writeto('/Users/rsimons/Desktop/clear/metal_pixels_all.fits', overwrite = True)

data = fits.open('/Users/rsimons/Desktop/clear/metal_pixels_all.fits')



cal_array = array(['r3', 'r2', 'r23', 'o32', 'all'])
cal_array = array(['r3', 'r2', 'r23', 'o32'])
fig, axes = plt.subplots(1,len(cal_array), figsize = (3*len(cal_array), 3))

for n, calib in enumerate(cal_array):
    r = data[calib.upper()].data['r']
    er = data[calib.upper()].data['er']
    z = data[calib.upper()].data['z']
    uez = data[calib.upper()].data['uez']
    lez = data[calib.upper()].data['lez']

    #axes[n].errorbar(r, z, xerr = er, yerr = [uez, lez], fmt = 'o')
    axes[n].plot(z, r, '.', markersize = 1)
    import metal_calibs
    reload(metal_calibs)
    OH_m = linspace(6, 10, 100)

    if calib == 'r2': clb = metal_calibs.OH_R2
    if calib == 'r3': clb = metal_calibs.OH_R3
    if calib == 'r23': clb = metal_calibs.OH_R23
    if calib == 'o32': clb = metal_calibs.OH_O32


    R_m = array([clb(OH) for OH in OH_m])

    axes[n].plot(OH_m, R_m,'r-')
    axes[n].set_xlabel(r'$\log$ (O/H) + 12')
    axes[n].set_ylabel(r'$\log$ ' + calib.upper())

    axes[n].set_xlim(7, 9.5)
    axes[n].set_ylim(-1, 1)


fig.tight_layout()
fig.savefig('/Users/rsimons/Desktop/clear/figures/metal_pixels_all.png', dpi = 300)












