import astropy
from astropy.io import fits
import glob
from glob import glob
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from numpy import *
plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['xtick.labelsize'] = 12





metal_dir = '/Volumes/pegasus/clear/metal_maps'


if False:
    fls = glob(metal_dir + '/*metals.fits')#[0:200]
    #fls = glob(metal_dir + '/ERSPRIME_40192_metals.fits')
    #fls = glob(metal_dir + '/ERSPRIME_40776_metals.fits')

    dicts =array([({}, 'r3'), ({}, 'r2'), ({}, 'r23'), ({}, 'o32'), ({}, 'all')])
    for cdict, _ in dicts:
        for key in ['z', 'uez', 'lez', 'r', 'er']: cdict[key] = []

    for f, fl in enumerate(fls):
        a = fits.open(fl)
        b = a.info(False)
        headernames = array([bb[1] for bb in b])
        for (cdict, calib) in array(dicts):
            if len(where(headernames == 'Z_' + calib.upper())[0] > 0):
                cdict['z'].extend(a['Z_' + calib.upper()].data[:,:,0].ravel())
                cdict['uez'].extend(a['Z_' + calib.upper()].data[:,:,1].ravel())
                cdict['lez'].extend(a['Z_' + calib.upper()].data[:,:,2].ravel())
                if calib != 'all':
                    cdict['r'].extend(a[calib.upper()].data[:,:].ravel())
                    cdict['er'].extend(a['E' + calib.upper()].data[:,:].ravel())
                else:
                    cdict['r'].extend(nan*a['Z_' + calib.upper()].data[:,:,0].ravel())
                    cdict['er'].extend(nan*a['Z_' + calib.upper()].data[:,:,0].ravel())

    master_hdulist = []
    prihdr = fits.Header()
    prihdu = fits.PrimaryHDU(header=prihdr)    
    master_hdulist.append(prihdu)
    for (cdict, calib) in dicts:
        cols = []
        gd = where(~isnan(cdict['z']))[0]
        gd = arange(len(cdict['z']))
        cols.append(fits.Column(name = 'z', array = np.array(cdict['z'])[gd], format = 'D'))
        cols.append(fits.Column(name = 'uez', array = np.array(cdict['uez'])[gd], format = 'D'))
        cols.append(fits.Column(name = 'lez', array = np.array(cdict['lez'])[gd], format = 'D'))    
        cols.append(fits.Column(name = 'r', array = np.array(cdict['r'])[gd], format = 'D'))
        cols.append(fits.Column(name = 'er', array = np.array(cdict['er'])[gd], format = 'D'))    

        cols = fits.ColDefs(cols)
        
        master_hdulist.append(fits.BinTableHDU.from_columns(cols, name = calib))
    thdulist = fits.HDUList(master_hdulist)
    thdulist.writeto('/Users/rsimons/Desktop/clear/metal_pixels.fits', overwrite = True)




if True:
    data = fits.open('/Users/rsimons/Desktop/clear/metal_pixels.fits')
    cal_array = array(['r3', 'r2', 'r23', 'o32', 'all'])
    cal_array = array(['r3', 'o32', 'all'])

    fig, axes = plt.subplots(1,len(cal_array), figsize = (3*len(cal_array), 3))


    for n, calib in enumerate(cal_array):
    
        if calib == 'all':
            r = data['r3'].data['r']
            er = data['r3'].data['er']/r/log(10)
            r = log10(r)
        else:
            r = data[calib.upper()].data['r']
            er = data[calib.upper()].data['er']/r/log(10)
            r = log10(r)
        z = data[calib.upper()].data['z']
        uez = data[calib.upper()].data['uez']
        lez = data[calib.upper()].data['lez']




        #axes[n].errorbar(r, z, xerr = er, yerr = [uez, lez], fmt = 'o')
        gd = where(er > 0.)#1./2./log(10))
        print nanmin(z[gd])
        axes[n].errorbar(z[gd], r[gd], xerr = [uez[gd], lez[gd]], yerr = er[gd], fmt = 'o', markersize = 0.10, linewidth = 0.0)
        import metal_calibs
        reload(metal_calibs)
        OH_m = linspace(6, 10, 100)

        if calib == 'r2': clb = metal_calibs.OH_R2
        if calib == 'r3': clb = metal_calibs.OH_R3
        if calib == 'r23': clb = metal_calibs.OH_R23
        if calib == 'o32': clb = metal_calibs.OH_O32
        if calib == 'all': clb = metal_calibs.OH_R23


        R_m = array([clb(OH) for OH in OH_m])

        axes[n].plot(OH_m, R_m,'r-')
        axes[n].set_xlabel(r'$\log$ (O/H)$_{\text{%s}}$ + 12'%calib.upper(), fontsize = 15)
        axes[n].set_ylabel(r'$\log$ ' + calib.upper(), fontsize = 15)
        if calib == 'all':  axes[n].set_ylabel(r'$\log$ R3', fontsize = 15)

        axes[n].set_xlim(6.8, 9.5)
        #axes[n].set_yscale('log')
       # axes[n].set_ylim(0.1, 40)
        axes[n].set_ylim(-1, 1.5)



    fig.tight_layout()
    fig.savefig('/Users/rsimons/Desktop/clear/figures/metal_pixels_all.png', dpi = 300)












