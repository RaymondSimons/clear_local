import numpy as np
import astropy
from astropy.io import fits

cat_dir = '../../clear/Catalogs'
cat_gdn = fits.open(cat_dir + '/GN_CLEAR.linefit.concat.v1.0.0.fits')
cat_gds = fits.open(cat_dir + '/GS_CLEAR.linefit.concat.v1.0.0.fits')





lines = ['OII', 'OIII', 'Hb']
fields = ['GDN', 'GDS']


if False:
    for sig in [3, 5,10]:
        for c, cat in enumerate([cat_gdn, cat_gds]):
            o2_flux = cat[1].data['OII_FLUX']
            o2_err  = cat[1].data['OII_FLUX_ERR']

            o3_flux = cat[1].data['OIII_FLUX']
            o3_err  = cat[1].data['OIII_FLUX_ERR']

            hb_flux = cat[1].data['Hb_FLUX']
            hb_err  = cat[1].data['Hb_FLUX_ERR']

            print fields[c], '(%i sigma detections)'%sig


            good_o2     = where((o2_flux > sig*o2_err) & (o2_flux != -99))[0]
            good_o3     = where((o3_flux > sig*o3_err) & (o3_flux != -99))[0]
            good_hb     = where((hb_flux > sig*hb_err) & (hb_flux != -99))[0]

            good_o2o3   = where((o2_flux > sig*o2_err) & (o2_flux != -99) & (o3_flux > sig*o3_err) & (o3_flux != -99))[0]
            good_o2o3hb = where((hb_flux > sig*hb_err) & (hb_flux != -99) & (o2_flux > sig*o2_err) & (o2_flux != -99) & (o3_flux > sig*o3_err) & (o3_flux != -99))[0]

            print '%30s: %i' %('[OII]',  len(good_o2))
            print '%30s: %i' %('[OIII]', len(good_o3))

            print '%30s: %i' %('Hb',     len(good_hb))


            print '%30s: %i' %('[OII] + [OIII]',  len(good_o2o3))
            if sig >= 3:
                for i in good_o2o3:
                    print '\t\t\t\t', cat[1].data['grism_id'][i], cat[1].data['z_max_grism'][i]
            print '%30s: %i' %('[OII] + [OIII] + [Hb]',  len(good_o2o3hb))

            if sig >= 3:
                for i in good_o2o3hb:
                    print '\t\t\t\t', cat[1].data['grism_id'][i], cat[1].data['z_max_grism'][i]

        print '\n\n'






cat_dir = '/Users/rsimons/Desktop/clear/Catalogs'
cat_GN1 = np.loadtxt(cat_dir + '/GN1_lines_grizli.cat')
cat_GN2 = np.loadtxt(cat_dir + '/GN2_lines_grizli.cat')
cat_GN3 = np.loadtxt(cat_dir + '/GN3_lines_grizli.cat')

fields = ['GN1', 'GN2', 'GN3']

if True:
    for sig in [3,5,10]:
        for c, cat in enumerate([cat_GN1, cat_GN2, cat_GN3]):
            o2_flux = cat[:,2]
            o2_err  = cat[:,3]

            o3_flux = cat[:,4]
            o3_err  = cat[:,5]

            hb_flux = cat[:,8]
            hb_err  = cat[:,9]

            print fields[c], '(%i sigma detections)'%sig


            good_o2     = where((o2_flux > sig*o2_err) & (o2_flux != -99))[0]
            good_o3     = where((o3_flux > sig*o3_err) & (o3_flux != -99))[0]
            good_hb     = where((hb_flux > sig*hb_err) & (hb_flux != -99))[0]

            good_o2o3   = where((o2_flux > sig*o2_err) & (o2_flux != -99) & (o3_flux > sig*o3_err) & (o3_flux != -99))[0]
            good_o2o3hb = where((hb_flux > sig*hb_err) & (hb_flux != -99) & (o2_flux > sig*o2_err) & (o2_flux != -99) & (o3_flux > sig*o3_err) & (o3_flux != -99))[0]

            print '%30s: %i' %('[OII]',  len(good_o2))
            print '%30s: %i' %('[OIII]', len(good_o3))

            print '%30s: %i' %('Hb',     len(good_hb))


            print '%30s: %i' %('[OII] + [OIII]',  len(good_o2o3))
            if sig >= 3:
                for i in good_o2o3:
                    print '\t\t\t\t', cat[i, 0], cat[i, 10]
            print '%30s: %i' %('[OII] + [OIII] + [Hb]',  len(good_o2o3hb))

            if sig >= 3:
                for i in good_o2o3hb:
                    print '\t\t\t\t', cat[i, 0], cat[i, 10]

        print '\n\n'































