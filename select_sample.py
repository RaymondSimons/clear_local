import astropy
from astropy.io import fits
import glob
from glob import glob


cat_fls = glob('/Users/rsimons/Dropbox/grizli_v2.0_cats/*_lines_grizli.fits')

f = open('/Users/rsimons/Dropbox/rcs_clear/z_r_sample.cat', 'w+')

for c, cat_fl in enumerate(cat_fls):
    cat = fits.open(cat_fl)
    O2_f = cat[1].data['OII_FLUX']
    O2_ef = cat[1].data['OII_FLUX_ERR']
    O3_f = cat[1].data['OIII_FLUX']
    O3_ef = cat[1].data['OIII_FLUX_ERR']
    good = where((O2_f/O2_ef > 5.) & ((O3_f/O3_ef > 5.)))[0]
    print len(good)
    field = cat_fl.split('/')[-1].strip('lines_grizli.fits')
    print field
    for g, gd in enumerate(good):
        f.write('%s\t%.5i\n'%(field, cat[1].data['ID'][gd]))


f.close()



