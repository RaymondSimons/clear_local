import astropy
from astropy.io import fits
import glob
from glob import glob

field = 'GN2'




#Create a catalog of redshifts for the GN2 field
field = 'GN2'

PATH_TO_PREP = '/user/rsimons/grizli_extractions/PREP'
PATH_TO_CATS= '/user/rsimons/grizli_extractions/Catalogs'


#fls = glob('../PREP/*.full.fits')
fls = glob(PATH_TO_PREP + '/*.full.fits')


#b = open('../Catalogs/zcat_%s.cat'%field, 'w+')
b = open(PATH_TO_CATS + '/zcat_%s_new.cat'%field, 'w+')

b.write('#(0) field\n')
b.write('#(1) ID\n')
b.write('#(2) z, 50th percentile\n')
b.write('#(3) z, 16th percentile\n')
b.write('#(4) z, 84th percentile\n')
b.write('#(5) z, 02nd percentile\n')
b.write('#(6) z, 97th percentile\n\n\n\n\n\n')



for f, fl in enumerate(fls):
	hdu = fits.open(fl)
	i = hdu[0].header['ID']
	hdr = hdu[1].header
	b.write('%s\t\t%i\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n'%(field, i, hdr['Z50'],  hdr['Z16'], hdr['Z84'], hdr['Z02'],hdr['Z97']))



b.close()
