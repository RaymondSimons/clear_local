import astropy
from astropy.io import fits
import glob
from glob import glob

field = 'GN2'
fls = glob('/user/rsimons/grizli_extractions/PREP/*.full.fits')
cat = open('/user/rsimons/grizli_extractions/Catalogs/%s_lines_grizli.cat'%field, 'w+')

lines = ['OII', 'OIII', 'Ha', 'Hb']


cat.write('#(0) Field\n')
cat.write('#(1) ID\n')
cat.write('#(2) n_lines\n')
cat.write('#(3) OII flux, 1e-17 erg/s/cm2\n')
cat.write('#(4) OII flux err, 1e-17 erg/s/cm2\n')
cat.write('#(5) OIII flux, 1e-17 erg/s/cm2\n')
cat.write('#(6) OIII flux err, 1e-17 erg/s/cm2\n')
cat.write('#(7) Ha flux, 1e-17 erg/s/cm2\n')
cat.write('#(8) Ha flux err, 1e-17 erg/s/cm2\n')
cat.write('#(9) Hb flux, 1e-17 erg/s/cm2\n')
cat.write('#(10) Hb flux err, 1e-17 erg/s/cm2\n\n\n\n\n')


fluxs = zeros((4,2, len(fls))) - 99.
for f, fl in enumerate(fls):
	a = fits.open(fl)
	ID = a[0].header['ID']
	lines_f = a[0].header['HASLINES'].split(' ')
	if (len(lines_f) > 0) & (lines_f[0] !=''):
		for l, ln in enumerate(lines_f):
			flux_ln =  a[0].header['FLUX%.3i'%(l+1)]
			eflux_ln =  a[0].header['ERR%.3i'%(l+1)]
			ln_name =  a[0].header['LINE%.3i'%(l+1)]
			if ln_name == 'OII': j = 0
			elif ln_name == 'OIII':  j = 1
			elif ln_name == 'Ha': j = 2
			elif ln_name == 'Hb': j = 3
			fluxs[j,0,f] = flux_ln * 1.e17
			fluxs[j,1,f] = eflux_ln  * 1.e17

	n_lines = len(lines_f)

	cat.write('%s\t\t%i\t\t%i\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n'%(field, ID, n_lines, fluxs[0,0,f], fluxs[0,1,f], fluxs[1,0,f], fluxs[1,1,f], fluxs[2,0,f], fluxs[2,1,f], fluxs[3,0,f], fluxs[3,1,f]))


cat.close()