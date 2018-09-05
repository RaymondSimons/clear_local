#Python 2.7

import os
import glob
from glob import glob
import astropy
from astropy.io import fits


overlapping_fields = {'GN1':['GDN20'],
                      'GN2':['GDN8', 'GDN12', 'GDN21', 'GDN25'],
                      'GN3':['GDN18', 'GDN19', 'GDN22', 'GDN23'],
                      'GN4':['GDN21', 'GDN22', 'GDN25', 'GDN26'],
                      'GN5':['GDN17', 'GDN18'],
                      'GN7':['GDN3', 'GDN6', 'GDN7', 'GDN11'],
                      'ERSPRIME':['WFC3-ERSII-G01']}





fls = glob('/user/rsimons/grizli_extractions/RAW/*flt.fits')

field_to_use = 'GN2'

for f, fl in enumerate(fls):
	data = fits.open(fl)
	field = data[0].header['TARGNAME']
	if field.upper() == field_to_use or field.upper() in overlapping_fields[field_to_use]:

		print field, '\t', fl
		associated_files = glob(fl.strip('_flt.fits')*)

		for fl_2 in associated_files:
			print '\t', fl_2


