import os
import glob
from glob import glob



overlapping_fields = {'GN1':['GDN20'],
                      'GN2':['GDN8', 'GDN12', 'GDN21', 'GDN25'],
                      'GN3':['GDN18', 'GDN19', 'GDN22', 'GDN23'],
                      'GN4':['GDN21', 'GDN22', 'GDN25', 'GDN26'],
                      'GN5':['GDN17', 'GDN18'],
                      'GN7':['GDN3', 'GDN6', 'GDN7', 'GDN11'],
                      'ERSPRIME':['WFC3-ERSII-G01']}





fls = glob('/user/rsimons/grizli_extractions/RAW/*flt.fits')

for f, fl in enumerate(fls):
	print fl




