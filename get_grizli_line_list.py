from astropy.io import fits
import glob
from glob import glob

full_dir = '/Volumes/pegasus/clear/grizli_extractions/*/*/Prep'
fls = glob(full_dir + '/*full.fits')

line_list_new = []
for f, fl in enumerate(fls):
    print (f, len(fls))
    a = fits.getheader(fl)
    lines_f = a['HASLINES'].split(' ')
    if (len(lines_f) > 0) & (lines_f[0] !=''):
        for l, ln in enumerate(lines_f):
            ln_name =  a['LINE%.3i'%(l+1)]
            line_list_new.append(ln_name)

line_list_new = array(line_list_new)
line_list_new = unique(line_list_new)