import astropy
from astropy.io import fits
import grizli
from grizli.pipeline import auto_script
from grizli.multifit import GroupFLT, MultiBeam, get_redshift_fit_defaults
import glob
from glob import glob




field = 'GN2'
prep_dir = '/user/rsimons/grizli_extractions/%s/*/Prep'%field
fls = glob(prep_dir + '/%s_*_*.full.fits'%field)





for fl in fls:
    print (fl)





