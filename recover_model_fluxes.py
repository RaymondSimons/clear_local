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


out_dir = '/user/rsimons/grizli_extractions/Catalogs/bestfit_model_fluxes/%s/'%field


for fl in fls[0:1]:
    out_file = out_dir + fl.split('/')[-1].replace('full.fits', 'fluxes.cat')
    data = fits.open(fl)
    print (out_file)





