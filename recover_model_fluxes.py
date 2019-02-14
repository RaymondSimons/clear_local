import astropy
from astropy.io import fits
import grizli
from grizli import utils

from grizli.pipeline import auto_script
from grizli.multifit import GroupFLT, MultiBeam, get_redshift_fit_defaults
import glob
from glob import glob




field = 'GN2'
prep_dir = '/user/rsimons/grizli_extractions/%s/j123652p6215/Prep'%field
fls = glob(prep_dir + '/%s_*_*.full.fits'%field)


out_dir = '/user/rsimons/grizli_extractions/Catalogs/bestfit_model_fluxes/%s/'%field


for fl in fls[0:1]:
    out_file = out_dir + fl.split('/')[-1].replace('full.fits', 'fluxes.cat')
    data = fits.open(fl)

    di = int(fl.split('/')[-1].split('_')[-1].strip('full.fits'))
    mb = MultiBeam(prep_dir + '/{0}_{1:05d}.beams.fits'.format(field, di), group_name=field)

    BOUNDED_DEFAULTS = {'method':'bvls', 'tol':1.e-8, 'verbose':0}
    templ1 = grizli.utils.load_templates(fwhm=1200, line_complexes=False, stars=False, 
                                     full_line_list=None, continuum_list=None, 
                                     fsps_templates=True)

    tfit = mb.template_at_z(z = data[1].header['Z_MAP'], templates = templ1, fit_background=True, fitter='nnls', bounded_kwargs=BOUNDED_DEFAULTS)

    A_phot = mb._interpolate_photometry(z=tfit['z'], templates=templ1)
    A_model = A_phot.T.dot(data[1].data['coeffs'])

    print (out_file)





