import matplotlib
matplotlib.use('agg')
import time
import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import drizzlepac
import grizli
import glob
from grizli import utils
import importlib
from grizli.prep import process_direct_grism_visit
#from hsaquery import query, overlaps
from grizli.pipeline import auto_script
from grizli.multifit import GroupFLT, MultiBeam, get_redshift_fit_defaults
import os, sys, argparse
from grizli.pipeline import photoz
from astropy.table import Table
import eazy
from joblib import Parallel, delayed
from glob import glob
from mastquery import query, overlaps



eazy.symlink_eazy_inputs(path=os.path.dirname(eazy.__file__)+'/data', 
             path_is_env=False)
templ0 = grizli.utils.load_templates(fwhm=1200, line_complexes=True, stars=False, 
                                     full_line_list=None,  continuum_list=None, 
                                     fsps_templates=True)

# Load individual line templates for fitting the line fluxes
templ1 = grizli.utils.load_templates(fwhm=1200, line_complexes=False, stars=False, 
                                     full_line_list=None, continuum_list=None, 
                                     fsps_templates=True)



t0, t1 = grizli.utils.load_quasar_templates(uv_line_complex = False, broad_fwhm = 2800, 
                                            narrow_fwhm = 1000, fixed_narrow_lines = True)




ID = 28121
phot_scale_order = 1

PATH_TO_CATS = '/user/rsimons/grizli_extractions/Catalogs'
translate_file = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Eazy/{0}_3dhst.{1}.translate'.format('goodsn', 'v4.1')

params = {}
params['CATALOG_FILE'] = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Catalog/{0}_3dhst.{1}.cat'.format('goodsn', 'v4.1')
params['Z_STEP'] = 0.002
params['Z_MAX'] = 4

params['MAIN_OUTPUT_FILE'] = '{0}_3dhst.{1}.eazypy'.format('goodsn', 'v4.1')
params['PRIOR_FILTER'] = 205


params['MW_EBV'] = {'aegis':0.0066, 'cosmos':0.0148, 'goodss':0.0069, 
                   'uds':0.0195, 'goodsn':0.0103}['goodsn']

params['TEMPLATES_FILE'] = 'templates/fsps_full/tweak_fsps_QSF_12_v3.param'




pline = {'kernel': 'point', 'pixfrac': 0.2, 'pixscale': 0.1, 'size': 8, 'wcs': None}
ez = eazy.photoz.PhotoZ(param_file=None, translate_file=translate_file, 
                        zeropoint_file=None, params=params, 
                        load_prior=True, load_products=False)

ep = photoz.EazyPhot(ez, zgrid=ez.zgrid)



field = 'GN4'
mb = grizli.multifit.MultiBeam('GN4_%.5i.beams.fits'%ID, fcontam=1.0, group_name=field)

tab = utils.GTable()
tab['ra'] = [mb.ra]
tab['dec'] = [mb.dec]

tab['id'] = id
phot, ii, dd = ep.get_phot_dict(tab['ra'][0], tab['dec'][0])


#gabe suggests mask_sn_limit



out = grizli.fitting.run_all(
    ID, 
    t0=t0, 
    t1=t1, 
    bad_pa_threshold = np.inf,
    fwhm=1200, 
    use_psf = True,
    fit_trace_shift = True,
    zr=[0.95, 1.15],       #zr = [0, 12.0]
    dz=[0.004, 0.0005], 
    fitter='nnls',
    group_name=field,
    fit_stacks=False,      #fit_stacks = False
    prior=None, 
    fcontam=0.,
    pline=pline, 
    mask_sn_limit=np.inf, 
    fit_only_beams=True, #fit_only_beams = True
    fit_beams=False,       #fit_beams = False
    root=field,
    phot=phot, 
    verbose=True, 
    scale_photometry=phot_scale_order, 
    show_beams=True)



