import astropy
from astropy.io import fits
import grizli
from grizli import utils

from grizli.pipeline import auto_script
from grizli.multifit import GroupFLT, MultiBeam, get_redshift_fit_defaults
import glob
from glob import glob
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


class Pointing():
    """ Generalization of GN1, GS1, ERSPRIME, etc

    To change field-dependent catalog, seg map, ref image, and padding
    only need to change them here.

    """

    def __init__(self, field, ref_filter):
        if 'N' in field.upper():
            self.pad = 200
            #self.radec_catalog = PATH_TO_CATS + '/goodsN_radec.cat'
            self.radec_catalog = PATH_TO_CATS + '/gdn_radec_f140_14_24.cat'
            
            self.seg_map =  PATH_TO_CATS + '/Goods_N_plus_seg.fits'
            self.catalog =  PATH_TO_CATS + '/goodsn-F105W-astrodrizzle-v4.4_drz_sub_plus.cat'

            self.ref_image =  PATH_TO_CATS + '/goodsn-F105W-astrodrizzle-v4.4_drz_sci.fits'

            #self.tempfilt, self.coeffs, self.temp_sed, self.pz = readEazyBinary(MAIN_OUTPUT_FILE='goodsn_3dhst.v4.4', OUTPUT_DIRECTORY=PATH_TO_CATS, CACHE_FILE='Same')


            self.params = {}
            #self.params['CATALOG_FILE'] = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Catalog/{0}_3dhst.{1}.cat'.format('goodsn', 'v4.3')
            self.params['CATALOG_FILE'] = PATH_TO_CATS + '/{0}_3dhst.{1}.cat'.format('goodsn', 'v4.4')

            self.params['Z_STEP'] = 0.002
            self.params['Z_MAX'] = 4

            self.params['MAIN_OUTPUT_FILE'] = '{0}_3dhst.{1}.eazypy'.format('goodsn', 'v4.4')
            self.params['PRIOR_FILTER'] = 205


            self.params['MW_EBV'] = {'aegis':0.0066, 'cosmos':0.0148, 'goodss':0.0069, 
                                    'uds':0.0195, 'goodsn':0.0103}['goodsn']

            self.params['TEMPLATES_FILE'] = 'templates/fsps_full/tweak_fsps_QSF_12_v3.param'
            #self.translate_file = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Eazy/{0}_3dhst.{1}.translate'.format('goodsn', 'v4.3')
            self.translate_file = PATH_TO_CATS + '/{0}_{1}.translate'.format('goodsn', 'v4.4')




        elif 'S' in field.upper():
            self.pad = 200 # grizli default
            #self.radec_catalog = '../Catalogs/goodsS_radec.cat'
            #self.radec_catalog = PATH_TO_CATS + '/goodsS_radec.cat'
            self.radec_catalog = PATH_TO_CATS + '/gds_radec_f140_14_24.cat'
            self.seg_map =  PATH_TO_CATS + '/Goods_S_plus_seg.fits'
            self.catalog =  PATH_TO_CATS + '/goodss-F105W-astrodrizzle-v4.3_drz_sub_plus.cat'
            self.ref_image =  PATH_TO_CATS + '/goodss-F105W-astrodrizzle-v4.3_drz_sci.fits' 

            #self.tempfilt, self.coeffs, self.temp_sed, self.pz = readEazyBinary(MAIN_OUTPUT_FILE='goodss_3dhst.v4.3', OUTPUT_DIRECTORY=PATH_TO_CATS, CACHE_FILE='Same')


            self.params = {}
            #self.params['CATALOG_FILE'] = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Catalog/{0}_3dhst.{1}.cat'.format('goodss', 'v4.3')
            self.params['CATALOG_FILE'] = PATH_TO_CATS + '/{0}_3dhst.{1}.cat'.format('goodss', 'v4.4')
            self.params['Z_STEP'] = 0.002
            self.params['Z_MAX'] = 4

            self.params['MAIN_OUTPUT_FILE'] = '{0}_3dhst.{1}.eazypy'.format('goodss', 'v4.4')
            self.params['PRIOR_FILTER'] = 205


            self.params['MW_EBV'] = {'aegis':0.0066, 'cosmos':0.0148, 'goodss':0.0069, 
                                    'uds':0.0195, 'goodsn':0.0103}['goodsn']

            self.params['TEMPLATES_FILE'] = 'templates/fsps_full/tweak_fsps_QSF_12_v3.param'
            #self.translate_file = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Eazy/{0}_3dhst.{1}.translate'.format('goodss', 'v4.3')
            self.translate_file = PATH_TO_CATS + '/{0}_{1}.translate'.format('goodss', 'v4.4')




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


    p = Pointing(field = field, ref_filter = 'F105W')


    ez = eazy.photoz.PhotoZ(param_file=None, translate_file=p.translate_file, 
                            zeropoint_file=None, params=p.params, 
                            load_prior=True, load_products=False)

    ep = photoz.EazyPhot(ez, grizli_templates=templ0, zgrid=ez.zgrid)

    tab = utils.GTable()
    tab['ra'], tab['dec'], tab['id']  = [mb.ra], [mb.dec], id
    phot, ii, dd = ep.get_phot_dict(tab['ra'][0], tab['dec'][0])

    mb.set_photometry(phot)


    A_phot = mb._interpolate_photometry(z=tfit['z'], templates=templ1)
    A_model = A_phot.T.dot(data[1].data['coeffs'])

    print (out_file)





