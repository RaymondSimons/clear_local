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




#field = 'GN2'
def per_field(field,templ0, BOUNDED_DEFAULTS):
    prep_dir = glob('/user/rsimons/grizli_extractions/%s/*/Prep'%field)[0]
    fls = glob(prep_dir + '/%s_*.full.fits'%field)


    os.chdir(prep_dir)

    out_dir = '/user/rsimons/grizli_extractions/Catalogs/bestfit_model_fluxes/%s/'%field
    try: os.mkdir(out_dir)
    except: pass
    



    p = Pointing(field = field, ref_filter = 'F105W')

    eazy.symlink_eazy_inputs(path=os.path.dirname(eazy.__file__)+'/data', path_is_env=False)

    ez = eazy.photoz.PhotoZ(param_file=None, translate_file=p.translate_file, 
                            zeropoint_file=None, params=p.params, 
                            load_prior=True, load_products=False)

    ep = photoz.EazyPhot(ez, grizli_templates=templ0, zgrid=ez.zgrid)


    tab = utils.GTable()

    for fl in fls:
        out_file = out_dir + fl.split('/')[-1].replace('full.fits', 'model_fluxes.npy')
        data = fits.open(fl)

        di = int(fl.split('/')[-1].split('_')[-1].strip('full.fits'))
        mb = MultiBeam(prep_dir + '/{0}_{1:05d}.beams.fits'.format(field, di), fcontam = 0.2, group_name=field, MW_EBV = 0., sys_err = 0.03, psf = False, min_mask=0.01, min_sens=0.08)
        mb.initialize_masked_arrays()
        tfit = mb.template_at_z(z = data[1].header['Z_MAP'], templates = templ0, fit_background=True, fitter='nnls', bounded_kwargs=BOUNDED_DEFAULTS)
        tab['ra'], tab['dec'], tab['id']  = [mb.ra], [mb.dec], di
        phot, ii, dd = ep.get_phot_dict(tab['ra'][0], tab['dec'][0])

        mb.set_photometry(**phot, min_err = 0.03)


        #zg = data[1].data['zgrid']
        #c = data[1].data['coeffs']
        #temp = data['TEMPL'].data


        #zg = tfit['z']
        #temp = tfit['templates']
        A_phot = mb._interpolate_photometry(z=tfit['z'], templates=tfit['templates'])
        A_model = A_phot.T.dot(tfit['coeffs'])


        print (out_file)
        to_save = np.array([mb.photom_pivot, mb.photom_flam, mb.photom_eflam, A_model])
        np.save(out_file, to_save)

    os.chdir('/home/rsimons/git/clear_local')



if __name__ == '__main__':
    global PATH_TO_CATS
    PATH_TO_CATS = '/user/rsimons/grizli_extractions/Catalogs'

    grizli.utils.symlink_templates(force=True)
    BOUNDED_DEFAULTS = {'method':'bvls', 'tol':1.e-8, 'verbose':0}
    templ0 = grizli.utils.load_templates(fwhm=1200, line_complexes=True, stars=False, 
                                         full_line_list=None,  continuum_list=None, 
                                         fsps_templates=True)

    fields = np.array(['GS1','GS2', 'GS3', 'GS5', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7'])
    Parallel(n_jobs = 4, backend = 'threading')(delayed(per_field)(field = field, templ0 = templ0, BOUNDED_DEFAULTS = BOUNDED_DEFAULTS) for field in fields)


























