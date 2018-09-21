import grizli
from grizli import utils
import eazy
import os
import grizli
from grizli.pipeline import photoz
import astropy.units as u

#PATH_TO_RAW = '/user/rsimons/grizli_extractions/RAW'
#PATH_TO_PREP = '/user/rsimons/grizli_extractions/PREP'
#PATH_TO_SCRIPTS = '/home/rsimons/git/clear_local'
#PATH_TO_CATS= '/user/rsimons/grizli_extractions/Catalogs'


PATH_TO_PREP = '/Users/rsimons/Desktop/clear/PREP'
PATH_TO_SCRIPTS = '/Users/rsimons/Desktop/clear/scripts'
PATH_TO_CATS= '/Users/rsimons/Desktop/clear/Catalogs'





os.chdir(PATH_TO_PREP)

eazy.symlink_eazy_inputs(path=os.path.dirname(eazy.__file__)+'/data', 
                         path_is_env=False)

field = 'goodsn'
version = 'v4.1'

params = {}
params['CATALOG_FILE'] = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Catalog/{0}_3dhst.{1}.cat'.format(field, version)
params['Z_STEP'] = 0.002
params['Z_MAX'] = 4

params['MAIN_OUTPUT_FILE'] = '{0}_3dhst.{1}.eazypy'.format(field, version)
params['PRIOR_FILTER'] = 205

# Galactic extinction
params['MW_EBV'] = {'aegis':0.0066, 'cosmos':0.0148, 'goodss':0.0069, 
                    'uds':0.0195, 'goodsn':0.0103}[field]

params['TEMPLATES_FILE'] = 'templates/fsps_full/tweak_fsps_QSF_12_v3.param'

translate_file = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Eazy/{0}_3dhst.{1}.translate'.format(field, version)

ez = eazy.photoz.PhotoZ(param_file=None, translate_file=translate_file, 
                        zeropoint_file=None, params=params, 
                        load_prior=True, load_products=False)

templ0 = grizli.utils.load_templates(fwhm=1200, line_complexes=True, stars=False, 
                                     full_line_list=None,  continuum_list=None, 
                                     fsps_templates=True)

tab = utils.GTable()
tab['ra'] = [189.20596196]
tab['dec'] = [62.22970839600013]

tab['id'] = 16736

idx3, dr3 = ez.cat.match_to_catalog_sky(tab)

ez.param['VERBOSITY'] = 1.
ez.fit_parallel(idx=idx3, verbose=False) 



ep = photoz.EazyPhot(ez, grizli_templates=templ0, zgrid=ez.zgrid)


























os.chdir(PATH_TO_SCRIPTS)