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

plt.ioff()
plt.close('all')

def parse():
    '''
    Parse command line arguments
    ''' 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='''CLEAR grizli extractions.''')
    parser.add_argument('-field',       '--field',          default='GS1', help='field to extract')
    parser.add_argument('-mag_lim',     '--mag_lim',        type = int, default=25, help='field to extract')
    parser.add_argument('-mag_max',     '--mag_max',        type = int, default= 0, help='field to extract')
    parser.add_argument('-do_files',    '--do_files',       default = True, help = 'bool to load files')
    parser.add_argument('-do_model',    '--do_model',       default = True, help = 'bool to model spectra')

    parser.add_argument('-fwop', '--fwop',    action = "store_true", default = False, help = 'fit with photometry')
    parser.add_argument('-do_retrieve', '--do_retrieve',    action = "store_true", default = False, help = 'bool to retrieve files from MAST')
    parser.add_argument('-do_prep',     '--do_prep',        action = "store_true", default = False, help = 'bool to PREP files with Grizli')
    parser.add_argument('-new_model',   '--new_model',      action = "store_true", default = False, help = 'bool to create new Grizli models')
    parser.add_argument('-do_fit',      '--do_fit',         action = "store_true", default = False, help = 'bool to fit modeled spectra')
    parser.add_argument('-do_beams',    '--do_beams',         action = "store_true", default = False, help = 'bool to write beams files')
    parser.add_argument('-use_psf',      '--use_psf',         action = "store_true", default = False, help = 'use psf extraction in fitting routine')

    parser.add_argument('-fit_min_id',  '--fit_min_id',     type = int, default = 0, help = 'ID to start on for the fit')
    parser.add_argument('-n_jobs',      '--n_jobs',         type = int, default = 2, help = 'number of threads')
    parser.add_argument('-id_choose',   '--id_choose',         type = int, default = None, help = 'ID to fit')
    parser.add_argument('-pso',         '--pso',         type = int, default = 1, help = 'phot_scale_order')
    parser.add_argument('-PATH_TO_RAW'    , '--PATH_TO_RAW'    , default = '/user/rsimons/grizli_extractions/RAW', help = 'path to RAW directory')
    parser.add_argument('-PATH_TO_PREP'   , '--PATH_TO_PREP'   , default = '/user/rsimons/grizli_extractions/PREP', help = 'path to prep directory')
    parser.add_argument('-PATH_TO_SCRIPTS', '--PATH_TO_SCRIPTS', default = '/home/rsimons/git/clear_local', help = 'path to scripts directory')
    parser.add_argument('-PATH_TO_CATS'   , '--PATH_TO_CATS'   , default = '/user/rsimons/grizli_extractions/Catalogs', help = 'path to catalog directory')
    parser.add_argument('-PATH_TO_HOME'   , '--PATH_TO_HOME'   , default = '/user/rsimons/grizli_extractions', help = 'path to home directory sans field')



    args = vars(parser.parse_args())
    return args

def readEazyBinary(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same'):

    """
    Author: Gabe Brammer
    This function has been clipped from eazyPy.py in thethreedhst git respository
    https://github.com/gbrammer/threedhst/tree/master/threedhst

    tempfilt, coeffs, temp_sed, pz = readEazyBinary(MAIN_OUTPUT_FILE='photz', \
                                                OUTPUT_DIRECTORY='./OUTPUT', \
                                                CACHE_FILE = 'Same')

    Read Eazy BINARY_OUTPUTS files into structure data.
    
    If the BINARY_OUTPUTS files are not in './OUTPUT', provide either a relative or absolute path
    in the OUTPUT_DIRECTORY keyword.
    
    By default assumes that CACHE_FILE is MAIN_OUTPUT_FILE+'.tempfilt'.
    Specify the full filename if otherwise. 
    """
    
    #root='COSMOS/OUTPUT/cat3.4_default_lines_zp33sspNoU'
    
    root = OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE
    
    ###### .tempfilt
    if CACHE_FILE == 'Same':
        CACHE_FILE = root+'.tempfilt'
    
    if os.path.exists(CACHE_FILE) is False:
        print(('File, %s, not found.' %(CACHE_FILE)))
        return -1,-1,-1,-1
    
    f = open(CACHE_FILE,'rb')
    
    s = np.fromfile(file=f,dtype=np.int32, count=4)
    NFILT=s[0]
    NTEMP=s[1]
    NZ=s[2]
    NOBJ=s[3]
    tempfilt = np.fromfile(file=f,dtype=np.double,count=NFILT*NTEMP*NZ).reshape((NZ,NTEMP,NFILT)).transpose()
    lc = np.fromfile(file=f,dtype=np.double,count=NFILT)
    zgrid = np.fromfile(file=f,dtype=np.double,count=NZ)
    fnu = np.fromfile(file=f,dtype=np.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
    efnu = np.fromfile(file=f,dtype=np.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
    
    f.close()
    
    tempfilt  = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
                 'tempfilt':tempfilt,'lc':lc,'zgrid':zgrid,'fnu':fnu,'efnu':efnu}
    
    ###### .coeff
    f = open(root+'.coeff','rb')
    
    s = np.fromfile(file=f,dtype=np.int32, count=4)
    NFILT=s[0]
    NTEMP=s[1]
    NZ=s[2]
    NOBJ=s[3]
    coeffs = np.fromfile(file=f,dtype=np.double,count=NTEMP*NOBJ).reshape((NOBJ,NTEMP)).transpose()
    izbest = np.fromfile(file=f,dtype=np.int32,count=NOBJ)
    tnorm = np.fromfile(file=f,dtype=np.double,count=NTEMP)
    
    f.close()
    
    coeffs = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
              'coeffs':coeffs,'izbest':izbest,'tnorm':tnorm}
              
    ###### .temp_sed
    f = open(root+'.temp_sed','rb')
    s = np.fromfile(file=f,dtype=np.int32, count=3)
    NTEMP=s[0]
    NTEMPL=s[1]
    NZ=s[2]
    templam = np.fromfile(file=f,dtype=np.double,count=NTEMPL)
    temp_seds = np.fromfile(file=f,dtype=np.double,count=NTEMPL*NTEMP).reshape((NTEMP,NTEMPL)).transpose()
    da = np.fromfile(file=f,dtype=np.double,count=NZ)
    db = np.fromfile(file=f,dtype=np.double,count=NZ)
    
    f.close()
    
    temp_sed = {'NTEMP':NTEMP,'NTEMPL':NTEMPL,'NZ':NZ,\
              'templam':templam,'temp_seds':temp_seds,'da':da,'db':db}
              
    ###### .pz
    if os.path.exists(root+'.pz'):
        f = open(root+'.pz','rb')
        s = np.fromfile(file=f,dtype=np.int32, count=2)
        NZ=s[0]
        NOBJ=s[1]
        chi2fit = np.fromfile(file=f,dtype=np.double,count=NZ*NOBJ).reshape((NOBJ,NZ)).transpose()

        ### This will break if APPLY_PRIOR No
        s = np.fromfile(file=f,dtype=np.int32, count=1)
        
        if len(s) > 0:
            NK = s[0]
            kbins = np.fromfile(file=f,dtype=np.double,count=NK)
            priorzk = np.fromfile(file=f, dtype=np.double, count=NZ*NK).reshape((NK,NZ)).transpose()
            kidx = np.fromfile(file=f,dtype=np.int32,count=NOBJ)
            pz = {'NZ':NZ,'NOBJ':NOBJ,'NK':NK, 'chi2fit':chi2fit, 'kbins':kbins, 'priorzk':priorzk,'kidx':kidx}
        else:
            pz = None
        
        f.close()
        
    else:
        pz = None
    
    if False:
        f = open(root+'.zbin','rb')
        s = np.fromfile(file=f,dtype=np.int32, count=1)
        NOBJ=s[0]
        z_a = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        z_p = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        z_m1 = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        z_m2 = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        z_peak = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        f.close()
        
    ###### Done.    
    return tempfilt, coeffs, temp_sed, pz

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
            self.params['CATALOG_FILE'] = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Catalog/{0}_3dhst.{1}.cat'.format('goodsn', 'v4.4', 'goodsn', 'v4.4')

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
            self.params['CATALOG_FILE'] = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Catalog/{0}_3dhst.{1}.cat'.format('goodss', 'v4.4', 'goodss', 'v4.4')
            self.params['Z_STEP'] = 0.002
            self.params['Z_MAX'] = 4

            self.params['MAIN_OUTPUT_FILE'] = '{0}_3dhst.{1}.eazypy'.format('goodss', 'v4.4')
            self.params['PRIOR_FILTER'] = 205


            self.params['MW_EBV'] = {'aegis':0.0066, 'cosmos':0.0148, 'goodss':0.0069, 
                                    'uds':0.0195, 'goodsn':0.0103}['goodsn']

            self.params['TEMPLATES_FILE'] = 'templates/fsps_full/tweak_fsps_QSF_12_v3.param'
            #self.translate_file = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Eazy/{0}_3dhst.{1}.translate'.format('goodss', 'v4.3')
            self.translate_file = PATH_TO_CATS + '/{0}_{1}.translate'.format('goodss', 'v4.4')





def grizli_getfiles(run = True):
    if run == False: return
    else: 'Running grizli_getfiles...'

    os.chdir(PATH_TO_PREP)
    files = glob('%s/*flt.fits'%PATH_TO_RAW)
    info = grizli.utils.get_flt_info(files)
    visits, filters = grizli.utils.parse_flt_files(info=info, uniquename=True)
    return visits, filters

def grizli_model(visits, field = '', ref_filter_1 = 'F105W', ref_grism_1 = 'G102', ref_filter_2 = 'F140W', ref_grism_2 = 'G141', run = True, new_model = False, mag_lim = 25):
    if run == False: return

    all_grism_files = []
    all_direct_files = []
    product_names = np.array([visit['product'] for visit in visits])
    filter_names = np.array([visit['product'].split('-')[-1] for visit in visits])
    basenames = np.array([visit['product'].split('.')[0]+'.0' for visit in visits])

    for v, visit in enumerate(visits):
        product = product_names[v]
        basename = basenames[v]
        filt1 = filter_names[v]
        if (ref_filter_1.lower() in filt1) or (ref_filter_2.lower() in filt1):
            all_direct_files.extend(visit['files'])
            grism_index_1 = np.where((basenames == basename) & (filter_names == ref_grism_1.lower()))[0]
            grism_index_2 = np.where((basenames == basename) & (filter_names == ref_grism_2.lower()))[0]
            if len(grism_index_1) > 0: all_grism_files.extend(visits[grism_index_1[0]]['files'])
            if len(grism_index_2) > 0: all_grism_files.extend(visits[grism_index_2[0]]['files'])
    p = Pointing(field=field, ref_filter=ref_filter_1)


    if not new_model: print('Loading contamination models...')
    else: print('Initializing contamination models...')
    
    grp = GroupFLT(
        grism_files=all_grism_files, 
        direct_files=[], 
        ref_file = p.ref_image,
        seg_file = p.seg_map,
        catalog  = p.catalog,
        pad=p.pad,
        cpu_count=4)
    
    if new_model:
        print('Computing contamination models with flat model...')
        grp.compute_full_model(mag_limit=25, cpu_count = 4)
    
        print('Refine continuum/contamination models with poly_order polynomial, subtracting off contamination..')
        grp.refine_list(poly_order=2, mag_limits=[16, 24], verbose=False)

        #poly_order = 3

        print('Saving contamination models')
        grp.save_full_data()
    
    return grp
   




if __name__ == '__main__':

    global PATH_TO_RAW, PATH_TO_PREP, PATH_TO_SCRIPTS, HOME_PATH, to_fits
    #to_fits = np.array([9116, 16736, 18108, 15610, 19451])
    args = parse()
    #to_fits = np.array([17829])
    #id_choose = 23116

    field               = args['field']
    mag_lim             = args['mag_lim']
    mag_max             = args['mag_max']
    files_bool          = args['do_files']
    retrieve_bool       = args['do_retrieve']
    prep_bool           = args['do_prep']
    model_bool          = args['do_model']
    new_model           = args['new_model']
    beams_model         = args['new_model']
    fit_bool            = args['do_fit']
    beams_bool          = args['do_beams']
    use_psf             = args['use_psf']
    fit_min_id          = args['fit_min_id']
    n_jobs              = args['n_jobs']
    id_choose           = args['id_choose']
    phot_scale_order    = args['pso']
    fit_without_phot    = args['fwop']
    PATH_TO_SCRIPTS     = args['PATH_TO_SCRIPTS'] 
    PATH_TO_CATS        = args['PATH_TO_CATS']    
    #PATH_TO_CATS = '/Users/rsimons/Desktop/clear/Catalogs'
    PATH_TO_HOME        = args['PATH_TO_HOME']
    HOME_PATH           = PATH_TO_HOME + '/' + field


    if fit_without_phot == True: phot_scale_order = -1


    print('\n\n\n\n###################\nParameters\n\n')
    print('field            ', field            )
    print('mag_lim          ', mag_lim          )
    print('mag_max          ', mag_max          )
    print('files_bool       ', files_bool       )
    print('retrieve_bool    ', retrieve_bool    )
    print('prep_bool        ', prep_bool        )
    print('model_bool       ', model_bool       )
    print('new_model        ', new_model        )
    print('beams_bool       ', beams_bool       )
    print('fit_bool         ', fit_bool         )
    print('use_psf          ', use_psf          )
    print('fit_min_id       ', fit_min_id       )
    print('n_jobs           ', n_jobs           )
    print('id_choose           ', id_choose           )
    print('phot_scale_order ', phot_scale_order )
    print('fit_without_phot ', fit_without_phot )
    print('PATH_TO_SCRIPTS  ', PATH_TO_SCRIPTS  )
    print('PATH_TO_CATS     ', PATH_TO_CATS     )
    print('PATH_TO_HOME     ', PATH_TO_HOME     )
    print('HOME_PATH        ', HOME_PATH        )
    print('\n\n\n\n####################\n\n\n\n')



    if not os.path.isdir(HOME_PATH): os.system('mkdir %s'%HOME_PATH)

    print ('Changing to %s'%HOME_PATH)
    os.chdir(HOME_PATH)

    PATH_TO_RAW         = glob(HOME_PATH + '/*/RAW')[0]
    PATH_TO_PREP        = glob(HOME_PATH + '/*/Prep')[0]

    print ('Changing to %s'%PATH_TO_PREP)
    os.chdir(PATH_TO_PREP)


    visits, filters = grizli_getfiles(run = files_bool)
    grp = grizli_model(visits, field = field, ref_filter_1 = 'F105W', ref_grism_1 = 'G102', ref_filter_2 = 'F140W', ref_grism_2 = 'G141',
                       run = model_bool, new_model = new_model, mag_lim = mag_lim)    


    to_save = np.array([grp.catalog['NUMBER'], grp.catalog['MAG_AUTO']])
    np.save('/user/rsimons/grizli_extractions/Catalogs/model_catalog/%s_catalog.npy'%field, to_save)

    print ('Changing to %s'%PATH_TO_SCRIPTS)
    os.chdir(PATH_TO_SCRIPTS)

























