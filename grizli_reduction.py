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
'''
try: 
    from mastquery import query, overlaps
    use_mquery = True
except: 
    from hsaquery import query, overlaps
    use_mquery = False
'''



plt.ioff()
plt.close('all')
### A lot of this is working from the Grizli cookbook described here:
### https://github.com/gbrammer/grizli/blob/master/examples/WFC3IR_Reduction.ipynb
### and Katie and Iva's grizli pipeline for CLEAR


'''
overlapping_fields = {'GN1':['GDN20'],
                      'GN2':['GDN8', 'GDN12', 'GDN21', 'GDN25'],
                      'GN3':['GDN18', 'GDN19', 'GDN22', 'GDN23'],
                      'GN4':['GDN21', 'GDN22', 'GDN25', 'GDN26'],
                      'GN5':['GDN17', 'GDN18'],
                      'GN7':['GDN3', 'GDN6', 'GDN7', 'GDN11'],
                      'ERSPRIME':['WFC3-ERSII-G01']}
'''



def parse():
    '''
    Parse command line arguments
    ''' 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='''CLEAR grizli extractions.''')
    parser.add_argument('-field',       '--field',          default='GS1', help='field to extract')
    parser.add_argument('-mag_lim',     '--mag_lim',        default=22, help='field to extract')
    parser.add_argument('-mag_max',     '--mag_max',        default= 0, help='field to extract')
    parser.add_argument('-do_files',    '--do_files',       default = True, help = 'bool to load files')
    parser.add_argument('-do_model',    '--do_model',       default = True, help = 'bool to model spectra')
    parser.add_argument('-do_retrieve', '--do_retrieve',    action = "store_true", default = False, help = 'bool to retrieve files from MAST')
    parser.add_argument('-do_prep',     '--do_prep',        action = "store_true", default = False, help = 'bool to PREP files with Grizli')
    parser.add_argument('-new_model',   '--new_model',      action = "store_true", default = False, help = 'bool to create new Grizli models')
    parser.add_argument('-do_fit',      '--do_fit',         action = "store_true", default = False, help = 'bool to fit modeled spectra')
    parser.add_argument('-fit_min_id',  '--fit_min_id',     default = 0, help = 'ID to start on for the fit')
    parser.add_argument('-n_jobs',      '--n_jobs',         default = 2, help = 'number of threads')


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
            self.pad = 500 # really only necessary for GDN
            #self.radec_catalog = PATH_TO_CATS + '/goodsN_radec.cat'
            self.radec_catalog = PATH_TO_CATS + '/gdn_radec_f140_14_24.cat'
            
            self.seg_map =  PATH_TO_CATS + '/Goods_N_plus_seg.fits'
            self.catalog =  PATH_TO_CATS + '/goodsn-F105W-astrodrizzle-v4.4_drz_sub_plus.cat'
            self.ref_image =  PATH_TO_CATS + '/goodsn-F105W-astrodrizzle-v4.4_drz_sci.fits'

            self.tempfilt, self.coeffs, self.temp_sed, self.pz = readEazyBinary(MAIN_OUTPUT_FILE='goodsn_3dhst.v4.1', OUTPUT_DIRECTORY=PATH_TO_CATS, CACHE_FILE='Same')


            self.params = {}
            self.params['CATALOG_FILE'] = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Catalog/{0}_3dhst.{1}.cat'.format('goodsn', 'v4.1')
            self.params['Z_STEP'] = 0.002
            self.params['Z_MAX'] = 4

            self.params['MAIN_OUTPUT_FILE'] = '{0}_3dhst.{1}.eazypy'.format('goodsn', 'v4.1')
            self.params['PRIOR_FILTER'] = 205


            self.params['MW_EBV'] = {'aegis':0.0066, 'cosmos':0.0148, 'goodss':0.0069, 
                                    'uds':0.0195, 'goodsn':0.0103}['goodsn']

            self.params['TEMPLATES_FILE'] = 'templates/fsps_full/tweak_fsps_QSF_12_v3.param'
            self.translate_file = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Eazy/{0}_3dhst.{1}.translate'.format('goodsn', 'v4.1')




        elif 'S' in field.upper():
            self.pad = 200 # grizli default
            #self.radec_catalog = '../Catalogs/goodsS_radec.cat'
            #self.radec_catalog = PATH_TO_CATS + '/goodsS_radec.cat'
            self.radec_catalog = PATH_TO_CATS + '/gds_radec_f140_14_24.cat'
            self.seg_map =  PATH_TO_CATS + '/Goods_S_plus_seg.fits'
            self.catalog =  PATH_TO_CATS + '/goodss-F105W-astrodrizzle-v4.3_drz_sub_plus.cat'
            self.ref_image =  PATH_TO_CATS + '/goodss-F105W-astrodrizzle-v4.3_drz_sci.fits' 

            self.tempfilt, self.coeffs, self.temp_sed, self.pz = readEazyBinary(MAIN_OUTPUT_FILE='goodss_3dhst.v4.1', OUTPUT_DIRECTORY=PATH_TO_CATS, CACHE_FILE='Same')


            self.params = {}
            self.params['CATALOG_FILE'] = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Catalog/{0}_3dhst.{1}.cat'.format('goodss', 'v4.1')
            self.params['Z_STEP'] = 0.002
            self.params['Z_MAX'] = 4

            self.params['MAIN_OUTPUT_FILE'] = '{0}_3dhst.{1}.eazypy'.format('goodss', 'v4.1')
            self.params['PRIOR_FILTER'] = 205


            self.params['MW_EBV'] = {'aegis':0.0066, 'cosmos':0.0148, 'goodss':0.0069, 
                                    'uds':0.0195, 'goodsn':0.0103}['goodsn']

            self.params['TEMPLATES_FILE'] = 'templates/fsps_full/tweak_fsps_QSF_12_v3.param'
            self.translate_file = PATH_TO_CATS + '/{0}_3dhst.{1}.cats/Eazy/{0}_3dhst.{1}.translate'.format('goodss', 'v4.1')





def grizli_getfiles(run = True):
    if run == False: return
    else: 'Running grizli_getfiles...'

    os.chdir(PATH_TO_PREP)
    files = glob('%s/*flt.fits'%PATH_TO_RAW)
    info = grizli.utils.get_flt_info(files)
    visits, filters = grizli.utils.parse_flt_files(info=info, uniquename=True)
    return visits, filters

def grizli_prep(visits, field = '', run = True):
    if run == False: return
    else: 'Running grizli_prep...'

    print ('\n\n\n\n\n\n\n')
    product_names = np.array([visit['product'] for visit in visits])
    filter_names = np.array([visit['product'].split('-')[-1] for visit in visits])
    basenames = np.array([visit['product'].split('.')[0]+'.0' for visit in visits])
    for ref_grism, ref_filter in [('G102', 'F105W'), ('G141', 'F140W')]:
        print ('Processing %s + %s visits'%(ref_grism, ref_filter))
        for v, visit in enumerate(visits):
            product = product_names[v]
            basename = basenames[v]
            filt1 = filter_names[v]
            #print (filt1.lower())
            field_in_contest = basename.split('-')[0]

            #print (field_in_contest)
            #if field_in_contest.upper() == field.upper() or field_in_contest.upper() in overlapping_fields[field]:
            if (ref_filter.lower() == filt1.lower()):
                #found a direct image, now search for grism counterpart
                if len(np.where((basenames == basename) & (filter_names == ref_grism.lower()))[0]) > 0:
                    grism_index= np.where((basenames == basename) & (filter_names == ref_grism.lower()))[0][0]
                    #print(grism_index)
                    p = Pointing(field = field, ref_filter = ref_filter)
                    radec_catalog = p.radec_catalog
                    print (field_in_contest, visits[grism_index], radec_catalog)
                    #radec_catalog = None
                    status = process_direct_grism_visit(direct = visit,
                                                        grism = visits[grism_index],
                                                        radec = radec_catalog, 
                                                        align_mag_limits = [14, 24])
            else:
                print ('no grism associated with direct image %s'%basename)
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
        cpu_count=8)
    
    if new_model:
        print('Computing contamination models...')
        grp.compute_full_model(mag_limit=mag_lim)
    
        print('Re-computing continuum models with higher-order polynomials and subtracting off contamination..')
        grp.refine_list(poly_order=2, mag_limits=[16, 24], verbose=False)

        print('Saving contamination models')
        grp.save_full_data()
    
    return grp
        
def grizli_fit(grp, id, min_id, mag, field = '', mag_lim = 35, mag_lim_lower = 35, run = True, id_choose = None, ref_filter = 'F105W', use_pz_prior = True, use_phot = True, scale_phot = True, templ0 = None, templ1 = None, ez = None, ep = None, pline = None):
    if fit_bool == False: return
    if (mag <= mag_lim) & (mag >=mag_lim_lower) & (id > min_id):
    #if id in to_fits:
    #if id == id_choose:
        print(id, mag)
        beams = grp.get_beams(id, size=80) #size??
        if beams != []:
            print("beams: ", beams)
            mb = grizli.multifit.MultiBeam(beams, fcontam=1.0, group_name=field)
            mb.write_master_fits()
            
            # Fit polynomial model for initial continuum subtraction
            wave = np.linspace(2000,2.5e4,100)
            poly_templates = grizli.utils.polynomial_templates(
                wave=wave, 
                order=7,
                line=False)

            pfit = mb.template_at_z(
                z=0, 
                templates=poly_templates, 
                fit_background=True, 
                fitter='lstsq', 
                fwhm=1400, 
                get_uncertainties=2)


            if pfit != None:
            #    pass
            # Drizzle grisms / PAs
                hdu, fig = mb.drizzle_grisms_and_PAs(
                    size=32, 
                    fcontam=0.2, 
                    flambda=False, 
                    scale=1, 
                    pixfrac=0.5, 
                    kernel='point', 
                    make_figure=True, 
                    usewcs=False, 
                    zfit=pfit,
                    diff=True)
                # Save drizzled ("stacked") 2D trace as PNG and FITS
                fig.savefig('{0}_{1:05d}.stack.png'.format(field, id))
                hdu.writeto('{0}_{1:05d}.stack.fits'.format(field, id), clobber=True)


                try:
                    if use_pz_prior:
                        #use redshift prior from z_phot
                        prior = np.zeros((2, len(p.tempfilt['zgrid'])))
                        prior[0] = p.tempfilt['zgrid']
                        prior[1] = p.pz['chi2fit'][:,id]
                    else:
                        prior = None 
                    phot_scale_order = 0



                    tab = utils.GTable()
                    tab['ra'] = [mb.ra]
                    tab['dec'] = [mb.dec]

                    tab['id'] = id
                    phot, ii, dd = ep.get_phot_dict(tab['ra'][0], tab['dec'][0])

                    out = grizli.fitting.run_all(
                        id, 
                        t0=templ0, 
                        t1=templ1, 
                        fwhm=1200, 
                        zr=[0.0, 3.5], 
                        dz=[0.004, 0.0005], 
                        fitter='nnls',
                        group_name=field,
                        fit_stacks=True, 
                        prior=None, 
                        fcontam=0.,
                        pline=pline, 
                        mask_sn_limit=7, 
                        fit_only_beams=False,
                        fit_beams=True, 
                        root=field,
                        fit_trace_shift=False, 
                        phot=phot, 
                        verbose=True, 
                        scale_photometry=phot_scale_order, 
                        show_beams=True)
                    mb, st, fit, tfit, line_hdu = out
                    fit_hdu = fits.open('{0}_{1:05d}.full.fits'.format(field, id)) 

                    fit_hdu.info()
                    # same as the fit table above, redshift fit to the stacked spectra
                    fit_stack = Table(fit_hdu['ZFIT_STACK'].data) 


                    # zoom in around the initial best-guess with the individual "beam" spectra
                    fit_beam = Table(fit_hdu['ZFIT_BEAM'].data)   

                    templ = Table(fit_hdu['TEMPL'].data)
                    print('{0} has lines [{1}]'.format(fit_hdu.filename(), fit_hdu[0].header['HASLINES']))

                    # Helper script for plotting them, not generated automatically
                    fig = grizli.fitting.show_drizzled_lines(fit_hdu, size_arcsec=1.6, cmap='plasma_r')
                    fig.savefig('{0}_{1:05d}.line.png'.format(field, id))
                    plt.close('all')
                except:
                    print ('Problem in fitting.run_all')

                    plt.close('all')
    else:
        return

def retrieve_archival_data(field, retrieve_bool = False):
    if retrieve_bool == False: return

    os.chdir(HOME_PATH)    
    
    parent = query.run_query(box = None, proposal_id = [14227], instruments=['WFC3/IR', 'ACS/WFC'], 
                         filters = ['G102'], target_name = field)


    tabs = overlaps.find_overlaps(parent, buffer_arcmin=0.01, 
                                  filters=['G102', 'G141'], 
                                  instruments=['WFC3/IR','WFC3/UVIS','ACS/WFC'], close=False)

    pids = list(np.unique(tabs[0]['proposal_id']))

    tabs = overlaps.find_overlaps(parent, buffer_arcmin=0.01, proposal_id = pids,
                                  filters=['G102', 'G141', 'F098M', 'F105W', 'F125W', 'F140W'], 
                                  instruments=['WFC3/IR','WFC3/UVIS','ACS/WFC'], close=False)
    footprint_fits_file = glob('*footprint.fits')[0]
    jtargname = footprint_fits_file.strip('_footprint.fits')

    auto_script.fetch_files(field_root=jtargname, HOME_PATH=HOME_PATH, remove_bad=True, reprocess_parallel=True)

    print (pids)


if __name__ == '__main__':

    global PATH_TO_RAW, PATH_TO_PREP, PATH_TO_SCRIPTS, HOME_PATH, to_fits
    #to_fits = np.array([9116, 16736, 18108, 15610, 19451])
    args = parse()
    #to_fits = np.array([17829])
    #id_choose = 23116

    field           = args['field']
    mag_lim         = args['mag_lim']
    mag_max         = args['mag_max']
    files_bool      = args['do_files']
    retrieve_bool   = args['do_retrieve']
    prep_bool       = args['do_prep']
    model_bool      = args['do_model']
    new_model       = args['new_model']
    fit_bool        = args['do_fit']
    fit_min_id      = int(args['fit_min_id'])
    n_jobs          = int(args['n_jobs'])

    PATH_TO_SCRIPTS     = args['PATH_TO_SCRIPTS'] 
    PATH_TO_CATS        = args['PATH_TO_CATS']    
    PATH_TO_HOME        = args['PATH_TO_HOME']

    HOME_PATH           = PATH_TO_HOME + '/' + field

    if not os.path.isdir(HOME_PATH): os.system('mkdir %s'%HOME_PATH)

    print ('Changing to %s'%HOME_PATH)
    os.chdir(HOME_PATH)



    extra = retrieve_archival_data(field = field, retrieve_bool = retrieve_bool)

    PATH_TO_RAW         = glob(HOME_PATH + '/*/RAW')[0]
    PATH_TO_PREP        = glob(HOME_PATH + '/*/Prep')[0]



    print ('Changing to %s'%PATH_TO_PREP)
    os.chdir(PATH_TO_PREP)

    visits, filters = grizli_getfiles(run = files_bool)

    grizli_prep(visits = visits, field = field, run = prep_bool)

    grp = grizli_model(visits, field = field, ref_filter_1 = 'F105W', ref_grism_1 = 'G102', ref_filter_2 = 'F140W', ref_grism_2 = 'G141',
                       run = model_bool, new_model = new_model, mag_lim = mag_lim)




    if fit_bool:
        eazy.symlink_eazy_inputs(path=os.path.dirname(eazy.__file__)+'/data', 
                     path_is_env=False)
        templ0 = grizli.utils.load_templates(fwhm=1200, line_complexes=True, stars=False, 
                                             full_line_list=None,  continuum_list=None, 
                                             fsps_templates=True)

        # Load individual line templates for fitting the line fluxes
        templ1 = grizli.utils.load_templates(fwhm=1200, line_complexes=False, stars=False, 
                                             full_line_list=None, continuum_list=None, 
                                             fsps_templates=True)


        p = Pointing(field = field, ref_filter = 'F105W')


        pline = {'kernel': 'point', 'pixfrac': 0.2, 'pixscale': 0.1, 'size': 8, 'wcs': None}
        ez = eazy.photoz.PhotoZ(param_file=None, translate_file=p.translate_file, 
                                zeropoint_file=None, params=p.params, 
                                load_prior=True, load_products=False)

        ep = photoz.EazyPhot(ez, grizli_templates=templ0, zgrid=ez.zgrid)
            

        Parallel(n_jobs = n_jobs, backend = 'threading')(delayed(grizli_fit)(grp, id = id, min_id = fit_min_id, mag = mag, field = field, mag_lim = mag_lim, mag_lim_lower = mag_max,
                                                                        run = fit_bool, id_choose = 22945, use_pz_prior = False, use_phot = True, scale_phot = True,
                                                                        templ0 = templ0, templ1 = templ1, ez = ez, ep = ep, pline = pline,) for id, mag in zip(np.array(grp.catalog['NUMBER']), np.array(grp.catalog['MAG_AUTO'])))






    print ('Changing to %s'%PATH_TO_SCRIPTS)
    os.chdir(PATH_TO_SCRIPTS)










    '''
    parent = query.run_query(box = None, proposal_id = [14227], instruments=['WFC3/IR', 'ACS/WFC'], 
                     filters = ['G102'], target_name = field)

    # Find all G102 and G141 observations overlapping the parent query in the archive
    tabs = overlaps.find_overlaps(parent, buffer_arcmin=0.01, 
                                  filters=['G102', 'G141'], 
                                  instruments=['WFC3/IR','WFC3/UVIS','ACS/WFC'], close=False)



    # rerun overlaps.find_overlaps
    #
    # pids = list(np.unique(tabs[0]['proposal_id']))
    # 
    # proposal_id = pids

    footprint_fits_file = glob('*footprint.fits')[0]
    jtargname = footprint_fits_file.strip('_footprint.fits')


    # A list of the target names
    fp_fits = fits.open(footprint_fits_file)
    overlapping_target_names = set(fp_fits[1].data['target'])


    # Move the footprint figure files to $HOME_PATH/query_results/ so that they are not overwritten
    os.system('cp %s/%s_footprint.fits %s/query_results/%s_footprint_%s.fits'%(HOME_PATH, jtargname, HOME_PATH, jtargname, 'all_G102_G141'))
    os.system('cp %s/%s_footprint.npy %s/query_results/%s_footprint_%s.npy'%(HOME_PATH, jtargname, HOME_PATH, jtargname,  'all_G102_G141'))
    os.system('cp %s/%s_footprint.pdf %s/query_results/%s_footprint_%s.pdf'%(HOME_PATH, jtargname, HOME_PATH, jtargname,  'all_G102_G141'))
    os.system('cp %s/%s_info.dat %s/query_results/%s_info_%s.dat'%(HOME_PATH, jtargname, HOME_PATH, jtargname,  'all_G102_G141'))

    # Loop targ_name by targ_name
    for t, targ_name in enumerate(overlapping_target_names):
        if use_mquery:
            extra = {'target_name':targ_name}
        else:
            extra = query.DEFAULT_EXTRA.copy()
            extra += ["TARGET.TARGET_NAME LIKE '%s'"%targ_name]
        
        # search the MAST archive again, this time looking for 
        # all grism and imaging observations with the given target name
        tabs = overlaps.find_overlaps(parent, buffer_arcmin=0.01, 
                                      filters=['G102', 'G141', 'F098M', 'F105W', 'F125W', 'F140W'], 
                                      instruments=['WFC3/IR','WFC3/UVIS','ACS/WFC'], 
                                      extra=extra, close=False)
        if True:
            # retrieve raw data from MAST
            s3_status = os.system('aws s3 ls s3://stpubdata --request-payer requester')
            auto_script.fetch_files(field_root=jtargname, HOME_PATH=HOME_PATH, remove_bad=True, 
                                    reprocess_parallel=True, s3_sync=(s3_status == 0))

        # Move the figure files to $HOME_PATH/query_results/ so that they are not overwritten
        os.system('mv %s/%s_footprint.fits %s/query_results/%s_footprint_%s.fits'%(HOME_PATH, jtargname, HOME_PATH, jtargname, targ_name))
        os.system('mv %s/%s_footprint.npy %s/query_results/%s_footprint_%s.npy'%(HOME_PATH, jtargname, HOME_PATH, jtargname, targ_name))
        os.system('mv %s/%s_footprint.pdf %s/query_results/%s_footprint_%s.pdf'%(HOME_PATH, jtargname, HOME_PATH, jtargname, targ_name))
        os.system('mv %s/%s_info.dat %s/query_results/%s_info_%s.dat'%(HOME_PATH, jtargname, HOME_PATH, jtargname, targ_name))

        os.chdir(HOME_PATH)
    extra = retrieve_archival_data(visits = visits, field = field, retrieve_bool = retrieve_bool)
    
    PATH_TO_RAW = glob.glob('/Volumes/wd/clear/%s/*/RAW'%field)[0]
    PATH_TO_PREP = glob.glob('/Volumes/wd/clear/%s/*/PREP'%field)[0]
    os.chdir(PATH_TO_PREP)
    #visits, filters = grizli_getfiles(run = files_bool)
    visits, filters = grizli_getfiles(run = files_bool)

    grizli_prep(visits = visits, ref_filter = 'F105W', ref_grism = 'G102', run = prep_bool)
    grizli_prep(visits = visits, ref_filter = 'F140W', ref_grism = 'G141', run = prep_bool)


    grp = grizli_model(visits, field = field, ref_filter_1 = 'F105W', ref_grism_1 = 'G102', ref_filter_2 = 'F140W', ref_grism_2 = 'G141',
                       run = model_bool, load_only = load_bool, mag_lim = mag_lim)

    if True:
        eazy.symlink_eazy_inputs(path=os.path.dirname(eazy.__file__)+'/data', 
                     path_is_env=False)
        templ0 = grizli.utils.load_templates(fwhm=1200, line_complexes=True, stars=False, 
                                             full_line_list=None,  continuum_list=None, 
                                             fsps_templates=True)

        # Load individual line templates for fitting the line fluxes
        templ1 = grizli.utils.load_templates(fwhm=1200, line_complexes=False, stars=False, 
                                             full_line_list=None, continuum_list=None, 
                                             fsps_templates=True)


        p = Pointing(field = field, ref_filter = 'F105W')


        pline = {'kernel': 'point', 'pixfrac': 0.2, 'pixscale': 0.1, 'size': 8, 'wcs': None}
        ez = eazy.photoz.PhotoZ(param_file=None, translate_file=p.translate_file, 
                                zeropoint_file=None, params=p.params, 
                                load_prior=True, load_products=False)

        ep = photoz.EazyPhot(ez, grizli_templates=templ0, zgrid=ez.zgrid)
        
    if True:
        Parallel(n_jobs = -1, backend = 'threading')(delayed(grizli_fit)(grp, id = id, min_id = 26130, mag = mag, field = field, mag_lim = mag_lim, mag_lim_lower = mag_lim_lower,
                                                                        run = fit_bool, id_choose = 22945, use_pz_prior = False, use_phot = True, scale_phot = True,
                                                                        templ0 = templ0, templ1 = templ1, ez = ez, ep = ep, pline = pline,) for id, mag in zip(np.array(grp.catalog['NUMBER']), np.array(grp.catalog['MAG_AUTO'])))

    '''

    #grizli_fit(grp, id = id, mag = mag, field = field, mag_lim = mag_lim, mag_lim_lower = mag_lim_lower, 
    #            run = fit_bool, id_choose = 22945, 
    #            use_pz_prior = False, use_phot = True, scale_phot = True, 
    #            templ0 = templ0, templ1 = templ1, ez = ez, ep = ep, pline = pline)




























