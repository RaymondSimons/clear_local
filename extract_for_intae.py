import grizli
from grizli import utils
from grizli.prep import process_direct_grism_visit
from grizli.pipeline import auto_script
from grizli.multifit import GroupFLT, MultiBeam, get_redshift_fit_defaults
from grizli.pipeline import photoz
import os
import glob
from glob import glob
import os, sys, argparse
from joblib import Parallel, delayed
from astropy.io import ascii, fits
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
plt.ioff()


PATH_TO_CATS = '/Volumes/pegasus/clear/grizli_extractions/Catalogs'




def parse():
    '''
    Parse command line arguments
    ''' 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='''CLEAR grizli extractions.''')
    parser.add_argument('-field',       '--field',          default='GS1', help='field to extract')
    args = vars(parser.parse_args())
    return args



class Pointing():
    def __init__(self, field, ref_filter):
        if 'N' in field.upper():
            self.pad = 200
            self.seg_map =  PATH_TO_CATS + '/Goods_N_plus_unmatched_seg.fits'
            self.catalog = PATH_TO_CATS + '/goodsn-v4.4-withunmatched.cat'
            self.ref_image =  PATH_TO_CATS + '/goodsn-F105W-astrodrizzle-v4.4_drz_sci.fits'

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
            self.pad = 200
            self.seg_map =  PATH_TO_CATS + '/Goods_S_plus_unmatched_seg.fits'
            self.catalog = PATH_TO_CATS + '/goodss-v4.4-withunmatched.cat'
            self.ref_image =  PATH_TO_CATS + '/goodss-F105W-astrodrizzle-v4.3_drz_sci.fits' 
            self.params = {}
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


            
            
            






def run_all(field):
    HOME_PATH = '/Volumes/pegasus/clear/grizli_extractions/%s'%field
    PATH_TO_RAW         = glob(HOME_PATH + '/*/RAW')[0]
    PATH_TO_PREP        = glob(HOME_PATH + '/*/Prep_z67')[0]
    PATH_TO_CATS = '/Volumes/pegasus/clear/grizli_extractions/Catalogs'
    if 'S' in field:
        dis_matched = ascii.read(PATH_TO_CATS + '/z67_in_CLEAR_GS_matched_Finkelstein.txt')
        dis_unmatched = ascii.read(PATH_TO_CATS + '/z67_in_CLEAR_GS_unmatched_Finkelstein.txt')

    else:
        dis_matched = ascii.read(PATH_TO_CATS + '/z67_in_CLEAR_GN_matched_Finkelstein.txt')
        dis_unmatched = ascii.read(PATH_TO_CATS + '/z67_in_CLEAR_GN_unmatched_Finkelstein.txt')

    os.chdir(PATH_TO_PREP)
    files = glob('%s/*flt.fits'%PATH_TO_RAW)
    info = grizli.utils.get_flt_info(files)
    visits, filters = grizli.utils.parse_flt_files(info=info, uniquename=True)


    ref_filter_1 = 'F105W'
    ref_grism_1 = 'G102'
    ref_filter_2 = 'F140W' 
    ref_grism_2 = 'G141'
    mag_lim = 25

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



    print('Loading contamination models...')
    p = Pointing(field=field, ref_filter=ref_filter_1)

    grp = GroupFLT(
        grism_files=all_grism_files, 
        direct_files=[], 
        ref_file = p.ref_image,
        seg_file = p.seg_map,
        catalog  = p.catalog,
        pad=p.pad,
        cpu_count=1)


    if False:
        #Rewrite GrismFLT files to blot new targets (we only need to do this once)
        def rewrite_flt(FLT, p):
            FLT.process_seg_file(p.seg_map)
            FLT.save_full_pickle()

        Parallel(n_jobs = -1, backend = 'threading')(delayed(rewrite_flt) (FLT = grp.FLTs[g], p = p) for g in arange(len(grp.FLTs)))


        grp = GroupFLT(
            grism_files=all_grism_files, 
            direct_files=[], 
            ref_file = p.ref_image,
            seg_file = p.seg_map,
            catalog  = p.catalog,
            pad=p.pad,
            cpu_count=1)







    for id in dis_unmatched['ID_Finkelstein']:
        if not os.path.isfile('%s_%s.beams.fits'%(field, id)):
            beams = grp.get_beams(id, size=80)

            # can separate beams extraction, save, load in without needing models
            fcontam = 0.2
            if beams != []:
                mb = grizli.multifit.MultiBeam(beams, fcontam=fcontam, group_name=field)
                mb.write_master_fits()            



    for id in dis_matched['ID_Finkelstein']:
        if not os.path.isfile('%s_%s.beams.fits'%(field, id)):
            beams = grp.get_beams(id, size=80)
            # can separate beams extraction, save, load in without needing models
            fcontam = 0.2
            if beams != []:
                mb = grizli.multifit.MultiBeam(beams, fcontam=fcontam, group_name=field)
                mb.write_master_fits()            





    import eazy
    eazy.symlink_eazy_inputs(path=os.path.dirname(eazy.__file__)+'/data')

    templ0 = grizli.utils.load_templates(fwhm=1200, line_complexes=True, stars=False, 
                                         full_line_list=None,  continuum_list=None, 
                                         fsps_templates=True)

    # Load individual line templates for fitting the line fluxes
    templ1 = grizli.utils.load_templates(fwhm=1200, line_complexes=False, stars=False, 
                                         full_line_list=None, continuum_list=None, 
                                         fsps_templates=True)





    def do_fit(id, field, templ0, templ1, fcontam = 0.2):
        print (field, id)
        print (os.path.isfile(field + '_' + '%.5i.beams.fits'%id))
        print ((not os.path.isfile('%s_%s.full.fits'%(field, id))))
        if (os.path.isfile(field + '_' + '%.5i.beams.fits'%id)) & (not os.path.isfile('%s_%s.full.fits'%(field, id))):
            pline = {'kernel': 'point', 'pixfrac': 0.2, 'pixscale': 0.1, 'size': 8, 'wcs': None}
            mb = grizli.multifit.MultiBeam(field + '_' + '%.5i.beams.fits'%id, fcontam=fcontam, group_name=field)
            wave = np.linspace(2000,2.5e4,100)

            print ('creating poly_templates...')
            poly_templates = grizli.utils.polynomial_templates(wave=wave, order=7,line=False)
            pfit = mb.template_at_z(z=0, templates=poly_templates, fit_background=True, fitter='lstsq', fwhm=1400, get_uncertainties=2)
            print ('drizzle_grisms_and_PAs...')

            if not os.path.isfile('{0}_{1:05d}.stack.fits'.format(field, id)):
                hdu, fig = mb.drizzle_grisms_and_PAs(size=32, fcontam=fcontam, flambda=False, scale=1, 
                                                    pixfrac=0.5, kernel='point', make_figure=True, usewcs=False, 
                                                    zfit=pfit,diff=True)
                # Save drizzled ("stacked") 2D trace as PNG and FITS
                fig.savefig('{0}_{1:05d}.stack.png'.format(field, id))
                hdu.writeto('{0}_{1:05d}.stack.fits'.format(field, id), clobber=True)


            try:
                out = grizli.fitting.run_all(
                    id, 
                    t0=templ0, 
                    t1=templ1, 
                    fwhm=1200, 
                    zr=[5., 12.0],              #zr=[0.0, 12.0],    #suggests zr = [0, 12.0] if we want to extend redshift fit
                    dz=[0.004, 0.0005], 
                    fitter='nnls',
                    group_name=field,# + '_%i'%phot_scale_order,
                    fit_stacks=False,          #suggests fit_stacks = False, fit to FLT files
                    prior=None, 
                    fcontam=fcontam,           #suggests fcontam = 0.2
                    pline=pline, 
                    mask_sn_limit=np.inf,      #suggests mask_sn_limit = np.inf
                    fit_only_beams=True,       #suggests fit_only_beams = True
                    fit_beams=False,           #suggests fit_beams = False
                    root=field,
                    fit_trace_shift=False,  
                    bad_pa_threshold = np.inf, #suggests bad_pa_threshold = np.inf
                    phot=None, 
                    verbose=True, 
                    scale_photometry=0, 
                    show_beams=True,
                    use_psf = True)          #default: False
            except:
                print ('exception in fit for %s %s'%(field, id))



    #for id in dis_unmatched['ID_Finkelstein']:
    #    do_fit(id, field, templ0, templ1)

    Parallel(n_jobs = -1)(delayed(do_fit)(id, field, templ0, templ1) for id in dis_unmatched['ID_Finkelstein'])
    Parallel(n_jobs = -1)(delayed(do_fit)(id, field, templ0, templ1) for id in dis_matched['ID_Finkelstein'])


    os.chdir('/Users/rsimons/Desktop/git/clear_local')



if __name__ == '__main__':
    fields = ['GN1',
            'GN2',
            'GN3',
            'GN4',
            'GN5',
            'GN7',
            'GS1',
            'GS2',
            'GS3',
            'GS4',
            'GS5',
            'ERSPRIME']

    fields = ['GN2']
    #fields = ['GN1']

    #Parallel(n_jobs = -1, backend = 'threading')(delayed(run_all) (field = field) for field in fields)
    args = parse()
    field               = args['field']
    run_all(field)
    
    #for field in fields: run_all(field)
    #for field in fields:
    #    print (field)
    #    run_all(field)

















