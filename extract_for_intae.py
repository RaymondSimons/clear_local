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
from astropy.table import Table

PATH_TO_CATS = '/Users/rsimons/Desktop/clear/catalogs'
#PATH_TO_CATS = '/Users/rsimons/Dropbox/clear/catalogs' #change


def parse():
    '''
    Parse command line arguments
    ''' 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='''CLEAR grizli extractions.''')
    parser.add_argument('-field',       '--field',          default='GS1', help='field to extract')
    parser.add_argument('-reblot',      '--reblot',         action = "store_true", default = False, help = 'use psf extraction in fitting routine')
    parser.add_argument('-make_beams',  '--make_beams',         action = "store_true", default = False, help = 'use psf extraction in fitting routine')
    parser.add_argument('-do_fit',      '--do_fit',         action = "store_true", default = False, help = 'use psf extraction in fitting routine')

    args = vars(parser.parse_args())
    return args



class Pointing():
    def __init__(self, field):
        if 'N' in field.upper():
            self.pad = 200
            self.seg_map =  PATH_TO_CATS + '/GND_3DHST_plus_v2_seg.fits'
            self.catalog = PATH_TO_CATS + '/GND_3DHST_plus_v2_seg.cat'
            self.ref_image =  PATH_TO_CATS + '/goodsn-F105W-astrodrizzle-v4.4_drz_sci.fits'
            self.intae_cat = PATH_TO_CATS + '/z68_selected_in_clear_GN.txt'

        elif 'S' in field.upper():
            self.pad = 200
            self.seg_map =  PATH_TO_CATS + '/GSD_3DHST_plus_v2_seg.fits'
            self.catalog = PATH_TO_CATS + '/GSD_3DHST_plus_v2_seg.cat'
            self.ref_image =  PATH_TO_CATS + '/goodss-F105W-astrodrizzle-v4.3_drz_sci.fits' 
            self.intae_cat = PATH_TO_CATS + '/z68_selected_in_clear_GS.txt'

def run_all(args):
    field = args['field']
    PATH_TO_PREP        = '/Users/rsimons/Desktop/clear/grizli_extractions_v3/%s/Prep_z67'%field
    #PATH_TO_PREP        = '/Users/rsimons/Dropbox/clear/grizli_extractions_v3.0/test_intae/%s/Prep_z67'%field #change
    os.chdir(PATH_TO_PREP)

    print('Loading contamination models...')
    p = Pointing(field=field)
    intae_cat = ascii.read(p.intae_cat)
    all_grism_files = [fl.replace('.01.GrismFLT.fits', '_flt.fits') for fl in glob('*.01.GrismFLT.fits')]
    grp = GroupFLT(
        grism_files=all_grism_files, 
        direct_files=[], 
        ref_file = p.ref_image,
        seg_file = p.seg_map,
        catalog  = p.catalog,
        pad=p.pad, cpu_count = -1)


    if args['reblot']:
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
            pad=p.pad, cpu_count = -1)
    
    sc_table = Table.read(p.catalog, format='ascii.commented_header')
    for FLT in grp.FLTs:
        FLT.process_seg_file(p.seg_map)
        FLT.catalog = FLT.blot_catalog(sc_table, sextractor = True)
    

    if args['make_beams']:
        for ob in intae_cat:
            di_3dhst = ob['ID_3DHST']
            if (di_3dhst == -1) | (di_3dhst > 5400000): di = 5400000 + ob['ID_SF'] 
            else: di = di_3dhst
            print (di_3dhst, di)
            beams = grp.get_beams(di, size=80)
            if beams != []:
                print ('Creating beams for:', field, di)
                mb = grizli.multifit.MultiBeam(beams, fcontam=0.2, group_name=field)
                mb.write_master_fits()            


    if args['do_fit']:
        templ0 = grizli.utils.load_templates(fwhm=1200, line_complexes=True, stars=False, 
                                             full_line_list=None,  continuum_list=None, 
                                             fsps_templates=True)

        # Load individual line templates for fitting the line fluxes
        templ1 = grizli.utils.load_templates(fwhm=1200, line_complexes=False, stars=False, 
                                             full_line_list=None, continuum_list=None, 
                                             fsps_templates=True)

        def do_fit(id, field, templ0, templ1,\
                   wave = np.linspace(2000,2.5e4,100),\
                   fcontam = 0.2, \
                   pline = {'kernel': 'point', 'pixfrac': 0.2, 'pixscale': 0.1, 'size': 16, 'wcs': None}):
            if os.path.exists(field + '_' + '%i.full.fits'%id): return
            print (field, id)            
            mb = grizli.multifit.MultiBeam(field + '_' + '%i.beams.fits'%id, fcontam=fcontam, group_name=field)

            print ('creating poly_templates...')
            poly_templates = grizli.utils.polynomial_templates(wave=wave, order=7,line=False)
            pfit = mb.template_at_z(z=0, templates=poly_templates, fit_background=True, fitter='lstsq', fwhm=1400, get_uncertainties=2)
            print ('drizzle_grisms_and_PAs...')


            hdu, fig = mb.drizzle_grisms_and_PAs(size=32, fcontam=fcontam, flambda=False, scale=1, 
                                                pixfrac=0.5, kernel='point', make_figure=True, usewcs=False, 
                                                zfit=pfit,diff=True)
            # Save drizzled ("stacked") 2D trace as PNG and FITS
            fig.savefig('{}_diff_{}.stack.png'.format(field, id))
            hdu.writeto('{}_diff_{}.stack.fits'.format(field, id), clobber=True)

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

        all_beams = [int(fl.split('_')[-1].replace('.beams.fits', '')) for fl in glob('*beams.fits')]
        #print (all_beams)
        #Parallel(n_jobs = -1)(delayed(do_fit)(di, field, templ0, templ1) for di in all_beams)
        for di in all_beams:
            do_fit(di, field, templ0, templ1)


    os.chdir('/Users/rsimons/Dropbox/git/clear_local')



if __name__ == '__main__':
    args = parse()
    run_all(args)


























