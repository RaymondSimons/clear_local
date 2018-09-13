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
from hsaquery import query, overlaps
from grizli.pipeline import auto_script
from grizli.multifit import GroupFLT, MultiBeam, get_redshift_fit_defaults
import os
from astropy.table import Table
plt.ioff()
### A lot of this is working from the Grizli cookbook described here:
### https://github.com/gbrammer/grizli/blob/master/examples/WFC3IR_Reduction.ipynb
### and Katie and Iva's grizli pipeline for CLEAR

overlapping_fields = {'GN1':['GDN20'],
                      'GN2':['GDN8', 'GDN12', 'GDN21', 'GDN25'],
                      'GN3':['GDN18', 'GDN19', 'GDN22', 'GDN23'],
                      'GN4':['GDN21', 'GDN22', 'GDN25', 'GDN26'],
                      'GN5':['GDN17', 'GDN18'],
                      'GN7':['GDN3', 'GDN6', 'GDN7', 'GDN11'],
                      'ERSPRIME':['WFC3-ERSII-G01']}






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
            self.radec_catalog = PATH_TO_CATS + '/old_radeccats/goodsN_radec.cat'
            #self.radec_catalog = '../Catalogs/new_radeccats/goodsn_radec.cat'
            self.seg_map =  PATH_TO_CATS + '/Goods_N_plus_seg.fits'
            self.catalog =  PATH_TO_CATS + '/goodsn-F105W-astrodrizzle-v4.4_drz_sub_plus.cat'
            self.ref_image =  PATH_TO_CATS + '/goodsn-F105W-astrodrizzle-v4.4_drz_sci.fits'

            self.tempfilt, self.coeffs, self.temp_sed, self.pz = readEazyBinary(MAIN_OUTPUT_FILE='goodsn_3dhst.v4.1', OUTPUT_DIRECTORY=PATH_TO_CATS, CACHE_FILE='Same')


            '''
            if '140' in ref_filter:
                self.catalog =  PATH_TO_CATS + '/GoodsN_plus_merged.cat'
                self.ref_image =  PATH_TO_CATS + '/goodsn_3dhst.v4.0.F125W_orig_sci.fits'
            elif '105' in ref_filter:
                #self.catalog = '../Catalogs/goodsn-F105W-astrodrizzle-v4.3_drz_sub_plus.cat'
                self.catalog =  PATH_TO_CATS + '/goodsn-F105W-astrodrizzle-v4.4_drz_sub_plus.cat'
                self.ref_image =  PATH_TO_CATS + '/goodsn-F105W-astrodrizzle-v4.4_drz_sci.fits'
                #self.ref_image = '../Catalogs/goodsn-F105W-astrodrizzle-v4.3_drz_sci.fits'
            '''           


        elif 'S' in field.upper():
            self.pad = 200 # grizli default
            #self.radec_catalog = '../Catalogs/goodsS_radec.cat'
            self.radec_catalog =  PATH_TO_CATS + '/goodss_3dhst.v4.1.radec.cat'
            self.seg_map =  PATH_TO_CATS + '/Goods_S_plus_seg.fits'
            if '140' in ref_filter:
                self.catalog =  PATH_TO_CATS + '/GoodsS_plus_merged.cat'
                self.ref_image =  PATH_TO_CATS + '/goodss_3dhst.v4.0.F125W_orig_sci.fits'  
            elif '105' in ref_filter:            
                self.catalog =  PATH_TO_CATS + '/goodss-F105W-astrodrizzle-v4.3_drz_sub_plus.cat'
                self.ref_image =  PATH_TO_CATS + '/goodss-F105W-astrodrizzle-v4.3_drz_sci.fits'
                #self.ref_image = '../Catalogs/goodss_3dhst.v4.0.F125W_orig_sci.fits'  
                #RCS need to retrieve the F015W image
            



def grizli_getfiles(run = True):
    if run == False: return
    else: 'Running grizli_getfiles...'

    os.chdir(PATH_TO_PREP)
    files = glob.glob('%s/*flt.fits'%PATH_TO_RAW)
    info = grizli.utils.get_flt_info(files)
    visits, filters = grizli.utils.parse_flt_files(info=info, uniquename=True)
    return visits, filters

def grizli_prep(visits, ref_filter = 'F105W', ref_grism = 'G102', field = 'GN2', run = True):
    if run == False: return
    else: 'Running grizli_prep...'

    print ('\n\n\n\n\n\n\n')
    product_names = np.array([visit['product'] for visit in visits])
    filter_names = np.array([visit['product'].split('-')[-1] for visit in visits])
    basenames = np.array([visit['product'].split('.')[0]+'.0' for visit in visits])

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
            grism_index= np.where((basenames == basename) & (filter_names == ref_grism.lower()))[0][0]
            #print(grism_index)
            p = Pointing(field = field, ref_filter = ref_filter)
            radec_catalog = p.radec_catalog
            print (field_in_contest, visits[grism_index])
            #print (visit[grism_index])
            #try:
            #status = process_direct_grism_visit(direct = visit,
            #                                    grism = visits[grism_index],
            #                                    radec = radec_catalog, 
            #                                    align_mag_limits = [14, 23])

    return visits, filters



def grizli_model(visits, field = 'GN2', ref_filter_1 = 'F105W', ref_grism_1 = 'G102', ref_filter_2 = 'F140W', ref_grism_2 = 'G141', run = True, load_only = True, mag_lim = 25):
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
        #filter_name = visit['product'].split('-')[-1]
        field_in_contest = visit['product'].split('-')[0].upper()
        if field_in_contest.upper() != 'GOODSN':
            #if field_in_contest == field or field_in_contest in overlapping_fields[field]:
            if (ref_filter_1.lower() in filt1) or (ref_filter_2.lower() in filt1):
                #Find grism files with a direct image
                all_direct_files.extend(visit['files'])
                grism_index_1= np.where((basenames == basename) & (filter_names == ref_grism_1.lower()))[0]
                grism_index_2= np.where((basenames == basename) & (filter_names == ref_grism_2.lower()))[0]

                if len(grism_index_1) > 0:
                    all_grism_files.extend(visits[grism_index_1[0]]['files'])
                    print(filter_names[grism_index_1[0]], visits[grism_index_1[0]]['product'])

                elif len(grism_index_2) > 0:
                    all_grism_files.extend(visits[grism_index_2[0]]['files'])
                    print(filter_names[grism_index_2[0]], visits[grism_index_2[0]]['product'])

            '''
            elif (ref_grism_1.lower() in filter_name) or (ref_grism_2.lower() in filter_name):
                all_grism_files.extend(visit['files'])
                print (field_in_contest, visit['files'], filter_name)
            '''

    #print (all_direct_files, all_grism_files)
    p = Pointing(field=field, ref_filter=ref_filter_1)
    if load_only:
        print('Loading contamination models...')

    grp = GroupFLT(
        grism_files=all_grism_files, 
        direct_files=[], 
        ref_file = p.ref_image,
        seg_file = p.seg_map,
        catalog  = p.catalog,
        pad=p.pad,
        cpu_count=8)

    if not load_only:
        print('Computing contamination models...')
        grp.compute_full_model(mag_limit=mag_lim)
        print('Refining List..')
        grp.refine_list(poly_order=2, mag_limits=[16, 24], verbose=False)
        print('Saving contamination models')
        grp.save_full_data()
    return grp
        
def grizli_fit(grp, field = '', mag_lim = 35, mag_lim_lower = 35, run = True, id_choose = None, ref_filter = 'F105W'):
    if fit_bool == False: return
    p = Pointing(field = field, ref_filter = ref_filter)

    templ0 = grizli.utils.load_templates(fwhm=1200, line_complexes=True, stars=False, 
                                         full_line_list=None,  continuum_list=None, 
                                         fsps_templates=True)

    # Load individual line templates for fitting the line fluxes
    templ1 = grizli.utils.load_templates(fwhm=1200, line_complexes=False, stars=False, 
                                         full_line_list=None, continuum_list=None, 
                                         fsps_templates=True)

    pline = {'kernel': 'point', 'pixfrac': 0.2, 'pixscale': 0.1, 'size': 8, 'wcs': None}
    for id, mag in zip(np.array(grp.catalog['NUMBER']), np.array(grp.catalog['MAG_AUTO'])):
        #if (mag <= mag_lim) & (mag >=mag_lim_lower) & (id > id_choose):
        if id == id_choose:
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


                #try:
                    if True:
                        #use redshift prior from z_phot
                        prior = np.zeros((2, len(p.tempfilt['zgrid'])))
                        prior[0] = p.tempfilt['zgrid']
                        prior[1] = p.pz['chi2fit'][:,id]
                    else:
                        prior = None 

                    out = grizli.fitting.run_all(
                        id, 
                        t0=templ0, 
                        t1=templ1, 
                        fwhm=1200, 
                        zr=[0.5, 2.3], 
                        dz=[0.004, 0.0005], 
                        fitter='nnls',
                        group_name=field,
                        fit_stacks=True, 
                        prior=prior, 
                        fcontam=0.,
                        pline=pline, 
                        mask_sn_limit=7, 
                        fit_only_beams=False,
                        fit_beams=True, 
                        root=field,
                        fit_trace_shift=False, 
                        phot=None, 
                        verbose=True, 
                        scale_photometry=False, 
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
                #except:
                    #print ('Problem in fitting.run_all')

                    plt.close('all')




def retrieve_archival_data(visits, field, retrieve_bool = False):
    if retrieve_bool == False: return
    os.chdir(PATH_TO_RAW)    

    product_names = np.array([visit['product'] for visit in visits])
    filter_names = np.array([visit['product'].split('-')[-1] for visit in visits])
    basenames = np.array([visit['product'].split('.')[0]+'.0' for visit in visits])
    print ('\n\n\n\n\n\n\n')
    for v, visit in enumerate(visits):
        product = product_names[v]
        basename = basenames[v]
        filt1 = filter_names[v]

        field_in_contest = basename.split('-')[0]
        filt1 = filter_names[v]

        print(field_in_contest)
        #if field_in_contest.upper() == field.upper() or field_in_contest.upper() in overlapping_fields[field]:





    ra_target = 189.22053
    dec_target = 62.24010833333
    radius_in_arcmin = 10


    #First run-through, ignore the imaging
    if False:
        parent = query.run_query(box=[ra_target, dec_target, radius_in_arcmin],instruments=['WFC3-IR', 'ACS-WFC'], 
                             extensions=['FLT'], filters=['G102', 'G141'], extra=[])
        tabs = overlaps.find_overlaps(parent, buffer_arcmin=0.01, filters=['G102', 'G141'], instruments=['WFC3-IR','WFC3-UVIS','ACS-WFC'], extra=[], close=False)
        s3_status = os.system('aws s3 ls s3://stpubdata --request-payer requester')
        HOME_PATH = os.getcwd()
        auto_script.fetch_files(field_root='j123625+621431', HOME_PATH=HOME_PATH, remove_bad=True, reprocess_parallel=True, s3_sync=(s3_status == 0))

    #Second run-through, retrieve the direct imaging
    if True:
        #Find targetnames
        fls_temp = glob.glob(PATH_TO_RAW+'/j123625+621431/RAW/*flt.fits')
        target_names = []
        for fl in fls_temp:
            data_temp = fits.open(fl)
            target_name = data_temp[0].header['TARGNAME']
            if target_name in target_names:
                pass
            else:
                target_names.append(target_name)

        print(target_names)
        for t, target_name in enumerate(target_names):
            parent = query.run_query(box=[ra_target, dec_target, radius_in_arcmin],instruments=['WFC3-IR', 'ACS-WFC'], extensions=['FLT'], filters=['F098M', 'F105W', 'F105W', 'F140W'], extra=[])
            extra = query.DEFAULT_EXTRA.copy()
            extra += ["TARGET.TARGET_NAME LIKE '%s'"%target_name]
            tabs = overlaps.find_overlaps(parent, buffer_arcmin=0.01, filters=['F098M', 'F105W', 'F125W', 'F140W'], instruments=['WFC3-IR','WFC3-UVIS','ACS-WFC'], extra=extra, close=False)
            s3_status = os.system('aws s3 ls s3://stpubdata --request-payer requester')
            HOME_PATH = os.getcwd()
            auto_script.fetch_files(field_root='j123625+621431', HOME_PATH=HOME_PATH, remove_bad=True, 
                                    reprocess_parallel=True, s3_sync=(s3_status == 0))



    os.chdir(PATH_TO_PREP)    




if __name__ == '__main__':

    global PATH_TO_RAW, PATH_TO_PREP, PATH_TO_SCRIPTS

    PATH_TO_RAW = '/user/rsimons/grizli_extractions/RAW'
    PATH_TO_PREP = '/user/rsimons/grizli_extractions/PREP'
    PATH_TO_SCRIPTS = '/home/rsimons/git/clear_local'
    PATH_TO_CATS= '/user/rsimons/grizli_extractions/Catalogs'

    os.chdir(PATH_TO_PREP)
    mag_lim = 25
    mag_lim_lower = 0

    id_choose = 23116
    if True:
        files_bool = True
        retrieve_bool = False
        prep_bool = False
        model_bool = True
        load_bool = True
        fit_bool = True

    #field = 'GN4'
    #for field in overlapping_fields.keys():

    for field in ['GN2']:
        visits, filters = grizli_getfiles(run = files_bool)
        extra = retrieve_archival_data(visits = visits, field = field, retrieve_bool = retrieve_bool)
        visits, filters = grizli_getfiles(run = files_bool)
        #grizli_prep(visits = visits, ref_filter = 'F140W', ref_grism = 'G141', run = prep_bool)
        grizli_prep(visits = visits, ref_filter = 'F105W', ref_grism = 'G102', run = prep_bool)
        grp = grizli_model(visits, field = field, ref_filter_1 = 'F105W', ref_grism_1 = 'G102', ref_filter_2 = 'F140W', ref_grism_2 = 'G141',
                           run = model_bool, load_only = load_bool, mag_lim = mag_lim)
        grizli_fit(grp, field = field, mag_lim = mag_lim, mag_lim_lower = mag_lim_lower, run = fit_bool, id_choose = 21706)
    os.chdir(PATH_TO_SCRIPTS)
























