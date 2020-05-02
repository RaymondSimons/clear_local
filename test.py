import grizli
from grizli.pipeline import auto_script
from grizli.pipeline.auto_script import get_yml_parameters
import glob
from grizli import utils, fitting, multifit, prep
import os

kwargs = get_yml_parameters()
ge_path = '/user/rsimons/grizli_extractions'



args = [('GS1','j033300m2742'),
        ('GS2','j033232m2742'),
        ('GS3','j033236m2744'),
        ('GS5','j033228m2743'),
        ('GN1','j123716p6222'),
        ('GN2','j123652p6215'),
        ('GN3','j123700p6219'),
        ('GN4','j123712p6216'),
        ('GN5','j123640p6219'),
        ('GN7','j123632p6213'),
        ('GS4','j033236m2748'),
        ('ERSPRIME','j033216m2743')]


for (field, root) in args:

    print (field, root)

    HOME_PATH = ge_path + '/' + field


    visit_prep_args = kwargs['visit_prep_args']
    path_to_prep = HOME_PATH + '/' + root + '/Prep'
    os.chdir(path_to_prep)


    ## IR filters
    if 'fix_stars' in visit_prep_args:
        fix_stars = visit_prep_args['fix_stars']
    else:
        fix_stars = False

    mosaic_args = kwargs['mosaic_args']
    mosaic_pixfrac = mosaic_args['mosaic_pixfrac']
    visits, all_groups, info = auto_script.parse_visits(field_root=root, 
                                                        HOME_PATH=HOME_PATH, 
                                                        **kwargs['parse_visits_args'])
    wcs_ref_file = '{0}_wcs-ref.fits'.format(root)
    auto_script.make_reference_wcs(info, output=wcs_ref_file, get_hdu=True, 
                           **mosaic_args['wcs_params'])
    combine_all_filters=True

    all_filters = mosaic_args['ir_filters'] + mosaic_args['optical_filters']
    auto_script.drizzle_overlaps(root, 
                             filters=all_filters, 
                             min_nexp=1, pixfrac=mosaic_pixfrac,
                             make_combined=True,
                             ref_image=wcs_ref_file,
                             drizzle_filters=False) 


    auto_script.drizzle_overlaps(root, filters=mosaic_args['ir_filters'], 
                                 min_nexp=1, pixfrac=mosaic_pixfrac,
                                 make_combined=(not combine_all_filters),
                                 ref_image=wcs_ref_file,
                                 include_saturated=fix_stars) 
    mask_spikes=True

    ir_mosaics = glob.glob('{0}-f*drz_sci.fits'.format(root))
    if (len(ir_mosaics) > 0) & (mask_spikes):
        cat = prep.make_SEP_catalog('{0}-ir'.format(root), threshold=4, 
                                    save_fits=False, 
                                    column_case=str.lower)

        selection = (cat['mag_auto'] < 17) & (cat['flux_radius'] < 4.5)
        for visit in visits:
            filt = visit['product'].split('-')[-1]
            if filt[:2] in ['f0','f1']:
                auto_script.mask_IR_psf_spikes(visit=visit, selection=selection,
                                   cat=cat, minR=5, dy=5)

        ## Remake mosaics
        auto_script.drizzle_overlaps(root, filters=mosaic_args['ir_filters'], 
                                     min_nexp=1, pixfrac=mosaic_pixfrac,
                                make_combined=(not combine_all_filters),
                                     ref_image=wcs_ref_file,
                                     include_saturated=True) 

    # Fill IR filter mosaics with scaled combined data so they can be used 
    # as grism reference
    fill_mosaics = mosaic_args['fill_mosaics']
    if fill_mosaics:
        if fill_mosaics == 'grism':
            # Only fill mosaics if grism filters exist
            has_grism = utils.column_string_operation(info['FILTER'], 
                                     ['G141','G102','G800L'],
                                     'count', 'or').sum() > 0
            if has_grism:
                auto_script.fill_filter_mosaics(root)                                             
        else:
            auto_script.fill_filter_mosaics(root)


    mosaics = glob.glob('{0}-ir_dr?_sci.fits'.format(root))
    wcs_ref_optical = wcs_ref_file

    preprocess_args = kwargs['preprocess_args']
    auto_script.drizzle_overlaps(root, 
            filters=mosaic_args['optical_filters'], 
            pixfrac=mosaic_pixfrac,
            make_combined=(len(mosaics) == 0), ref_image=wcs_ref_optical,
            min_nexp=1+preprocess_args['skip_single_optical_visits']*1)


    multiband_catalog_args=kwargs['multiband_catalog_args']
    tab = auto_script.multiband_catalog(field_root=root,
                                    **multiband_catalog_args)







