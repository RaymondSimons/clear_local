def load_paper_catalog(add_Z = True, cat_version = 'v4.5', cat_name = 'simons_sample', cat_dir = '/Users/rsimons/Dropbox/clear/catalogs'):
    print ('Loading catalog from paper sample...')
    import astropy.io
    from astropy.io import ascii
    from astropy.table import join
    cat_sample = ascii.read(cat_dir + '/' + cat_name + '_' + cat_version + '.cat')
    if add_Z:
        cat_Z = ascii.read(cat_dir + '/' + 'metal_highZbranch_profile_fits_v3.cat')
        cat = join(cat_sample, cat_Z)
    else: cat = cat_sample

    return cat


def load_eazy_catalog(field_name = 's', cat_version = 'v4.5', cat_dir = '/Users/rsimons/Dropbox/clear/catalogs'):
    print ('Loading Eazy %s catalog for GOODS-%s...'%(cat_version, field_name.upper()))
    import astropy.io
    from astropy.io import fits
    cat_eazy_zout = fits.open(cat_dir + '/goods%s_3dhst.%s.cats/Eazy/goods%s_3dhst.%s.zout.fits'%(field_name, cat_version, field_name, cat_version))
    return cat_eazy_zout

def load_grizli_catalog(field_name = 's', cat_version = 'v3.0', cat_dir = '/Users/rsimons/Dropbox/clear'):
    print ('Loading Grizli %s catalog for GOODS-%s...'%(cat_version, field_name.upper()))
    import astropy.io
    from astropy.io import fits
    
    cat_grizli_zout = fits.open(cat_dir + '/grizli_extractions_%s/grizli_%s_cats/GD%s_lines_grizli_master.fits'%(cat_version, cat_version, field_name.upper()))
    return cat_grizli_zout