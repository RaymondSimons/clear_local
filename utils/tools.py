def load_paper_catalog(add_Z = True, cat_version = 'v4.5', cat_name = 'simons_sample', cat_dir = '/Users/rsimons/Dropbox/clear/catalogs'):
    print ('Loading catalog from paper sample...')
    import astropy.io
    from astropy.io import ascii
    from astropy.table import join
    cat_sample = ascii.read(cat_dir + '/' + cat_name + '_' + cat_version + '.cat')
    if add_Z:
        cat_Z = ascii.read(cat_dir + '/' + 'metal_highZbranch_profile_fits.cat')
        cat = join(cat_sample, cat_Z)
    else: cat = cat_sample

    return cat


def load_eazy_catalog(field_name = 'goodss', cat_version = 'v4.5', cat_dir = '/Users/rsimons/Desktop/clear/catalogs'):
    print ('Loading Eazy v%s catalog for %s...'%(cat_version, field_name))
    import astropy.io
    from astropy.io import fits
    cat_eazy_zout = fits.open(cat_dir + '/%s_3dhst.%s.cats/Eazy/%s_3dhst.%s.zout.fits'%(field_name, cat_version, field_name, cat_version))
    return cat_eazy_zout