def load_paper_catalog(cat_name = 'simons_sample.cat', cat_dir = '/Users/rsimons/Desktop/clear/catalogs'):
    print ('Loading catalog from paper sample...')
    import astropy.io
    from astropy.io import ascii
    cat = ascii.read(cat_dir + '/' + cat_name)
    return cat


def load_eazy_catalog(field_name = 'goodss', cat_version = '4.4', cat_dir = '/Users/rsimons/Desktop/clear/catalogs'):
    print ('Loading Eazy v%s catalog for %s...'%(cat_version, field_name))
    import astropy.io
    from astropy.io import fits
    cat_eazy_zout = fits.open(cat_dir + '/%s_3dhst.v%s.cats/Eazy/%s_3dhst.v%s.zout.fits'%(field_name, cat_version, field_name, cat_version))
    return cat_eazy_zout