def load_paper_catalog(cat_name = 'simons_sample.cat', cat_dir = '/Users/rsimons/Desktop/clear/catalogs'):
    print ('Loading catalog from paper sample...')
    import astropy.io
    from astropy.io import ascii
    cat = ascii.read(cat_dir + '/' + cat_name)
    return cat