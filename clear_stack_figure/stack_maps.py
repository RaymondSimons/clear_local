import glob
from glob import glob
import astropy
from astropy.io import fits



if __name__ == '__main__':
    prep_dir = '/Volumes/gdrive/clear/grizli_extractions/*/*/Prep'
    full_fls = glob(prep_dir + '/*full.fits')
    if False:
        all_full_Ha = []
        all_full_OII = []
        all_full_OIII = []
        all_full_Hb = []

        for f, full_fl in enumerate(full_fls):
            print f, len(full_fls)
            full = fits.open(full_fl)
            '''
            if 'Ha' in full[0].header['HASLINES']:
                all_full_Ha.append(full['LINE', 'Ha'].data)

            if 'Hb' in full[0].header['HASLINES']:
                all_full_Hb.append(full['LINE', 'Hb'].data)
            '''
            if 'OII ' in full[0].header['HASLINES']:
                all_full_OII.append(full['LINE', 'OII'].data)
            elif "OII'" in full[0].header['HASLINES']:
                all_full_OII.append(full['LINE', 'OII'].data)

            if 'OIII ' in full[0].header['HASLINES']:
                all_full_OIII.append(full['LINE', 'OIII'].data)
            elif "OIII'" in full[0].header['HASLINES']:
                all_full_OIII.append(full['LINE', 'OIII'].data)

            full.close()

    if True:
        #all_full_OIII = array(all_full_OIII)
        #all_full_OIII[all_full_OIII == 0.] = nan
        #ha_median = nanmedian(all_full_OIII, axis = 0)
        clf()
        imshow(ha_median, interpolation = 'nearest', cmap = 'Greys_r')