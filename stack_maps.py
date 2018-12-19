import glob
from glob import glob
import astropy
from astropy.io import fits



if __name__ == '__main__':
    prep_dir = '/Volumes/gdrive/clear/grizli_extractions/*/*/Prep'
    full_fls = glob(prep_dir + '/*full.fits')
    if False:
        all_full = []
        for f, full_fl in enumerate(full_fls):
            print f, len(full_fls)
            full = fits.open(full_fl)
            if 'Ha' in full[0].header['HASLINES']:
                all_full.append(full['LINE', 'Ha'].data)
    if True:
        all_full = array(all_full)
        all_full[all_full == 0.] = nan
        ha_median = nanmedian(all_full, axis = 0)
        clf()
        imshow(ha_median, interpolation = 'nearest', cmap = 'Greys_r')