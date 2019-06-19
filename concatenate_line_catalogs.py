import astropy
from astropy.io import fits
import glob
from glob import glob



fls = glob('/Volumes/pegasus/clear/catalogs/grizli_v2.1_cats/*fits')

master_hdulist = []
prihdr = fits.Header()
prihdu = fits.PrimaryHDU(header=prihdr)    
master_hdulist.append(prihdu)


for fl in fls:
    field = fl.split('/')[-1].split('_')[0]

    a = fits.open(fl)

    hdu = fits.BinTableHDU(data = a[1].data, header = a[1].header, name = field)
    master_hdulist.append(hdu)


thdulist = fits.HDUList(master_hdulist)
thdulist.writeto('/Volumes/pegasus/clear/catalogs/grizli_v2.1_cats/master_cat.fits', overwrite = True)
