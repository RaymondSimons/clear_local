import astropy
from astropy.io import fits
import os
import glob
from glob import glob
from numpy import *
os.system('mkdir /user/rsimons/metal_maps_reduced')


fls = glob('/user/rsimons/metal_maps/*')

for f, fl in enumerate(fls):
    out_file = fl.replace('/metal_maps/', '/metal_maps_reduced/')
    a = fits.open(fl)
    new_hdulist = []
    for i in arange(len(a)):
        if '_FULL' not in a[i].name: new_hdulist.append(a[i])

    thdulist = fits.HDUList(new_hdulist)
    print ('Saving to {}'.format(out_file))
    thdulist.writeto(out_file, overwrite = True)
