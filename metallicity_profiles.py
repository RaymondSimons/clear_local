import astropy
from astropy.io import fits
import argparse
import glob
from glob import glob



def make_metal_profile(fl):
    print (fl)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--field', default = 'GS1')
    args = vars(parser.parse_args())

    field = args['field']

    fls = glob('/user/rsimons/metal_maps/%s_*_metals.fits')
    print (fls)
    
    for fl in fls:
        make_metal_profile(fl = fl)
