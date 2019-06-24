import astropy
from astropy.io import fits
import argparse
import glob
from glob import glob





if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--field', default = 'GS1')
    args = vars(parser.parse_args())

    field = args['field']

    fls = glob('/user/rsimons/metal_maps/%s_*_metals.fits')
    
    for fl in fls:
        make_metal_profile(fl = fl)
