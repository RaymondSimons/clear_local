import numpy as np
from numpy import *
import astropy
from astropy.io import fits
gals = ['VELA%.2i'%i for i in arange(36)]


for g, gal in enumerate(gals):
    print g, gal