#!/home/rsimons/miniconda2/envs/grizli/bin/python
import time
from multiprocessing import Pool
import os
import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
import glob
import importlib
import os
from astropy.table import Table
from matplotlib.colors import LogNorm
from IPython.display import Image
from numpy import *
import photutils
from astropy.cosmology import Planck15 as cosmo
from matplotlib.backends.backend_pdf import PdfPages
import glob
from glob import glob
from astropy.convolution import Gaussian2DKernel, convolve_fft, Box2DKernel
from scipy.interpolate import interp1d
import joblib
from joblib import Parallel, delayed
from astropy.stats import sigma_clip
import emcee
import scipy.optimize as op
import time
import emcee

mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 14
from sys import argv


plt.ioff()

np.random.seed(1)
prt = False


def OH_R23(OH, use = 'M08'):
    #taken from table in Patricio+
    if use == 'M08': 
        x = OH - 8.69
        if (OH > 9.2) | (OH < 7.0):
            if prt: print ('R23 outside of Maiolino+ 08 calibration range...')
        c = [0.7462, -0.7149, -0.9401, -0.6154, -0.2524]
    if use == 'J15': 
        x = OH
        if (OH > 9.0) | (OH < 7.6):
            if prt: print ('R23 outside of Jones+ 15 calibration range...')
        c = [-54.1003, 13.9083, -0.8782, 0.0, 0.0]
    if use == 'C17': 
        x = OH - 8.69
        if (OH > 8.85) | (OH < 8.4):
            if prt: print ('R23 outside of Curti+ 17 calibration range...')
        c = [0.527, -1.569, -1.652, -0.421, 0.0]
    if use == 'B18': 
        x = OH
        if (OH > 8.4) | (OH < 7.8):
            if prt: print ('R23 outside of Bian+ 18 calibration range...')
        c = [138.0430, -54.8284, 7.2954, -0.32293, 0.0]

    result = c[0] * x**0. + c[1] * x**1. + c[2] * x**2. + c[3] * x**3. + c[4] * x**4.
    return result

def OH_R2(OH, use = 'C17'):
    #taken from table in Patricio+
    if use == 'C17': 
        x = OH - 8.69
        if (OH > 8.3) | (OH < 7.6):
            if prt: print ('R2 outside of Curti+ 17 calibration range...')
        c = [0.418, -0.961, -3.505, -1.949, 0.0]
    result = c[0] * x**0. + c[1] * x**1. + c[2] * x**2. + c[3] * x**3. + c[4] * x**4.
    return result

def OH_R3(OH, use = 'C17'):
    #taken from table in Patricio+
    if use == 'C17': 
        x = OH - 8.69
        if (OH > 8.85) | (OH < 8.3):
            if prt: print ('R2 outside of Curti+ 17 calibration range...')
        c = [-0.277 -3.549 -3.593 -0.981, 0.0]
    result = c[0] * x**0. + c[1] * x**1. + c[2] * x**2. + c[3] * x**3. + c[4] * x**4.
    return result

def OH_O3(OH, use = 'M08'):
    #taken from table in Patricio+
    if use == 'M08': 
        x = OH - 8.69
        if (OH > 9.2) | (OH < 7.0):
            if prt: print ('O3 outside of Maiolino+ 08 calibration range...')
        c = [0.1549, -1.5031, -0.9790, -0.0297, 0.0]
    if use == 'J15': 
        x = OH
        if (OH > 9.0) | (OH < 7.6):
            if prt: print ('O3 outside of Jones+ 15 calibration range...')
        c = [-88.4378, 22.7529, -1.4501, 0.0, 0.0]
    if use == 'C17': 
        x = OH - 8.69
        if (OH > 8.85) | (OH < 8.3):
            if prt: print ('O3 outside of Curti+ 17 calibration range...')
        c = [-0.277, -3.549, -3.593, -0.981, 0.0]
    if use == 'B18': 
        x = OH
        if (OH > 8.4) | (OH < 7.8):
            if prt: print ('O3 outside of Bian+ 18 calibration range...')
        c = [43.9836, -21.6211, 3.4277, -0.1747, 0.0]

    result = c[0] * x**0. + c[1] * x**1. + c[2] * x**2. + c[3] * x**3. + c[4] * x**4.
    return result

def OH_O2(OH, use = 'M08'):
    #taken from table in Patricio+
    if use == 'M08': 
        x = OH - 8.69
        if (OH > 9.2) | (OH < 7.0):
            if prt: print ('O2 outside of Maiolino+ 08 calibration range...')
        c = [0.5603, 0.0450, -1.8017, -1.8434, -0.6549]
    if use == 'J15': 
        x = OH
        if (OH > 9.0) | (OH < 7.6):
            if prt: print ('O2 outside of Jones+ 15 calibration range...')
        c = [-154.9571, 36.9128, -2.1921, 0.0, 0.0]
    if use == 'C17': 
        x = OH - 8.69
        if (OH > 8.3) | (OH < 7.6):
            if prt: print ('O2 outside of Curti+ 17 calibration range...')
        c = [0.418, -0.961, -3.505, -1.949, 0.0]

    result = c[0] * x**0. + c[1] * x**1. + c[2] * x**2. + c[3] * x**3. + c[4] * x**4.
    return result

def OH_O32(OH, use = 'M08'):
    #taken from table in Patricio+
    if use == 'M08': 
        x = OH - 8.69
        if (OH > 9.2) | (OH < 7.0):
            if prt: print ('O32 outside of Maiolino+ 08 calibration range...')
        c = [-0.2839, -1.3881, -0.3172, 0., 0.]
    if use == 'J15': 
        x = OH
        if (OH > 9.0) | (OH < 7.6):
            if prt: print ('O32 outside of Jones+ 15 calibration range...')
        c = [17.9828, -2.1552, 0.0, 0.0, 0.0]
    if use == 'C17': 
        x = OH - 8.69
        if (OH > 8.85) | (OH < 7.6):
            if prt: print ('O32 outside of Curti+ 17 calibration range...')
        c = [-0.691, -2.944, -1.308, 0.0, 0.0]

    result = c[0] * x**0. + c[1] * x**1. + c[2] * x**2. + c[3] * x**3. + c[4] * x**4.
    return result

def OH_S2(OH, use = 'C19'):
    #Curti+ 19 in prep
    return nan

def OH_O3S2(OH, use = 'C19'):
    #Curti+ 19 in prep
    return nan

def determine_R(OH, diagnostic):
    if diagnostic == 'R23':  return OH_R23(OH)
    if diagnostic == 'R2':   return OH_R2(OH)
    if diagnostic == 'R3':   return OH_R3(OH)
    if diagnostic == 'O3':   return OH_O3(OH)
    if diagnostic == 'O2':   return OH_O2(OH)
    if diagnostic == 'O32':  return OH_O32(OH)
    if diagnostic == 'S2':   return OH_S2(OH)
    if diagnostic == 'O3S2': return OH_O3S2(OH)

def lnlike(OH, R, Rerr, diagnostics):
    model = nan * zeros(len(R))
    for RR in arange(len(R)):
        model[RR] = determine_R(OH, diagnostics[RR])
    inv_sigma2 = 1.0/Rerr**2
    return -0.5*(np.sum((R-model)**2*inv_sigma2))

def lnprior(OH):
    if 7.0 < OH < 9.5: return 0.0
    return -np.inf

def lnprob(OH, R, Rerr, diagnostics):
    lp = lnprior(OH)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(OH, R, Rerr, diagnostics)


def write_fits():

    return




if __name__ == '__main__':

    print (argv[1], argv[2])

    
    np.random.seed()
    OH_true = 8.6
    Re1 = 0.3
    Re2 = 0.3
    R = array([OH_R23(OH_true) + np.random.normal(0, Re1), OH_O32(OH_true)+ np.random.normal(0, Re2)])
    Rerr = array([Re1, Re2])
    diagnostics = array(['R23', 'O32'])

    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, [OH_true], args=(R, Rerr, diagnostics))
    OH_ml = result["x"]

    ndim, nwalkers = 1, 100

    pos = [result["x"] + 1e-4*np.random.randn(1) for i in range(nwalkers)]

    a = time.time()
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(R, Rerr, diagnostics), pool = pool)
        sampler.run_mcmc(pos, 1000)       
        samples = sampler.chain[:, 100:, :].reshape((-1, ndim))
        OH_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))    
        print (list(OH_mcmc))
    b = time.time()
    print (b-a)
    
    a = time.time()
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(R, Rerr, diagnostics))
    sampler.run_mcmc(pos, 1000)       
    samples = sampler.chain[:, 100:, :].reshape((-1, ndim))

    OH_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))    
    print (list(OH_mcmc))
    b = time.time()
    print (b-a)






































