#!/home/rsimons/miniconda2/envs/grizli/bin/python
import time
from multiprocessing import Pool
import os
import numpy as np
import matplotlib.pyplot as plt
import astropy
from astropy.io import fits
import matplotlib as mpl
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
from math import *
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



def run_mcmc(pos, R, eR, diagnostics, Nsteps = 300, ndim = 1, nwalkers = 100):
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(R, eR, diagnostics))
    a = time.time()
    sampler.run_mcmc(pos, Nsteps)
    b = time.time()
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

    OH_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))    
    print (diagnostics, list(OH_mcmc))


if __name__ == '__main__':
    kern = Box2DKernel(2)

    np.random.seed()
    field, di = argv[1], argv[2]
    fl = glob('/user/rsimons/grizli_extractions/%s/j*/Prep/*%s.full.fits'%(field, di))[0]

    if os.path.isfile(fl):
        a = fits.open(fl)
        O3  = a['LINE', 'OIII'].data
        eO3 = 1./np.sqrt(a['LINEWHT', 'OIII'].data)
        O2  = a['LINE', 'OII'].data
        eO2 = 1./np.sqrt(a['LINEWHT', 'OII'].data)

        Hb  = a['LINE', 'Hb'].data
        eHb = 1./np.sqrt(a['LINEWHT', 'Hb'].data)


        O3 = convolve_fft(O3, kern)
        O2 = convolve_fft(O2, kern)
        Hb = convolve_fft(Hb, kern)

        R_O32 = O3/O2
        eR_O32 = R_O32 * np.sqrt((eO3/O3)**2. + (eO2/O2)**2.)

        R_R23 = (O3 + O2)/Hb
        eR_R23 = R_R23 * np.sqrt((eO3**2. + eO2**2.)/(O3 + O2)**2. + (eHb/Hb)**2.)


        Rs  = [R_O32, R_R23]
        eRs = [eR_O32, eR_R23]

        nwalkers = 100


        for i in arange(shape(O3)[0]):
            for j in arange(shape(O3)[1]):
                for d, diagnostic in enumerate(array([['O32'], ['O32', 'R23']])):
                    if d == 0: 
                        Rs = array([R_O32[i,j]])
                        eRs = array([eR_O32[i,j]])
                    if d == 1: 
                        Rs = array([R_O32[i,j], R_R23[i,j]])
                        eRs = array([eR_O32[i,j], eR_R23[i,j]])

                    if (O3[i,j]/eO3[i,j] > 0.5) & (O2[i,j]/eO2[i,j] > 0.5) & (Hb[i,j]/eHb[i,j] > 0.5) :
                        nll = lambda *args: -lnlike(*args)

                        result = op.minimize(nll, [8.5], args=(Rs, eRs, diagnostic))
                        OH_ml = result["x"]
                        pos = [result["x"] + 1e-4*np.random.randn(1) for nn in range(nwalkers)]
                        run_mcmc(pos = pos, R = Rs, eR = eRs, 
                                 diagnostics = diagnostic, nwalkers = nwalkers)



    #OH_true = 8.6
    #Re1, Re2 = 0.3, 0.3
    #R = array([OH_R23(OH_true) + np.random.normal(0, Re1), OH_O32(OH_true)+ np.random.normal(0, Re2)])    
    #Rerr = array([Re1, Re2])
    #diagnostics = array(['R23', 'O32'])
    '''
    nll = lambda *args: -lnlike(*args)
    result = op.minimize(nll, [OH_true], args=(R, Rerr, diagnostics))
    OH_ml = result["x"]
    nwalkers = 100
    pos = [result["x"] + 1e-4*np.random.randn(1) for i in range(nwalkers)]
    Ntotal = 10
    print ('Running Parallel')
    print ('_________________')
    c = time.time()
    Parallel(n_jobs = -1)(delayed(run_mcmc)(pos = pos, R = R, Rerr = Rerr, diagnostics = diagnostics, nwalkers = nwalkers) for i in arange(Ntotal))
    d = time.time()
    print ('Total parallel: %.2f s\n\n'%(d - c))

    c = time.time()




    print ('Running Serial')
    print ('_________________')
    for i in arange(Ntotal):
        run_mcmc(pos = pos, R = R, Rerr = Rerr, diagnostics = diagnostics)
    d = time.time()
    print ('Total serial: %.2f s '%(d - c))
    '''
























