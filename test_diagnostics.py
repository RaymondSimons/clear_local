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
import metal_calibs as calib
from joblib import Parallel, delayed
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


def generate_R(OH, diagnostic):
    if diagnostic == 'R23':  
        return calib.OH_R23(OH)
    if diagnostic == 'O32':  
        return calib.OH_O32(OH)


def determine_R(OH, diagnostic):
    if diagnostic == 'R23':  return calib.OH_R23(OH)
    if diagnostic == 'R2':   return calib.OH_R2(OH)
    if diagnostic == 'R3':   return calib.OH_R3(OH)
    if diagnostic == 'O3':   return calib.OH_O3(OH)
    if diagnostic == 'O2':   return calib.OH_O2(OH)
    if diagnostic == 'O32':  return calib.OH_O32(OH)
    if diagnostic == 'S2':   return calib.OH_S2(OH)
    if diagnostic == 'O3S2': return calib.OH_O3S2(OH)

def lnlike(OH, R, Rerr, diagnostics):
    model = nan * zeros(len(R))
    for RR in arange(len(R)):
        model[RR] = determine_R(OH, diagnostics[RR])
    inv_sigma2 = 1.0/Rerr**2
    return -0.5*(np.sum((R-model)**2*inv_sigma2))

def lnprior(OH, OHmin = 7.0, OHmax = 9.5):
    if OHmin < OH < OHmax: return 0.0
    return -np.inf

def lnGaussianprior(OH, OHmin = 4., OHmax = 12.):
    # Gaussian prior on m
    # uniform prior on c
    # OHmin lower range of prior
    # OHmax upper range of prior
    # set prior to 1 (log prior to 0) if in the range and zero (-inf) outside the range 
    lp = 0. if OHmin < OH < OHmax else -np.inf
    OHmu = 8.0     # mean of the Gaussian prior
    OHsigma = 2.0 # standard deviation of the Gaussian prior
    lp -= 0.5*((OH - OHmu)/OHsigma)**2   
    return lp

def lnprob(OH, R, Rerr, diagnostics, use_prior = 'gaussian'):
    if use_prior == 'gaussian': lp = lnGaussianprior(OH)
    if use_prior == 'top': lp = lnprior(OH)

    if not np.isfinite(lp): return -np.inf
    return lp + lnlike(OH, R, Rerr, diagnostics)


def run_mcmc(pos, R, eR, diagnostics, Nsteps = 300, Nburn = 50, Ndim = 1, Nwalkers = 100, use_prior = 'gaussian'):
    sampler = emcee.EnsembleSampler(Nwalkers, Ndim, lnprob, args=(R, eR, diagnostics, use_prior))
    sampler.run_mcmc(pos, Nsteps)
    samples = sampler.chain[:, Nburn:, :].reshape((-1, Ndim))
    OH_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))    
    return list(OH_mcmc), samples



def test_diag(OH_true):
    Nwalkers = 200
    Nsteps = 2000.
    Nburn = 50
    Ndim = 1
    SN = 100.
    np.random.seed()
    diagnostics = array([['O32'], 
                         ['R2'],
                         ['R3'],                         
                         ['R23'],
                         ['O32', 'R23'],
                         ['O32', 'R2', 'R3', 'R23']
                         ])

    Z = nan * zeros((len(diagnostics), 3))

    plt.close('all')
    fig = plt.figure(figsize = (len(diagnostics)*5, 5))

    R_dict = {}
    R_dict['eO32'] = 0.4
    R_dict['eR23'] = 0.1
    R_dict['eR2'] = 0.1
    R_dict['eR3'] = 0.1
    

    R_dict['O32'] = generate_R(OH_true, diagnostic = 'O32')# + np.random.normal(0, R_dict['eO32'])
    R_dict['R2'] = generate_R(OH_true, diagnostic = 'R2')# + np.random.normal(0, R_dict['eO32'])
    R_dict['R3'] = generate_R(OH_true, diagnostic = 'R3')# + np.random.normal(0, R_dict['eO32'])
    R_dict['R23'] = generate_R(OH_true, diagnostic = 'R23')# + np.random.normal(0, R_dict['eR23'])


    for d, diagnostic in enumerate(diagnostics):
        R = array([R_dict[diag] for diag in diagnostic])
        eR = array([R_dict['e'+diag] for diag in diagnostic])

        print ('calculating metallicity using ', diagnostic)
        #Z = nan * zeros((shape(Rs)[1], shape(Rs)[2], 3))
        nll = lambda *args: -lnlike(*args)
        result = op.minimize(nll, [8.0], args=(R, eR, diagnostic))
        OH_ml = result["x"]
        pos = [result["x"] + 1e-4*np.random.randn(1) for nn in range(Nwalkers)]

        OH_result, samples = run_mcmc(pos = pos, R = R, eR = eR, 
                             diagnostics = diagnostic, Nsteps = Nsteps, 
                             Nburn = Nburn, Ndim = Ndim, Nwalkers = Nwalkers)

        ax = fig.add_subplot(1,len(diagnostics),d+1)

        ax.hist(samples[50:,0], bins = linspace(6.0, 10.0, 200), histtype = 'step')

        ax.axvline(OH_true, color = 'black')
        Z[d,0]  = OH_result[0][0]
        Z[d,1]  = OH_result[0][1]
        Z[d,2]  = OH_result[0][2]

        print  (OH_result[0][0],  OH_result[0][1], OH_result[0][2])
        ax.set_title(diagnostic, fontsize = 15)
        ax.set_xlabel(r'12 + $\log$(O/H)', fontsize = 20)
        fit_str = ('%s'%diagnostic).strip('[').strip(']').strip("'").replace("', '", '_')
        fits.writeto('/Users/rsimons/Desktop/clear/diagnostics/out/%.3f_%s_samples.fits'%(OH_true, fit_str), 
                      samples, 
                      overwrite = True)
    fig.tight_layout()
    fig.savefig('/Users/rsimons/Desktop/clear/diagnostics/figures/%.3f_diagnostics.png'%(OH_true))






if __name__ == '__main__':
    Parallel(n_jobs = -1)(delayed(test_diag)(OH_true) for OH_true in arange(7, 9.5, 0.1))























































