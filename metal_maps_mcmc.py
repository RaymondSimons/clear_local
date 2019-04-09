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
        cf = [-0.277, -3.549, -3.593, -0.981, 0.0]
    result = cf[0] * x**0. + cf[1] * x**1. + cf[2] * x**2. + cf[3] * x**3. + cf[4] * x**4.
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






def run_mcmc(pos, R, eR, diagnostics, Nsteps = 300, ndim = 1, nwalkers = 100):
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(R, eR, diagnostics))
    sampler.run_mcmc(pos, Nsteps)
    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
    OH_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0)))    
    return list(OH_mcmc)










if __name__ == '__main__':
    kern = Box2DKernel(3)
    nwalkers = 100
    outdir = '/user/rsimons/metal_maps'
    SN_limit = 3
    np.random.seed()
    field, di = argv[1], argv[2]
    fl = glob('/user/rsimons/grizli_extractions/%s/j*/Prep/*%s.full.fits'%(field, di))[0]

    if os.path.isfile(fl):


        master_hdulist = []
        prihdr = fits.Header()
        prihdr['COMMENT'] = "Storing the metallicity maps in this FITS file."
        prihdr['field']   = field
        prihdr['ID']      = di
        prihdu = fits.PrimaryHDU(header=prihdr)    
        master_hdulist.append(prihdu)

        colhdr = fits.Header()

        full = fits.open(fl)
        Rs  = []
        eRs = []
        diagnostics = []
        haslines = full[0].header['haslines']

        #do we have OII?
        if 'OII' in haslines:
            O2  = full['LINE', 'OII'].data
            eO2 = 1./np.sqrt(full['LINEWHT', 'OII'].data)

            O2 = convolve_fft(O2, kern)
            eO2 /= sqrt(3.)

            master_hdulist.append(fits.ImageHDU(data = O2, header = colhdr, name = 'OII_BC'))
            master_hdulist.append(fits.ImageHDU(data = eO2, header = colhdr, name = 'eOII'))

        #do we have OIII?
        if 'OIII' in haslines:
            O3  = full['LINE', 'OIII'].data
            eO3 = 1./np.sqrt(full['LINEWHT', 'OIII'].data)
           
            O3 = convolve_fft(O3, kern)
            eO3 /= sqrt(3.)
            master_hdulist.append(fits.ImageHDU(data = O3, header = colhdr, name = 'OIII_BC'))
            master_hdulist.append(fits.ImageHDU(data = eO3, header = colhdr, name = 'eOIII'))

        #do we have Hb?
        if 'Hb' in haslines:
            Hb  = full['LINE', 'Hb'].data
            eHb = 1./np.sqrt(full['LINEWHT', 'Hb'].data)

            Hb = convolve_fft(Hb, kern)
            eHb /= sqrt(3.)

            master_hdulist.append(fits.ImageHDU(data = Hb, header = colhdr, name = 'Hb_BC'))
            master_hdulist.append(fits.ImageHDU(data = eHb, header = colhdr, name = 'eHb'))

        #do we have O32?
        if ('OII' in haslines) & ('OIII' in haslines):
            R_O32 = O3/O2
            eR_O32 = R_O32 * np.sqrt((eO3/O3)**2. + (eO2/O2)**2.)
            diagnostics.append(['O32'])
            Rs.append(R_O32)
            eRs.append(eR_O32)

            master_hdulist.append(fits.ImageHDU(data = R_O32, header = colhdr, name = 'O32'))
            master_hdulist.append(fits.ImageHDU(data = eR_O32, header = colhdr, name = 'eO32'))



        #do we have R2?
        if ('OII' in haslines) & ('Hb' in haslines):
            R_R2 = O2/Hb
            eR_R2 = R_R2 * np.sqrt((eO2/O2)**2. + (eHb/Hb)**2.)
            diagnostics.append(['R2'])
            Rs.append(R_R2)
            eRs.append(eR_R2)
            master_hdulist.append(fits.ImageHDU(data = R_R2, header = colhdr, name = 'R2'))
            master_hdulist.append(fits.ImageHDU(data = eR_R2, header = colhdr, name = 'eR2'))


        #do we have R3?
        if ('OIII' in haslines) & ('Hb' in haslines):
            R_R3 = O3/Hb
            eR_R3 = R_R3 * np.sqrt((eO3/O3)**2. + (eHb/Hb)**2.)
            diagnostics.append(['R3'])
            Rs.append(R_R3)
            eRs.append(eR_R3)
            master_hdulist.append(fits.ImageHDU(data = R_R3, header = colhdr, name = 'R3'))
            master_hdulist.append(fits.ImageHDU(data = eR_R3, header = colhdr, name = 'eR3'))

        #do we have R23?
        if ('OII' in haslines) & ('OIII' in haslines) & ('Hb' in haslines):
            R_R23 = (O2 + O3)/Hb
            eR_R23 = R_R23 * np.sqrt((eO3**2. + eO2**2.)/(O3 + O2)**2. + (eHb/Hb)**2.)
            diagnostics.append(['R23'])
            Rs.append(R_R23)
            eRs.append(eR_R23)
            master_hdulist.append(fits.ImageHDU(data = R_R23, header = colhdr, name = 'R23'))
            master_hdulist.append(fits.ImageHDU(data = eR_R23, header = colhdr, name = 'eR23'))

        all_diags = []
        all_Rs = []
        all_eRs = []

        for dd, d in enumerate(array(diagnostics)): all_diags.append(d[0])
        diagnostics.append(all_diags)



        diagnostics = array(diagnostics)
        Rs = array(Rs)
        eRs = array(eRs)



        all_Rs = np.empty((shape(Rs)[1], shape(Rs)[2]), dtype = 'object')
        all_eRs = np.empty((shape(Rs)[1], shape(Rs)[2]), dtype = 'object')

        for d, diagnostic in enumerate(diagnostics[0:-1]):
            print (diagnostic)
            Z = nan * zeros((shape(Rs)[1], shape(Rs)[2], 3))

            for i in arange(shape(Rs)[1]):
                for j in arange(shape(Rs)[2]):
                    print (i,j)
                    Rs_ij = array([Rs[d][i,j]])
                    eRs_ij = array([eRs[d][i,j]])

                    if all_Rs[i,j] == None: all_Rs[i,j] = [Rs_ij]
                    else: all_Rs[i,j].append(Rs_ij)
                    if all_eRs[i,j] == None: all_eRs[i,j] = [eRs_ij]
                    else: all_eRs[i,j].append(eRs_ij)

                    nll = lambda *args: -lnlike(*args)
                    if Rs_ij[0]/eRs_ij[0] > SN_limit:
                        result = op.minimize(nll, [8.5], args=(Rs_ij, eRs_ij, diagnostic))
                        OH_ml = result["x"]
                        pos = [result["x"] + 1e-4*np.random.randn(1) for nn in range(nwalkers)]
                        OH_result = run_mcmc(pos = pos, R = Rs_ij, eR = eRs_ij, 
                                             diagnostics = diagnostic, nwalkers = nwalkers)
                        Z[i,j,0]  = OH_result[0][0]
                        Z[i,j,1]  = OH_result[0][1]
                        Z[i,j,2]  = OH_result[0][2]

            master_hdulist.append(fits.ImageHDU(data = Z, header = colhdr, name = 'Z_%s'%diagnostic[0]))


        Z = nan * zeros((shape(Rs)[1], shape(Rs)[2], 3))
        for i in arange(shape(Rs)[1]):
            for j in arange(shape(Rs)[2]):
                Ndet = len(where(array(all_Rs[i,j])/array(all_eRs[i,j]) > SN_limit)[0])
                Ntot = len(all_Rs[i,j])
                if Ndet > Ntot - 1:
                    nll = lambda *args: -lnlike(*args)
                    result = op.minimize(nll, [8.5], args=(all_Rs[i,j], all_eRs[i,j], diagnostics[-1]))
                    OH_ml = result["x"]
                    pos = [result["x"] + 1e-4*np.random.randn(1) for nn in range(nwalkers)]
                    OH_result = run_mcmc(pos = pos, R = all_Rs[i,j], eR = all_eRs[i,j], 
                                         diagnostics = diagnostics[-1], nwalkers = nwalkers)
                    Z[i,j,0]  = OH_result[0][0]
                    Z[i,j,1]  = OH_result[0][1]
                    Z[i,j,2]  = OH_result[0][2]

        master_hdulist.append(fits.ImageHDU(data = Z, header = colhdr, name = 'Z_all'%diagnostic[0]))


        fits_name = out_dir + '/%s_%s_metals.fits'%(field, di)
        print ('\tSaving to ' + fits_name)
        thdulist = fits.HDUList(master_hdulist)
        thdulist.writeto(fits_name, clobber = True)

















