#!/home/rsimons/miniconda2/envs/grizli/bin/python
import time
from multiprocessing import Pool
import os
import numpy as np
from numpy import *
import astropy
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import Planck15 as cosmo
from astropy.convolution import Gaussian2DKernel, convolve_fft, Box2DKernel
from astropy.stats import sigma_clip
import importlib
import photutils
import glob
from glob import glob
from scipy.interpolate import interp1d
import joblib
from joblib import Parallel, delayed
import emcee
import metal_calibs as calib
import scipy.optimize as op
import time
from math import *
from sys import argv


np.random.seed(1)

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
    OHmu = 8.7     # mean of the Gaussian prior
    OHsigma = 0.5 # standard deviation of the Gaussian prior
    lp -= 0.5*((OH - OHmu)/OHsigma)**2   
    return lp

def lnprob(OH, R, Rerr, diagnostics, use_prior = 'top'):
    if use_prior == 'top': lp = lnprior(OH)
    if use_prior == 'gaussian': lp = lnGaussianprior(OH)

    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(OH, R, Rerr, diagnostics)


def run_mcmc(pos, R, eR, diagnostics, Nsteps = 300, Nburn = 50, Ndim = 1, Nwalkers = 100, use_prior = 'gaussian'):
    sampler = emcee.EnsembleSampler(Nwalkers, Ndim, lnprob, args=(R, eR, diagnostics, use_prior))
    sampler.run_mcmc(pos, Nsteps)
    samples = sampler.chain[:, Nburn:, :].reshape((-1, Ndim))
    OH_mcmc = map(lambda v: (v[2], v[3]-v[2], v[2]-v[1], v[4]-v[2], v[2]-v[0]), zip(*np.percentile(samples, [2.5, 16, 50, 84, 97.5], axis=0)))    
    return list(OH_mcmc)[0], samples




if __name__ == '__main__':
    boxcar_size = 3
    kern = Box2DKernel(boxcar_size)
    Nwalkers = 100
    Nsteps = 300
    Nburn = 50
    Ndim = 1
    Nchain_saved = 2000

    SN_limit = 1./log(10) #we are going to look at pixels with ratios of S/N > 1 or better, in the absolute sense
    np.random.seed()
    field, di = argv[1], argv[2]

    on_science = True
    if on_science:
        out_dir = '/user/rsimons/metal_maps'
        full_dir = '/user/rsimons/grizli_extractions'
    else:
        out_dir  = '/Volumes/pegasus/clear/metal_maps/local_testing'
        full_dir = '/Volumes/pegasus/clear/grizli_extractions'
    

    print ('%s/%s/j*/Prep/*%s.full.fits'%(full_dir, field, di))
    fl = glob('%s/%s/j*/Prep/*%s.full.fits'%(full_dir, field, di))[0]



    if os.path.isfile(fl):


        master_hdulist = []
        prihdr = fits.Header()
        prihdr['COMMENT'] = "Storing the metallicity maps in this FITS file."
        prihdr['field']   = field
        prihdr['ID']      = di
        prihdu = fits.PrimaryHDU(header=prihdr)    
        master_hdulist.append(prihdu)

        colhdr = fits.Header()
        Zcolhdr = fits.Header()

        colhdr['boxcar_size']=boxcar_size

        Zcolhdr['SN_limit']=SN_limit
        Zcolhdr['Nwalkers']=Nwalkers
        Zcolhdr['Nburn']=Nburn
        Zcolhdr['Ndim']=Ndim
        Zcolhdr['Nsteps']=Nsteps




        full = fits.open(fl)
        Rs  = []
        eRs = []
        diagnostics = []
        haslines = full[0].header['haslines']
        if True:
            #do we have OII?
            if 'OII ' in haslines:
                O2  = full['LINE', 'OII'].data
                eO2 = 1./np.sqrt(full['LINEWHT', 'OII'].data)

                O2 = convolve_fft(O2, kern)
                eO2 /= sqrt(boxcar_size**2.)

                master_hdulist.append(fits.ImageHDU(data = O2, header = colhdr, name =  'OII'))
                master_hdulist.append(fits.ImageHDU(data = eO2, header = colhdr, name = 'eOII'))

            #do we have OIII?
            if 'OIII ' in haslines:
                O3  = full['LINE', 'OIII'].data
                eO3 = 1./np.sqrt(full['LINEWHT', 'OIII'].data)
               
                O3 = convolve_fft(O3, kern)
                eO3 /= sqrt(boxcar_size**2.)
                master_hdulist.append(fits.ImageHDU(data = O3, header = colhdr, name = 'OIII'))
                master_hdulist.append(fits.ImageHDU(data = eO3, header = colhdr, name = 'eOIII'))

            #do we have Hb?
            if 'Hb' in haslines:
                Hb  = full['LINE', 'Hb'].data
                eHb = 1./np.sqrt(full['LINEWHT', 'Hb'].data)

                Hb = convolve_fft(Hb, kern)
                eHb /= sqrt(boxcar_size**2.)

                master_hdulist.append(fits.ImageHDU(data = Hb, header = colhdr, name = 'Hb'))
                master_hdulist.append(fits.ImageHDU(data = eHb, header = colhdr, name = 'eHb'))

            #do we have O32?
            if ('OII ' in haslines) & ('OIII ' in haslines):
                R_O32 = O3/O2
                eR_O32 = R_O32 * np.sqrt((eO3/O3)**2. + (eO2/O2)**2.)
                diagnostics.append(['O32'])
                Rs.append(R_O32)
                eRs.append(eR_O32)

                master_hdulist.append(fits.ImageHDU(data = R_O32, header = colhdr, name = 'O32'))
                master_hdulist.append(fits.ImageHDU(data = eR_O32, header = colhdr, name = 'eO32'))



            #do we have R2?
            if ('OII ' in haslines) & ('Hb' in haslines):
                R_R2 = O2/Hb
                eR_R2 = R_R2 * np.sqrt((eO2/O2)**2. + (eHb/Hb)**2.)
                diagnostics.append(['R2'])
                Rs.append(R_R2)
                eRs.append(eR_R2)
                master_hdulist.append(fits.ImageHDU(data = R_R2, header = colhdr, name = 'R2'))
                master_hdulist.append(fits.ImageHDU(data = eR_R2, header = colhdr, name = 'eR2'))


            #do we have R3?
            if ('OIII ' in haslines) & ('Hb' in haslines):
                R_R3 = O3/Hb
                eR_R3 = R_R3 * np.sqrt((eO3/O3)**2. + (eHb/Hb)**2.)
                diagnostics.append(['R3'])
                Rs.append(R_R3)
                eRs.append(eR_R3)
                master_hdulist.append(fits.ImageHDU(data = R_R3, header = colhdr, name = 'R3'))
                master_hdulist.append(fits.ImageHDU(data = eR_R3, header = colhdr, name = 'eR3'))

            #do we have R23?
            if ('OII ' in haslines) & ('OIII ' in haslines) & ('Hb' in haslines):
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


        minx = int(shape(Rs)[1]/2 - 20)
        maxx = int(shape(Rs)[1]/2 + 20)

        full_hdulist = []
        for d, diagnostic in enumerate(diagnostics[0:-1]):
            print ('calculating metallicity using ', diagnostic)
            Z = nan * zeros((shape(Rs)[1], shape(Rs)[2], 5))
            Z_full = nan * zeros((shape(Rs)[1], shape(Rs)[2],  Nchain_saved))
            for i in arange(minx, maxx):
                for j in arange(minx, maxx):
                    if Rs[d][i,j] > 0:
                        Rs_ij = array([log10(Rs[d][i,j])])
                        eRs_ij = array([eRs[d][i,j]/Rs[d][i,j]/log(10)])

                        if all_Rs[i,j] == None: all_Rs[i,j] = [Rs_ij]
                        else: all_Rs[i,j].append(Rs_ij)
                        if all_eRs[i,j] == None: all_eRs[i,j] = [eRs_ij]
                        else: all_eRs[i,j].append(eRs_ij)

                        nll = lambda *args: -lnlike(*args)
                        #if Rs_ij[0]/eRs_ij[0] > SN_limit:
                        if eRs_ij[0] < SN_limit:
                            result = op.minimize(nll, [8.5], args=(Rs_ij, eRs_ij, diagnostic))
                            OH_ml = result["x"]
                            pos = [result["x"] + 1e-4*np.random.randn(1) for nn in range(Nwalkers)]

                            OH_result, samples = run_mcmc(pos = pos, R = Rs_ij, eR = eRs_ij, 
                                                          diagnostics = diagnostic, Nsteps = Nsteps,
                                                          Nburn = Nburn, Ndim = Ndim, Nwalkers = Nwalkers)

                            Z[i,j]  = OH_result
                            Z_full[i,j]  = samples[np.random.randint(0, len(samples[:,0]), Nchain_saved),0]


            if diagnostic[0] == 'R23': use = 'M08'
            if diagnostic[0] == 'R2':  use = 'M08'
            if diagnostic[0] == 'R3':  use = 'M08'
            if diagnostic[0] == 'O32': use = 'M08'
            Zcolhdr['calibration'] = use
            master_hdulist.append(fits.ImageHDU(data = Z, header = Zcolhdr, name = 'Z_%s'%diagnostic[0]))
            full_hdulist.append(fits.ImageHDU(data = Z_full, header = Zcolhdr, name = 'Z_%s_full'%diagnostic[0]))


        Z = nan * zeros((shape(Rs)[1], shape(Rs)[2], 5))
        Z_full = nan * zeros((shape(Rs)[1], shape(Rs)[2], Nchain_saved))
        #Z_full = np.empty((shape(Rs)[1], shape(Rs)[2]), dtype = 'object')
        print ('calculating metallicity using all available diagnostics: ', diagnostics[-1])

        for i in arange(minx, maxx):
            for j in arange(minx, maxx):
                if all_Rs[i,j] != None:
                    #Ndet = len(where(array(all_Rs[i,j])/array(all_eRs[i,j]) > SN_limit)[0])
                    Ndet = len(where(array(all_eRs[i,j]) < SN_limit)[0])
                    if Ndet > 0:
                        nll = lambda *args: -lnlike(*args)
                        result = op.minimize(nll, [8.5], args=(array(all_Rs[i,j]), array(all_eRs[i,j]), diagnostics[-1]))
                        OH_ml = result["x"]
                        pos = [result["x"] + 1e-4*np.random.randn(1) for nn in range(Nwalkers)]
                        OH_result, samples = run_mcmc(pos = pos, R = array(all_Rs[i,j]), eR = array(all_eRs[i,j]), 
                                             diagnostics = diagnostics[-1], Nsteps = Nsteps, 
                                             Nburn = Nburn, Ndim = Ndim, Nwalkers = Nwalkers)
                        Z[i,j] = OH_result
                        Z_full[i,j]  = samples[np.random.randint(0, len(samples[:,0]), Nchain_saved),0]

        master_hdulist.append(fits.ImageHDU(data = Z, header = Zcolhdr, name = 'Z_all'))
        full_hdulist.append(fits.ImageHDU(data = Z_full, header = Zcolhdr, name = 'Z_all_full'))

        master_hdulist.extend(full_hdulist)
        fits_name = out_dir + '/%s_%s_metals.fits'%(field, di)
        print ('\tSaving to ' + fits_name)
        thdulist = fits.HDUList(master_hdulist)
        thdulist.writeto(fits_name, overwrite = True)

















