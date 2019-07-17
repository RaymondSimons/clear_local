#!/home/rsimons/miniconda2/envs/grizli/bin/python
import pidly
from calzetti import k as calk
import numpy as np
import os
import time
from multiprocessing import Pool
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
import metal_calibs as calib
import scipy.optimize as op
import time
from math import *
from sys import argv
from hri import hri

def izi(fluxes, errors, lines, idl=None, dosave=False, savfile='res.sav', 
            grid=os.path.join(os.environ['IZI_DIR'],'grids','l09_high_csf_n1e2_6.0Myr.fits')) :

            #idl = pidly.IDL()
            idl('fluxes = {0}'.format(np.array2string(fluxes, separator=',',max_line_width=1000)))
            idl('errors = {0}'.format(np.array2string(errors, separator=',',max_line_width=1000)))
            idl('lines = {0}'.format(np.array2string(lines, separator=',',max_line_width=1000)))

            idl('res=izi(fluxes, errors, lines, NZ=100, gridfile="{0}")'.format(grid))
            if dosave :
                print ('saving')
                idl('save, file="{0}", res'.format(savfile))
            res = idl.ev('res', use_cache=True)
            return(res)




if __name__ == '__main__':
    np.random.seed(1)
    boxcar_size = 3
    kern = Box2DKernel(boxcar_size)
    field, di = argv[1], argv[2]

    out_dir = '/user/rsimons/metal_maps'
    full_dir = '/user/rsimons/grizli_extractions'

    print ('%s/%s/j*/Prep/*%s.full.fits'%(full_dir, field, di))
    fl = glob('%s/%s/j*/Prep/*%s.full.fits'%(full_dir, field, di))[0]


    wdth = 15
    xmd = 40

    xmn = xmd - wdth
    xmx = xmd + wdth

    ymn = xmd - wdth
    ymx = xmd + wdth

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
        full = fits.open(fl)
        haslines = full[0].header['haslines']


        lines = [('OII', 'oii3726;oii3729'),
                 ('OIII', 'oiii4959;oiii5007'),
                 ('Hb', 'hbeta'),
                 ('Ha', 'nii6548;halpha;nii6584'),
                 ('SII', 'sii6717;sii6731')
                ]


        lines_use = []
        for l, (line, izi_line) in enumerate(lines):
            if line in haslines:
                lines_use.append((line, izi_line))

                lmap  = full['LINE', line].data[xmn:xmx, ymn:ymx]
                elmap = 1./np.sqrt(full['LINEWHT', line].data[xmn:xmx, ymn:ymx])

                lmap_smoothed = convolve_fft(lmap, kern)
                elmap_smoothed = elmap/np.sqrt(boxcar_size**2.)

                master_hdulist.append(fits.ImageHDU(data = lmap, header = colhdr, name =  '%s'%line))
                master_hdulist.append(fits.ImageHDU(data = elmap, header = colhdr, name = 'e%s'%line))

                master_hdulist.append(fits.ImageHDU(data = lmap_smoothed, header = colhdr, name =  '%s_s'%line))
                master_hdulist.append(fits.ImageHDU(data = elmap_smoothed, header = colhdr, name = 'e%s_s'%line))


        thdulist_temp = fits.HDUList(master_hdulist)
        Z = nan * zeros((shape(lmap)[0], shape(lmap)[1], 4))


        idl_path = '/grp/software/Linux/itt/idl/idl84/idl/bin/idl'
        idl = pidly.IDL(idl_path)



        for i in arange(shape(Z)[0]):
            print (i)
            for j in arange(shape(Z)[0]):
                savfile = out_dir + '/%s_%s_%i_%i.sav'%(field, di, i, j)
                fluxes_for_izi = []
                errors_for_izi = []
                lines_for_izi = []
                for l, (line, izi_line) in enumerate(lines_use):
                    fluxes_for_izi.append(thdulist_temp['%s_s'%line].data[i,j])
                    errors_for_izi.append(thdulist_temp['e%s_s'%line].data[i,j])
                    lines_for_izi.append(izi_line)
                fluxes_for_izi = np.array(fluxes_for_izi)
                errors_for_izi = np.array(errors_for_izi)
                lines_for_izi  = np.array(lines_for_izi)

                n_detected = len(np.where(fluxes_for_izi/errors_for_izi > 1.)[0])
                if n_detected > 1:
                    res = izi(fluxes_for_izi, errors_for_izi, lines_for_izi, idl=idl, dosave=False, savfile=savfile,
                                  grid=os.environ['IZI_DIR']+'/grids/d13_kappa20.fits')
                    (tZmod, tZlo, tZhi, tnpeaks) = hri( res['zarr'][0], res['zpdfmar'][0])

                    Z[i,j,0] = tZmod
                    Z[i,j,1] = tZlo
                    Z[i,j,2] = tZhi
                    Z[i,j,3] = tnpeaks


        master_hdulist.append(fits.ImageHDU(data = Z, header = Zcolhdr, name = 'Z'))
        fits_name = out_dir + '/%s_%s_metals.fits'%(field, di)
        print ('\tSaving to ' + fits_name)
        thdulist = fits.HDUList(master_hdulist)
        thdulist.writeto(fits_name, overwrite = True)

















'''


def izi(fluxes, errors, lines, idl=None, dosave=False, savfile='res.sav', 
            grid=os.path.join(os.environ['IZI_DIR'],'grids','l09_high_csf_n1e2_6.0Myr.fits')) :

            #idl = pidly.IDL()
            idl('fluxes = {0}'.format(np.array2string(fluxes, separator=',',max_line_width=1000)))
            idl('errors = {0}'.format(np.array2string(errors, separator=',',max_line_width=1000)))
            idl('lines = {0}'.format(np.array2string(lines, separator=',',max_line_width=1000)))
            #idl('forprint, fluxes, errors, lines')
            #print(grid, os.path.isfile(grid))
            #print('gridfile={0})'.format(grid))
            idl('res=izi(fluxes, errors, lines, NZ=100, gridfile="{0}")'.format(grid))
            if dosave :
                idl('save, file="{0}", res'.format(savfile))
            res = idl.ev('res', use_cache=True)
            return(res)


idl_path = '/grp/software/Linux/itt/idl/idl84/idl/bin/idl'
idl = pidly.IDL(idl_path)
Vlam = 5470. # from Johnson Cousins_V 
calkV = calk(Vlam)
tlam = np.array([3727., 5007.,  4863.])

lines = np.array(['oii3726;oii3729','oiii4959;oiii5007',  'hbeta'])
fluxes = np.array([0.4e-17, 0.8e-17, 1.e-17])
errors = np.array([0.1e-17, 0.2e-17, 1.e-18])

savfile = 'test.sav'

res = izi(fluxes, errors, lines, idl=idl, dosave=True, savfile=savfile,
              grid=os.environ['IZI_DIR']+'/grids/d13_kappa20.fits')

'''
# take AV from nZ:
'''
Av = nZ['Av'].iloc[i]

for l in range(len(tlam)) :
    Alam = 1*Av * calk(tlam[l]) / calkV
    fluxes[l] = fluxes[l] * np.power(10, 0.4*Alam)
    errors[l] = errors[l] * np.power(10, 0.4*Alam)
    #print(Av, Alam, np.power(10,0.4*Alam))


errors = errors / fluxes[2]
fluxes = fluxes / fluxes[2]
'''



