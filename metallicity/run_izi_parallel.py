import pidly
from clear_local.caseys_izi_tools.calzetti import k as calk
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
from astropy.convolution import Gaussian1DKernel, Gaussian2DKernel, convolve_fft, Box2DKernel
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
import os
from joblib import Parallel, delayed



def izi(fluxes, errors, lines, logzprior = None, idl=None, dosave=False, savfile='res.sav', 
            grid=os.path.join(os.environ['IZI_DIR'],'grids','d13_kappa20.fits')) :

            #idl = pidly.IDL()

            idl('fluxes = {0}'.format(np.array2string(fluxes, separator=',',max_line_width=1000)))
            idl('errors = {0}'.format(np.array2string(errors, separator=',',max_line_width=1000)))
            idl('lines = {0}'.format(np.array2string(lines, separator=',',max_line_width=1000)))
            idl('logzprior = {0}'.format(np.array2string(logzprior, separator=',',max_line_width=1000)).replace('\n', ''))
            idl('res=izi(fluxes, errors, lines, LOGZPRIOR=logzprior, NZ=100, gridfile="{0}")'.format(grid))
            if dosave :
                idl('save, file="{0}", res'.format(savfile))
            res = idl.ev('res', use_cache=True)
            return(res)


def run_izi(Z, Z_pdf, idl, thdulist_temp, lines_use, Av = None, do_extinction = True, smooth = True):
    start_prior = 8.5
    end_prior = 9.5
    dZ = 0.05
    z_arr = np.arange(start_prior, end_prior, dZ)
    prior = zeros(len(z_arr)) + 1.
    #prior[int((8.5 - 7)/dZ)] = 1.
    #gauss_kernel = Gaussian1DKernel(0.5/dZ)
    #prior = convolve_fft(prior, gauss_kernel)
    logzprior = vstack((z_arr, prior))


    for i in arange(shape(Z)[0]):
        for j in arange(shape(Z)[0]):
            fluxes_for_izi = []
            errors_for_izi = []
            lines_for_izi = []
            for l, (line, izi_line) in enumerate(lines_use):
                fl_str = line
                er_str = 'e'+line
                if smooth: 
                    fl_str+='_s'
                    er_str+='_s'
                if do_extinction: 
                    fl_str+='_ec'
                    er_str+='_ec'

                fluxes_for_izi.append(thdulist_temp[fl_str].data[i,j])
                errors_for_izi.append(thdulist_temp[er_str].data[i,j])
                lines_for_izi.append(izi_line)


            fluxes_for_izi = np.array(fluxes_for_izi)
            errors_for_izi = np.array(errors_for_izi)
            lines_for_izi  = np.array(lines_for_izi)

            gd = where((np.isfinite(fluxes_for_izi)) & (np.isfinite(errors_for_izi)) & (fluxes_for_izi > 0.))[0]
            fluxes_for_izi = fluxes_for_izi[gd]
            errors_for_izi = errors_for_izi[gd]
            lines_for_izi = lines_for_izi[gd]




            n_detected = len(np.where(fluxes_for_izi/errors_for_izi > 1.)[0])
            if n_detected > 1:
                res = izi(fluxes_for_izi, errors_for_izi, lines_for_izi, logzprior=logzprior, idl=idl, dosave=False, savfile=None,
                              grid=os.environ['IZI_DIR']+'/grids/d13_kappa20.fits')
                (tZmod, tZlo, tZhi, tnpeaks) = hri( res['zarr'][0], res['zpdfmar'][0])

                Z_pdf[i,j,:,0] = res['zarr'][0]
                Z_pdf[i,j,:,1] = res['zpdfmar'][0]
                Z[i,j,0] = tZmod
                Z[i,j,1] = tZlo
                Z[i,j,2] = tZhi
                Z[i,j,3] = tnpeaks

    return Z, Z_pdf

def run_all(field, di, out_dir = '/user/rsimons/metal_maps', full_dir = '/user/rsimons/grizli_extractions'):
    
    np.random.seed(1)
    boxcar_size = 3
    kern = Box2DKernel(boxcar_size)

    if 'S' in field: fld = 'goodss'
    if 'N' in field: fld = 'goodsn'


    eazy_fits = fits.open('/user/rsimons/grizli_extractions/Catalogs/%s_3dhst.v4.4.cats/Eazy/%s_3dhst.v4.4.zout.fits'%(fld, fld))
    gd = where(eazy_fits[1].data['id'].astype('int') == int(di))[0][0]
    Av = eazy_fits[1].data['Av'][gd]


    out_dir = '/user/rsimons/metal_maps_v3'
    full_dir = '/user/rsimons/grizli_extractions_v3'

    #print ('%s/%s/j*/Prep/*%s.full.fits'%(full_dir, field, di))
    #fl = glob('%s/%s/j*/Prep/*%s.full.fits'%(full_dir, field, di))[0]
    fl = glob('%s/%s/*%s.full.fits'%(full_dir, field, di))[0]

    wdth = 20
    xmd = 80


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
        prihdr['Av']      = Av
        prihdu = fits.PrimaryHDU(header=prihdr)    
        master_hdulist.append(prihdu)


        colhdr = fits.Header()
        Zcolhdr = fits.Header()

        colhdr['bc_kern']=boxcar_size
        full = fits.open(fl)
        master_hdulist.append(full['DSCI'])
        master_hdulist.append(full['DWHT'])
        master_hdulist.append(full['SEG'])


        haslines = full[0].header['haslines']


        lines = [('OII', 'oii3726;oii3729', 3727.),
                 ('OIII', 'oiii4959;oiii5007', 5007.),
                 ('Hb', 'hbeta', 4863.),
                 ('Ha', 'nii6548;halpha;nii6584', 6563.),
                 ('SII', 'sii6717;sii6731', 6725.)
                ]
        
        #lines = [('OIII', 'oiii4959;oiii5007', 5007.),
        #         ('Hb', 'hbeta', 4863.)
        #        ]



        '''
        lines = [('OII', 'oii3726;oii3729', 3727.),
                 ('OIII', 'oiii4959;oiii5007', 5007.),
                 ('Hb', 'hbeta', 4863.)
                ]
        '''

        Vlam = 5470. # from Johnson Cousins_V 
        calkV = calk(Vlam)


        lines_use = []
        haslines = np.array(haslines.split(' '))
        for l, (line, izi_line, line_wav) in enumerate(lines):
            gd = where(line == haslines)[0]
            if len(gd) > 0:
                lines_use.append((line, izi_line))

                lmap  = full['LINE', line].data
                elmap = 1./np.sqrt(full['LINEWHT', line].data)

                lmap_smoothed = convolve_fft(lmap, kern)
                elmap_smoothed = elmap/np.sqrt(boxcar_size**2.)

                Alam = 1*Av * calk(line_wav) / calkV
                ec = np.power(10, 0.4*Alam)
                '''
                master_hdulist.append(fits.ImageHDU(data = lmap[xmn:xmx, ymn:ymx], header = colhdr, name =  '%s'%line))
                master_hdulist.append(fits.ImageHDU(data = elmap[xmn:xmx, ymn:ymx], header = colhdr, name = 'e%s'%line))

                master_hdulist.append(fits.ImageHDU(data = lmap_smoothed[xmn:xmx, ymn:ymx], header = colhdr, name =  '%s_s'%line))
                master_hdulist.append(fits.ImageHDU(data = elmap_smoothed[xmn:xmx, ymn:ymx], header = colhdr, name = 'e%s_s'%line))


                master_hdulist.append(fits.ImageHDU(data = lmap[xmn:xmx, ymn:ymx] * ec, header = colhdr, name =  '%s_ec'%line))
                master_hdulist.append(fits.ImageHDU(data = elmap[xmn:xmx, ymn:ymx] * ec, header = colhdr, name = 'e%s_ec'%line))
                '''
                master_hdulist.append(fits.ImageHDU(data = lmap_smoothed[xmn:xmx, ymn:ymx] * ec, header = colhdr, name =  '%s_s_ec'%line))
                master_hdulist.append(fits.ImageHDU(data = elmap_smoothed[xmn:xmx, ymn:ymx] * ec, header = colhdr, name = 'e%s_s_ec'%line))



        thdulist_temp = fits.HDUList(master_hdulist)
        Z_empty     = nan * zeros((wdth*2, wdth*2, 4))
        Z_pdf_empty = nan * zeros((wdth*2, wdth*2, 100, 2))

        idl_path = '/Applications/harris/idl87/bin/idl'#'/Applications/harris/idl87'#'/grp/software/Linux/itt/idl/idl84/idl/bin/idl'
        idl = pidly.IDL(idl_path)


        #Z, Z_pdf     = run_izi(Z = Z_empty, Z_pdf = Z_pdf_empty, idl = idl, thdulist_temp = thdulist_temp,  lines_use = lines_use, do_extinction = False, smooth = False)
        #Z_s, Z_pdf_s = run_izi(Z = Z_empty, Z_pdf = Z_pdf_empty, idl = idl, thdulist_temp = thdulist_temp, lines_use = lines_use,  do_extinction = False,smooth = True)
        #Z_ec, Z_pdf_ec     = run_izi(Z = Z_empty, Z_pdf = Z_pdf_empty, idl = idl, thdulist_temp = thdulist_temp,  lines_use = lines_use, do_extinction = True, smooth = False)
        Z_s_ec, Z_pdf_s_ec = run_izi(Z = Z_empty, Z_pdf = Z_pdf_empty, idl = idl, thdulist_temp = thdulist_temp, lines_use = lines_use,  do_extinction = True,smooth = True)


        #master_hdulist.append(fits.ImageHDU(data = Z, header = Zcolhdr, name = 'Z'))
        #master_hdulist.append(fits.ImageHDU(data = Z_s, header = Zcolhdr, name = 'Z_s'))
        #master_hdulist.append(fits.ImageHDU(data = Z_ec, header = Zcolhdr, name = 'Z_ec'))
        master_hdulist.append(fits.ImageHDU(data = Z_s_ec, header = Zcolhdr, name = 'Z_s_ec'))


        #master_hdulist.append(fits.ImageHDU(data = Z_pdf, header = Zcolhdr, name = 'Z_pdf'))
        #master_hdulist.append(fits.ImageHDU(data = Z_pdf_s, header = Zcolhdr, name = 'Z_pdf_s'))
        #master_hdulist.append(fits.ImageHDU(data = Z_pdf_ec, header = Zcolhdr, name = 'Z_pdf_ec'))
        master_hdulist.append(fits.ImageHDU(data = Z_pdf_s_ec, header = Zcolhdr, name = 'Z_pdf_s_ec'))



        fits_name = out_dir + '/%s_%s_metals_highZbranch.fits'%(field, di)
        print ('\tSaving to ' + fits_name)
        thdulist = fits.HDUList(master_hdulist)
        thdulist.writeto(fits_name, overwrite = True)
if __name__ == '__main__':
    izi_cat = ascii.read('/user/rsimons/good_izi.cat', header_start = 0)
    #Parallel(n_jobs = -1)(delayed(run_all)(fld, di) for f, (fld, di) in enumerate(zip(izi_cat['field'], izi_cat['id'])))        


    for f, (fld, di) in enumerate(zip(izi_cat['field'], izi_cat['id'])):
        run_all(fld, di)
    

