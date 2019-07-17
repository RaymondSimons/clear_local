#!/home/rsimons/miniconda2/bin/python
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
                idl('save, file="{0}", res'.format(savfile))
            res = idl.ev('res', use_cache=True)
            return(res)


def run_izi(Z, Z_pdf, idl, thdulist_temp, lines_use, Av = None, do_extinction = True, smooth = True):
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

            gd = where((np.isfinite(fluxes_for_izi)) & (np.isfinite(errors_for_izi)))[0]
            fluxes_for_izi = fluxes_for_izi[gd]
            errors_for_izi = errors_for_izi[gd]
            lines_for_izi = lines_for_izi[gd]


            n_detected = len(np.where(fluxes_for_izi/errors_for_izi > 1.)[0])
            if n_detected > 1:
                res = izi(fluxes_for_izi, errors_for_izi, lines_for_izi, idl=idl, dosave=False, savfile=None,
                              grid=os.environ['IZI_DIR']+'/grids/d13_kappa20.fits')
                (tZmod, tZlo, tZhi, tnpeaks) = hri( res['zarr'][0], res['zpdfmar'][0])

                Z_pdf[i,j,:,0] = res['zarr'][0]
                Z_pdf[i,j,:,1] = res['zpdfmar'][0]
                Z[i,j,0] = tZmod
                Z[i,j,1] = tZlo
                Z[i,j,2] = tZhi
                Z[i,j,3] = tnpeaks

    return Z, Z_pdf

if __name__ == '__main__':
    np.random.seed(1)
    boxcar_size = 3
    kern = Box2DKernel(boxcar_size)
    field, di = argv[1], argv[2]

    if 'S' in field: fld = 'goodss'
    if 'N' in field: fld = 'goodsn'


    eazy_fits = fits.open('/user/rsimons/grizli_extractions/Catalogs/%s_3dhst.v4.4.cats/Eazy/%s_3dhst.v4.4.zout.fits'%(fld, fld))
    gd = where(eazy_fits[1].data['id'].astype('int') == int(di))[0][0]
    Av = eazy_fits[1].data['Av'][gd]


    out_dir = '/user/rsimons/metal_maps'
    full_dir = '/user/rsimons/grizli_extractions'

    print ('%s/%s/j*/Prep/*%s.full.fits'%(full_dir, field, di))
    fl = glob('%s/%s/j*/Prep/*%s.full.fits'%(full_dir, field, di))[0]


    wdth = 20
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

                master_hdulist.append(fits.ImageHDU(data = lmap[xmn:xmx, ymn:ymx], header = colhdr, name =  '%s'%line))
                master_hdulist.append(fits.ImageHDU(data = elmap[xmn:xmx, ymn:ymx], header = colhdr, name = 'e%s'%line))

                master_hdulist.append(fits.ImageHDU(data = lmap_smoothed[xmn:xmx, ymn:ymx], header = colhdr, name =  '%s_s'%line))
                master_hdulist.append(fits.ImageHDU(data = elmap_smoothed[xmn:xmx, ymn:ymx], header = colhdr, name = 'e%s_s'%line))

                Alam = 1*Av * calk(line_wav) / calkV

                ec = np.power(10, 0.4*Alam)
                master_hdulist.append(fits.ImageHDU(data = lmap[xmn:xmx, ymn:ymx] * ec, header = colhdr, name =  '%s_ec'%line))
                master_hdulist.append(fits.ImageHDU(data = elmap[xmn:xmx, ymn:ymx] * ec, header = colhdr, name = 'e%s_ec'%line))

                master_hdulist.append(fits.ImageHDU(data = lmap_smoothed[xmn:xmx, ymn:ymx] * ec, header = colhdr, name =  '%s_s_ec'%line))
                master_hdulist.append(fits.ImageHDU(data = elmap_smoothed[xmn:xmx, ymn:ymx] * ec, header = colhdr, name = 'e%s_s_ec'%line))



        thdulist_temp = fits.HDUList(master_hdulist)
        Z = nan * zeros((wdth*2, wdth*2, 4))
        Z_s = nan * zeros((wdth*2, wdth*2, 4))

        Z_ec = nan * zeros((wdth*2, wdth*2, 4))
        Z_s_ec = nan * zeros((wdth*2, wdth*2, 4))



        Z_pdf = nan * zeros((wdth*2, wdth*2, 100, 2))
        Z_pdf_s = nan * zeros((wdth*2, wdth*2, 100, 2))

        Z_pdf_ec = nan * zeros((wdth*2, wdth*2, 100, 2))
        Z_pdf_s_ec = nan * zeros((wdth*2, wdth*2, 100, 2))



        idl_path = '/grp/software/Linux/itt/idl/idl84/idl/bin/idl'
        idl = pidly.IDL(idl_path)


        Z, Z_pdf     = run_izi(Z = Z, Z_pdf = Z_pdf, idl = idl, thdulist_temp = thdulist_temp,  lines_use = lines_use, do_extinction = False, smooth = False)
        Z_s, Z_pdf_s = run_izi(Z = Z_s, Z_pdf = Z_pdf_s, idl = idl, thdulist_temp = thdulist_temp, lines_use = lines_use,  do_extinction = False,smooth = True)

        Z_ec, Z_pdf_ec     = run_izi(Z = Z_ec, Z_pdf = Z_pdf_ec, idl = idl, thdulist_temp = thdulist_temp,  lines_use = lines_use, do_extinction = True, smooth = False)
        Z_s_ec, Z_pdf_s_ec = run_izi(Z = Z_s_ec, Z_pdf = Z_pdf_s_ec, idl = idl, thdulist_temp = thdulist_temp, lines_use = lines_use,  do_extinction = True,smooth = True)



        master_hdulist.append(fits.ImageHDU(data = Z, header = Zcolhdr, name = 'Z'))
        master_hdulist.append(fits.ImageHDU(data = Z_s, header = Zcolhdr, name = 'Z_s'))
        master_hdulist.append(fits.ImageHDU(data = Z_ec, header = Zcolhdr, name = 'Z_ec'))
        master_hdulist.append(fits.ImageHDU(data = Z_s_ec, header = Zcolhdr, name = 'Z_s_ec'))


        master_hdulist.append(fits.ImageHDU(data = Z_pdf, header = Zcolhdr, name = 'Z_pdf'))
        master_hdulist.append(fits.ImageHDU(data = Z_pdf_s, header = Zcolhdr, name = 'Z_pdf_s'))
        master_hdulist.append(fits.ImageHDU(data = Z_pdf_ec, header = Zcolhdr, name = 'Z_pdf_ec'))
        master_hdulist.append(fits.ImageHDU(data = Z_pdf_s_ec, header = Zcolhdr, name = 'Z_pdf_s_ec'))



        fits_name = out_dir + '/%s_%s_metals.fits'%(field, di)
        print ('\tSaving to ' + fits_name)
        thdulist = fits.HDUList(master_hdulist)
        thdulist.writeto(fits_name, overwrite = True)







