import astropy
from astropy.io import fits
import glob
from glob import glob
import numpy as np
from numpy import *

for field in ['GS1']:#, 'GS2', 'GS3', 'GS5', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7']:

    fls = glob('/user/rsimons/grizli_extractions/%s/*/Prep/*.full.fits'%field)
    fits_name = '/user/rsimons/grizli_extractions/Catalogs/%s_lines_grizli.fits'%field

    lines = ['Lya',
             'CIV',
             'MgII',
             'OII',
             'Hd',
             'Hg',
             'OIIIx',
             'HeII',
             'Hb',
             'OIII',
             'Ha',
             'SII',
             'SIII',
             'HeI',
             'HeIb',
             'NeIII',
             'NeV',
             'NeVI',
             'OI']

    fluxs = zeros((len(lines),2, len(fls))) - 99.
    exptime = zeros((2, len(fls)))

    zs = zeros((5, len(fls))) - 99.
    IDs = []
    ras = []
    decs = []
    nlines = []
    for f, fl in enumerate(fls):
        print f, len(fls)
        a = fits.open(fl)
        IDs.append(a[0].header['ID'])
        ras.append(a[0].header['ra'])
        decs.append(a[0].header['dec'])
        nlines.append(a[0].header['NUMLINES'])
        exptime[0,f] = float(a[0].header['T_G102'])
        exptime[1,f] = float(a[0].header['T_G141'])
        lines_f = a[0].header['HASLINES'].split(' ')
        print lines_f
        if (len(lines_f) > 0) & (lines_f[0] !=''):
            zs[0, f] = a[1].header['Z50']
            zs[1, f] = a[1].header['Z02']
            zs[2, f] = a[1].header['Z16']
            zs[3, f] = a[1].header['Z84']
            zs[4, f] = a[1].header['Z97']

            for l, ln in enumerate(lines_f):
                flux_ln =  a[0].header['FLUX%.3i'%(l+1)]
                eflux_ln =  a[0].header['ERR%.3i'%(l+1)]
                ln_name =  a[0].header['LINE%.3i'%(l+1)]
                for ll, line in lines:	
                    if ln_name == line:
                        fluxs[j,0,f] = flux_ln * 1.e17
                        fluxs[j,1,f] = eflux_ln  * 1.e17

    master_hdulist = []
    prihdr = fits.Header()

    prihdu = fits.PrimaryHDU(header=prihdr)    
    master_hdulist.append(prihdu)

    colhdr = fits.Header()

    col_list = [fits.Column(name='ID', format = 'D', array=array(IDs)),
    fits.Column(name='RA', format = 'D', array=array(ras)),
    fits.Column(name='DEC',format = 'D', array=array(decs)),
    fits.Column(name='nlines',format = 'D', array=array(nlines)),
    fits.Column(name='z_50',format = 'D', array=zs[0,:]),
    fits.Column(name='z_02',format = 'D', array=zs[1,:]),
    fits.Column(name='z_16',format = 'D', array=zs[2,:]),
    fits.Column(name='z_84',format = 'D', array=zs[3,:]),
    fits.Column(name='z_97',format = 'D', array=zs[4,:])]
    for l, line in lines:
        col_list.append(fits.Column(name='%s_FLUX'%line,format = 'D', array=fluxs[l,0,:]))
        col_list.append(fits.Column(name='%s_FLUX_ERR'%line,format = 'D', array=fluxs[l,1,:]))
    fits.Column(name='T_G102',format = 'D', array=exptime[0,:])
    fits.Column(name='T_G141',format = 'D', array=exptime[1,:])

    coldefs = fits.ColDefs(col_list)
    table_hdu = fits.BinTableHDU.from_columns(coldefs)
    master_hdulist.append(table_hdu)
    thdulist = fits.HDUList(master_hdulist)
    thdulist.writeto(fits_name, overwrite = True)








