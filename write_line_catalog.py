import astropy
from astropy.io import fits
import glob
from glob import glob
import numpy as np
from numpy import *

for field in ['GS1']:#, 'GS2', 'GS3', 'GS5', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7']:



    fls = glob('/user/rsimons/grizli_extractions/%s/*/Prep/*.full.fits'%field)
    cat = open('/user/rsimons/grizli_extractions/Catalogs/%s_lines_grizli.cat'%field, 'w+')

    fits_name = '/user/rsimons/grizli_extractions/Catalogs/%s_lines_grizli.fits'%field


    lines = ['OII', 'OIII', 'Ha', 'Hb']


    cat.write('#(0) ID\n')
    cat.write('#(1) ra\n')
    cat.write('#(2) dec\n')
    cat.write('#(3) number of lines\n')
    cat.write('#(4) OII flux, 1e-17 erg/s/cm2\n')
    cat.write('#(5) OII flux err, 1e-17 erg/s/cm2\n')
    cat.write('#(6) OIII flux, 1e-17 erg/s/cm2\n')
    cat.write('#(7) OIII flux err, 1e-17 erg/s/cm2\n')
    cat.write('#(8) Ha flux, 1e-17 erg/s/cm2\n')
    cat.write('#(9) Ha flux err, 1e-17 erg/s/cm2\n')
    cat.write('#(10) Hb flux, 1e-17 erg/s/cm2\n')
    cat.write('#(11) Hb flux err, 1e-17 erg/s/cm2\n')
    cat.write('#(12) z_50\n')
    cat.write('#(13) z_02\n')
    cat.write('#(14) z_16\n')
    cat.write('#(15) z_84\n')
    cat.write('#(16) z_97\n')

    fluxs = zeros((4,2, len(fls))) - 99.
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



        lines_f = a[0].header['HASLINES'].split(' ')
        if (len(lines_f) > 0) & (lines_f[0] !=''):
            j = -99

            zs[0, f] = a[1].header['Z50']
            zs[1, f] = a[1].header['Z02']
            zs[2, f] = a[1].header['Z16']
            zs[3, f] = a[1].header['Z84']
            zs[4, f] = a[1].header['Z97']

            for l, ln in enumerate(lines_f):
                flux_ln =  a[0].header['FLUX%.3i'%(l+1)]
                eflux_ln =  a[0].header['ERR%.3i'%(l+1)]
                ln_name =  a[0].header['LINE%.3i'%(l+1)]			
                if ln_name == 'OII': 		j = 0
                elif ln_name == 'OIII':  	j = 1
                elif ln_name == 'Ha': 		j = 2
                elif ln_name == 'Hb': 		j = 3

                if j != -99:
                    fluxs[j,0,f] = flux_ln * 1.e17
                    fluxs[j,1,f] = eflux_ln  * 1.e17
                j = -99


        cat.write('%i\t\t%.8f\t\t%.8f\t\t%i\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%.2f\n'%
            (IDs[f],  ras[f], decs[f], nlines[f], fluxs[0,0,f], fluxs[0,1,f], fluxs[1,0,f], fluxs[1,1,f], fluxs[2,0,f], fluxs[2,1,f], fluxs[3,0,f], fluxs[3,1,f], zs[0, f], zs[1, f], zs[2, f], zs[3, f], zs[4, f]))

    master_hdulist = []
    prihdr = fits.Header()

    prihdu = fits.PrimaryHDU(header=prihdr)    
    master_hdulist.append(prihdu)

    colhdr = fits.Header()

    col_list = [fits.Column(name='ID', format = 'D', array=array(IDs)),
    fits.Column(name='RA', format = 'D', array=array(ras)),
    fits.Column(name='DEC',format = 'D', array=array(decs)),
    fits.Column(name='nlines',format = 'D', array=array(nlines)),
    fits.Column(name='OII_f',format = 'D', array=fluxs[0,0,:]),
    fits.Column(name='OII_e',format = 'D', array=fluxs[0,1,:]),
    fits.Column(name='OIII_f',format = 'D', array=fluxs[1,0,:]),
    fits.Column(name='OIII_e',format = 'D', array=fluxs[1,1,:]),
    fits.Column(name='Ha_f',format = 'D', array=fluxs[2,0,:]),
    fits.Column(name='Ha_e',format = 'D', array=fluxs[2,1,:]),
    fits.Column(name='Hb_f',format = 'D', array=fluxs[3,0,:]),
    fits.Column(name='Hb_e',format = 'D', array=fluxs[3,1,:]),
    fits.Column(name='z_50',format = 'D', array=zs[0,:]),
    fits.Column(name='z_02',format = 'D', array=zs[1,:]),
    fits.Column(name='z_16',format = 'D', array=zs[2,:]),
    fits.Column(name='z_84',format = 'D', array=zs[3,:]),
    fits.Column(name='z_97',format = 'D', array=zs[4,:])]

    coldefs = fits.ColDefs(col_list)
    table_hdu = fits.BinTableHDU.from_columns(coldefs)
    master_hdulist.append(table_hdu)
    thdulist = fits.HDUList(master_hdulist)
    thdulist.writeto(fits_name, overwrite = True)








    cat.close()











