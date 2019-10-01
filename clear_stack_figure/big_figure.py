import astropy
from astropy.io import fits, ascii
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
import os
import astropy
from astropy.convolution import Box1DKernel, convolve_fft
plt.ioff()
plt.close('all')


kern = Box1DKernel(1)

zmin = 0.
zmax = 4.
zbins = 10.**linspace(np.log10(1 + zmin), np.log10(1 + zmax), 250.) - 1.



G102_waves = fits.open('/Users/rsimons/Desktop/clear/grizli_v2.1/1Dspec/GS4/GS4_22868.1D.fits')['G102'].data['wave']
G141_waves = fits.open('/Users/rsimons/Desktop/clear/grizli_v2.1/1Dspec/GS4/GS4_22868.1D.fits')['G141'].data['wave']

G102_minx = 14
G102_maxx = -12

G141_minx = 24
G141_maxx = -11



G102_waves = G102_waves[G102_minx:G102_maxx]
G141_waves = G141_waves[G141_minx:G141_maxx]

nG102 = len(G102_waves)
nG141 = len(G141_waves)


full_wave = concatenate((G102_waves, G141_waves))


#fields = ['ERSPRIME', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7', 'GS1', 'GS2', 'GS3', 'GS4', 'GS5']

fields = ['ERSPRIME', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7', 'GS1', 'GS2', 'GS3', 'GS4', 'GS5']

nrows = len(zbins)
ncols = len(full_wave)

if True:
    cat_gds = ascii.read('/Users/rsimons/Desktop/clear/Catalogs/goodss_3dhst.v4.4.cat')
    cat_gdn = ascii.read('/Users/rsimons/Desktop/clear/Catalogs/goodsn_3dhst.v4.4.cat')


mag_gds = 25.0-2.5*log10(cat_gds['f_H'])
mag_gdn = 25.0-2.5*log10(cat_gdn['f_H'])
id_gds = cat_gds['id']
id_gdn = cat_gdn['id']

lines = ['ArIII-7138',
'CIII-1908',
'CIV-1549',
'H8',
'H9',
'Ha',
'Hb',
'Hd',
'HeI-1083',
'HeI-5877',
'HeII-1640',
'Hg',
'Lya',
'MgII',
'NIII-1750',
'NIV-1487',
'NV-1240',
'NeIII-3867',
'NeV-3346',
'NeVI-3426',
'OI-6302',
'OII',
'OII-7325',
'OIII',
'OIII-1663',
'OIII-4363',
'PaB',
'SII',
'SIII']


if True:
    line_cat_dir = '/Users/rsimons/Desktop/clear/Catalogs/grizli_v2.1_cats'
    big_fig = np.empty((nrows, ncols), dtype = 'object')
    for i in arange(nrows):
        for j in arange(ncols):
            big_fig[i,j] = []


    for f, field in enumerate(fields[0:1]):
        print (field)
        spec_fits_dir = '/Users/rsimons/Desktop/clear/grizli_v2.1/1Dspec/%s'%field
        cat = fits.open(line_cat_dir + '/' + field+ '_lines_grizli.fits')
        ndet = zeros(len(cat[1].data['ID']))
        sns = zeros((len(cat[1].data['ID']), len(lines)))

        '''
        for l, line in enumerate(lines):
            sns[:,l] = cat[1].data['%s_FLUX'%line]/cat[1].data['%s_FLUX_ERR'%line]
            gd = where(sns[:,l] > 5.)[0]
            ndet[gd] +=1


        ''' 
        if 'N' in field: 
            mag_3d = mag_gdn
            id_3d = id_gdn
        if 'S' in field: 
            mag_3d = mag_gds
            id_3d = id_gds


        for cc, c in enumerate(cat[1].data):
            di = c['ID']
            z = c['z_50']


            cat_match = where(id_3d == di)[0]
            mag_cat = mag_3d[cat_match]

            spec_fits_file = spec_fits_dir + '/%s_%.5i.1D.fits'%(field, di)
            if (z > zmin) & (os.path.isfile(spec_fits_file)) & (z < zmax) & (ndet[cc] > -1.) & (mag_cat < 25):
                spec_fits = fits.open(spec_fits_file)

                row = argmin(abs(zbins - z))
                fl_info = spec_fits.info(False)
                names = [hdu[1] for hdu in fl_info]

                if 'G102' in names: 
                    flx  = spec_fits['G102'].data['flux'][G102_minx:G102_maxx]
                    err  = spec_fits['G102'].data['err'][G102_minx:G102_maxx]
                    cont = spec_fits['G102'].data['cont'][G102_minx:G102_maxx]

                    crit1 = np.isfinite(flx)
                    crit2 = ~isnan(flx)
                    crit3 = flx != 0.0
                    crit4 = err != 0.0
                    

                    gd = where(crit1 & crit2 & crit3 & crit4)[0]
                    chi2 = sum(abs(flx[gd] - cont[gd])/err[gd])/len(flx[gd])
                    if chi2 < 2:
                        for g in gd: big_fig[row, :nG102][g].append((flx[g] - cont[g])/cont[g])

                if 'G141' in names: 
                    flx  = spec_fits['G141'].data['flux'][G141_minx:G141_maxx]
                    err  = spec_fits['G141'].data['err'][G141_minx:G141_maxx]
                    cont = spec_fits['G141'].data['cont'][G141_minx:G141_maxx]


                    crit1 = np.isfinite(flx)
                    crit2 = ~isnan(flx)
                    crit3 = flx != 0.0
                    crit4 = err != 0.0
                    

                    gd = where(crit1 & crit2 & crit3 & crit4)[0]
                    chi2 = sum(abs(flx[gd] - cont[gd])/err[gd])/len(flx[gd])
                    if chi2 < 2:
                        for g in gd: big_fig[row, nG102:][g].append((flx[g] - cont[g])/cont[g])


if True:
    big_fig_new = np.zeros((nrows, ncols))

    for i in arange(nrows):
        for j in arange(ncols):
            big_fig_new[i,j] = np.nanmedian(np.array(big_fig[i,j]))

#np.save('big_fig.npy', big_fig)
#np.save('big_fig_new.npy', big_fig_new)
#big_fig = np.load('big_fig.npy', allow_pickle = True)[()]
#big_fig_new = np.load('big_fig_new.npy', allow_pickle = True)[()]

cmap = plt.cm.Greys_r
cmap.set_bad('k', 1.)






fig, ax = plt.subplots(1,1, figsize = (ncols/20., nrows/20.))


rvl_arr = big_fig_new.ravel()
rvl_arr = rvl_arr[~isnan(rvl_arr)]
sort_rvl = sort(rvl_arr)
vmin = sort_rvl[int(0.10 * len(sort_rvl))]
vmax = sort_rvl[int(0.95 * len(sort_rvl))]


ax.imshow(big_fig_new, cmap = cmap, origin = 'lower', vmin = 0., vmax = 0.5)#vmin*0.0000001, vmax = vmax)
ax.axis('off')
#ax.set_ylim(0, 100)



figdir = '/Users/rsimons/Desktop/clear/figures'
fig.subplots_adjust(left = 0.0, right = 1.0, top = 1.0, bottom = 0.0)
fig.savefig(figdir + '/big_fig_new.png')
#plt.close('all')












