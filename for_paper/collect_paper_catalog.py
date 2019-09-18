import astropy
from astropy.io import fits, ascii
import numpy as np
from numpy import *



grizli_cat_dir = '/Users/rsimons/Desktop/clear/catalogs/grizli_v2.1_cats'

#if False:

gfit_gds_cat_125 = ascii.read('/Users/rsimons/Desktop/clear/catalogs/galfit_allfields/goodss/goodss_3dhst.v4.1_f125w.galfit')
gfit_gds_cat_160 = ascii.read('/Users/rsimons/Desktop/clear/catalogs/galfit_allfields/goodss/goodss_3dhst.v4.1_f160w.galfit')

gfit_gdn_cat_125 = ascii.read('/Users/rsimons/Desktop/clear/catalogs/galfit_allfields/goodsn/goodsn_3dhst.v4.1_f125w.galfit')
gfit_gdn_cat_160 = ascii.read('/Users/rsimons/Desktop/clear/catalogs/galfit_allfields/goodsn/goodsn_3dhst.v4.1_f160w.galfit')


eazy_gdn_zout = fits.open('/Users/rsimons/Desktop/clear/catalogs/goodsn_3dhst.v4.4.cats/Eazy/goodsn_3dhst.v4.4.zout.fits')
eazy_gdn_data = fits.open('/Users/rsimons/Desktop/clear/catalogs/goodsn_3dhst.v4.4.cats/Eazy/goodsn_3dhst.v4.4.data.fits')


eazy_gds_zout = fits.open('/Users/rsimons/Desktop/clear/catalogs/goodss_3dhst.v4.4.cats/Eazy/goodss_3dhst.v4.4.zout.fits')
eazy_gds_data = fits.open('/Users/rsimons/Desktop/clear/catalogs/goodss_3dhst.v4.4.cats/Eazy/goodss_3dhst.v4.4.data.fits')




gds_xray = ascii.read('/Users/rsimons/Desktop/clear/catalogs/clear_gdsxray.dat')
gdn_xray = ascii.read('/Users/rsimons/Desktop/clear/catalogs/clear_gdnxray.dat')

sample_cat = np.loadtxt('/Users/rsimons/Desktop/clear/catalogs/any_sample.cat', dtype = 'str')

whitaker_gds = np.loadtxt('/Users/rsimons/Desktop/clear/catalogs/sfr_3dhst.v4.1/goodss_3dhst.v4.1.sfr')
whitaker_gdn = np.loadtxt('/Users/rsimons/Desktop/clear/catalogs/sfr_3dhst.v4.1/goodsn_3dhst.v4.1.sfr')




f = open('/Users/rsimons/Desktop/clear/catalogs/simons_sample.cat', 'w+')


f.write('# field id ra dec z_02 z_16 z_50 z_84 z_97 z_map z_risk sfr_w sfr_w_ir sfr_w_uv sfr lmass z025_eazy z160_eazy z500_eazy z840_eazy z975_eazy pa_125 dpa_125 q_125 dq_125 n_125 dn_125 re_125 dre_125 pa_160 dpa_160 q_160 dq_160 n_160 dn_160 re_160 dre_160 xclass\n')


for o, obj in enumerate(sample_cat):
    fld = obj[0]
    di = int(obj[1])
    print (o, fld, di)



    grizli_cat = fits.open('/Users/rsimons/Desktop/clear/catalogs/grizli_v2.1_cats/%s_lines_grizli.fits'%fld)
    grizli_cat = grizli_cat[1].data[grizli_cat[1].data['ID'] == int(di)]


    if 'S' in fld:
        gfit_125 = gfit_gds_cat_125[gfit_gds_cat_125['NUMBER'] == di]
        gfit_160 = gfit_gds_cat_160[gfit_gds_cat_160['NUMBER'] == di]
        xray_cat =         gds_xray[        gds_xray['ID'] == di]
        eazy_zout=    eazy_gds_zout[1].data[eazy_gds_zout[1].data['id'] == di]
        whitaker =     whitaker_gds[    whitaker_gds[:,0] == di]

    else:
        gfit_125 = gfit_gdn_cat_125[gfit_gdn_cat_125['NUMBER'] == di]
        gfit_160 = gfit_gdn_cat_160[gfit_gdn_cat_160['NUMBER'] == di]
        xray_cat =         gdn_xray[        gdn_xray['ID'] == di]
        eazy_zout=    eazy_gdn_zout[1].data[eazy_gdn_zout[1].data['id'] == di]
        whitaker =     whitaker_gdn[    whitaker_gdn[:,0] == di]

    ra  = float(grizli_cat['RA'])
    dec = float(grizli_cat['DEC'])
    if len(whitaker) > 0:
        wsfr     = float(whitaker[0][1])
        wsfr_ir  = float(whitaker[0][2])
        wsfr_uv  = float(whitaker[0][3])

    else:
        wsfr    = nan
        wsfr_ir = nan
        wsfr_uv = nan


    if len(eazy_zout) > 0:
        esfr   = float(eazy_zout['sfr'])
        elmass = float(log10(eazy_zout['mass']))

        z025   = float(eazy_zout['z025'])
        z160   = float(eazy_zout['z160'])
        z500   = float(eazy_zout['z500'])
        z840   = float(eazy_zout['z840'])
        z975   = float(eazy_zout['z975'])
    else:
        esfr    = nan
        elmass  = nan
        z025    = nan
        z160    = nan
        z500    = nan
        z840    = nan
        z975    = nan

    z_02   = float(grizli_cat['z_02'])
    z_16   = float(grizli_cat['z_16'])
    z_50   = float(grizli_cat['z_50'])
    z_84   = float(grizli_cat['z_84'])
    z_97   = float(grizli_cat['z_97'])


    z_map   = float(grizli_cat['z_MAP'])
    z_risk   = float(grizli_cat['z_RISK'])

    if len(gfit_125) > 0:
        pa_125  =   float(gfit_125['pa'].data)
        dpa_125 =   float(gfit_125['dpa'].data)
        q_125   =   float(gfit_125['q'].data)
        dq_125  =   float(gfit_125['dq'].data)
        n_125   =   float(gfit_125['n'].data)
        dn_125  =   float(gfit_125['dn'].data)
        re_125  =   float(gfit_125['re'].data)
        dre_125 =   float(gfit_125['dre'].data)
    else:
        pa_125  =  nan
        dpa_125 =  nan
        q_125   =  nan
        dq_125  =  nan
        n_125   =  nan
        dn_125  =  nan
        re_125  =  nan
        dre_125 =  nan

    if len(gfit_160) > 0:
        pa_160  =   float(gfit_160['pa'].data)
        dpa_160 =   float(gfit_160['dpa'].data)
        q_160   =   float(gfit_160['q'].data)
        dq_160  =   float(gfit_160['dq'].data)
        n_160   =   float(gfit_160['n'].data)
        dn_160  =   float(gfit_160['dn'].data)
        re_160  =   float(gfit_160['re'].data)
        dre_160 =   float(gfit_160['dre'].data)
    else:
        pa_160  =   nan
        dpa_160 =   nan
        q_160   =   nan
        dq_160  =   nan
        n_160   =   nan
        dn_160  =   nan
        re_160  =   nan
        dre_160 =   nan

    if len(xray_cat) > 0:
        xclass = xray_cat[0]['Xclass']
    else:
        xclass = 'NA'

    f.write('%8s  %i  %.7f  %.7f  \
             %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  %.5f  \
             %.5f  %.5f  %.5f  \
             %.5f  %.5f  \
             %.5f  %.5f  %.5f  %.5f  %.5f  \
             %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f  %.3f \
             %.3f  %.3f  %.3f  %.3f  %s\n'\
             %(fld, di, ra, dec, \
               z_02, z_16, z_50, z_84, z_97, z_map, z_risk, \
               wsfr, wsfr_ir, wsfr_uv, \
               esfr, elmass, \
               z025,z160,z500,z840,z975,\
               pa_125, dpa_125, q_125, dq_125, n_125, dn_125, re_125, dre_125 , pa_160, dpa_160, q_160, \
               dq_160, n_160, dn_160, re_160, dre_160, xclass))

f.close()



















