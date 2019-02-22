import numpy as np
import astropy
from astropy.io import fits
plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True



fields = np.array(['GS1','GS2', 'GS3', 'GS5', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7'])

gdn_zout = fits.open('/Users/rsimons/Dropbox/rcs_clear/catalogs/goodsn_3dhst.v4.4.cats/Eazy/goodsn_3dhst.v4.4.zout.fits')
gds_zout = fits.open('/Users/rsimons/Dropbox/rcs_clear/catalogs/goodss_3dhst.v4.4.cats/Eazy/goodss_3dhst.v4.4.zout.fits')




fig, ax = plt.subplots(1,3, figsize = (12,4))

ax[0].set_xlim(0,4)
ax[1].set_xlim(0,4)
ax[0].set_ylim(0,4)
ax[1].set_ylim(-1,1)
ax[2].set_ylim(-1,1)

ax[0].plot([0,10], [0,10], '-', color = 'grey')

fh = []
for f, field in enumerate(fields):
    if   'N' in field: zout = gdn_zout
    elif 'S' in field: zout = gds_zout
    zs = zout[1].data['z_spec']
    grizli_cat = fits.open('/Users/rsimons/Desktop/clear/Catalogs/grizli_v2.1_cats/%s_lines_grizli.fits'%field)
    #zgrism = grizli_cat[1].

    nlines = zeros(len(grizli_cat[1].data['ID']))
    names = ['Lya_FLUX',
             'CIV_FLUX',
             'MgII_FLUX',
             'OII_FLUX',
             'Hd_FLUX',
             'Hg_FLUX',
             'OIIIx_FLUX',
             'HeII_FLUX',
             'Hb_FLUX',
             'OIII_FLUX',
             'Ha_FLUX',
             'SII_FLUX',
             'SIII_FLUX',
             'HeI_FLUX',
             'HeIb_FLUX',
             'NeIII_FLUX',
             'NeV_FLUX',
             'NeVI_FLUX']


    for name in names:
        gd = where((grizli_cat[1].data[name] > 5 * grizli_cat[1].data[name + '_ERR']) & (grizli_cat[1].data[name] > 0.))[0]
        nlines[gd]+=1


    di =  grizli_cat[1].data['ID']
    zg = grizli_cat[1].data['z_50']
    zg_ue =  grizli_cat[1].data['z_84'] - grizli_cat[1].data['z_50']
    zg_le = grizli_cat[1].data['z_50'] - grizli_cat[1].data['z_16']

    for i in arange(len(zg)):
        gd = where(zout[1].data['id'] == di[i])[0][0]
        if zs[gd] > 0:
            ax[0].errorbar(zs[gd], zg[i], yerr = zg_le[i], fmt = 'o', color = 'black', markersize = 3)
            ax[1].plot(zs[gd], zg[i] - zs[gd],'k.')
            if False:
                if nlines[i] == 0:
                    ax[1].plot(zs[gd], zg[i] - zs[gd],'k.')
                elif nlines[i] == 1:
                    ax[1].plot(zs[gd], zg[i] - zs[gd],'g.')
                elif nlines[i] > 1:
                    ax[1].plot(zs[gd], zg[i] - zs[gd],'r.')

            fh.append( zg[i] - zs[gd])

ax[2].hist(fh, bins = linspace(-1,1, 100), color = 'black', orientation = 'horizontal')
ax[2].axis('off')
ax[0].set_xlabel(r'z$_{spec}$', fontsize = 18)
ax[1].set_xlabel(r'z$_{spec}$', fontsize = 18)
ax[0].set_ylabel(r'z$_{grism}$', fontsize = 18)
ax[1].set_ylabel(r'z$_{grism}$ - z$_{spec}$', fontsize = 18)

fig.tight_layout()
fig.savefig('/Users/rsimons/Desktop/clear/figures/zspec_zgrism_v2.1.png', dpi = 300)