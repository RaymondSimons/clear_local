import numpy as np
import astropy
from astropy.io import fits
from matplotlib.ticker import NullFormatter
plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True



fields = np.array(['GS1','GS2', 'GS3', 'GS4', 'GS5', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7'])

gdn_zout = fits.open('/Users/rsimons/Dropbox/rcs_clear/catalogs/goodsn_3dhst.v4.4.cats/Eazy/goodsn_3dhst.v4.4.zout.fits')
gds_zout = fits.open('/Users/rsimons/Dropbox/rcs_clear/catalogs/goodss_3dhst.v4.4.cats/Eazy/goodss_3dhst.v4.4.zout.fits')




fig1, ax1 = plt.subplots(1,3, figsize = (12,4))
fig2, ax2 = plt.subplots(1,3, figsize = (12,4))

for ax in [ax1, ax2]:
    ax[0].plot([0,10], [0,10], '-', color = 'blue')
    ax[1].plot([0,10], [0,0], '-', color = 'blue')




fh1, fh2_0, fh2_1 = [], [], []

for f, field in enumerate(fields):
    if   'N' in field: zout = gdn_zout
    elif 'S' in field: zout = gds_zout
    zs = zout[1].data['z_spec']
    zp = zout[1].data['z_phot']

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
            if False:
                if nlines[i]   == 0: clr = 'black'
                elif nlines[i] == 1: clr = 'green'
                elif nlines[i] >  1: clr = 'red'
            else: clr = 'black'
            #ax[0].errorbar(1 + zs[gd], 1+zg[i], yerr = zg_le[i], fmt = 'o', color = clr, markersize = 2.)
            ax1[0].plot(1 + zs[gd], 1+zg[i], '.', color = clr, markersize = 3, alpha = 0.2)
            ax1[1].plot(1 + zs[gd], (zg[i] - zs[gd])/(1+zs[gd]),'.', color = clr, markersize = 3, alpha = 0.2)

            fh1.append((zg[i] - zs[gd])/(1+zs[gd]))
        #clr = 'black'
        if nlines[i]   == 0: 
            clr = 'black'
            fh2_0.append((zg[i] - zp[gd])/(1+zg[i]))

        elif nlines[i] > 0: 
            clr = 'red'
            fh2_1.append((zg[i] - zp[gd])/(1+zg[i]))

        ax2[0].plot(1 + zg[i], 1+zp[gd], '.', color = clr, markersize = 2, alpha = 0.2)
        ax2[1].plot(1 + zg[i], (zg[i] - zp[gd])/(1+zg[i]),'.', color = clr, markersize = 2, alpha = 0.2)









ax1[2].hist(fh1, bins = linspace(-1,1, 100), color = 'black', orientation = 'horizontal')
ax1[0].set_xlabel(r'z$_{spec}$', fontsize = 18)
ax1[1].set_xlabel(r'z$_{spec}$', fontsize = 18)
ax1[1].set_ylabel(r'(z$_{grism}$ - z$_{spec}$)/(1 + z$_{spec}$)', fontsize = 18)


ax2[2].hist(fh2_0, bins = linspace(-1,1, 100), color = 'black', orientation = 'horizontal', histtype = 'step')
ax2[2].hist(fh2_1, bins = linspace(-1,1, 100), color = 'red', orientation = 'horizontal', histtype = 'step')
ax2[0].set_xlabel(r'z$_{grism}$', fontsize = 18)
ax2[1].set_xlabel(r'z$_{grism}$', fontsize = 18)
ax2[1].set_ylabel(r'(z$_{grism}$ - z$_{phot}$)/(1 + z$_{grism}$)', fontsize = 18)



ax1[0].set_ylabel(r'z$_{grism}$', fontsize = 18)
ax2[0].set_ylabel(r'z$_{phot}$', fontsize = 18)



#Ha confused with OIII

for ax in [ax1, ax2]:
    ax[2].axis('off')
    ax[1].set_xlim(1,6)
    ax[1].set_xscale('log')
    ax[1].set_ylim(-1,1)
    ax[2].set_ylim(-1,1)



    ax[0].set_xscale('log')
    ax[0].set_yscale('log')


    for a in [ax[0], ax[1]]:
        #a.yaxis.set_major_formatter(NullFormatter())
        a.yaxis.set_minor_formatter(NullFormatter())
        #a.xaxis.set_major_formatter(NullFormatter())
        a.xaxis.set_minor_formatter(NullFormatter())
        a.set_xlim(1, 6)
        a.set_xticks(arange(1, 7))
        a.set_xticklabels('%i'%(i-1) for i in arange(1, 7))

    ax[0].set_ylim(1, 6)
    ax[0].set_yticks(arange(1, 7))


    ax[0].set_yticklabels('%i'%(i-1) for i in arange(1, 7))


z_temp = arange(8)
ax1[0].plot(1 + z_temp, (1+z_temp)*6563./(5007), 'k--', alpha = 0.2)
ax1[0].plot(1 + z_temp, (1+z_temp)*5007./(6563), 'k--', alpha = 0.2)
ax1[0].plot(1 + z_temp, (1+z_temp)*3726./(5007), 'k--', alpha = 0.2)
ax1[0].plot(1 + z_temp, (1+z_temp)*5007./(3726), 'k--', alpha = 0.2)
ax1[0].plot(1 + z_temp, (1+z_temp)*3726./(6563), 'k--', alpha = 0.2)
ax1[0].plot(1 + z_temp, (1+z_temp)*6563./(3726), 'k--', alpha = 0.2)
ax1[0].annotate(r'H$\alpha$ : [OIII]', (0.70, 0.95) , xycoords = 'axes fraction', color = 'grey', rotation = 45)
ax1[0].annotate(r'[OII] : [OIII]',(0.51, 0.88) , xycoords = 'axes fraction', color = 'grey', rotation = 45)
ax1[0].annotate(r'H$\alpha$ : [OII]', (0.37, 0.86), xycoords = 'axes fraction', color = 'grey', rotation = 45)




fig1.tight_layout()
fig1.savefig('/Users/rsimons/Desktop/clear/figures/zspec_zgrism_v2.1.png', dpi = 300)


ax2[0].annotate('Nlines = 0', (0.1, 0.9), xycoords = 'axes fraction', fontsize = 20, color = 'black')
ax2[0].annotate('Nlines $>$ 0', (0.1, 0.8), xycoords = 'axes fraction', fontsize = 20, color = 'red')
fig2.tight_layout()
fig2.savefig('/Users/rsimons/Desktop/clear/figures/zphot_zgrism_v2.1.png', dpi = 300)














