import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord 
import astropy.units as u
plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True


cat = np.loadtxt('/Users/rsimons/Desktop/clear/Catalogs/z_r_O32.cat',dtype = 'str')





ma_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/ma17.cat')
wang_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/wang17.cat')
wang18_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/wang18.cat')
jones_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/jones+13.cat')
swinbank_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/swinbank12.cat')

#GN1_physcat = np.loadtxt('/Users/rsimons/Desktop/clear/Catalogs/GN1_physcat.cat')
#GN2_physcat = np.loadtxt('/Users/rsimons/Desktop/clear/Catalogs/GN2_physcat.cat')
#GN3_physcat = np.loadtxt('/Users/rsimons/Desktop/clear/Catalogs/GN3_physcat.cat')


gs_fout = fits.open('/Users/rsimons/Desktop/clear/Catalogs/goodss_3dhst.v4.1.fout.FITS')
gn_fout = fits.open('/Users/rsimons/Desktop/clear/Catalogs/goodsn_3dhst.v4.1.fout.FITS')


gn_cat = fits.open('/Users/rsimons/Desktop/clear/Catalogs/goodsn_3dhst.v4.1.cats/Catalog/goodsn_3dhst.v4.1.cat.FITS')
gs_cat = fits.open('/Users/rsimons/Desktop/clear/Catalogs/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS')


x_gds = ascii.read('/Users/rsimons/Desktop/clear/Catalogs/xray_GDS.txt')
x_gdn = ascii.read('/Users/rsimons/Desktop/clear/Catalogs/xray_GDN.txt')


if True:
    with PdfPages('/Users/rsimons/Dropbox/rcs_clear/z_radius_plots/z_radius_new2.pdf') as pdf:

        fig, ax = plt.subplots(1,1, figsize = (10, 5))
        ax.errorbar(wang_cat[:,5], wang_cat[:,3], yerr =wang_cat[:,4], fmt = 'o', color = 'blue', label = 'Wang+ 17', zorder = 1)
        ax.errorbar(wang18_cat[:,3], wang18_cat[:,1], yerr =wang18_cat[:,2], fmt = 'o', color = 'darkblue', label = 'Wang+ 18', zorder = 1)
        ax.errorbar(jones_cat[:,0], jones_cat[:,1], yerr =jones_cat[:,2], fmt = 'o', color = 'skyblue', label = 'Jones+ 13', zorder = 1)
        ax.errorbar(swinbank_cat[:,0], swinbank_cat[:,1], yerr =swinbank_cat[:,2], fmt = 'o', color = 'darkgreen', label = 'Swinbank+ 12', zorder = 1)
        ax.errorbar(log10(ma_cat[:,1]), ma_cat[:,7], yerr =ma_cat[:,8], fmt = 's', markeredgecolor = 'black', markerfacecolor = "None", color = 'black', label = 'Ma+ 17; simulations', zorder = 1)

        ax.axhline(y = 0, linestyle = '-', color = 'grey', alpha = 0.4, linewidth = 2, zorder = 0)

        wuyts_x = linspace(10.0, 11.5, 100)
        wuyts_y = -0.017*(wuyts_x - 10) + 0.0

        ax.plot(wuyts_x, wuyts_y, '--', linewidth = 3, color = 'midnightblue', label = 'Wuyts+ 16', zorder = 1)
        to_pl = []
        for c in cat:
            fld = c[0]
            di = c[1]
            if 'N'   in fld: ct = gn_cat
            elif 'S' in fld: ct = gs_cat

            match_cat = where(int(di) == ct[1].data['id'])[0]
            ra_c = ct[1].data['ra'][match_cat]
            dec_c = ct[1].data['dec'][match_cat]

            if 'N'   in fld: 
                fout = gn_fout
                xcat = x_gdn
                deg_hms = SkyCoord(ra = ra_c * u.degree, dec = dec_c * u.degree)
                dist = sqrt(((xcat['RAm'] * 60. + xcat['RAs']) - (deg_hms.ra.hms.m * 60. + deg_hms.ra.hms.s))**2. + ((xcat['DEm'] * 60. + xcat['DEs']) - (deg_hms.dec.dms.m * 60. + deg_hms.dec.dms.s))**2.)
                tp = 'Type'

            elif 'S' in fld: 
                fout = gs_fout
                xcat = x_gds
                dist = sqrt((xcat['RAdeg'] - ra_c)**2. + (xcat['DEdeg'] - dec_c)**2.) * 3600.
                tp = 'OType'


            gd_xcat = argmin(dist)
            
            print fld, dist[gd_xcat]
            if dist[gd_xcat] > 2.: xobj = 'none'
            else: xobj = xcat[tp][gd_xcat]

            if xobj == 'AGN': mrker = '*'
            elif xobj == 'none': mrker = 's'
            elif xobj == 'Galaxy': mrker = 'o'
            elif xobj == 'Star': mrker = 'D'

            mstar = fout[1].data['lmass'][fout[1].data['id'] == int(di)]

            ax.errorbar(mstar, float(c[4]), yerr = float(c[5]), fmt = mrker, color = 'red', fillstyle = 'none', markeredgecolor = 'red', ms = 8)


        ax.errorbar(-99, -1, yerr = 0.01, fmt = 's', color = 'red', fillstyle = 'none', markeredgecolor = 'red', label = 'CLEAR, (N = 112)', zorder = 10)
        ax.errorbar(9.25, 0.0246, xerr = 0.25, yerr = 0.003, fmt = 'o', color = 'red', markeredgecolor = 'black', ms = 10, label = 'CLEAR, STACK', zorder = 10)
        ax.errorbar(9.75, 0.0163, xerr = 0.25, yerr = 0.003, fmt = 'o', color = 'red', markeredgecolor = 'black', ms = 10, zorder = 10)
        ax.errorbar(10.25, 0.0121, xerr = 0.25, yerr = 0.004,   fmt = 'o', color = 'red', markeredgecolor = 'black', ms = 10,  zorder = 10)
        ax.errorbar(10.75, 0.0055, xerr = 0.25, yerr = 0.008,  fmt = 'o', color = 'red', markeredgecolor = 'black', ms = 10, zorder = 10)








        ax.annotate(r'$0.7 < z < 1.5$', (0.63, 0.08), xycoords = 'axes fraction', fontsize = 25, fontweight = 'bold')
        ax.set_xlabel(r'$\log$ M$_{*}$ (M$_{\odot}$)', fontsize = 20)
        ax.set_ylabel(r'$\frac{\Delta \log(O/H)}{\Delta R}$ (dex kpc$^{-1}$)', rotation = 90, fontsize = 20)
        ax.legend(bbox_to_anchor=(1.0, 1.05), frameon = False, fontsize = 18)

        ax.set_ylim(-0.33, 0.50)
        ax.set_xlim(8.0, 11.5)

        fig.subplots_adjust(bottom = 0.20, left = 0.15, right = 0.70, top = 0.95)
        pdf.savefig()


if False:
    with PdfPages('/Users/rsimons/Dropbox/rcs_clear/z_radius_plots/mstar_sfr.pdf') as pdf:

            fig, ax = plt.subplots(1,1, figsize = (6.5, 3))


            for c in cat:
                fld = c[0]
                if fld == 'GN1': physcat = GN1_physcat
                if fld == 'GN2': physcat = GN2_physcat
                if fld == 'GN3': physcat = GN3_physcat

                #mstar = physcat[where(physcat[:,0] == int(c[1]))[0][0],1]
                mstar = float(c[-5])
                print mstar
                sfr = float(c[4])
                ax.plot(mstar, sfr, color = 'red', marker = 'o', zorder = 10, markeredgecolor = 'black', markersize = 10)


            ax.plot(-99, -99, color = 'red', marker = 'o', zorder = 10, markeredgecolor = 'black', label = 'CLEAR G102\n+ archival G141\n(2 of 10 pointings)')
            whitaker_all = np.loadtxt('/Users/rsimons/Downloads/goodsn_3dhst_v4.1.5_catalogs/goodsn_3dhst.v4.1.5.zbest.sfr')
            fast_all = fits.open('/Users/rsimons/Desktop/clear/Catalogs/goodsn_3dhst.v4.1.cats/Fast/goodsn_3dhst.v4.1.fout.FITS')




            whit = []

            for i in arange(len(whitaker_all)):
                good = where(fast_all[1].data['id'] == whitaker_all[i,0])[0][0]
                whit.append([fast_all[1].data['lmass'][good], whitaker_all[i,1], fast_all[1].data['z'][good]])




            whit_sfr = array([[8.8 ,  -0.03, 0.13],
                            [9.1 ,   0.17, 0.08],
                            [9.3 ,   0.38, 0.05],
                            [9.5 ,   0.64, 0.02],
                            [9.7 ,   0.81, 0.02],
                            [9.9 ,   1.02, 0.03],
                            [10.1,   1.18, 0.04],
                            [10.3,   1.35, 0.04],
                            [10.5,   1.47, 0.05],
                            [10.7,   1.58, 0.05],
                            [10.9,   1.69, 0.08],
                            [11.1,   1.74, 0.13],
                            [11.3,   1.81, 0.11]])

            whit = array(whit)
            ax.errorbar(whit_sfr[:,0], 10**whit_sfr[:,1], yerr = 10.**whit_sfr[:,2], fmt = 'o', color = 'grey', zorder = 1, label = 'GDN (medians)')
            gz = where((whit[:,2] > 1.0) & (whit[:,2] < 1.5))[0]
            ax.plot(whit[gz,0], whit[gz,1], marker = '.', linewidth = 0., color = 'grey', markersize = 1, zorder = 0, alpha = 0.4, label = 'GDN (all)')

            M_whit = linspace(9, 12, 1000)
            sfr = -24 + 4.17*M_whit - 0.16*M_whit**2.



            ax.set_xlabel(r'$\log$ M$_{*}$ (M$_{\odot}$) [FAST]', fontsize = 15)
            ax.set_ylabel(r'SFR$_{UV+IR}$ (M$_{\odot}$ yr$^{-1}$)', rotation = 90, fontsize = 15)
            ax.legend(bbox_to_anchor=(1.0, 1.05), frameon = False)
            ax.set_xlim(8.0, 11.5)
            ax.set_ylim(0.1, 300)
            ax.set_yscale('log')
            ax.legend(bbox_to_anchor=(0.99, 0.95), frameon = False)

            fig.subplots_adjust(bottom = 0.20, left = 0.18, right = 0.70, top = 0.95)
            pdf.savefig()





