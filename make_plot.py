import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits
plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True


cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/z_radius.cat',dtype = 'str')
ma_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/ma17.cat')
wang_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/wang17.cat')
wang18_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/wang18.cat')
jones_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/jones+13.cat')
swinbank_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/swinbank12.cat')

GN1_physcat = np.loadtxt('/Users/rsimons/Desktop/clear/Catalogs/GN1_physcat.cat')
GN2_physcat = np.loadtxt('/Users/rsimons/Desktop/clear/Catalogs/GN2_physcat.cat')
GN3_physcat = np.loadtxt('/Users/rsimons/Desktop/clear/Catalogs/GN3_physcat.cat')





if False:
    with PdfPages('/Users/rsimons/Dropbox/rcs_clear/z_radius_plots/z_radius.pdf') as pdf:

        fig, ax = plt.subplots(1,1, figsize = (6.0, 3))
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
            if fld == 'GN1': physcat = GN1_physcat
            if fld == 'GN2': physcat = GN2_physcat
            if fld == 'GN3': physcat = GN3_physcat

            #mstar = physcat[where(physcat[:,0] == int(c[1]))[0][0],1]
            mstar = float(c[-5])
            print mstar
            ax.errorbar(mstar, float(c[2]), yerr = float(c[3]), fmt = 'o', color = 'red', markeredgecolor = 'black', ms = 10)
        ax.errorbar(-99, float(c[2]), yerr = float(c[3]), fmt = 'o', color = 'red', markeredgecolor = 'black', ms = 10, label = 'CLEAR G102\n+ archival G141\n(2 of 10 pointings)', zorder = 10)

        ax.annotate(r'$z\,\sim\,1.5$', (0.05, 0.85), xycoords = 'axes fraction', fontsize = 20, fontweight = 'bold')
        ax.set_xlabel(r'$\log$ M$_{*}$ (M$_{\odot}$)', fontsize = 15)
        ax.set_ylabel(r'$\frac{\Delta \log(O/H)}{\Delta R}$ (dex kpc$^{-1}$)', rotation = 90, fontsize = 15)
        ax.legend(bbox_to_anchor=(1.0, 1.05), frameon = False)

        ax.set_ylim(-0.33, 0.3)
        ax.set_xlim(8.0, 11.5)

        fig.subplots_adjust(bottom = 0.20, left = 0.15, right = 0.70, top = 0.95)
        pdf.savefig()


if True:
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





