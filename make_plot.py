import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord 
import astropy.units as u
import matplotlib as mpl
from astropy.cosmology import Planck15 as cosmo
import matplotlib.pyplot as plt
from numpy import *
import glob
from glob import glob
plt.rcParams['xtick.labelsize']=14
plt.rcParams['ytick.labelsize']=14
plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

cat = np.loadtxt('/Users/rsimons/Desktop/clear/Catalogs/z_r_O32.cat',dtype = 'str')
fls = glob('/Users/rsimons/Desktop/clear/izi_metal_profiles/fits/*npy')




ma_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/ma17.cat')
wang_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/wang17.cat')
wang18_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/wang18.cat')
jones_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/jones+13.cat')
swinbank_cat = np.loadtxt('/Users/rsimons/Dropbox/rcs_clear/catalogs/swinbank12.cat')

#GN1_physcat = np.loadtxt('/Users/rsimons/Desktop/clear/Catalogs/GN1_physcat.cat')
#GN2_physcat = np.loadtxt('/Users/rsimons/Desktop/clear/Catalogs/GN2_physcat.cat')
#GN3_physcat = np.loadtxt('/Users/rsimons/Desktop/clear/Catalogs/GN3_physcat.cat')


#gs_fout = fits.open('/Users/rsimons/Desktop/clear/Catalogs/goodss_3dhst.v4.4.fout.FITS')
#gn_fout = fits.open('/Users/rsimons/Desktop/clear/Catalogs/goodsn_3dhst.v4.4.fout.FITS')


#gn_cat = fits.open('/Users/rsimons/Desktop/clear/Catalogs/goodsn_3dhst.v4.4.cats/Catalog/goodsn_3dhst.v4.4.cat.FITS')
#gs_cat = fits.open('/Users/rsimons/Desktop/clear/Catalogs/goodss_3dhst.v4.4.cats/Catalog/goodss_3dhst.v4.4.cat.FITS')
count_all = 0
count_flat = 0
count_falling = 0
count_rising = 0

if False:
    gn_cat = fits.open('/Users/rsimons/Desktop/clear/Catalogs/goodsn_3dhst.v4.4.zout.fits')
    gs_cat = fits.open('/Users/rsimons/Desktop/clear/Catalogs/goodss_3dhst.v4.4.zout.fits')

    x_gds = ascii.read('/Users/rsimons/Desktop/clear/Catalogs/clear_gdsxray.dat')
    x_gdn = ascii.read('/Users/rsimons/Desktop/clear/Catalogs/clear_gdnxray.dat')

    gf_gds = ascii.read('/Users/rsimons/Desktop/clear/Catalogs/allfields/goodss/goodss_3dhst.v4.1_f160w.galfit')
    gf_gdn = ascii.read('/Users/rsimons/Desktop/clear/Catalogs/allfields/goodsn/goodsn_3dhst.v4.1_f160w.galfit')


re_arr = []
ms_arr = []

fit_types = array(['', '_S', '_EC', '_S_EC'])


if True:
    for ft, fit_type in enumerate(fit_types):
        with PdfPages('/Users/rsimons/Desktop/clear/figures/izi_z_radius%s.pdf'%fit_type) as pdf:
            ms = 5
            fig, ax = plt.subplots(1,1, figsize = (9, 4))
            if False:
                ax.errorbar(wang_cat[:,5], wang_cat[:,3], yerr =wang_cat[:,4], fmt = 'o', ms = ms,color = 'blue', label = 'Wang+ 17', zorder = 1)
                ax.errorbar(wang18_cat[:,3], wang18_cat[:,1], yerr =wang18_cat[:,2], fmt = 'o',  ms = ms,color = 'darkblue', label = 'Wang+ 18', zorder = 1)
                ax.errorbar(jones_cat[:,0], jones_cat[:,1], yerr =jones_cat[:,2], fmt = 'o',  ms = ms,color = 'skyblue', label = 'Jones+ 13', zorder = 1)
                ax.errorbar(swinbank_cat[:,0], swinbank_cat[:,1], yerr =swinbank_cat[:,2],  ms = ms,fmt = 'o', color = 'darkgreen', label = 'Swinbank+ 12', zorder = 1)
                ax.errorbar(log10(ma_cat[:,1]), ma_cat[:,7], yerr =ma_cat[:,8], fmt = 's',  ms = ms,markeredgecolor = 'black', markerfacecolor = "None", color = 'black', label = 'Ma+ 17; simulations', zorder = 1)
                wuyts_x = linspace(10.0, 11.5, 100)
                wuyts_y = -0.017*(wuyts_x - 10) + 0.0
                ax.errorbar(-99, -1, yerr = 1., fmt = mrker, color = 'red',  markeredgecolor = 'black', ms = 6, alpha = 1.0, label = 'Simons+ in prep')

                ax.plot(wuyts_x, wuyts_y, '--', linewidth = 3, color = 'midnightblue', label = 'Wuyts+ 16', zorder = 1)

            ax.axhline(y = 0, linestyle = '-', color = 'grey', alpha = 0.4, linewidth = 2, zorder = 0)

            to_pl = []
            for fl in fls:

                fld = fl.split('/')[-1].split('_')[0]
                di = fl.split('/')[-1].split('_')[1]
                if 'N'   in fld: 
                    ct = gn_cat
                    xcat = x_gdn
                    gcat = gf_gdn
                elif 'S' in fld: 
                    ct = gs_cat
                    xcat = x_gds
                    gcat = gf_gds
                match_cat = where(int(di) == ct[1].data['id'])[0][0]
                match_xcat = where(int(di) == xcat['ID'])[0]
                match_gcat = where(int(di) == gcat['NUMBER'])[0]

                if len(match_xcat) > 0: match_xcat = match_xcat[0]
                xclass = xcat['Xclass'][match_xcat]
                ra_c = ct[1].data['ra'][match_cat]
                dec_c = ct[1].data['dec'][match_cat]


                z = ct[1].data['z500'][match_cat]
                re = gcat['re'][match_cat]

                re_kpc = re/cosmo.arcsec_per_kpc_proper(z).value 
                mrker = 'o'
                mstar = np.log10(ct[1].data['mass'][match_cat])


                ft = np.load(fl, allow_pickle = True)[()]
                if xclass != 'AGN':
                    kpc_per_pix = 0.1 / cosmo.arcsec_per_kpc_proper(z).value 
                    try: 
                        crit1 = ft['p%s'%fit_type][0]/kpc_per_pix + 2*max(float(sqrt(ft['V%s'%fit_type][0,0])/kpc_per_pix), 0.03) < 0.0
                        crit2 = ft['p%s'%fit_type][0]/kpc_per_pix - 2*max(float(sqrt(ft['V%s'%fit_type][0,0])/kpc_per_pix), 0.03) > 0.0
                        if crit1 | crit2:
                            pass
                        else:
                            count_flat+=1 

                        alp = 1.0
                        if False:
                            if crit1:
                                count_falling+=1
                                alp = 1.0   
                        if False:             
                            if crit2:
                                count_rising+=1
                                alp = 1.0                    

                        ax.errorbar(mstar, float(ft['p%s'%fit_type][0]/kpc_per_pix), yerr = float(sqrt(ft['V%s'%fit_type][0,0])/kpc_per_pix), fmt = mrker, color = 'red',  markeredgecolor = 'black', ms = 5, alpha = alp)
                        #ax.errorbar(mstar, float(ft['p_S_EC'][0]/kpc_per_pix*re_kpc), yerr = 0., fmt = mrker, color = 'red',  markeredgecolor = 'black', ms = 5, alpha = alp)
                        re_arr.append(ft['p%s'%fit_type][0]/kpc_per_pix*re_kpc)
                        ms_arr.append(mstar)
                        count_all+=1


                    except: pass

            #ax.errorbar(-99, -1, yerr = 0.01, fmt = 's', color = 'red', fillstyle = 'none', markeredgecolor = 'red', label = 'CLEAR, (N = 112)', zorder = 10)
            #ax.errorbar(9.25, 0.0246, xerr = 0.25, yerr = 0.003, fmt = 'o', color = 'red', markeredgecolor = 'black', ms = 10, label = 'CLEAR, STACK', zorder = 10)
            #ax.errorbar(9.75, 0.0163, xerr = 0.25, yerr = 0.003, fmt = 'o', color = 'red', markeredgecolor = 'black', ms = 10, zorder = 10)
            #ax.errorbar(10.25, 0.0121, xerr = 0.25, yerr = 0.004,   fmt = 'o', color = 'red', markeredgecolor = 'black', ms = 10,  zorder = 10)
            #ax.errorbar(10.75, 0.0055, xerr = 0.25, yerr = 0.008,  fmt = 'o', color = 'red', markeredgecolor = 'black', ms = 10, zorder = 10)








            ax.annotate(r'$0.7 < z < 2.0$', (0.55, 0.85), xycoords = 'axes fraction', fontsize = 25, fontweight = 'bold')
            ax.set_xlabel(r'$\log$ M$_{*}$ (M$_{\odot}$)', fontsize = 20)
            ax.set_ylabel(r'$\frac{\Delta \log(O/H)}{\Delta R}$ (dex kpc$^{-1}$)', rotation = 90, fontsize = 20)
            ax.legend(bbox_to_anchor=(1.0, 1.05), frameon = False, fontsize = 18)

            ax.set_ylim(-0.33, 0.5)

            ax.set_xlim(8.2, 11.5)

            fig.subplots_adjust(bottom = 0.20, left = 0.15, right = 0.65, top = 0.95)
            pdf.savefig()




res = np.load('test.npy', allow_pickle = True)[()]
re_arr = res['re_arr']
ms_arr = res['mstar']





if False:
    with PdfPages('/Users/rsimons/Desktop/clear/figures/izi_z_effradius.pdf') as pdf:
        ms = 5

        fig, ax = plt.subplots(1,1, figsize = (9, 4))

        ax.axhline(y = 0, linestyle = '-', color = 'grey', alpha = 0.4, linewidth = 2, zorder = 0)


        bins = np.arange(9, 10.50, 0.25)

        for b, bn in enumerate(bins):
            gd = where((ms_arr > bn) & (ms_arr < bn+0.25) & (abs(re_arr) < 0.5))[0]
            re_bin = re_arr[gd]

            med = np.median(re_bin)
            d = np.abs(re_bin - med)
            mdev = np.median(d)

            mean_mstar = bn + 0.25/2.

            mrker = 'D'
            alp = 1.0

            med = np.mean(re_bin)
            #mdev = np.std(re_bin)/np.sqrt(len(re_bin))

            print (med)

            ax.errorbar(mean_mstar, med, yerr = mdev, fmt = mrker, color = 'red',  markeredgecolor = 'black', ms = 10., alpha = alp)

        #eyeballed from figure 11
        belfiore_re  = np.array([0.020, -0.04, -0.045, -0.08, -0.13, -0.14])
        belfiore_ere = np.array([0.022, 0.01, 0.05, 0.05, 0.02, 0.02])
        belfiore_mass = np.arange(9, 10.5, 0.25)

        ax.errorbar(belfiore_mass, belfiore_re, yerr = belfiore_ere, fmt = mrker, color = 'black',  markeredgecolor = 'black', ms = 10., alpha = alp)

        #ax.annotate(r'$0.7 < z < 2.0$', (0.55, 0.85), xycoords = 'axes fraction', fontsize = 25, fontweight = 'bold')
        ax.set_xlabel(r'$\log$ M$_{*}$ (M$_{\odot}$)', fontsize = 20)
        ax.set_ylabel(r'$\frac{\Delta \log(O/H)}{\Delta R}$ (dex R$_{\text{eff}}$$^{-1}$)', rotation = 90, fontsize = 20)
        ax.legend(bbox_to_anchor=(1.0, 1.05), frameon = False, fontsize = 18)

        ax.set_ylim(-0.2, 0.15)
        ax.set_xlim(8.7, 10.8)

        ax.set_xticks(arange(9, 11, 0.5))
        ax.set_yticks(arange(-0.2, 0.2, 0.1))

        fig.subplots_adjust(bottom = 0.20, left = 0.15, right = 0.65, top = 0.95)
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
                print (mstar)
                sfr = float(c[4])
                ax.plot(mstar, sfr, color = 'red', marker = 'o', zorder = 10, markeredgecolor = 'black', markersize = 10)


            #ax.plot(-99, -99, color = 'red', marker = 'o', zorder = 10, markeredgecolor = 'black', label = 'CLEAR G102\n+ archival G141\n(2 of 10 pointings)')
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





