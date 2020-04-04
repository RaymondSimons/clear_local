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
from clear_local.utils.tools import *
plt.rcParams['xtick.labelsize']=12
plt.rcParams['ytick.labelsize']=12
plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']


cat = load_paper_catalog()


add_lit_points = True
add_wuyts_line = False
lit_cats = ['jones13', 'wang17', 'wang18', 'wang19', 'swinbank12', 'curti20']
lit_cat_dir = '/Users/rsimons/Dropbox/clear/catalogs/literature_gradients'

stack_dic = np.load('/Users/rsimons/Dropbox/clear/catalogs/stack_profiles.npy', allow_pickle = True)[()]

count_all = 0
count_flat = 0
count_falling = 0
count_rising = 0

fit_types = array(['', '_S', '_EC', '_S_EC'])

fit_types = array([''])
mrker = 'o'

if True:
    for ff, fit_type in enumerate(fit_types):
        with PdfPages('/Users/rsimons/Dropbox/clear/figures/for_paper/izi_z_radius%s_highZbranch.pdf'%fit_type) as pdf:
            ms = 5
            fig, (ax1, ax2) = plt.subplots(1,2, figsize = (10, 4), sharey = True)
            axes = [ax1, ax2]
            if add_lit_points:
                for lit_cat in lit_cats:
                    lcat = ascii.read(lit_cat_dir + '/' + lit_cat + '.cat')

                    mass  = lcat['lmass']
                    zgrad = lcat['zgrad']
                    uezgrad = lcat['uezgrad']
                    lezgrad = lcat['lezgrad']


                    ax1.errorbar(mass, zgrad,  \
                                yerr = [lezgrad, uezgrad], ms = 3, \
                                fmt = 'o', color = 'grey', zorder = 1)
                if add_wuyts_line:
                    wuyts_x = linspace(10.0, 11.5, 100)
                    wuyts_y = -0.017*(wuyts_x - 10) + 0.0
                    ax1.plot(wuyts_x, wuyts_y, '--', linewidth = 3, color = 'midnightblue', label = 'Wuyts+ 16', zorder = 1)

            to_pl = []


            mstars = []
            ms = []
            ems = []
            for f, ft in enumerate(cat):
                xclass = ft['xclass']


                z       = ft['z_map']
                re      = ft['re_125']
                mstar   = ft['mass_eazy']
                re_kpc  = re/cosmo.arcsec_per_kpc_proper(z).value 

                if xclass != 'AGN':
                    kpc_per_pix = 0.1 / cosmo.arcsec_per_kpc_proper(z).value 
                    m = ft['m%s_S_EC'%fit_type]/kpc_per_pix
                    m_err = ft['m%s_S_EC_err'%fit_type]/kpc_per_pix
                    ax1.errorbar(mstar, m, yerr = m_err, fmt = mrker, \
                                color = 'red',  markeredgecolor = 'black', \
                                ms = 5, zorder = 10)
                    mstars.append(mstar)
                    ms.append(m)
                    ems.append(m_err)
            if True:
                for i in stack_dic:
                    mass_mid = stack_dic[i]['mass_mid']
                    p = stack_dic[i]['r_kpc']['p']
                    V = stack_dic[i]['r_kpc']['V']
                    ax2.errorbar(mass_mid, p[0], yerr = np.sqrt(V[0,0]), fmt = 's', \
                                color = 'red',  markeredgecolor = 'black', \
                                ms = 10, zorder = 1)

            else:
                mstars  = np.array(mstars)
                ms      = np.array(ms    )
                ems     = np.array(ems   )

                dm = 0.25
                mass_bin_arr = np.arange(8.5, 11.0 + dm, dm)
                mass_bins = [(mass_bin_arr[i], mass_bin_arr[i+1]) for i in arange(len(mass_bin_arr) - 1)]

                for i in arange(len(mass_bins)):

                    gd_i = where((mstars > mass_bins[i][0]) & (mstars < mass_bins[i][1]))[0]
                    mass_mid = np.mean(mass_bins[i])
                    print (len(gd_i))
                    ms_i  = ms[gd_i]
                    ems_i = ems[gd_i]
                    bn = 250
                    from astropy.stats import bootstrap

                    dz_bs = bootstrap(ms_i, bootnum = bn)
                    median_Z_array = zeros(bn)      
                    for b in arange(bn):
                        median_Z_array[b] = np.nanmean(dz_bs[b] + np.random.normal(0., ems_i, len(dz_bs[b])))
                    emean_dZ_dr = np.std(median_Z_array)
                    mean_dZ_dr =  np.mean(median_Z_array)
                    print (emean_dZ_dr, mean_dZ_dr)

                    ax2.errorbar(mass_mid, mean_dZ_dr, yerr = emean_dZ_dr, fmt = 's', \
                    color = 'red',  markeredgecolor = 'black', \
                    ms = 10, zorder = 1)




            for a, ax in enumerate(axes):
                ax.annotate(r'$0.6 < z < 2.6$', (0.98, 0.96), ha = 'right', va = 'top',  xycoords = 'axes fraction', fontsize = 20, fontweight = 'bold')
                ax.set_xlabel(r'$\log$ M$_{*}$ (M$_{\odot}$)', fontsize = 20)
                ax.set_ylim(-0.33, 0.5)
                ax.set_xlim(8.4, 11.1)
                if a == 0:
                    ax.axhline(y = 0, xmin = 0.0, xmax = 1.0, linestyle = '-', color = 'black', alpha = 1.0, linewidth = 2, zorder = 0)
                if a == 1:
                    ax.annotate("no\ngradient", (8.45, 0.05), xycoords = 'data', va = 'top', ha = 'left', color = 'black', fontsize = 20)
                    ax.axhline(y = 0, xmin = 0.15, xmax = 1.0, linestyle = '-', color = 'black', alpha = 1.0, linewidth = 2, zorder = 0)
            ax1.set_ylabel(r'$\frac{\Delta \log(O/H)}{\Delta R}$ (dex kpc$^{-1}$)', rotation = 90, fontsize = 20)
            #ax1.annotate('individual galaxies', (0.05, 0.96), xycoords = 'axes fraction',
            #             ha = 'left', va = 'top', fontsize = 20, color = 'black')
            ax2.annotate('population stack'   , (0.05, 0.96), xycoords = 'axes fraction',
                         ha = 'left', va = 'top', fontsize = 20,  color = 'black')

            ax1.annotate("CLEAR", (0.98, 0.15), ha = 'right', va = 'center', xycoords = 'axes fraction', color = 'red', fontsize = 20, fontweight = 'bold')
            ax2.annotate("CLEAR", (0.98, 0.15), ha = 'right', va = 'center', xycoords = 'axes fraction', color = 'red', fontsize = 20, fontweight = 'bold')
            ax1.annotate("literature", (0.98, 0.05), ha = 'right', va = 'center', xycoords = 'axes fraction', color = 'grey', fontsize = 20, fontweight = 'bold')


            fig.subplots_adjust(wspace = 0.05, bottom = 0.20, left = 0.1, right = 0.98, top = 0.95)
            pdf.savefig()





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





