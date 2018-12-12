import numpy as np
import astropy
from astropy.io import fits
from matplotlib.colors import LogNorm
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['axes.linewidth'] = 1

plt.ioff()
plt.close('all')

cat_dir = '/Volumes/gdrive/clear/catalogs'
fig_dir = '/Volumes/gdrive/clear/figures'

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

if True:
    exp_102 = array([])
    exp_141 = array([])

    fields = ['GS1','GS2', 'GS3', 'GS5', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7']
    
    #flux_f = zeros((2, len(lines), len(fields)))
    flux_Ha = []
    flux_Hb = []
    flux_O2 = []
    flux_O3 = []

    eflux_Ha = []
    eflux_Hb = []
    eflux_O2 = []
    eflux_O3 = []

    tdhstN_lf = fits.open(cat_dir + '/GN_CLEAR.linefit.concat.v1.0.0.fits')
    tdhstS_lf = fits.open(cat_dir + '/GN_CLEAR.linefit.concat.v1.0.0.fits')
    tdhstN_zf = fits.open(cat_dir + '/GN_CLEAR.zfit.concat.v1.0.0.fits')
    tdhstS_zf = fits.open(cat_dir + '/GS_CLEAR.zfit.concat.v1.0.0.fits')

    gds_fast = fits.open(cat_dir + '/goodss_3dhst.v4.1.cats/Fast/goodss_3dhst.v4.1.fout.FITS')
    gdn_fast = fits.open(cat_dir + '/goodsn_3dhst.v4.1.cats/Fast/goodsn_3dhst.v4.1.fout.FITS')

    gds_eazy = fits.open(cat_dir + '/goodss_3dhst.v4.1.cats/Fast/goodss_3dhst.v4.1.fout.FITS')
    gdn_eazy = fits.open(cat_dir + '/goodsn_3dhst.v4.1.cats/Fast/goodsn_3dhst.v4.1.fout.FITS')

    gds_cat = fits.open(cat_dir + '/goodss_3dhst.v4.1.cats/Catalog/goodss_3dhst.v4.1.cat.FITS')
    gdn_cat = fits.open(cat_dir + '/goodsn_3dhst.v4.1.cats/Catalog/goodsn_3dhst.v4.1.cat.FITS')










    z = {}
    z['ID'] = []
    z['jh_mag'] = []
    z['z_max_grism'] = []
    z['z_peak_grism'] = []
    z['l95'] = []
    z['l68'] = []
    z['u68'] = []
    z['u95'] = []

    z['z_peak_phot'] = []
    z['z_phot_l95'] = []
    z['z_phot_l68'] = []
    z['z_phot_u68'] = []
    z['z_phot_u95'] = []

    z['nlines_g'] = []
    z['z_50'] = []
    z['z_02'] = []
    z['z_16'] = []
    z['z_84'] = []
    z['z_97'] = []


    z['Ha_FLUX_t'] = []
    z['Hb_FLUX_t'] = []
    z['OII_FLUX_t'] = []
    z['OIII_FLUX_t'] = []

    z['Ha_FLUX_ERR_t'] = []
    z['Hb_FLUX_ERR_t'] = []
    z['OII_FLUX_ERR_t'] = []
    z['OIII_FLUX_ERR_t'] = []


    z['Ha_FLUX_g'] = []
    z['Hb_FLUX_g'] = []
    z['OII_FLUX_g'] = []
    z['OIII_FLUX_g'] = []

    z['Ha_FLUX_ERR_g'] = []
    z['Hb_FLUX_ERR_g'] = []
    z['OII_FLUX_ERR_g'] = []
    z['OIII_FLUX_ERR_g'] = []

    fl = {}
    fl['nlines'] = []
    fl['Ha_FLUX'] = []
    fl['Hb_FLUX'] = []
    fl['O2_FLUX'] = []
    fl['O3_FLUX'] = []
    fl['Ha_FLUX_ERR'] = []
    fl['Hb_FLUX_ERR'] = []
    fl['O2_FLUX_ERR'] = []
    fl['O3_FLUX_ERR'] = []


    fl['lmass'] = []
    fl['lsfr'] = []
    fl['m_f160w'] = []
    fl['m_f140w'] = []
    fl['m_f125w'] = []









    for f, field in enumerate(fields):
        if 'N' in field: 
            lfcat = tdhstN_lf
            zfcat = tdhstN_zf
            fast_cat = gdn_fast
            eazy_cat = gdn_eazy
            mag_cat  = gdn_cat 



        if 'S' in field: 
            lfcat = tdhstS_lf
            zfcat = tdhstS_zf
            fast_cat = gds_fast
            eazy_cat = gds_eazy
            mag_cat = gds_cat 





        g_cat = fits.open(cat_dir + '/%s_lines_grizli.fits'%field)


        for i in arange(len(g_cat[1].data['ID'])):

            gd_z = where(g_cat[1].data['ID'][i] == zfcat[1].data['phot_id'])[0]
            gd_l = where(g_cat[1].data['ID'][i] == lfcat[1].data['phot_id'])[0]

            gd_f = where(g_cat[1].data['ID'][i] == fast_cat[1].data['id'])[0]
            gd_e = where(g_cat[1].data['ID'][i] == eazy_cat[1].data['id'])[0]
            gd_c = where(g_cat[1].data['ID'][i] == mag_cat[1].data['id'])[0]

            if ((len(gd_f) > 0) & (len(gd_e) > 0) & (len(gd_c) > 0)):
                if len(gd_f) > 1: gd_f = gd_f[0]
                if len(gd_e) > 1: gd_e = gd_e[0]
                if len(gd_c) > 1: gd_c = gd_c[0]
                gd_f = int(gd_f)
                gd_e = int(gd_e)
                gd_c = int(gd_c)

                fl['lmass'].append(fast_cat[1].data['lmass'][gd_f])
                fl['lsfr'].append(fast_cat[1].data['lsfr'][gd_f])
                fl['m_f160w'].append(25.0 - 2.5*log10(mag_cat[1].data['f_f160w'][gd_c]))
                fl['m_f140w'].append(25.0 - 2.5*log10(mag_cat[1].data['f_f140w'][gd_c]))
                fl['m_f125w'].append(25.0 - 2.5*log10(mag_cat[1].data['f_f125w'][gd_c]))

                fl['nlines'].append(g_cat[1].data['nlines'][i])

                fl['Ha_FLUX'].append(g_cat[1].data['Ha_FLUX'][i])
                fl['Hb_FLUX'].append(g_cat[1].data['Hb_FLUX'][i])
                fl['O2_FLUX'].append(g_cat[1].data['OII_FLUX'][i])
                fl['O3_FLUX'].append(g_cat[1].data['OIII_Flux'][i])
                fl['Ha_FLUX_ERR'].append(g_cat[1].data['Ha_FLUX_ERR'][i])
                fl['Hb_FLUX_ERR'].append(g_cat[1].data['Hb_FLUX_ERR'][i])
                fl['O2_FLUX_ERR'].append(g_cat[1].data['OII_FLUX_ERR'][i])
                fl['O3_FLUX_ERR'].append(g_cat[1].data['OIII_FLUX_ERR'][i])





            if (len(gd_z) > 0) & (len(gd_l) > 0):
                if len(gd_z) > 1: gd_z = gd_z[0]
                if len(gd_l) > 1: gd_l = gd_l[0]

                gd_z = int(gd_z)
                gd_l = int(gd_l)
                z['ID'].append(g_cat[1].data['ID'][i])
                z['jh_mag'].append(lfcat[1].data['jh_mag'][gd_l])       
                z['z_max_grism'].append(zfcat[1].data['z_max_grism'][gd_z])
                z['z_peak_grism'].append(zfcat[1].data['z_peak_grism'][gd_z])
                z['l95'].append(zfcat[1].data['l95'][gd_z])
                z['l68'].append(zfcat[1].data['l68'][gd_z])
                z['u68'].append(zfcat[1].data['u68'][gd_z])
                z['u95'].append(zfcat[1].data['u95'][gd_z])

                z['z_peak_phot'].append(zfcat[1].data['z_peak_phot'][gd_z])
                z['z_phot_l95'].append(zfcat[1].data['z_phot_l95'][gd_z])
                z['z_phot_l68'].append(zfcat[1].data['z_phot_l68'][gd_z])
                z['z_phot_u68'].append(zfcat[1].data['z_phot_u68'][gd_z])
                z['z_phot_u95'].append(zfcat[1].data['z_phot_u95'][gd_z])


                z['nlines_g'].append(g_cat[1].data['nlines'][i])
                z['z_50'].append(g_cat[1].data['z_50'][i])
                z['z_02'].append(g_cat[1].data['z_02'][i])
                z['z_16'].append(g_cat[1].data['z_16'][i])
                z['z_84'].append(g_cat[1].data['z_84'][i])
                z['z_97'].append(g_cat[1].data['z_97'][i])





                z['Ha_FLUX_t'].append(lfcat[1].data['Ha_FLUX'][gd_l])
                z['Hb_FLUX_t'].append(lfcat[1].data['Hb_FLUX'][gd_l])
                z['OII_FLUX_t'].append(lfcat[1].data['OII_FLUX'][gd_l])
                z['OIII_FLUX_t'].append(lfcat[1].data['OIII_FLUX'][gd_l])
                z['Ha_FLUX_ERR_t'].append(lfcat[1].data['Ha_FLUX_ERR'][gd_l])
                z['Hb_FLUX_ERR_t'].append(lfcat[1].data['Hb_FLUX_ERR'][gd_l])
                z['OII_FLUX_ERR_t'].append(lfcat[1].data['OII_FLUX_ERR'][gd_l])
                z['OIII_FLUX_ERR_t'].append(lfcat[1].data['OIII_FLUX_ERR'][gd_l])


                z['Ha_FLUX_g'].append(g_cat[1].data['Ha_FLUX'][i])
                z['Hb_FLUX_g'].append(g_cat[1].data['Hb_FLUX'][i])
                z['OII_FLUX_g'].append(g_cat[1].data['OII_FLUX'][i])
                z['OIII_FLUX_g'].append(g_cat[1].data['OIII_FLUX'][i])
                z['Ha_FLUX_ERR_g'].append(g_cat[1].data['Ha_FLUX_ERR'][i])
                z['Hb_FLUX_ERR_g'].append(g_cat[1].data['Hb_FLUX_ERR'][i])
                z['OII_FLUX_ERR_g'].append(g_cat[1].data['OII_FLUX_ERR'][i])
                z['OIII_FLUX_ERR_g'].append(g_cat[1].data['OIII_FLUX_ERR'][i])










        exp_102 = concatenate((exp_102, g_cat[1].data['T_G102']/60./60.))
        exp_141 = concatenate((exp_141, g_cat[1].data['T_G141']/60./60.))


if False:
    fig, axes = plt.subplots(1,2, figsize = (7,3.5))
    axes[0].hist(exp_102, color = 'darkblue', linewidth = 2, bins = linspace(0, 10, 20))
    axes[1].hist(exp_141, color = 'darkred', linewidth = 2, bins = linspace(0, 10, 20))

    fs = 10
    for ax in axes:
        ax.set_yscale('log')
        ax.set_xlabel('Total Expsoure Time (hr)', fontsize = fs)
        ax.set_ylim(5, 6000)

    axes[0].set_ylabel('Number of Objects', fontsize = fs)
    axes[0].set_title('G102', fontweight = 'bold', fontsize = fs*2.0, color = 'darkblue')
    axes[1].set_title('G141', fontweight = 'bold', fontsize = fs*2.0, color = 'darkred')

    fig.subplots_adjust(bottom = 0.15)
    fig.savefig(fig_dir + '/exposure_times.png', dpi = 400)



#Redshift Comparison with TDHST fits
if False:
    matplotlib.rcParams['axes.linewidth'] = 2

    fig, axes = plt.subplots(2,1, figsize = (4.6, 8))
    fig2, ax1 = plt.subplots(1,1, figsize = (5, 5))



    ax1.hist(array(z['z_peak_phot']) - array(z['z_50']), color = 'red', bins = linspace(-0.3, 0.3, 100), histtype = 'step', linewidth = 1)
    ax1.hist(array(z['z_peak_grism']) - array(z['z_50']), color = 'blue', bins = linspace(-0.3, 0.3, 100), histtype = 'step', linewidth = 1)

    ax1.set_xlabel('z$_{X}$ - z$_{grizli}$', fontsize = 15)
    ax1.set_ylabel('Number of objects', fontsize = 15)


    fs = 25
    ax1.annotate('$z_{phot}$', (0.04, 0.85), color = 'red', fontweight = 'bold', xycoords = 'axes fraction', fontsize = fs)
    ax1.annotate('$z_{grism, 3DHST}$',(0.04, 0.76),  color = 'blue',  fontweight = 'bold',xycoords = 'axes fraction', fontsize = fs)



    axes[0].plot(z['z_peak_phot'], z['z_50'], '.', markersize = 1, color = 'red', alpha = 0.3)

    axes[1].plot(z['z_peak_grism'], z['z_50'], '.', markersize = 1, color = 'blue', alpha = 0.3)
    #axes[1, 1].hist2d(array(z['z_peak_grism']), array(z['z_50']), bins = [linspace(0, 3.5, 50), linspace(0, 3.5, 50)], vmin = 0, vmax = 5)

    for ax in axes:
        ax.plot([0, 5], [0,5], 'k-', alpha = 0.2, zorder = 1)
    



    for ax in axes.ravel():
        ax.set_xlim(0,3.5)
        ax.set_ylim(0,3.5)


    axes[0].set_xlabel('', fontsize = 25)
    axes[1].set_xlabel('$z_{grizli}$', fontsize = 25)

    axes[0].set_ylabel('$z_{phot}$', fontsize = 25)
    axes[1].set_ylabel('$z_{3DHST, pipeline}$', fontsize = 25)




    fig.subplots_adjust(left = 0.15,bottom = 0.10, top = 0.95)
    fig.savefig(fig_dir + '/z_comparisons.png', dpi = 400)
    fig2.savefig(fig_dir + '/z_comparisons_hist.png', dpi = 400)


    fig, axes = plt.subplots(1,2, figsize = (10, 5))

    ax1 = axes[0]
    ax2 = axes[1]

    ez_grizli = abs(array(z['z_50']) - array(z['z_16']))
    ez_td = abs(array(z['z_peak_grism']) - array(z['l68']))
    ez_phot = abs(array(z['z_peak_phot'])  -  array(z['z_phot_l68']))
    delz_phot = array(z['z_peak_phot']) - array(z['z_50'])
    edelz_phot = sqrt(ez_grizli**2. + ez_phot**2.)

    delz_td = array(z['z_peak_grism']) - array(z['z_50'])
    edelz_td = sqrt(ez_grizli**2. + ez_td**2.)


    ax2.hist(delz_td/edelz_td, color = 'blue', bins = linspace(-5, 5., 50), histtype = 'step', linewidth = 1)
    ax1.hist(delz_phot/edelz_phot, color = 'red', bins = linspace(-5, 5., 50), histtype = 'step', linewidth = 1)

    gdt = where((abs(delz_td/edelz_td) < 5))[0]
    gdp = where((abs(delz_phot/edelz_phot) < 5))[0]

    stdt = std(delz_td[gdt]/edelz_td[gdt])
    stdp = std(delz_phot[gdp]/edelz_phot[gdp])

    mnt = mean(delz_td[gdt]/edelz_td[gdt])
    mnp = mean(delz_phot[gdp]/edelz_phot[gdp])


    ax2.axvline(mnt, color = 'blue', linestyle = '-')
    ax2.axvline(mnt+stdt, color = 'blue', linestyle = '--')
    ax2.axvline(mnt+-stdt, color = 'blue', linestyle = '--')


    ax1.axvline(mnp, color = 'red', linestyle = '-')
    ax1.axvline(mnp+stdp, color = 'red', linestyle = '-.')
    ax1.axvline(mnp+-stdp, color = 'red', linestyle = '-.')

    for ax in axes:
        ax.set_xticks(arange(-4,4))
        ax.set_xlim(-4.5,4.5)
        ax.set_xlabel(r'$\Delta z/\sigma_{\Delta z}$', fontsize = 20)
    #ax1.set_xlabel(r'($F_{Grizli} - F_{3DHST})/(\sigma_{Grizli}^{2} + \sigma_{F3DHST}^{2})^{1/2}$')
    ax1.set_ylabel('Number of objects', fontsize = 15)

    #ax1.annotate('$z_{phot}$', (0.04, 0.85), color = 'red', fontweight = 'bold', xycoords = 'axes fraction', fontsize = fs)
    #ax1.annotate('$z_{grism, 3DHST}$',(0.04, 0.76),  color = 'blue',  fontweight = 'bold',xycoords = 'axes fraction', fontsize = fs)

    fig.subplots_adjust(left = 0.15, right = 0.95)
    fig.savefig(fig_dir + '/redshift_diff_hist.png', dpi = 400)
    print 'Saved ' + fig_dir + '/%s_redshift_diff_comparison.png'%line



#Flux S/N versus line strength
if False:
    matplotlib.rcParams['axes.linewidth'] = 2

    for line in ['Ha', 'Hb', 'OII', 'OIII']:
        fig, ax = plt.subplots(1,1, figsize = (5,5))
        f1 = 1.e-17*array(z['%s_FLUX_t'%line])
        e1 = 1.e-17*array(z['%s_FLUX_ERR_t'%line])



        ax.plot(f1, f1/e1, 'k.', markersize = 4, zorder = 2)

        line_str = line
        if line == 'Ha':
            ax.annotate(r'H$\alpha$', (0.75, 0.10), xycoords = 'axes fraction', fontweight = 'bold', color = 'black', fontsize = 40)
        elif line == 'Hb':
            ax.annotate(r'H$\beta$', (0.75, 0.10), xycoords = 'axes fraction', fontweight = 'bold', color = 'black', fontsize = 40)
        elif line == 'OII':
            ax.annotate('[OII]', (0.75, 0.10), xycoords = 'axes fraction', fontweight = 'bold', color = 'black', fontsize = 40)
        elif line == 'OIII':
            ax.annotate('[OIII]', (0.70, 0.10), xycoords = 'axes fraction', fontweight = 'bold', color = 'black', fontsize = 40)

        ax.set_xlabel(r'$F_{Grizli}$ (erg s$^{-1}$ cm$^{-2}$)', fontsize = 15)
        ax.set_ylabel('line flux signal-to-noise', fontsize = 15, rotation = 90)

        ax.set_xscale('log')
        ax.set_yscale('log')

        ax.set_xlim(5.e-19, 5.e-15)
        ax.set_ylim(0.1, 50)

        ax.set_yticks([0.1, 1, 3, 5, 10, 25])
        ax.set_yticklabels(['0.1', '1', '3', '5', '10', '25'])

        ax.axhline(y = 5, alpha = 0.2, color = 'grey')
        ax.axvline(x = 5e-17, alpha = 0.2, color = 'grey')


        fig.subplots_adjust(left = 0.15, right = 0.95)
        fig.savefig(fig_dir + '/%s_flux_SN.png'%line, dpi = 400)
        print 'Saved ' + fig_dir + '/%s_flux_SN.png'%line







#Flux Comparison with TDHST fits
if False:
    matplotlib.rcParams['axes.linewidth'] = 2

    f1s = array([])
    f2s = array([])
    e1s = array([])
    e2s = array([])



    for line in ['Ha', 'Hb', 'OII', 'OIII']:
        fig, ax = plt.subplots(1,1, figsize = (5,5))
        f1 = 1.e-17*array(z['%s_FLUX_t'%line])
        f2 = 1.e-17*array(z['%s_FLUX_g'%line])
        e1 = 1.e-17*array(z['%s_FLUX_ERR_t'%line])
        e2 = 1.e-17*array(z['%s_FLUX_ERR_g'%line])

        f1s = concatenate((f1s,f1))
        f2s = concatenate((f2s,f2))
        e1s = concatenate((e1s,e1))
        e2s = concatenate((e2s,e2))




        ax.plot(f1, f2, 'k.', markersize = 4, zorder = 2)
        ax.errorbar(f1, f2, xerr = e1, yerr = e2,  fmt = 'o', color = 'grey', markersize = 0., linewidth = 0.3, alpha = 0.6, zorder = 2)

        line_str = line
        if line == 'Ha':
            ax.annotate(r'H$\alpha$', (0.75, 0.10), xycoords = 'axes fraction', fontweight = 'bold', color = 'black', fontsize = 40)
        elif line == 'Hb':
            ax.annotate(r'H$\beta$', (0.75, 0.10), xycoords = 'axes fraction', fontweight = 'bold', color = 'black', fontsize = 40)
        elif line == 'OII':
            ax.annotate('[OII]', (0.75, 0.10), xycoords = 'axes fraction', fontweight = 'bold', color = 'black', fontsize = 40)
        elif line == 'OIII':
            ax.annotate('[OIII]', (0.70, 0.10), xycoords = 'axes fraction', fontweight = 'bold', color = 'black', fontsize = 40)



        ax.plot([0, 1], [0,1], 'k-', alpha = 0.2, zorder = 1)

        ax.set_xlabel(r'$F_{Grizli}$ (erg s$^{-1}$ cm$^{-2}$)', fontsize = 15)
        ax.set_ylabel(r'$F_{3DHST pipeline}$ (erg s$^{-1}$ cm$^{-2}$)', fontsize = 15)

        ax.set_yscale('log')
        ax.set_xscale('log')

        ax.set_xlim(5.e-19, 5.e-15)
        ax.set_ylim(5.e-19, 5.e-15)



        fig.subplots_adjust(left = 0.15, right = 0.95)
        fig.savefig(fig_dir + '/%s_flux_comparison.png'%line, dpi = 400)
        print 'Saved ' + fig_dir + '/%s_flux_comparison.png'%line


    fig, ax1 = plt.subplots(1,1, figsize = (5, 5))


    gd = where((f1s - f2s != 0) & (f1s > 0 ) & (f2s > 0 ))[0]

    f1s = f1s[gd]
    f2s = f2s[gd]
    e1s = e1s[gd]
    e2s = e2s[gd]

    ax1.hist((f1s - f2s)/sqrt(e1s**2. + e2s**2.), color = 'black', bins = linspace(-5, 5., 50), histtype = 'step', linewidth = 1)


    gd = where((abs((f1s - f2s)/sqrt(e1s**2. + e2s**2.)) < 5))[0]
    mn = mean((f1s[gd] - f2s[gd])/sqrt(e1s[gd]**2. + e2s[gd]**2.))
    st = std((f1s[gd] - f2s[gd])/sqrt(e1s[gd]**2. + e2s[gd]**2.))
    ax1.axvline(mn, color = 'grey', linestyle = '-')
    ax1.axvline(mn+st, color = 'grey', linestyle = '--')
    ax1.axvline(mn-st, color = 'grey', linestyle = '--')

    ax1.set_xticks(arange(-4,4))
    ax1.set_xlim(-4.5,4.5)
    ax1.set_xlabel(r'$\Delta F/\sigma_{\Delta F}$', fontsize = 20)
    #ax1.set_xlabel(r'($F_{Grizli} - F_{3DHST})/(\sigma_{Grizli}^{2} + \sigma_{F3DHST}^{2})^{1/2}$')
    ax1.set_ylabel('Number of Objects', fontsize = 20)


    fig.subplots_adjust(left = 0.15, right = 0.95)
    fig.savefig(fig_dir + '/flux_diff_hist.png', dpi = 400)
    print 'Saved ' + fig_dir + '/%s_flux_comparison.png'%line


#Mass hist
if False:
    matplotlib.rcParams['axes.linewidth'] = 2

    fig, ax = plt.subplots(1,1, figsize = (5,5))
    gd = where(~isnan(array(fl['lmass'])))[0]
    ax.hist(array(fl['lmass'])[gd], bins = linspace(6, 12, 100), color = 'black')
    ax.set_xlabel(r'$\log$ M$_{*}$/M$_{\odot}$', fontsize = 20)
    ax.set_ylabel('Number of objects', fontsize = 20)
    fig.subplots_adjust(left = 0.15, right = 0.95)
    fig.savefig(fig_dir + '/mass_histogram.png', dpi = 300)
    print 'Saved ' + fig_dir + '/mass_histogram.png'




    fig, ax = plt.subplots(1,1, figsize = (5,5))
    gd = where((~isnan(array(fl['m_f125w']))) & (isfinite(array(fl['m_f125w']))))[0]
    ax.hist(array(fl['m_f125w'])[gd], bins = linspace(16, 28, 100), color = 'black')
    ax.set_xlabel(r'm$_{125, AB}$', fontsize = 20)
    ax.set_ylabel('Number of objects', fontsize = 20)
    fig.subplots_adjust(left = 0.15, right = 0.95)
    fig.savefig(fig_dir + '/mf125w_histogram.png', dpi = 300)
    print 'Saved ' + fig_dir + '/mf125w_histogram.png'







#Linemap Figure
if False:
    matplotlib.rcParams['axes.linewidth'] = 2

    n = 48918
    full = fits.open('/Volumes/gdrive/clear/grizli_extractions/GS1/j033300m2742/Prep/GS1_%.5i.full.fits'%n)

    fig, axes = plt.subplots(2,3)

    vmn = 0.2
    vmx = 2.5
    axes[0,0].imshow(full['DSCI', 'F105W'].data - full['DSCI', 'F105W'].data.min() + 0.01, cmap = 'Greys', interpolation = 'nearest', norm = LogNorm(vmin = 0.01, vmax = 1))
    axes[0,1].imshow(full['Line', 'OII'].data - full['Line', 'OII'].data.min() + 0.001,  cmap = 'viridis', interpolation = 'nearest', norm = LogNorm(vmin = vmn, vmax = vmx))
    axes[0,2].imshow(full['Line', 'OIII'].data - full['Line', 'OIII'].data.min() + 0.001, cmap = 'viridis', interpolation = 'nearest', norm = LogNorm(vmin = vmn, vmax = vmx))
    axes[1,0].imshow(full['Line', 'Ha'].data - full['Line', 'Ha'].data.min() + 0.001,   cmap = 'viridis', interpolation = 'nearest', norm = LogNorm(vmin = vmn, vmax = vmx))
    axes[1,1].imshow(full['Line', 'Hb'].data - full['Line', 'Hb'].data.min() + 0.001,   cmap = 'viridis', interpolation = 'nearest', norm = LogNorm(vmin = vmn, vmax = vmx))
    axes[1,2].imshow(full['Line', 'Hg'].data - full['Line', 'Hg'].data.min() + 0.001,   cmap = 'viridis', interpolation = 'nearest', norm = LogNorm(vmin = vmn, vmax = vmx))



    for ax in axes.ravel():
        ax.set_xlim(25,55)
        ax.set_ylim(25,55)
        ax.axis('off')
    fig.subplots_adjust(left = 0.15, right = 0.95)


    figname = fig_dir + '/GS1_%.5i_fit.png'%n
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0, left = 0.0, right = 1.0, top = 1.0, bottom = 0.0)
    fig.savefig(figname, dpi = 300)
    print 'Saved ' + figname


#Make direct image
if False:
    fig, axes = subplots(1,2, figsize = (10,7.6))
    gfl = fits.open('/Volumes/gdrive/clear/grizli_extractions/GS1/j033300m2742/Prep/gs1-cxt-08-249.0-g102_drz_sci.fits')
    dfl = fits.open('/Volumes/gdrive/clear/grizli_extractions/GS1/j033300m2742/Prep/gs1-cxt-08-249.0-f105w_drz_sci.fits')

    dd = dfl[0].data
    gd = gfl[0].data

    cmap_v = plt.get_cmap('viridis')
    #cmap_v.set_bad(color = 'k', alpha = 1.)

    cmap_b = plt.get_cmap('Greys_r')
    #cmap_b.set_bad(color = 'k', alpha = 1.0)



    axes[0].imshow(log10(dd), alpha = 1.0, cmap = cmap_b, vmin = -1.5, vmax = 1., interpolation = 'nearest')
    axes[1].imshow(log10(gd), alpha = 1.0, cmap = cmap_v, vmin = -2.0, vmax = 0.2, interpolation = 'nearest')
    




    for ax in axes:
        ax.axis('off')
        ax.set_xlim(450,1091)




    fig.subplots_adjust(left = 0., right = 1.0, top = 1.0, bottom = 0.0)





    figname = fig_dir + '/GS1_249_g102_f105w.png'
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0, left = 0.0, right = 1.0, top = 1.0, bottom = 0.0)
    fig.savefig(figname, dpi = 300)
    print 'Saved ' + figname



#Detector space figure
if False:
    fig, axes = subplots(2,2, figsize = (10,10))
    a = fits.open('/Volumes/gdrive/clear/grizli_extractions/GS1/j033300m2742/Prep/icxt12r2q.01.GrismFLT.fits')
    di = a[4].data
    dd = a[5].data
    gd = a[10].data

    '''
    mn = min((dd.min(), gd.min()))

    dd = dd - mn + 0.001
    gd = gd - mn + 0.001
    '''
    cmap_v = plt.get_cmap('viridis')
    cmap_b = plt.get_cmap('Greys')

    #axes[0].imshow(log10(dd), alpha = 1.0, cmap = cmap_v, vmin = -0.5, vmax = 1.0, interpolation = 'nearest')
    #axes[1].imshow(log10(gd), alpha = 1.0, cmap = cmap_v, vmin = -0.5, vmax = 1.0, interpolation = 'nearest')

    grv = dd[dd!=0].ravel()
    mn_g = nanmean(grv)
    st_g = nanstd(grv)
    vmn = mn_g
    vmx = mn_g + 0.05*st_g


    lrv = log10(di[di!=0]).ravel()
    mn = nanmean(lrv)
    st = nanstd(lrv)


    axes[0,0].imshow(log10(di), alpha = 1.0, cmap = cmap_b, vmin = mn-2*st, vmax = mn + 5*st, interpolation = 'nearest')
    axes[0,1].imshow(dd, alpha = 1.0, cmap = cmap_v, vmin = vmn, vmax = vmx, interpolation = 'nearest')
    axes[1,0].imshow(gd, alpha = 1.0, cmap = cmap_v, vmin = vmn, vmax = vmx, interpolation = 'nearest')
    axes[1,1].imshow(dd - gd, alpha = 1.0, cmap = cmap_v, vmin = vmn, vmax = vmx, interpolation = 'nearest')




    for ax in axes.ravel():
        ax.axis('off')
        ax.set_xlim(200,900)
        ax.set_ylim(200,900)




    fig.subplots_adjust(left = 0., right = 1.0, top = 1.0, bottom = 0.0)





    figname = fig_dir + '/GS1_detector.png'
    fig.subplots_adjust(hspace = 0.0, wspace = 0.0, left = 0.0, right = 1.0, top = 1.0, bottom = 0.0)
    fig.savefig(figname, dpi = 300)
    print 'Saved ' + figname
































