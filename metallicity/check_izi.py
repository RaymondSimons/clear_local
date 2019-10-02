import numpy as np
from astropy.io import fits, ascii
from numpy import *
import matplotlib.pyplot as plt




izi_cat = ascii.read('/Users/rsimons/Dropbox/clear/catalogs/good_izi.cat', header_start = 0)
izi_cat = izi_cat[:50]
fit_types = array(['', '_S', '_EC', '_S_EC'])

fit_types = array(['_S'])

for fit_type in fit_types:
    if False:
        to_save = {}
        to_save['HB'] = []
        to_save['eHB'] = []
        to_save['O2'] = []
        to_save['eO2'] = []
        to_save['O3'] = []
        to_save['eO3'] = []
        to_save['Z'] = []
        to_save['euZ'] = []
        to_save['elZ'] = []

        to_save['npeaks'] = []


        to_save['HB']     = array(to_save['HB']    )
        to_save['eHB']    = array(to_save['eHB']   )
        to_save['O2']     = array(to_save['O2']    )
        to_save['eO2']    = array(to_save['eO2']   )
        to_save['O3']     = array(to_save['O3']    )
        to_save['eO3']    = array(to_save['eO3']   )
        to_save['Z']      = array(to_save['Z']     )
        to_save['euZ']    = array(to_save['euZ']   )
        to_save['elZ']    = array(to_save['elZ']   )
        to_save['npeaks'] = array(to_save['npeaks'])








        for f, (fld, di) in enumerate(zip(izi_cat['field'], izi_cat['id'])):
            indir = '/Users/rsimons/Dropbox/clear/products/metals/metal_maps'
            fl = indir + '/{}_{}_metals.fits'.format(fld, di)
            mm_fits = fits.open(fl)

            names = []

            for hdu in mm_fits:
                names.append(hdu.name)

            if ('HB%s'%fit_type in names) & ('OII%s'%fit_type in names) & ('OIII%s'%fit_type in names):
                flux_Hb =   mm_fits['%s%s'%('HB', fit_type)].data.ravel()
                flux_eHb =  mm_fits['%s%s'%('EHB', fit_type)].data.ravel()
                flux_Oii =   mm_fits['%s%s'%('OII', fit_type)].data.ravel()
                flux_eOii =  mm_fits['%s%s'%('EOII', fit_type)].data.ravel()
                flux_Oiii =   mm_fits['%s%s'%('OIII', fit_type)].data.ravel()
                flux_eOiii =  mm_fits['%s%s'%('EOIII', fit_type)].data.ravel()
                Z = mm_fits['Z%s'%fit_type].data
                n_peaks = Z[:,:,3].ravel()
                elZ = abs(Z[:,:,1] - Z[:,:,0]).ravel()
                euZ = abs(Z[:,:,2] - Z[:,:,0]).ravel()
                Z = Z[:,:,0].ravel()
                

                gd = where(~isnan(Z))

                to_save['HB']     = concatenate(( to_save['HB']    , flux_Hb[gd]))
                to_save['eHB']    = concatenate(( to_save['eHB']   , flux_eHb[gd]))
                to_save['O2']     = concatenate(( to_save['O2']    , flux_Oii[gd]))
                to_save['eO2']    = concatenate(( to_save['eO2']   , flux_eOii[gd]))
                to_save['O3']     = concatenate(( to_save['O3']    , flux_Oiii[gd]))
                to_save['eO3']    = concatenate(( to_save['eO3']   , flux_eOiii[gd]))
                to_save['Z']      = concatenate(( to_save['Z']     , Z[gd]))
                to_save['euZ']    = concatenate(( to_save['euZ']   , euZ[gd]))
                to_save['elZ']    = concatenate(( to_save['elZ']   , elZ[gd]))
                to_save['npeaks'] = concatenate(( to_save['npeaks'], n_peaks[gd]))



        np.save('/Users/rsimons/Dropbox/clear/catalogs/temp_izicheck.npy', to_save)

    if True:
        to_save = np.load('/Users/rsimons/Dropbox/clear/catalogs/temp_izicheck.npy', allow_pickle = True)[()]


        for sn_cut in [0, 1, 3, 5]:
            plt.close('all')

            fig1, ax1 = plt.subplots(1,1,figsize = (6,5.))
            fig2, ax2 = plt.subplots(1,1, figsize = (6,5.))
            O3 = to_save['O3']
            O2 = to_save['O2']
            HB = to_save['HB']


            eO3 = to_save['eO3']
            eO2 = to_save['eO2']
            eHB = to_save['eHB']


            npeaks = to_save['npeaks']
            O32 = O3/O2
            R23 = (O3 + O2)/HB
            Z = to_save['Z']
            gd = where((O3/eO3 > sn_cut) & (O2/eO2 > sn_cut) & (HB/eHB > sn_cut))[0]


            cmin = 6.5
            cmax = 9.5
            Z_for_col = array([min(max(ZZ - cmin, 0.), cmax - cmin) for ZZ in Z])/(cmax - cmin)
            Z_col = plt.cm.viridis(Z_for_col)

            sc1 = ax1.scatter(R23[gd], O32[gd], color = Z_col[gd], marker = 'o', cmap = plt.cm.viridis, s = 0.3, zorder = 10)
            cbar = fig1.colorbar(sc1, ax = ax1)
            cbar.set_label(r'12 + $\log$ (O/H)')
            cbar_ticks = np.arange(cmin, cmax+0.5, 0.5)
            cbar_ticks_str = ['%.1f'%c for c in cbar_ticks]
            cbar.set_ticks((cbar_ticks - cmin)/(cmax - cmin))
            cbar.set_ticklabels(cbar_ticks_str)


            cmin = 1
            cmax = 2
            P_for_col = array([min(max(ZZ - cmin, 0.), cmax - cmin) for ZZ in npeaks])/(cmax - cmin)
            P_col = plt.cm.viridis(P_for_col)

            ax2_sc = ax2.scatter(R23[gd], O32[gd], color = P_col[gd], marker = 'o', cmap = plt.cm.viridis, s = 0.3, zorder = 10)
            cbar2 = fig2.colorbar(ax2_sc, ax = ax2)
            cbar2.set_label(r'Number of peaks in posterior')
            cbar_ticks = np.linspace(cmin, cmax,2)
            cbar_ticks_str = ['%i'%c for c in cbar_ticks]
            cbar2.set_ticks((cbar_ticks - cmin)/(cmax - cmin))
            cbar2.set_ticklabels(cbar_ticks_str)


            for ax in [ax1, ax2]:
                ax.set_xlabel(r'R23 = (OIII + OII)/H$\beta$')
                ax.set_ylabel(r'O32 = OIII/OII')        
                ax.set_xscale('log')
                ax.set_yscale('log')
                ax.set_xlim(1.e-2, 1.e4)
                ax.set_ylim(5.e-4, 3.e3)
            


            for fig in [fig1, fig2]:
                fig.tight_layout()

            fig1.savefig('/Users/rsimons/Dropbox/clear/figures/izicheck_Z_%i.png'%sn_cut, dpi = 300)
            fig2.savefig('/Users/rsimons/Dropbox/clear/figures/izicheck_npeaks_%i.png'%sn_cut, dpi = 300)






        #ax.set_xlim(-10, 10)
        #ax.set_ylim(-10, 10)
        '''
        np.random.seed(1)

        x = np.random.normal(0, 1, 10000)
        y = np.random.normal(0, 1, 10000)
        z = np.random.normal(0, 1, 10000)

        gd = where((x>0) & (y>0) & (z>0))
        x = x[gd]
        y = y[gd]
        z = z[gd]

        sim_o32 = x/z
        sim_r23 = (x+y)/z

        noise_o32 = np.mean(sim_o32[sim_o32>0])
        noise_r23 = np.mean(sim_r23[sim_r23>0])


        ax.axvline(x = noise_r23, color = 'grey', zorder = 0)
        ax.axhline(y = noise_o32, color = 'grey', zorder = 0)
        '''





























