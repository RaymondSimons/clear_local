import astropy
from astropy.io import fits
from astropy.cosmology import Planck15 as cosmo
import argparse
import glob
from glob import glob
import numpy as np
from numpy import *
import photutils
from photutils import detect_sources
import matplotlib
import matplotlib.pyplot as plt




plt.ioff()
plt.close('all')
def make_metal_profile(fl):
    a = np.load(fl, allow_pickle = True)[()]

    fit_types = array(['', '_S', '_EC', '_S_EC'])

    fig = plt.figure(figsize = (8, 8))
    nrows = 4
    ncols = 3
    res = {}
    for ft, fit_type in enumerate(fit_types):
        res['p{}'.format(fit_type)]  = np.nan 
        res['V{}'.format(fit_type)]  = np.nan
        res['r{}'.format(fit_type)]  = np.nan
        res['Z{}'.format(fit_type)]  = np.nan
        res['eZ{}'.format(fit_type)] = np.nan

        ax_Z =  plt.subplot2grid((nrows, ncols), (ft, 0))
        ax_p =  plt.subplot2grid((nrows, ncols), (ft, 1), colspan = 1)
        ax_pf =  plt.subplot2grid((nrows, ncols), (ft, 2), colspan = 1)

        cmap = plt.cm.viridis
        cmap.set_bad('k')

        vmin = 7.8
        vmax = 9.2

        Zmap = a['Z{}'.format(fit_type)]
        Npeaks = a['Npeaks{}'.format(fit_type)].ravel()

        ax_Z.imshow(Zmap, cmap = cmap, vmin = vmin,  vmax = vmax)

        Zmap_mask = Zmap.mask

        elZmap = a['elZ{}'.format(fit_type)]
        euZmap = a['euZ{}'.format(fit_type)]

        crit1 = (~Zmap_mask.ravel()) & (~isnan(Zmap.data.ravel()))
        crit2 = Npeaks == 1

        r = a['r{}'.format(fit_type)].ravel()[crit1]
        Z = Zmap.data.ravel()[crit1]
        elZ = elZmap.data.ravel()[crit1]
        euZ = euZmap.data.ravel()[crit1]


        ax_p.errorbar(r, Z, yerr = [elZ, euZ], linestyle = 'None', fmt = 'x', color = 'grey', alpha = 0.2)




        r = a['r{}'.format(fit_type)].ravel()[(crit1) & (crit2)]
        Z = Zmap.data.ravel()[(crit1) & (crit2)]
        elZ = elZmap.data.ravel()[(crit1) & (crit2)]
        euZ = euZmap.data.ravel()[(crit1) & (crit2)]




        elZ = array([max(eZ, 0.2) for eZ in elZ])
        seuZ = array([max(eZ, 0.2) for eZ in euZ])
        ax_p.errorbar(r, Z, yerr = [elZ, euZ], linestyle = 'None', fmt = 'x', color = 'blue', alpha = 0.3)

        ax_p.set_ylim(6.9, 9.6)




        ax_Z.plot(a['midx'], a['midy'], 'x', color = 'white')
        if len(r) == 0: 
            ax_p.axis('off')
            ax_pf.axis('off')

        if len(Z) > 5:
            outl = nanstd(Z)


            def reject_outliers(data_for_mdev, data, m = 2.):
                d = np.abs(data_for_mdev - np.median(data_for_mdev))
                mdev = np.median(d)
                d2 = np.abs(data - np.median(data_for_mdev))

                s = d2/mdev if mdev else 0.

                return s<m

            Z_outlier = Z[r < 4]
            gd_tofit = reject_outliers(Z_outlier, Z)
            r   = r  [gd_tofit]
            Z   = Z  [gd_tofit]
            elZ = elZ[gd_tofit]
            euZ = euZ[gd_tofit]
            if len(Z) > 5.:
                ax_p.errorbar(r, Z, yerr = [elZ, euZ], linestyle = 'None', fmt = 'o', color = 'blue', alpha = 1.0)
                ax_pf.errorbar(r, Z, yerr = [elZ, euZ], linestyle = 'None', fmt = 'o', color = 'blue', alpha = 1.0)

                ax_pf.set_xlim(ax_p.get_xlim())
                if (max(r) - min(r) > 2) & (len(where(r < 4.)[0]) > 3.):


                    eZ = np.mean((elZ, euZ), axis= 0 )
                    try:
                        p, V = np.polyfit(r, Z, deg = 1., w = 1./(eZ**2.), cov = True)
                        draws = np.random.multivariate_normal(p, V, size = 100)
                        x = linspace(0, max(r) + 1., 100)
                        res['p{}'.format(fit_type)]  = p
                        res['V{}'.format(fit_type)]  = V
                        res['r{}'.format(fit_type)]  = r
                        res['Z{}'.format(fit_type)]  = Z
                        res['eZ{}'.format(fit_type)] = eZ


                        for d in draws:
                            ax_pf.plot(x, x*d[0] + d[1], color = 'blue', alpha = 0.1)
                    except:
                        pass


                    ax_Z.axis('off')
                    ylm_min = max(min(Z-elZ) - 0.2, 7.0)
                    ylm_max = min(max(Z+euZ) + 0.2, 9.5)

                    ax_pf.set_ylim(ylm_min, ylm_max)
                    ax_p.axhline(y = ylm_min, xmin = 0.8, xmax = 1.0, linestyle = '--', color = 'grey')
                    ax_p.axhline(y = ylm_max, xmin = 0.8, xmax = 1.0, linestyle = '--', color = 'grey')





    figdir = '/Users/rsimons/Desktop/clear/figures/izi_metal_maps/metal_profiles'



    npsavename = fl.split('/')[-1].replace('.npy', '_profile_fit.npy')
    np.save('/Users/rsimons/Desktop/clear/izi_metal_profiles/fits/{}'.format(npsavename), res)
    figname = figdir + '/' + fl.split('/')[-1].replace('npy', 'png')


    print ('saving %s'%figname)
    fig.tight_layout()
    fig.savefig(figname)

    plt.close('all')





if __name__ == '__main__':
    fls = glob('/Users/rsimons/Desktop/clear/izi_metal_profiles/*npy')

    for f, fl in enumerate(fls):
        #if '26739' in fl:
        make_metal_profile(fl)






