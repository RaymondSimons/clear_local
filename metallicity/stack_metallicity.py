from astropy.io import fits, ascii
from astropy.cosmology import Planck15 as cosmo
from clear_local.utils.tools import *
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
from astropy.stats import bootstrap
plt.close('all')

fitdir = '/Users/rsimons/Dropbox/clear/products/metals/metal_maps_cleaned'
cat = load_paper_catalog()


dm = 0.5
mass_bin_arr = np.arange(8.5, 11.0 + dm, dm)
mass_bins = [(mass_bin_arr[i], mass_bin_arr[i+1]) for i in arange(len(mass_bin_arr) - 1)]



full_dic = {}
for i in arange(len(mass_bins)):
    full_dic[i]         = {}    
    full_dic[i]['mass_min'] = mass_bins[i][0]
    full_dic[i]['mass_max'] = mass_bins[i][1]
    full_dic[i]['mass_mid'] = np.mean(mass_bins[i])




r_types = ['r_kpc', 'r_eff']
r_types = ['r_kpc']#, 'r_eff']


for r_type in r_types:
    print (r_type)
    fig, axes = plt.subplots(len(mass_bins), 1, figsize = (4,3*len(mass_bins)))
    for i in arange(len(mass_bins)):
        print ('\t', mass_bins[i])
        full_dic[i][r_type] = {}    

        full_dic[i][r_type]['r'] = np.array([])
        full_dic[i][r_type]['z'] = np.array([])
        full_dic[i][r_type]['r_binned'] = []
        full_dic[i][r_type]['z_binned'] = []
        full_dic[i][r_type]['ez_binned'] = []



    for o, obj in enumerate(cat):
        z       = obj['z_map']
        mstar   = obj['mass_eazy']
        scale  = cosmo.arcsec_per_kpc_proper(z).value 

        for i in arange(len(mass_bins)):
            if (mstar > mass_bins[i][0]) & (mstar < mass_bins[i][1]):
                break
        mm = fits.open(fitdir + '/%s_%s_metals_highZbranch_cleaned.fits'%(obj['field'], obj['id']))

        if r_type == 'r_kpc':
            r_ravel = mm['R'].data.ravel() * 0.1/scale
            dr = 2.

        elif r_type == 'r_eff':
            r_ravel = mm['R'].data.ravel() * 0.1/obj['re_125']
            dr = 0.25

        z_ravel = mm['Z_clean'].data.ravel()

        full_dic[i][r_type]['r'] = np.concatenate((full_dic[i][r_type]['r'], r_ravel))
        full_dic[i][r_type]['z'] = np.concatenate((full_dic[i][r_type]['z'], z_ravel))


        axes[i].plot(r_ravel, z_ravel, 'k.', markersize = 2)

    for i, ax in enumerate(axes):

        rmin = 0.
        rmax = nanmax(full_dic[i][r_type]['r'])

        r_arr = np.arange(rmin, rmax, dr)

        for r in r_arr:
            rmn = r
            rmx = r+dr
            gd = where((full_dic[i][r_type]['r'] > rmn) & \
                       (full_dic[i][r_type]['r'] < rmx) & \
                       (~np.isnan(full_dic[i][r_type]['z'])) &\
                       (np.isfinite(full_dic[i][r_type]['z']))
                       )[0]
            if len(gd) > 20:
                zs = full_dic[i][r_type]['z'][gd]
                r_mid  = np.mean([rmn, rmx])
                z_mean =  np.nanmean(zs)
                dx = dr/2.
                if False:
                    dz =  np.nanstd(zs)/len(gd)**0.5
                else:
                    #bootstrap the uncertainties
                    bn = 250
                    dz_bs = bootstrap(zs, bootnum = bn)
                    median_Z_array = zeros(bn)      
                    for b in arange(bn):
                        median_Z_array[b] = np.nanmedian(dz_bs[b] + np.random.normal(0., 0.2, len(dz_bs[b])))
                    dz = np.std(median_Z_array)

                #dz = max(0.25, dz)
                ax.errorbar(r_mid, z_mean, xerr = dr/2., yerr = dz, fmt = 'o', color = 'red', markersize = 10)
                full_dic[i][r_type]['r_binned'].append(r_mid)
                full_dic[i][r_type]['z_binned'].append(z_mean)
                full_dic[i][r_type]['ez_binned'].append(dz)


        r  = np.array(full_dic[i][r_type]['r_binned'])
        Z  = np.array(full_dic[i][r_type]['z_binned'])
        eZ = np.array(full_dic[i][r_type]['ez_binned'])

        if len(Z) > 2:
            p, V = np.polyfit(r, Z, deg = 1., w = 1./(eZ**2.), cov = True)
            draws = np.random.multivariate_normal(p, V, size = 100)
            x = linspace(0, max(r) + 1., 100)
            for d in draws:
                ax.plot(x, x*d[0] + d[1], color = 'red', alpha = 0.05, zorder = 2)

        else:
            p = np.array([np.nan, np.nan])
            V = np.array([[np.nan, np.nan], [np.nan, np.nan]])




        ax.annotate('m = %.3f +- %.3f dex/kpc'%(p[0], np.sqrt(V[0,0])), (0.95, 0.90), \
                    xycoords = 'axes fraction', ha = 'right', va = 'top')

        ax.annotate('b = %.3f +- %.3f dex'%(p[1], np.sqrt(V[1,1])), (0.95, 0.85),\
                    xycoords = 'axes fraction', ha = 'right', va = 'top')




        full_dic[i][r_type]['p'] = p
        full_dic[i][r_type]['V'] = V

        ax.set_ylim(8.3, 9.8)
        ax.set_xlim(0., nanmax(r)*1.1)
        ax.set_ylabel('12 + log(O/H)')
        ax.annotate('%.1f < '%mass_bins[i][0]+r'M$_*$/M$_{\odot}$' + ' < %.1f'%(mass_bins[i][1]), (0.95, 0.95),
                    xycoords = 'axes fraction', ha = 'right', va = 'top')



    if r_type == 'r_kpc':
        axes[-1].set_xlabel('distance from center (kpc)')
    else:
        axes[-1].set_xlabel('distance from center (r/R$_{eff}$)')
    fig.tight_layout()
    fig.savefig('/Users/rsimons/Dropbox/clear/figures/stack_metals_%s.png'%r_type, dpi = 300)


np.save('/Users/rsimons/Dropbox/clear/catalogs/stack_profiles.npy', full_dic)





















