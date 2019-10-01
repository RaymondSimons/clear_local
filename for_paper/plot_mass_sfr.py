import importlib
from clear_local.utils import tools
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
plt.ioff()
plt.close('all')



if __name__ == '__main__':
    cat      = tools.load_paper_catalog()
    eazy_cat = tools.load_eazy_catalog()


    gd = np.where(abs(cat['z_50'] - cat['z500_eazy'])/cat['z_50'] < 0.3)

    lmass = cat['mass_eazy'][gd]
    
    #use whitaker's UV+IR SFR
    sfr   = cat['sfr_w'][gd]

    z = cat['z_map'][gd]

    #where not available, use Eazy v4.4 derived SFR
    sfr[[sfr < 0]] = cat['sfr_eazy'][gd][sfr < 0]

    lmass_all = np.log10(eazy_cat[1].data['mass'])
    sfr_all = eazy_cat[1].data['sfr']



    #0.6 - 1.0 is 2 Gyr
    #
    z_bins = [(0.6, 1.0), (1.0, 1.5), (1.5, 2.6)]

    fig, axes = plt.subplots(1,3, figsize = (9.5, 10/3.))




    z_05_10 = [-27.40, 5.02, -0.22,8.5, 11.2 ]
    z_10_15 = [-26.03, 4.62, -0.19, 9.1, 11.3]
    z_15_20 = [-24.04, 4.17, -0.16, 9.2, 11.5]
    z_20_25 = [-19.99, 3.44, -0.13, 9.5, 11.5]


    coeffs_all = [z_05_10, z_10_15, (z_15_20, z_20_25)]

    for a, ax in enumerate(axes):
        zmin = z_bins[a][0]
        zmax = z_bins[a][1]        
        in_bin = where((z > zmin) & (z < zmax))
        ax.plot(lmass[in_bin], sfr[in_bin],'o', color = 'red', markersize = 5, markeredgecolor = 'black')

        coeffs = coeffs_all[a]

        whit_lw = 1.5

        ls_exp = '--'
        ax.annotate(r'%.1f\,$<$\,z\,$<$\,%.1f'%(zmin, zmax), (0.98, 0.03), va = 'bottom', ha = 'right', xycoords = 'axes fraction', fontsize = 20)
        for fac in [1, 10**-0.5, 10**0.5]:
            if fac == 1: lw_fac = whit_lw
            else: lw_fac = whit_lw/3.

            if a < 2.:
                xmin_xmax_arr = [(coeffs[3],coeffs[4], '-'), \
                                 (8,coeffs[3], '--'), \
                                 (coeffs[4],12,'--')]
            else:

                xmin_xmax_arr = [(coeffs[0][3],coeffs[0][4], '-'), \
                                 (8,coeffs[0][3], '--'), \
                                 (coeffs[0][4],12,'--')]

            for xx, (xmin, xmax, ls) in enumerate(xmin_xmax_arr):
                x_plot = np.linspace(xmin, xmax, 100)
                if a < 2:
                    y_plot = coeffs[0] + coeffs[1] * x_plot + coeffs[2]*x_plot**2.

                else:
                    y_plot_1 = coeffs[0][0] + coeffs[0][1] * x_plot + coeffs[0][2]*x_plot**2.
                    y_plot_2 = coeffs[1][0] + coeffs[1][1] * x_plot + coeffs[1][2]*x_plot**2.
                    y_plot = np.mean([y_plot_1, y_plot_2], axis = 0)


                ax.plot(x_plot, 10**y_plot*fac, ls, color = 'black', linewidth = lw_fac, zorder = 0)


        ax.set_xlabel(r'$\log$ M$_*$ (M$_{\odot}$)')
        if a == 0: ax.set_ylabel(r'star-formation rate (M$_{\odot}$ yr$^{-1}$)')





        ax.set_yscale('log')
        ax.set_xlim(8, 12)
        ax.set_ylim(1.e-3, 1.e3)

        if a > 0: 
            ax.set_yticklabels([''])

    #axes[0].annotate('main\nsequence', (11.95, 12), ha = 'right', fontsize = 10)

    fig.tight_layout()
    fig.savefig('/Users/rsimons/Desktop/clear/figures/for_paper/m_sfr.png', dpi = 300)




    plt.close('all')
