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

    fig, axes = plt.subplots(1,3, figsize = (15, 5))




    z_05_10 = [-27.40, 5.02, -0.22,8.5, 11.2 , 'black']
    z_10_15 = [-26.03, 4.62, -0.19, 9.1, 11.3,  'black']
    z_15_20 = [-24.04, 4.17, -0.16, 9.2, 11.5,  'black']
    z_20_25 = [-19.99, 3.44, -0.13, 9.5, 11.5, 'black']


    coeffs_all = [z_05_10, z_10_15, z_20_25, z_15_20 ]

    for a, ax in enumerate(axes):
        zmin = z_bins[a][0]
        zmax = z_bins[a][1]        
        in_bin = where((z > zmin) & (z < zmax))
        ax.plot(lmass[in_bin], sfr[in_bin],'o', color = 'red', markersize = 7, markeredgecolor = 'black')

        coeffs = coeffs_all[a]

        whit_lw = 4

        ax.annotate(r'%.1f\,$<$\,z\,$<$\,%.1f'%(zmin, zmax), (0.98, 0.03), va = 'bottom', ha = 'right', xycoords = 'axes fraction', fontsize = 30)
        for fac in [1, 3.3]:
            x_plot = np.linspace(coeffs[3],coeffs[4], 100)
            y_plot = coeffs[0] + coeffs[1] * x_plot + coeffs[2]*x_plot**2.
            ax.plot(x_plot, 10**y_plot/fac, '-', color = coeffs[-1], linewidth = whit_lw/fac, zorder = 0)




            x_plot_exp_lower = np.linspace(8,coeffs[3], 100)
            y_plot = coeffs[0] +coeffs[1] * x_plot_exp_lower + coeffs[2]*x_plot_exp_lower**2.
            ax.plot(x_plot_exp_lower, 10**y_plot/fac, '-', color = coeffs[-1], linewidth = whit_lw/fac, zorder = 0)


            x_plot_exp_upper = np.linspace(coeffs[4],12, 100)
            y_plot = coeffs[0] +coeffs[1] * x_plot_exp_upper + coeffs[2]*x_plot_exp_upper**2.
            ax.plot(x_plot_exp_upper, 10**y_plot/fac, '-', color = coeffs[-1], linewidth = whit_lw/fac, zorder = 0)




            x_plot = np.linspace(coeffs[3],coeffs[4], 100)
            y_plot = coeffs[0] + coeffs[1] * x_plot + coeffs[2]*x_plot**2.
            ax.plot(x_plot, 10**y_plot*fac, '-', color = coeffs[-1], linewidth = whit_lw/fac, zorder = 0)


            x_plot_exp_lower = np.linspace(8,coeffs[3], 100)
            y_plot = coeffs[0] +coeffs[1] * x_plot_exp_lower + coeffs[2]*x_plot_exp_lower**2.
            ax.plot(x_plot_exp_lower, 10**y_plot*fac, '-', color = coeffs[-1], linewidth = whit_lw/fac, zorder = 0)


            x_plot_exp_upper = np.linspace(coeffs[4],12, 100)
            y_plot = coeffs[0] +coeffs[1] * x_plot_exp_upper + coeffs[2]*x_plot_exp_upper**2.
            ax.plot(x_plot_exp_upper, 10**y_plot*fac, '-', color = coeffs[-1], linewidth = whit_lw/fac, zorder = 0)



        ax.set_xlabel(r'$\log$ M$_*$ (M$_{\odot}$)')
        if a == 0: ax.set_ylabel(r'star-formation rate (M$_{\odot}$ yr$^{-1}$)')

        ax.set_yscale('log')
        ax.set_xlim(8, 12)
        ax.set_ylim(1.e-3, 1.e3)
    fig.tight_layout()
    fig.savefig('/Users/rsimons/Desktop/clear/figures/for_paper/m_sfr.png', dpi = 300)




    plt.close('all')
