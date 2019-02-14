import glob
from glob import glob
import astropy
from astropy.io import fits
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
plt.close('all')
plt.ioff()


full_dir = '/Volumes/gdrive/clear/scalephot_testing'
full_dir = '/user/rsimons/grizli_extractions/GN2/j123652p6215/Prep'

fls = glob(full_dir + '/GN2_2_*full.fits')

ids = [fl.split('/')[-1].split('_')[-1].strip('full.fits') for fl in fls]



scale_orders = arange(-1,4)
fig, axes = plt.subplots(1,2, figsize = (14, 5))
fig2, axes2 = plt.subplots(2,2, figsize = (10, 10))
fig3, axes3 = plt.subplots(1,2, figsize = (14, 5))

for i, di in enumerate(ids):
    bics = []
    chi = []
    zs = []
    for s, scale_order in enumerate(scale_orders):
        try:
            fit_d = fits.open(full_dir + '/GN2_%i_%s.full.fits'%(scale_order, di))
            bic = fit_d[1].header['BIC_TEMP']
            bics.append(bic)
            chi.append(fit_d[1].header['CHIMIN']/fit_d[1].header['DOF'])
            zs.append([fit_d[1].header['Z50'], fit_d[1].header['Z16'], fit_d[1].header['Z84']])
        except:
            bics.append(nan)
            chi.append(nan)
            zs.append([nan, nan, nan])

    axes[0].plot(scale_orders, bics, '-', marker = 'o', alpha = 0.3, linewidth = 1.)
    axes[1].plot(scale_orders, chi, '-', marker = 'o', alpha = 0.3, linewidth = 1.)


    axes3[0].plot(scale_orders + np.random.normal(0, 0.1, len(scale_orders)), (array(bics) - bics[1])/array(bics[1]), 'k.', marker = 'o', markersize = 3, alpha = 0.3, linewidth = 1.)
    axes3[1].plot(scale_orders + np.random.normal(0, 0.1, len(scale_orders)), (array(chi )- chi[1])/array(chi[1]), 'k.', marker = 'o', markersize = 3,alpha = 0.3, linewidth = 1.)



    axes2[0, 0].plot(zs[1][0], zs[0][0],'k.', alpha = 0.3, linewidth = 1.)
    axes2[0, 1].plot(zs[1][0], zs[2][0],'k.', alpha = 0.3, linewidth = 1.)
    axes2[1, 0].plot(zs[1][0], zs[3][0],'k.', alpha = 0.3, linewidth = 1.)
    axes2[1, 1].plot(zs[1][0], zs[4][0],'k.', alpha = 0.3, linewidth = 1.)



axes2[0, 0].set_ylabel('z(-1)', fontsize = 20)
axes2[0, 1].set_ylabel('z(+1)', fontsize = 20)
axes2[1, 0].set_ylabel('z(+2)', fontsize = 20)
axes2[1, 1].set_ylabel('z(+3)', fontsize = 20)


for ax in axes2.ravel():
    ax.set_xlim(0,3)
    ax.set_ylim(0,3)
    ax.set_xlabel('z(0)', fontsize = 20)

for ax in axes:
    ax.set_xlim(-1.1,3.3)
    ax.set_xlabel('Spectrum/Photometry Scaling Order')
    ax.set_xticks([-1, 0, 1, 2, 3])
    ax.set_xticklabels(['-1', '0', '+1', '+2', '+3'])

for ax in axes3:
    ax.set_xlim(-1.1,3.3)
    ax.set_xlabel('Spectrum/Photometry Scaling Order')
    ax.set_xticks([-1, 0, 1, 2, 3])
    ax.set_xticklabels(['-1', '0', '+1', '+2', '+3'])



axes[0].annotate('-1 = no photometry, 0 = no scaling, +1 = constant scaling, +2 = first-order scaling, +3 = second-order scaling', (0.5, 0.9), ha = 'center', xycoords = 'figure fraction')    
axes2[0,0].annotate('-1 = no photometry, 0 = no scaling, +1 = constant scaling, +2 = first-order scaling, +3 = second-order scaling', (0.5, 0.9), ha = 'center', xycoords = 'figure fraction')    
axes[0].set_ylabel('Bayesian Information Criterion, Template Fit \n("BIC_TEMP")')
axes[1].set_ylabel(r'Minimium reduced $\chi^2$, Template Fit \n("CHIMIN"/"DOF")')


axes3[0].set_ylabel(r'$\Delta$BIC/BIC[0]')
axes3[1].set_ylabel(r'$\Delta\chi_{min}^2$/$\chi_{min}^2$[0]')


axes3[0].axhline(y = 0, color = 'grey', linestyle = '--')
axes3[1].axhline(y = 0, color = 'grey', linestyle = '--')


axes[1].set_ylim(0.5,100)
axes[1].set_yscale('log')

axes3[0].set_ylim(-1, 1)
axes3[1].set_ylim(-1, 1)

fig.savefig('/home/rsimons/git/clear_local/bics_scale_phot.png', dpi = 300)
fig2.savefig('/home/rsimons/git/clear_local/z_comparison.png', dpi = 300)
fig3.savefig('/home/rsimons/git/clear_local/scale_phot.png', dpi = 300)







