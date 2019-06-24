import astropy
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
plt.rcParams['text.usetex'] = True
plt.rcParams['axes.linewidth'] = 2
plt.ioff()
plt.close('all')

fig, axes = plt.subplots(2,2, figsize = (8,8))

for a, ax in enumerate(axes.ravel()):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(0.1e-17, 400e-17)
    ax.set_ylim(0.1e-17, 400e-17)
    ax.plot([0, 1], [0, 1], color = 'black', linewidth = 1, linestyle = '--')
    ax.set_xlabel('F$_{Grizli}$ (erg s$^{-1}$ cm$^{-2}$)', fontsize = 16)
    ax.set_ylabel('F$_{3DHST pipeline}$ (erg s$^{-1}$ cm$^{-2}$)', fontsize = 16)

flds = ['GN1',
        'GN2',
        'GN3',
        'GN4',
        'GN5',
        'GN7',
        'GS1',
        'GS2',
        'GS3',
        'GS4',
        'GS5']
#flds = ['GN1']

if False:
    f_over_sigf = []
    for f, fld in enumerate(flds):
        grizli_cat = fits.open('/Users/rsimons/Desktop/clear/Catalogs/grizli_v2.1_cats/%s_lines_grizli.fits'%fld)
        if 'S' in fld: tdhst_fld = 'GS'
        if 'N' in fld: tdhst_fld = 'GN'
        tdhst_cat = fits.open('/Users/rsimons/Desktop/clear/Catalogs/%s_CLEAR.fits'%tdhst_fld)


        lines = ['OII', 'OIII', 'Hb', 'Ha']
        line_str = ['[OII]', '[OIII]', r'H$\beta$', r'H$\alpha$']

        for c, di in enumerate(grizli_cat[1].data['ID']):
            good = where(tdhst_cat[1].data['phot_id'] == int(di))[0]
            if len(good) > 0:
                good = good[0]
                for l, line in enumerate(lines):
                    ax = axes.ravel()[l]
                    ax.annotate('%s'%line_str[l], (0.05, 0.85), xycoords = 'axes fraction', fontsize = 30, fontweight = 'bold')

                    flux_1 = grizli_cat[1].data[line+'_FLUX'][c]
                    eflux_1 = grizli_cat[1].data[line+'_FLUX_ERR'][c]
                    flux_2 = tdhst_cat[1].data['%s_FLUX'%line][good]
                    eflux_2 = tdhst_cat[1].data['%s_FLUX_ERR'%line][good]
                    clr = 'black'
                    mkr = 'o'
                    ms = 1.5
                    #if (flux_1 != -99) & (flux_2 != -99):
                    if (flux_1/eflux_1 > 2.) & (flux_2/eflux_2 > 2.):
                        ax.errorbar(1.e-17 * flux_1, 1.e-17 * flux_2, xerr = 1.e-17 * eflux_1, yerr = 1.e-17 * eflux_2, color = clr, fmt = mkr, markersize = ms, alpha = 0.3, lw = 0.)
                        f_over_sigf.append((flux_2 - flux_1)/(sqrt(eflux_1**2. + eflux_2**2.)))                    


f_over_sigf = array(f_over_sigf)

fig2, ax = plt.subplots(1,1, figsize = (6,6))
ax.set_xlabel(r'$\Delta F$/$\sigma_{\delta F}$')
ax.set_ylabel(r'Number of Objects')

ax.hist(f_over_sigf, bins = arange(-5, 5, 0.2))

fig2.tight_layout()
fig2.savefig('/Users/rsimons/Desktop/clear/figures/f_over_sigf.png', dpi = 300)


fig.tight_layout()
fig.savefig('/Users/rsimons/Desktop/clear/figures/lines_comparison.png', dpi = 300)

plt.close('all')










