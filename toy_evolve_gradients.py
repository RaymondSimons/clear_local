import astropy
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.io import ascii
import numpy as np
from numpy import *
from astropy import units as u
from mpl_toolkits.axes_grid.inset_locator import inset_axes

np.random.seed(1)

plt.close('all')
plt.ioff()
nelson = ascii.read('/Users/rsimons/Desktop/clear/catalogs/nelson16.cat')

def mzr(mass, calibration = 'O3N2'):
    """Taken from Table 2 of R Sanders+ 18"""
    if calibration == 'O3N2':
        m = 0.30
        b = 5.30
    return m*mass+b


fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2, figsize = (8, 8))

clrs = ['blue', 'green', 'orange', 'red']
re = [] #from van der wel 14
ls = ['-', '-.']


ax1.set_ylabel(r'$\log$ $\Sigma_{sfr}$ (M$_{\odot}$ yr$^{-1}$ kpc$^{-2}$)')


#taken from review talk, Kennicutt upcoming
n_med_low, n_sig_low, KS_norm_low = 1.39, 0.07, -3.82
n_med_high, n_sig_high, KS_norm_high = 1.06, 0.06, -2.42


y_med, y_sig = 0.03, 0.005
mass_bins = [[9.0, 9.5], [9.5, 10.0], [10.0, 10.5], [10.5, 11.0]]

n_draws = 100

KS_slope_low  = np.random.normal(n_med_low, n_sig_low, n_draws)
KS_slope_high = np.random.normal(n_med_high, n_sig_high, n_draws)







for c in arange(4):
    fs = 15
    da = 0.07
    xann = 0.95
    yann = 0.95 - da*c
    if c == 0: ax1.annotate(r'$\log$ M$_{*}$/M$_{\odot}$ = %.1f - %.1f'%(mass_bins[c][0], mass_bins[c][0]), \
                 xy = (xann, yann), xycoords = 'axes fraction', ha = 'right', va = 'top', color = clrs[c], fontsize = fs)
    else: ax1.annotate('%.1f - %.1f'%(mass_bins[c][0], mass_bins[c][0]), \
                 xy = (xann, yann), xycoords = 'axes fraction', ha = 'right', va = 'top', color = clrs[c], fontsize = fs)

ax1.annotate('Nelson+ 16\n$z\sim1$', xy = (0.05, 0.05), xycoords = 'axes fraction', ha = 'left', va = 'bottom', color = 'black', fontsize = 15)


ax2.set_ylim(-5, 3)
ax2.set_xlim(-1, 5)
ax2.set_xlabel(r'$\log$ $\Sigma_{gas}$ (M$_{\odot}$ pc$^{-2}$)')
ax2.set_ylabel(r'$\log$ $\Sigma_{sfr}$ (M$_{\odot}$ yr$^{-1}$ kpc$^{-2}$)')

ax3.set_ylabel('dZ/dt (Myr$^{-1}$)')
ax4.set_ylabel('dZ/dR (dex kpc$^{-1}$)')
ax4.set_xlabel('time (Myr)')


for ax in [ax1, ax3]: ax.set_xlabel('R (kpc)')







y = np.random.normal(y_med, y_sig, n_draws)

y[y<0] = 0.001

#dt = 1.e5
js = 10**np.arange(-1, 3, 0.05) * 1.e6


dzdrs = np.zeros((len(js), n_draws))*np.nan









for i in arange(4):
    cat = nelson[nelson['group'] == i]
    sfr = cat['sfr'] * u.Msun/u.yr/u.kpc**2
    r = cat['r'] * u.kpc
    ax1.plot(r, log10(sfr.value), 'o-', color = clrs[i])
    Z0 = mzr(mean(mass_bins[i])) - 12

    for nn, n in enumerate([n_med_low, n_med_high]):
        A = 10**(KS_norm_low) * u.pc**(2.*n)/u.yr/u.kpc**2./u.Msun**(n-1.)
        const = y_med * A**(1./n)
        dzdt = ((const * sfr**(1-1./n)).value*u.pc**2./u.kpc**2/u.yr).to('1/yr')
        if i == 3: lbl='$n$ = %.2f'%(n)
        else:  lbl = '_nolegend_'
        ax3.plot(cat['r'], dzdt * 1.e6, linestyle = ls[nn], color = clrs[i], label = lbl)



    for nn, (n_low, n_high, yy) in enumerate(zip(KS_slope_low, KS_slope_high, y)):
        if nn%10 == 0: print (nn, len(y))
        n = np.nan * np.zeros(len(sfr))
        KS_norm = np.nan * np.zeros(len(sfr))

        low_sequence  = np.where(log10(sfr.value) < -1)[0]
        high_sequence = np.where(log10(sfr.value) > -1)[0]

        n[low_sequence] = n_low
        n[high_sequence] = n_high
        KS_norm[low_sequence] = KS_norm_low
        KS_norm[high_sequence] = KS_norm_high

        dzdt = []
        for kk in arange(len(sfr)):


            A = 10**(KS_norm[kk]) * u.pc**(2.*n[kk])/u.yr/u.kpc**2./u.Msun**(n[kk]-1.)
            const = yy * A**(1./n[kk])
            dzdt.append(((const * sfr[kk]**(1-1./n[kk])).value*u.pc**2./u.kpc**2/u.yr).to('1/yr').value)

        dzdt = np.array(dzdt) /u.yr


        rinner = r[0]
        router = mean(r[3:8])
        dzdtinner = dzdt[0]
        dzdtouter = mean(dzdt[3:8])
        t_prev = 0.0
        zinner, zouter = Z0, Z0
        for jj, j in enumerate(js):

            zinner=np.log10(dzdtinner.value*(j - t_prev) + 10**zinner)
            zouter=np.log10(dzdtouter.value*(j - t_prev) + 10**zouter)
            dzdr = (zouter - zinner)/(router - rinner)
            dzdrs[jj, nn] = dzdr.value  
            t_prev = j

    p16, p84 = np.percentile(dzdrs, [16, 84], axis = 1) 
    ax4.fill_between(js * 1.e-6, y1 = p16, y2 = p84, color = clrs[i], alpha = 0.3)
    #ax3.plot(cat['r'], dzdt, linestyle = ls[nn], color = clrs[i])

#ax4.set_xscale('symlog', linxthresh = 10.)
ax4.set_xscale('log')
ax4.set_xlim(0.5, max(js) * 1.e-6)





lgas_low = np.linspace(-1, 2, 100)
lgas_high = np.linspace(2, 10, 100)

for nn, (n_low, n_high) in enumerate(zip(KS_slope_low, KS_slope_high)):
    lsfr = n_low * lgas_low + KS_norm_low
    ax2.plot(lgas_low, lsfr, color = 'grey', linestyle = '-', alpha = 0.2)

    lsfr = n_high * lgas_high + KS_norm_high
    ax2.plot(lgas_high, lsfr, color = 'black', linestyle = '-', alpha = 0.2)




ax2_mini = inset_axes(ax2,width=1.1, height=0.8,
                    bbox_to_anchor=(0.60, 0.18),
                    bbox_transform=ax2.transAxes, loc=3, borderpad=0)


ax2_mini.hist(KS_slope_low, color = 'grey')
ax2_mini.hist(KS_slope_high, color = 'black')
ax2_mini.set_xlabel('$n$')
ax2_mini.set_xlim(0.75, 1.75)
ax2_mini.set_xticks([0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7])
ax2_mini.set_xticklabels(['0.8', '', '1.0', '', '1.2', '', '1.4', '', '1.6', ''])









#ax2_mini.spines['left'].set_visible(False)
#ax2_mini.spines['top'].set_visible(False)
#ax2_mini.spines['right'].set_visible(False)

ax2_mini.set_yticks([])





ax4.axhline(y = 0.0, xmin = 0.0, xmax = 1.0, color = 'black')
ax4.annotate("flat", (0.95, 0.93), xycoords = 'axes fraction',\
            va = 'top', ha = 'right', fontsize = 15, color = 'black')

ax3.get_yaxis().get_major_formatter().set_scientific(True)
ax4.set_ylim(ax4.get_ylim()[0], 0.05)

ax2.annotate('Kennicutt-Schmidt relation\n$\Sigma_{SFR}\propto(\Sigma_{gas})^{n}$', xy = (0.05, 0.95), xycoords = 'axes fraction', ha = 'left', va = 'top', color = 'black', fontsize = 15)

ax2.annotate('low-SFR\nsequence',  xy = (0.3, -4.5),ha = 'left', va = 'bottom',  xycoords = 'data', color = 'grey')
ax2.annotate('high-SFR\nsequence', xy = (3, -0.5), ha = 'left', va = 'bottom', xycoords = 'data',  color = 'black')

ax3.legend(loc = 1)

fig.tight_layout()
fig.savefig('/Users/rsimons/Dropbox/clear/figures/for_paper/toy_models/ks.png', dpi = 400)


















