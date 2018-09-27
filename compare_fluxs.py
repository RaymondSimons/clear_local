import astropy
from astropy.io import fits
plt.ioff()
plt.close('all')


fig, axes = plt.subplots(1,4, figsize = (20,5))

for a, ax in enumerate(axes):
	ax.set_xscale('log')
	ax.set_yscale('log')
	ax.set_xlim(0.1, 400)
	ax.set_ylim(0.1, 400)
	ax.plot([0, 10000], [0, 10000], color = 'black', linewidth = 1, linestyle = '--')
	ax.set_xlabel('Grizli pipeline-derived line flux (1e-17 erg s$^{-1}$ cm$^{-2}$)')
	if a == 0: ax.set_ylabel('3DHST pipeline-derived line flux (1e-17 erg s$^{-1}$ cm$^{-2}$)')

cat1 = np.loadtxt('../Catalogs/lines_grizli.cat', dtype = 'str')
cat2 = fits.open('../Catalogs/GN_CLEAR.linefit.concat.v1.0.0.fits')


lines = ['OII', 'OIII', 'Ha', 'Hb']

for c, cat in enumerate(cat1):
	good = where(cat2[1].data['phot_id'] == int(cat[1]))[0]
	if len(good) > 0:
		good = good[0]
		for l, line in enumerate(lines):
			ax = axes[l]
			ax.annotate('%s'%line, (0.05, 0.9), xycoords = 'axes fraction', fontsize = 20, fontweight = 'bold')

			flux_1 = float(cat[3+l*2])
			eflux_1 = float(cat[4+l*2])
			flux_2 = cat2[1].data['%s_FLUX'%line][good]
			eflux_2 = cat2[1].data['%s_FLUX_ERR'%line][good]
			clr = 'black'
			mkr = 'o'
			if flux_1 == -99:
				flux_1 = 0.2
				eflux_2 = 0.0
				eflux_1 = 0.0
				clr = 'grey'
				mkr = '<'
			if flux_2 == -99:
				flux_2 = 0.2
				eflux_2 = 0.0
				eflux_1 = 0.0

				clr = 'grey'
				mkr = 'v'


			ax.errorbar(flux_1, flux_2, xerr = eflux_1, yerr = eflux_2, color = clr, fmt = mkr)





fig.savefig('../figures/lines_comparison.png', dpi = 300)

plt.close('all')