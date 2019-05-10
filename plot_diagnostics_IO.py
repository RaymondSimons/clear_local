import astropy
from astropy.io import fits
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['xtick.labelsize'] = 14

plt.ioff()

diagnostics = array([['O32'], 
                     ['R23'],
                     ['O32', 'R23']])

x = arange(7, 9.5, 0.025)
a = zeros((3, len(x), len(x)))



plt.close('all')
fig = plt.figure(figsize = (len(diagnostics)*5, 5))
for d, diagnostic in enumerate(diagnostics):
    for o, OH_true in enumerate(x):
        fit_str = ('%s'%diagnostic).strip('[').strip(']').strip("'").replace("', '", '_')
        out = fits.getdata('/Users/rsimons/Desktop/clear/diagnostics/out/%.3f_%s_samples.fits'%(OH_true, fit_str))[:,0]
        p, xedge = histogram(out, bins = concatenate((x, array([9.5]))))
        a[d, :, o] = p/(1.*nanmax(p))

    ax = fig.add_subplot(1,3,d+1)

    cm = mpl.cm.viridis
    cm.set_bad('k', 1)
    ax.imshow(a[d], cmap = cm, interpolation = 'nearest', origin = 'upper', vmin = 0.15, vmax = 0.98)

    ax.set_xlabel(r'12 + $\log$(O/H), true', fontsize = 20)
    if d == 0: ax.set_ylabel(r'12 + $\log$(O/H), out', fontsize = 20)
    ax.set_title(diagnostic, fontsize = 15)


    xticks = arange(7, 9.5, 0.5)
    gd_xticks = array([where(abs(x-xx) < 0.001)[0][0] for xx in xticks])
    st_xticks = array(['%.1f'%xx for xx in xticks])

    ax.set_xticks(gd_xticks)
    ax.set_xticklabels(st_xticks)

    ax.set_yticks(gd_xticks)
    ax.set_yticklabels(st_xticks)
    ax.plot([0, 150],[0, 150], color = 'white', linestyle = '-')
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)


fig.savefig('/Users/rsimons/Desktop/clear/diagnostics/figures/2D_hist.png', dpi = 300)