import astropy
from astropy.io import fits, ascii
from astropy.convolution import Gaussian2DKernel, convolve_fft, Box2DKernel
from astropy.cosmology import Planck15 as cosmo
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 


seed(3)



fig = figure(figsize = (8, 6))

ax = fig.add_subplot(111)


z_bin = [0.7, 1.3, 2.0]

new_cat_dir = '/Users/rsimons/Desktop/clear/Catalogs'


X,Y = np.meshgrid(np.arange(-40, 40, 1)*0.1, np.arange(-40, 40, 1)*0.1)
r = sqrt(X**2. + Y**2.)
r+= np.random.normal(0, 0.1, shape(r))

clrs = ['blue', 'red']
for mm, m_bin in enumerate(array([[9, 9.5,10], [10, 10.5, 11]])):
    rmx = 1.4
    bins_median = arange(0., rmx, 0.05)

    data_1 = fits.getdata(new_cat_dir + '/o32_%.1f_%.1f_%.1f_%.1f.fits'%(z_bin[0], z_bin[1], m_bin[0], m_bin[1]))
    data_2 = fits.getdata(new_cat_dir + '/o32_%.1f_%.1f_%.1f_%.1f.fits'%(z_bin[0], z_bin[1], m_bin[1], m_bin[2]))
    data_3 = fits.getdata(new_cat_dir + '/o32_%.1f_%.1f_%.1f_%.1f.fits'%(z_bin[1], z_bin[2], m_bin[0], m_bin[1]))
    data_4 = fits.getdata(new_cat_dir + '/o32_%.1f_%.1f_%.1f_%.1f.fits'%(z_bin[1], z_bin[2], m_bin[1], m_bin[2]))

    rs = array([])
    fs = array([])
    zs = array([])

    SN = 1.5
    ms = 2.
    alp = 0.15

    for data in [data_1, data_2, data_3, data_4]:
        for i in arange(len(data)):
            O3 = data[i,:,:,2]
            O2 = data[i,:,:,0]
            eO3 = data[i,:,:,3]
            eO2 = data[i,:,:,1]
            f = log10(O3/O2)
            gd = where((O3 > SN * eO3) & (O2 > SN * eO2))
            z = 8.54 - 0.59 * f[gd]
            #z, eZ = OH(O3[gd].ravel(), O2[gd].ravel(), eO3[gd].ravel(), eO2[gd].ravel())
            #re = r/cosmo.arcsec_per_kpc_proper(mean([z_bin[0], z_bin[1]])).value                        
            ax.plot(r[gd].ravel(), z, 'k.', 
                    markersize = ms, alpha = alp, color = clrs[mm])
            
            rs = concatenate([rs, r[gd].ravel()])
            fs = concatenate([fs, f[gd].ravel()])
            zs = concatenate([zs, z])




    to_fitr = []
    to_fitz = []
    to_fitez = []



    for rr, rfit in enumerate(bins_median[1:len(bins_median) - 1]):
        gi = where((rs > bins_median[rr]) & (rs < bins_median[rr+1]))[0]
        mid_x = mean([bins_median[rr], bins_median[rr+1]])

        stf = std(log10(fs[gi]))
        stf_u = 10**(log10(median(fs[gi])) + stf) - median(fs[gi])
        stf_l = median(fs[gi]) - 10**(log10(median(fs[gi])) + stf) 
        
        stz = std(zs[gi])/sqrt(len(zs[gi]))
        
        if stf > 1.:
            clr = 'blue'
            alp = 1.0
            
        else:
            clr = 'blue'
            alp = 1.0
            
        ax.errorbar(mid_x, median(zs[gi]), yerr = stz, marker = 'D', markeredgecolor = 'black', color = clrs[mm])
        
        to_fitr.append(mid_x)
        to_fitz.append(median(zs[gi]))
        to_fitez.append(stz)

    to_fitr = array(to_fitr)
    to_fitz = array(to_fitz)
    to_fitez = array(to_fitez)
    
    to_fitr = to_fitr[~isnan(to_fitz)]
    to_fitz = to_fitz[~isnan(to_fitz)]
    to_fitez = to_fitez[~isnan(to_fitez)]
    
    to_fitez[to_fitez == 0.0] = 0.1
    p, V = np.polyfit(to_fitr, to_fitz, deg = 1., w = 1./to_fitez, cov = True)
    draws = np.random.multivariate_normal(p, V, size = 100)
    x = np.linspace(0, rmx, 1000)
    for d in draws:
        ax.plot(x, x*d[0] + d[1], alpha = 0.1, color = clrs[mm])


ax.set_xlim(0, 1.5)
ax.set_ylim(7.8,9.2)

bbox_props = dict(boxstyle="square", fc="w", ec='w', alpha=0.4)
ax.annotate(r'0.7 $\boldsymbol{<}$ z $\boldsymbol{<}$ 2.0 (N = 130)', (0.982, 0.918), bbox = bbox_props, ha = 'right', xycoords = 'axes fraction', color = 'black', fontweight = 'bold', fontsize = 25)

ax.annotate(r'$\mathbf{9\,\,<\boldsymbol{\log}\, M_*/M_{\odot}< 10}$', (0.70, 0.15), ha = 'center', xycoords = 'axes fraction', color = 'blue', fontweight = 'bold', fontsize = 25)
ax.annotate(r'$\mathbf{10< \boldsymbol{\log}\, M_*/M_{\odot} < 11}$', (0.70, 0.05), ha = 'center', xycoords = 'axes fraction', color = 'red', fontweight = 'bold', fontsize = 25)


ax.set_ylabel(r'12 + log(O/H)', fontsize = 18)
ax.set_xlabel('distance from center (arcsec)', fontsize = 18)





fig.savefig('/Users/rsimons/Desktop/random/Z_radius_HST.png', dpi = 600)

















