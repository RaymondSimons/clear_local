import astropy
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import os
plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True



with PdfPages('/Users/rsimons/Desktop/clear/figures/exposure_times.pdf') as pdf:

    fields = ['GS1','GS2', 'GS3', 'GS4', 'GS5', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7', 'ERSPRIME']

    for field in fields:
        cat = fits.open('/Users/rsimons/Desktop/clear/Catalogs/grizli_v2.1_cats/%s_lines_grizli.fits'%field)
        t102 = cat[1].data['T_G102']
        t141 = cat[1].data['T_G141']


        fig, axes = plt.subplots(1,2, figsize = (9, 4.5))

        axes[0].hist(t102/3600.,  label = field, bins = linspace(0, 30, 40), alpha = 1.0, color = 'darkblue')
        axes[1].hist(t141/3600.,  label = field, bins = linspace(0, 30, 40), alpha = 1.0, color = 'darkred')

        #axes[0].annotate('%s'%field, (0.4, 1.0), xycoords = 'figure fraction', fontweight = 'bold', fontsize = 30)
        axes[0].set_title('%s'%field, fontweight = 'bold', fontsize = 30)

        axes[0].set_xlabel('G102 Exposure Time (hr)', fontsize = 12)
        axes[1].set_xlabel('G141 Exposure Time (hr)', fontsize = 12)

        '''
        gd = where(t102/3600. > 25)[0]
        for g in gd: 
            print cat[1].data['ID'][g]
            os.system('cp /Users/rsimons/all_png/GS4_%.5i.full.png /Users/rsimons/Desktop/clear/figures/g102_gt_25hrs/'%cat[1].data['ID'][g])
            os.system('cp /Users/rsimons/all_png/GS4_%.5i.stack.png /Users/rsimons/Desktop/clear/figures/g102_gt_25hrs/'%cat[1].data['ID'][g])
            os.system('cp /Users/rsimons/all_png/GS4_%.5i.line.png /Users/rsimons/Desktop/clear/figures/g102_gt_25hrs/'%cat[1].data['ID'][g])
            os.system('cp /Users/rsimons/all_png/GS4_%.5i.sed.png /Users/rsimons/Desktop/clear/figures/g102_gt_25hrs/'%cat[1].data['ID'][g])

        '''


        for ax in axes:

            ax.set_yscale('log')
        axes[0].set_ylabel('Number of Objects', fontsize = 12)


        pdf.savefig()