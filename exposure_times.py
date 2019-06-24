import astropy
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from matplotlib import *
plt.ioff()
plt.close('all')
mpl.rcParams['text.usetex'] = True
mpl.rcParams['axes.linewidth'] = 2
mpl.rcParams['ytick.minor.size'] = 2
mpl.rcParams['ytick.major.size'] = 4 
mpl.rcParams['xtick.minor.size'] = 2 
mpl.rcParams['xtick.major.size'] = 4 


with PdfPages('/Users/rsimons/Desktop/clear/figures/exposure_times.pdf') as pdf:

    fields = ['GS1','GS2', 'GS3', 'GS4', 'GS5', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7', 'ERSPRIME']
    t102_all = np.array([])
    t141_all = np.array([])
    NG102 = 0
    NG102_G141 = 0

    for field in fields:
        cat = fits.open('/Users/rsimons/Desktop/clear/Catalogs/grizli_v2.1_cats/%s_lines_grizli.fits'%field)
        t102 = cat[1].data['T_G102']
        t141 = cat[1].data['T_G141']

        NG102 += len(where(t102 > 0)[0])
        NG102_G141 += len(where((t102 > 0) & (t141 > 0))[0])

        t102_all = concatenate([t102_all, t102])
        t141_all = concatenate([t141_all, t141])

        fig, axes = plt.subplots(1,2, figsize = (8.5, 4.5))

        axes[0].hist(t102/3600.,  label = field, bins = linspace(0, 30, 40), alpha = 1.0, color = 'darkblue')
        axes[1].hist(t141/3600.,  label = field, bins = linspace(0, 30, 40), alpha = 1.0, color = 'darkred')

        #axes[0].annotate('%s'%field, (0.4, 1.0), xycoords = 'figure fraction', fontweight = 'bold', fontsize = 30)
        axes[0].set_title('%s'%field, fontweight = 'bold', fontsize = 30)

        axes[0].set_xlabel('G102 Exposure Time (hr)', fontsize = 16)
        axes[1].set_xlabel('G141 Exposure Time (hr)', fontsize = 16)


        for ax in axes: 
            ax.set_yscale('log')
            ax.tick_params(axis="x", labelsize=15)
            ax.tick_params(axis="y", labelsize=15)
            ax.set_ylim(0.7, 1000)
        axes[0].set_ylabel('Number of Objects', fontsize = 16)

        fig.tight_layout()
        pdf.savefig()


    fig, axes = plt.subplots(1,2, figsize = (8.5, 4.5))

    axes[0].hist(t102_all/3600.,  label = field, bins = linspace(0, 30, 40), alpha = 1.0, color = 'darkblue')
    axes[1].hist(t141_all/3600.,  label = field, bins = linspace(0, 30, 40), alpha = 1.0, color = 'darkred')

    #axes[0].annotate('%s'%field, (0.4, 1.0), xycoords = 'figure fraction', fontweight = 'bold', fontsize = 30)
    axes[0].set_title('G102', color = 'darkblue', fontweight = 'bold', fontsize = 30)
    axes[1].set_title('G141', color = 'darkred', fontweight = 'bold', fontsize = 30)

    for ax in axes: 
        ax.set_yscale('log')
        ax.tick_params(axis="x", labelsize=15)
        ax.tick_params(axis="y", labelsize=15)
        ax.set_ylim(0.7, 6000)

    axes[0].set_xlabel('G102 Exposure Time (hr)', fontsize = 16)
    axes[1].set_xlabel('G141 Exposure Time (hr)', fontsize = 16)
    axes[0].set_ylabel('Number of Objects', fontsize = 16)

    fig.tight_layout()
    pdf.savefig()


























