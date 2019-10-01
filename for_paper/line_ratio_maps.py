import astropy
from astropy.io import fits, ascii
from astropy.stats import bootstrap
import glob
from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import joblib
from numpy import *
from joblib import Parallel, delayed
from scipy.stats import median_absolute_deviation
from clear_local.utils import tools


plt.close('all')
plt.ioff()



def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    # Fast and numerically precise:
    variance = np.average((values-average)**2, weights=weights)
    return (average)#, math.sqrt(variance))


def make_big_plot(lr):
    fls = glob('/Users/rsimons/Desktop/clear/maps/line_ratio_maps/*line_ratios.fits')

    fig, axes = plt.subplots(27, 20, figsize = (20, 27))


    for f, fl in enumerate(fls):
        ax = axes.ravel()[f]
        lr_fits = fits.open(fl)
        hdu_names = [hdu.name for hdu in lr_fits]
        if '%s_S_EC'%lr in hdu_names:
            ax.imshow(np.log10(lr_fits['%s_S_EC'%lr].data[0]), cmap = 'viridis', vmin = -1, vmax = 1)


    for ax in axes.ravel():
        ax.set_xticks([])
        ax.set_yticks([])
        ax.axis('on')




    fig.subplots_adjust(left = 0.0, right = 1.0, top = 1.0, bottom = 0.0, hspace = 0.0, wspace = 0.0)

    fig.savefig('/Users/rsimons/Desktop/clear/figures/line_ratio_maps/%s_all.pdf'%lr)

    plt.close('all')




def save_stack(lr, mass_min = 7, mass_max = 13, do_all = False, typ = ''):
    if do_all:
        fls = glob('/Users/rsimons/Desktop/clear/maps/line_ratio_maps/*line_ratios.fits')

        sample_cat = tools.load_paper_catalog()
        lmass = sample_cat['lmass'].data
        flds  = sample_cat['field'].data
        dis   = sample_cat['id'].data
        xclass = sample_cat['xclass'].data

        stack_header = fits.Header()
        cnt = 0
        big_arr = np.nan*np.zeros((2, 40, 40, len(fls)))
        for f, fl in enumerate(fls):
            fld = fl.split('/')[-1].split('_')[0]
            di = int(fl.split('/')[-1].split('_')[1])
            where_cat = (flds == fld) & (dis == di)
            lms = lmass[where_cat]
            xcls = xclass[where_cat]
            if (lms > mass_min) & (lms < mass_max):
                lr_fits = fits.open(fl)
                hdu_names = [hdu.name for hdu in lr_fits]
                if '%s%s'%(lr, typ) in hdu_names:
                    non_nan = where(~np.isnan(lr_fits['%s%s'%(lr, typ)].data[:,:,:].ravel()))[0]
                    if (len(non_nan) > 0) & ('AGN' not in xcls):
                        big_arr[:,:,:,f] = lr_fits['%s%s'%(lr, typ)].data[:,:,:]
                        cnt+=1
        stack_header['count'] = cnt

        fits.writeto('/Users/rsimons/Desktop/clear/maps/line_ratio_maps/all/%s_%.2f_%.2f%s.fits'%(lr, mass_min, mass_max, typ), big_arr, header = stack_header, overwrite = True)

    big_arr = fits.open('/Users/rsimons/Desktop/clear/maps/line_ratio_maps/all/%s_%.2f_%.2f%s.fits'%(lr, mass_min, mass_max, typ))

    nx = big_arr.shape[1]
    ny = big_arr.shape[2]
    stacked_map = np.nan * np.zeros((2, nx, ny))
    for i in arange(nx):
        for j in arange(ny):
            good = where((~np.isnan(big_arr[0,i,j])) & (big_arr[0,i,j]/big_arr[1,i,j] > 0.5))[0]
            if len(good) > 5.:
                values  = big_arr[0,i,j, good]
                boostrap_values = bootstrap(values, bootnum = 100)
                mds = [median(bs_sample) for bs_sample in boostrap_values]

                stacked_map[0, i,j] = np.mean(mds)#, weights)
                stacked_map[1, i,j] = np.std(mds)#, weights)

    fits.writeto('/Users/rsimons/Desktop/clear/maps/line_ratio_maps/all/%s_%.2f_%.2f%s_stack.fits'%(lr, mass_min, mass_max, typ), stacked_map, header = stack_header, overwrite = True)


def make_stacked_figure(lrs,mass_min = 7, mass_max = 13, typ = '', mm = 0, mm_max = 1.):
    fig, axes = plt.subplots(1,4, figsize = (20,5))
    cmap = plt.cm.viridis
    cmap.set_bad('k')
    for l, lr in enumerate(lrs):
        stacked_map = fits.getdata('/Users/rsimons/Desktop/clear/maps/line_ratio_maps/all/%s_%.2f_%.2f%s_stack.fits'%(lr, mass_min, mass_max, typ))
        stacked_header = fits.getheader('/Users/rsimons/Desktop/clear/maps/line_ratio_maps/all/%s_%.2f_%.2f%s_stack.fits'%(lr, mass_min, mass_max, typ))


        sorted_ravel = sort(stacked_map[0].ravel())
        sorted_ravel = sorted_ravel[~np.isnan(sorted_ravel)]
        mad = median_absolute_deviation(sorted_ravel)
        sorted_ravel = sorted_ravel[abs(sorted_ravel - median(sorted_ravel))/mad < 5.]


        if len(sorted_ravel) > 5:
            vmn = sorted_ravel[int(0.10 * len(sorted_ravel))] - 0.5
            vmx = sorted_ravel[int(0.98 * len(sorted_ravel))] + 0.5

        else:
            vmn = 0.0
            vmx = 10.0


        axes.ravel()[l].imshow(stacked_map[0], interpolation = 'nearest', cmap = cmap, vmin = vmn, vmax = vmx)







        if mm == 0:
            axes.ravel()[l].annotate(lr, (0.05,0.98), va = 'top', ha = 'left',  xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = 50)
        if l == 0:
            if mm == 0:
                axes.ravel()[l].annotate(r'$\log$ M$_*$', (0.05,0.16), va = 'bottom', ha = 'left',  xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = 40)
                axes.ravel()[l].annotate("N={}\nstack".format(stacked_header['count']), (0.95, 0.98), ha = 'right', va = 'top', xycoords = 'axes fraction', color = 'grey', fontsize = 30)
                axes.ravel()[l].annotate("1''", (36, 15), ha = 'right', va = 'top', xycoords = 'data', color = 'white', fontsize = 30, fontweight = 'bold')
                axes.ravel()[l].plot([37,37], [15, 25], 'w-', linewidth = 4)

            else:
                axes.ravel()[l].annotate("{}".format(stacked_header['count']), (0.95, 0.98), ha = 'right', va = 'top', xycoords = 'axes fraction', color = 'grey', fontsize = 30)


            axes.ravel()[l].annotate(r'%.1f - %.1f M$_{\odot}$'%(mass_min, mass_max), (0.05,0.035), va = 'bottom', ha = 'left',  xycoords = 'axes fraction', color = 'white', fontweight = 'bold', fontsize = 40)


        else:
            axes.ravel()[l].annotate("{}".format(stacked_header['count']), (0.95, 0.98), ha = 'right', va = 'top', xycoords = 'axes fraction', color = 'grey', fontsize = 30)


    for ax in axes.ravel(): ax.axis('off')

    #axes[0,0].annotate('test', xy = (35, 35), xytext = (35, 25), xycoords = 'data', textcoords = 'data', arrowprops = dict(color = 'white', arrowstyle = '-'))

    #for ax in axes.ravel(): ax.plot([37,37], [28, 38], 'w-', linewidth = 4)

    fig.tight_layout()
    if 'EC' in typ:
        fig.savefig('/Users/rsimons/Desktop/clear/figures/for_paper/stacked_maps/EC/stacked_maps_%.4i_%.4i%s.pdf'%(100*mass_min, 100*mass_max, typ))
    else:
        fig.savefig('/Users/rsimons/Desktop/clear/figures/for_paper/stacked_maps/nonEC/stacked_maps_%.4i_%.4i%s.pdf'%(100*mass_min, 100*mass_max, typ))

    plt.close('all')



def stack_profile_plot(lr):
    fig, ax = plt.subplots(1,1, figsize = (8,6))
    stacked_map = fits.getdata('/Users/rsimons/Desktop/clear/maps/line_ratio_maps/all/%s_stack.fits'%lr)

    nx = stacked_map.shape[1]
    ny = stacked_map.shape[2]

    xs = arange(-nx/2., nx/2.)*0.1
    ys = arange(-nx/2., nx/2.)*0.1

    X, Y = np.meshgrid(xs, ys)

    R = np.sqrt(X**2. + Y**2.)




    ax.errorbar(R.ravel(), stacked_map[0].ravel(), yerr = stacked_map[1].ravel(), linestyle = 'None', fmt = 'o')





    fig.savefig('/Users/rsimons/Desktop/clear/figures/line_ratio_maps/stack/%s_stack_profiles.pdf'%lr)
    plt.close('all')




dx = 0.50
mass_arrs = [(xmn, xmn+dx) for xmn in arange(9, 11.5, dx)]


lrs = ['O32', 'R2', 'R3', 'R23']
typs =  ['', '_S', '_EC', '_S_EC']
typs =  ['_EC', '_S_EC']
typs =  ['_S_EC']


#Parallel(n_jobs = -1)(delayed(do_all_types)(typ) for typ in typs)
#Parallel(n_jobs = -1)(delayed(stack_profile_plot)(lr) for lr in lrs)



for typ in typs:
    Parallel(n_jobs = -1)(delayed(make_stacked_figure)(lrs, mass_min, mass_max, typ, mm, mm_max = len(mass_arrs)) for mm, (mass_min, mass_max) in enumerate(mass_arrs))




































