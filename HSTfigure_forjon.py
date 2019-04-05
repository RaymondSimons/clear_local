import astropy
from astropy.io import fits, ascii
from astropy.convolution import Gaussian2DKernel, convolve_fft, Box2DKernel
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}'] 

plt.ioff()
plt.close('all')

def O32_OH(O3, O2, eO3, eO2):
    O32_arr = zeros(O3.shape) * nan
    O32_e_arr = zeros(O3.shape) * nan

    OH_arr = zeros(O3.shape) * nan
    OH_e_arr = zeros(O3.shape) * nan

    for i in arange(O3.shape[0]):
        for j in arange(O3.shape[1]):
            O2_arr_temp = np.random.normal(O2[i,j], eO2[i,j], 200)
            O3_arr_temp = np.random.normal(O3[i,j], eO3[i,j], 200)
            O32_arr[i,j] = O3[i,j]/O2[i,j]
            OH_arr[i,j]   = np.nanmean(8.54 - 0.59 * O3_arr_temp/O2_arr_temp)
    
    return O32_arr, O32_e_arr, OH_arr, OH_e_arr#np.std(OH_z_arr)

fig = figure(figsize = (8, 5.2))




objects = array([('GS5', 42159), ('GS2', 46938)])




vmn_O2 = [-0.05, -0.05]
vmn_Hb = [-0.05, -0.05]
vmn_O3 = [-0.05, -0.05]


vmx_O2 = [0.12, 0.33]
vmx_Hb = [0.12, 0.30]
vmx_O3 = [0.12, 0.6]

vmn_z = [8.35, 8.2]
vmx_z = [8.75, 8.6]

sn_lim = [0.4, 0.75]

cm = mpl.cm.cool
cm.set_bad('black', 1.)
kern = Box2DKernel(3.)

axes_coords = [[0.04, 0.08, 0.18, 0.02], [0.55, 0.08, 0.18, 0.02]]

for o, obj in enumerate(objects):
    axes_O2   = plt.subplot2grid((4, 19), (0,0 + o*10), rowspan = 1, colspan = 3, fig = fig)
    axes_Hb   = plt.subplot2grid((4, 19), (0,3 + o*10), rowspan = 1, colspan = 3, fig = fig)
    axes_O3   = plt.subplot2grid((4, 19), (0,6 + o*10), rowspan = 1, colspan = 3, fig = fig)
    axes_Z    = plt.subplot2grid((4, 19), (1,0 + o*10), rowspan = 3, colspan = 9, fig = fig)
    axes_list = [axes_O2, axes_O3, axes_Hb, axes_Z]

    full = fits.open('/Users/rsimons/Desktop/random/%s_%s.full.fits'%(obj[0], obj[1]))

    O2, eO2 = full['LINE', 'OII'].data, 1./np.sqrt(full['LINEWHT', 'OII'].data)
    O3, eO3 = full['LINE', 'OIII'].data, 1./np.sqrt(full['LINEWHT', 'OIII'].data)    
    Hb, eHb = full['LINE', 'Hb'].data, 1./np.sqrt(full['LINEWHT', 'Hb'].data)    

    O2 = convolve_fft(O2, kern)
    Hb = convolve_fft(Hb, kern)
    O3 = convolve_fft(O3, kern)


    axes_O2.imshow(flipud(O2), vmin = vmn_O2[o], vmax = vmx_O2[o])
    axes_O3.imshow(flipud(O3), vmin = vmn_O3[o], vmax = vmx_O3[o])
    axes_Hb.imshow(flipud(Hb), vmin = vmn_Hb[o], vmax = vmx_Hb[o])

    Z = 8.54 - 0.59 * log10(O3/O2)
    if o == 0: eO2[:,43:] = np.inf

    Z_masked =np.ma.masked_where((O2/eO2 < sn_lim[o]) | (O3/eO3 < sn_lim[o]), Z)

    mp = axes_Z.imshow(flipud(Z_masked), vmin = vmn_z[o], vmax = vmx_z[o], cmap = cm)

    cbaxes = fig.add_axes(axes_coords[o]) 
    cbr = plt.colorbar(mp, cax = cbaxes, orientation = 'horizontal')
    cbr.ax.axes.tick_params(color = 'white', labelcolor = 'white')
    cbr.set_ticks(arange(7, 10, 0.2))

    for ax in axes_list:
        ax.set_xticks([])
        ax.set_yticks([])


    #axes_Z.annotate(r'12 + $\log$(O/H)', (0.52, 0.13),fontsize = 18,  xycoords = 'axes fraction', fontweight = 'bold', color = 'white')
    if o == 0: 
        axes_Z.annotate(r'12 + $\log$(O/H)', (0.05, 0.13),fontsize = 16,  xycoords = 'axes fraction', fontweight = 'bold', color = 'white')






    if o == 0:
        axes_Z.annotate('gas-phase\nmetallicity', (0.97, 0.03), ha = 'right', fontsize = 20,  xycoords = 'axes fraction', fontweight = 'bold', color = 'white')        
        txt2 = axes_Z.annotate('metal-rich\ncenter', (0.31, 0.54) , ha = 'right', color = cm(np.inf), fontweight = 'bold', xycoords = 'axes fraction', fontsize = 18)
        txt1 = axes_Z.annotate('metal-poor\noutskirts',    (0.36, 0.25), ha = 'right',  color = cm(0),fontweight = 'bold', xycoords = 'axes fraction', fontsize = 18)

        axes_Z.annotate('', xy=(0.44, 0.51), xytext=(0.28, 0.52), xycoords = 'axes fraction', arrowprops=dict(facecolor=cm(np.inf), connectionstyle = "arc3, rad = 0.2", shrink=0.05))
        axes_Z.annotate('', xy=(0.50, 0.26), xytext=(0.365, 0.24), xycoords = 'axes fraction', arrowprops=dict(facecolor=cm(0), connectionstyle = "arc3, rad = 0.2", shrink=0.05))


    if o == 1:
        txt1 = axes_Z.annotate('metal-poor\ncenter',(0.33, 0.54), ha = 'right',  color = cm(0),fontweight = 'bold', xycoords = 'axes fraction', fontsize = 18)
        txt2 = axes_Z.annotate('metal-rich\noutskirts', (0.29, 0.2), ha = 'right', color = cm(np.inf), fontweight = 'bold', xycoords = 'axes fraction', fontsize = 18)
        axes_Z.annotate('', xy=(0.36, 0.48), xytext=(0.25, 0.52), xycoords = 'axes fraction', arrowprops=dict(facecolor=cm(0), connectionstyle = "arc3, rad = 0.2", shrink=0.05))
        axes_Z.annotate('', xy=(0.38, 0.20), xytext=(0.21, 0.18), xycoords = 'axes fraction', arrowprops=dict(facecolor=cm(np.inf), connectionstyle = "arc3, rad = 0.2", shrink=0.05))


    #import matplotlib.patheffects as PathEffects
    #for txt in [txt1, txt2]:
    #    txt.set_path_effects([PathEffects.withStroke(linewidth=1.0, foreground='grey')])





    axes_O2.annotate('[OII]', (0.97, 0.03), va = 'bottom', ha = 'right', fontsize = 20, xycoords = 'axes fraction', fontweight = 'bold', color = 'white')
    axes_Hb.annotate(r'H$\beta$', (0.97, 0.03), va = 'bottom', ha = 'right',fontsize = 20,  xycoords = 'axes fraction', fontweight = 'bold', color = 'white')
    axes_O3.annotate('[OIII]', (0.97, 0.03), va = 'bottom', ha = 'right',fontsize = 20,  xycoords = 'axes fraction', fontweight = 'bold', color = 'white')



    for ax in [axes_O2, axes_O3, axes_Hb]:
        ax.set_xticks([])
        ax.set_yticks([])
        dx = -2
        dy = -1
        ax.set_xlim(32 + dx, 48 + dx)
        ax.set_ylim(32 + dy, 48 + dy)

    for ax in [axes_Z]: 
        ax.set_xticks([])
        ax.set_yticks([])
        dx = -3
        dy = 0
        ax.set_xlim(32 + dx, 48 + dx)
        ax.set_ylim(32 + dy, 48 + dy)


    #axes_Z.plot([44, 44], [45, 47], '-', color = 'white', linewidth = 4)
    #if o == 0: axes_Z.annotate("0.2''", (42, 46), color = 'white', fontsize = 15)




fig.subplots_adjust(hspace = 0.0, wspace = 0.0, left = 0.02, right = 0.98, top = 0.98, bottom = 0.02)
#fig.tight_layout()
fig.savefig('/Users/rsimons/Desktop/random/Z_objects_HST.png', dpi = 600)

plt.close('all')



    #full = fits.open('/Users/rsimons/Desktop/random/%s_%i.full.fits'%(field, di))

'''

lineo2 = full['LINE', 'OII'].data
lineo3 = full['LINE', 'OIII'].data

elineo2 = 1./np.sqrt(full['LINEWHT', 'OII'].data)
elineo3 = 1./np.sqrt(full['LINEWHT', 'OIII'].data)


kern = Box2DKernel(3.)
#kern = Gaussian2DKernel(0.8)

O2 = convolve_fft(lineo2, kern)
O3 = convolve_fft(lineo3, kern)





Z = 8.54 - 0.59 * O3/O2


#elineo2[:,:] = np.inf



sn_lim = 0.40
Z_masked =np.ma.masked_where((O2/elineo2 < sn_lim) | (O3/elineo3 < sn_lim), Z)
cm = mpl.cm.cool
cm.set_bad('black', 1.)

axes_O2.imshow(flipud(O2), vmin = -0.05, vmax = 0.15)
axes_O3.imshow(flipud(O3), vmin = -0.01, vmax = 0.10)

zmn = 6.8 
zmx = 8.1
zmn = 7.4
zmx = 8.5

mp = axes_Z.imshow(flipud(Z_masked), vmin = zmn, vmax = zmx, cmap = cm)



#txt1 = axes_Z.annotate('metal-poor\ncenter',(0.32, 0.54), ha = 'right',  color = cm(0),fontweight = 'bold', xycoords = 'axes fraction', fontsize = 15)
#txt2 = axes_Z.annotate('metal-rich\noutskirts', (0.28, 0.2), ha = 'right', color = cm(np.inf), fontweight = 'bold', xycoords = 'axes fraction', fontsize = 15)
#axes_Z.annotate('', xy=(0.36, 0.48), xytext=(0.25, 0.52), xycoords = 'axes fraction', arrowprops=dict(facecolor=cm(0), connectionstyle = "arc3, rad = 0.2", shrink=0.05))
#axes_Z.annotate('', xy=(0.38, 0.20), xytext=(0.21, 0.18), xycoords = 'axes fraction', arrowprops=dict(facecolor=cm(np.inf), connectionstyle = "arc3, rad = 0.2", shrink=0.05))



txt2 = axes_Z.annotate('metal-rich\ncenter', (0.31, 0.54) , ha = 'right', color = cm(np.inf), fontweight = 'bold', xycoords = 'axes fraction', fontsize = 15)
txt1 = axes_Z.annotate('metal-poor\noutskirts',    (0.36, 0.25), ha = 'right',  color = cm(0),fontweight = 'bold', xycoords = 'axes fraction', fontsize = 15)

axes_Z.annotate('', xy=(0.40, 0.54), xytext=(0.28, 0.52), xycoords = 'axes fraction', arrowprops=dict(facecolor=cm(np.inf), connectionstyle = "arc3, rad = 0.2", shrink=0.05))
axes_Z.annotate('', xy=(0.50, 0.26), xytext=(0.365, 0.24), xycoords = 'axes fraction', arrowprops=dict(facecolor=cm(0), connectionstyle = "arc3, rad = 0.2", shrink=0.05))


import matplotlib.patheffects as PathEffects
#for txt in [txt1, txt2]:
#    txt.set_path_effects([PathEffects.withStroke(linewidth=0.8, foreground='grey')])





for ax in [axes_O2, axes_O3]:
    ax.set_xticks([])
    ax.set_yticks([])
    dx = -2
    dy = -1
    ax.set_xlim(32 + dx, 49 + dx)
    ax.set_ylim(32 + dy, 49 + dy)

for ax in [axes_Z]: 
    ax.set_xticks([])
    ax.set_yticks([])
    dx = -3
    dy = -1
    ax.set_xlim(32 + dx, 49 + dx)
    ax.set_ylim(32 + dy, 49 + dy)


axes_O2.annotate('[OII]', (0.64, 0.06), fontsize = 20, xycoords = 'axes fraction', fontweight = 'bold', color = 'white')
axes_O3.annotate('[OIII]', (0.64, 0.06),fontsize = 20,  xycoords = 'axes fraction', fontweight = 'bold', color = 'white')
#axes_Z.annotate(r'12 + $\log$(O/H)', (0.52, 0.13),fontsize = 18,  xycoords = 'axes fraction', fontweight = 'bold', color = 'white')
axes_Z.annotate(r'12 + $\log$(O/H)', (0.50, 0.13),fontsize = 18,  xycoords = 'axes fraction', fontweight = 'bold', color = 'white')
#axes_Z.annotate('gas-phase\nmetallicity', (0.98, 0.83), ha = 'right', fontsize = 20,  xycoords = 'axes fraction', fontweight = 'bold', color = 'white')



axes_p.set_ylabel(r'12 + log(O/H)', fontsize = 18)
axes_p.set_xlabel('distance from center (arcsec)', fontsize = 18)




cbaxes = fig.add_axes([0.265, 0.22, 0.18, 0.02]) 
cbr = plt.colorbar(mp, cax = cbaxes, orientation = 'horizontal')
cbr.ax.axes.tick_params(color = 'white', labelcolor = 'white')
cbr.set_ticks([7, 7.5, 8, 8.5])

cat = ascii.read('/Users/rsimons/Dropbox/rcs_clear/catalogs/z_r_O32.cat')






x = linspace(0, 2., 10)
for c in cat[0:50]:
    if (c[3] > 7) & (c[3] < 9):
        axes_p.plot(x, x * c[2] + c[3], color = 'black', alpha = 0.2)

axes_p.set_xlim(0, 1.5)
axes_p.set_ylim(6.5,9.5)



'''

#for c in cat

#fig.tight_layout()
