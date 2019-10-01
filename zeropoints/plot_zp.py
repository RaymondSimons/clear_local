import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
plt.close('all')
plt.ioff()


fig, axes = plt.subplots(1,2, figsize = (8,4))

for f, field in enumerate(array(['gds', 'gdn'])):
    ax = axes[f]
    a = np.loadtxt('zp_v44_%s'%field, usecols = (1,2,3))
    ax.axhline(y = 1.0, linestyle = 'dashed', color = 'black')
    ax.plot(a[:,2], a[:,1], 'kD', label = 'Skelton+ v4.1 zp corrections')
    ax.plot(a[:,2], a[:,0], 'bo', label = 'new zp corrections')
    ax.set_xscale('log')
    ax.set_xlabel('central wavelength (um)')
    if f ==0 : ax.set_ylabel('zeropoint flux correction')
    ax.set_ylim(0.70, 1.25)
    if f == 1:
        ax2 = ax.twinx()

        #mag_ticks = array([-2.5*log10(fl) for fl in ax2.get_yticks()])
        mag_ticks = arange(-50, 50, 10)/100.

        ax2.set_yticks(10**(mag_ticks/-2.5))
        temp_a = ['', '+', '-']
        ax2yticks_str = array([temp_a[int(sign(mag))] + '%.1f'%mag for mag in mag_ticks])


        ax2.set_ylim(0.70, 1.25)


        ax2.set_yticklabels(ax2yticks_str)
        ax2.set_ylabel('zeropoint magnitude offset')
        ax.set_yticklabels([])
    ax.set_xticks([0.5, 1.0, 2, 3, 4, 5])
    ax.set_xticklabels(['0.5', '1', '2', '3', '4', '5'])
    ax.set_title('%s'%field.upper(), fontsize = 15)
    if f == 1: ax.legend(loc = 1)
    
fig.subplots_adjust(wspace = 0.05)
fig.savefig('/Users/rsimons/Desktop/clear/figures/zp_corrections.png', dpi = 300)