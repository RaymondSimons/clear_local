import numpy as np
import astropy
from astropy.io import fits, ascii
plt.ioff()
plt.close('all')

fields = ['s', 'n']
fields = ['n']
cat_dir = '/Users/rsimons/Desktop/clear/Catalogs'
v44cat_dir = '/Users/rsimons/Dropbox/v44cats'
fig_dir = '/Users/rsimons/Desktop/clear/figures'

for f, field in enumerate(fields):
    if True:
        v41    = ascii.read(cat_dir + '/goods%s_3dhst.v4.1.cats/Catalog/goods%s_3dhst.v4.1.cat'%(field, field))
        #v43    = ascii.read(cat_dir + '/goods%s_3dhst.v4.3.cat'%(field))
        v44    = ascii.read(v44cat_dir + '/goods%s_3dhst.v4.4.cats/Catalog/goods%s_3dhst.v4.4.cat'%(field, field))
        #v43nzp = ascii.read(cat_dir + '/goods%s_3dhst.v4.3.nzpcat'%(field))


    fig, axes = subplots(1,4, figsize = (16,4))

    v44 = v44[0:len(v41)]
    m_Z_v41     = v41['f_B']
    m_850_v41   = v44['f_B']

    #m_Z_v43     = 25 - 2.5 * log10(array(v43['f_Z']))
    #m_850_v43   = 25 - 2.5 * log10(array(v43['f_F850LP']))

    #m_Z_v43nzp   = 25 - 2.5 * log10(array(v43nzp['f_Z']))
    #m_850_v43nzp = 25 - 2.5 * log10(array(v43nzp['f_F850LP']))

    axes[0].plot(m_Z_v41, m_850_v41/m_Z_v41, '.', markersize = 0.2)
    #axes[1].plot(m_Z_v43, m_Z_v41, '.', markersize = 0.2)
    #axes[2].plot(m_Z_v43nzp, m_Z_v41, '.', markersize = 0.2)
    #axes[3].plot(m_Z_v43, m_Z_v43nzp, '.', markersize = 0.2)

    for ax in axes:
        #ax.plot([15, 35], [15,35], 'k-')
        #ax.set_xlim(15, 30)
        ax.set_ylim(0.8, 1.2)

    fig.savefig(fig_dir + '/test.png')