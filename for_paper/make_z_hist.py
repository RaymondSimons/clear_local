import astropy
from astropy.io import ascii
import matplotlib
import matplotlib.pyplot as plt
from numpy import *
plt.ioff()
plt.close('all')


simons_cat = ascii.read('/Users/rsimons/Desktop/clear/catalogs/simons_sample.cat')




z = simons_cat['z_map'].data

lz = log10(1+z)


fig = plt.figure(figsize = (4,3.))
ax = fig.add_subplot(1,1,1)



ax.hist(lz, bins = linspace(log10(1+0.5),log10(1+3.5), 20),color = 'black')


ax.set_xlabel('')
ax.set_ylabel('number of objects', fontsize = 16)

ax.annotate('sample', (0.05, 0.90), va = 'top', xycoords = 'axes fraction', fontsize = 20)

xticks = arange(0, 4, 1)
xticklabels = ['%i'%i for i in xticks]
lxticks = log10(1+xticks)


ax.set_xticks(lxticks)
#ax.set_xticklabels(xticklabels)
ax.set_xticklabels([''])

ax.set_xlim(log10(1), log10(4.9))

fig.tight_layout()
fig.savefig('/Users/rsimons/Desktop/clear/figures/for_paper/z_hist.png', dpi = 800)
