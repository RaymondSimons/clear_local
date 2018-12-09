import astropy
from astropy.io import fits
import glob
from glob import glob
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch
from astropy.coordinates import SkyCoord


plt.close('all')
plt.ioff()





fp_dir = '/Volumes/gdrive/clear/grizli_extractions/footprints'


fits_fls = glob(fp_dir + '/*fits')

figN = figure()
figS = figure()
axN = figN.add_subplot(111)
axS = figS.add_subplot(111)





xS0 = 53.07
dSx = 0.25
yS0 = -27.85
dSy = 0.25


axS.set_xlim(xS0, xS0+dSx)
axS.set_ylim(yS0, yS0+dSy)
axS.set_aspect(dSy/dSx)
axS.set_title('GOODS-S')




xN0 = -171
dNx = 0.4
yN0 = 62.1
dNy = 0.4


axN.set_xlim(xN0, xN0+dNx)
axN.set_ylim(yN0, yN0+dNy)
axN.set_aspect(dNy/dNx)
axN.set_title('GOODS-N')



for f in fits_fls:
    b = fits.open(f)

    footprints = b[1].data['footprint']

    BLUE = '#6699cc'
    for i, fp in enumerate(footprints):
        flt = b[1].data['filter'][i]
        if flt.startswith('G'):
            print flt
            exptime = b[1].data['exptime'][i]/60./60.
            print exptime/60./60.
            fp_s = fp.strip('POLYGON')
            l = []
            fp_split = fp_s.split(' ')
            try:
                fp_split= fp_split[1:len(fp_split)]
                coords = [(float(fp_split[2*i]), float(fp_split[2*i+1])) for i in arange(len(fp_split)/2)]
            except:
                fp_split= fp_split[2:len(fp_split)]
                coords = [(float(fp_split[2*i]), float(fp_split[2*i+1])) for i in arange(len(fp_split)/2)]

            print coords
            p = Polygon(coords)
            if coords[0][1] >0:

                if flt == 'G102': clr = 'blue'
                elif flt == 'G141': clr = 'red'
                patch = PolygonPatch(p, fc=clr, ec='black', alpha=min(exptime/28., 1.), zorder=2)
                axN.add_patch(patch)

            if coords[0][1] < 0:
                if flt == 'G102': clr = 'blue'
                elif flt == 'G141': clr = 'red'
                patch = PolygonPatch(p, fc=clr, ec='black', alpha=min(exptime/28., 1.), zorder=2)
                axS.add_patch(patch)



        #x, y = p.exterior.xy
        #ax.plot(x, y, color='#6699cc', alpha=0.7,
        #    linewidth=3, solid_capstyle='round', zorder=2)

figN.savefig('/Users/rsimons/Dropbox/rcs_clear/figures/coverage_goodsN.png', dpi = 300)
figS.savefig('/Users/rsimons/Dropbox/rcs_clear/figures/coverage_goodsS.png', dpi = 300)




'''
'''