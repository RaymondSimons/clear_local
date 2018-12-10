import astropy
from astropy.io import fits
import glob
from glob import glob
from shapely.geometry.polygon import Polygon
from descartes import PolygonPatch
from astropy.coordinates import SkyCoord


matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['axes.linewidth'] = 3

plt.close('all')
plt.ioff()





fp_dir = '/Volumes/gdrive/clear/grizli_extractions/footprints'


fits_fls = glob(fp_dir + '/*fits')

figN = figure()
figS = figure()
axN = figN.add_subplot(111)
axS = figS.add_subplot(111)



c0 = SkyCoord("3h33m00s -27d58m00s")
c1 = SkyCoord("3h32m00s -27d40m00s")





gds_dticks = ['-27d40m00s', '-27d45m00s', '-27d50m00s', '-27d50m00s', '-27d55m00s']
gds_aticks = ['3h32m00s', '3h32m30s', '3h33m00s', '3h33m30s']

gdn_dticks = ['+62d05m00s', '+62d10m00s', '+62d15m00s', '+62d20m00s', '+62d25m00s']
gdn_aticks = ['12h35m00s', '12h36m00s', '12h37m00s', '12h38m00s']




for ax, aticks, dticks, tit in [(axN,  gdn_aticks, gdn_dticks, 'GOODS-N'), (axS,  gds_aticks, gds_dticks,'GOODS-S')]:

    dt = [SkyCoord("3h33m00s %s"%g).dec.value for g in dticks]

    if tit == "GOODS-N": 
        at = [SkyCoord("%s -27d40m00s"%g).ra.value - 360. for g in aticks]
    else:
        at = [SkyCoord("%s -27d40m00s"%g).ra.value for g in aticks]





    dticklbls = ["$%s$"%dtick.replace('d', '^\circ').replace('m', '^m').replace('s','^s') for dtick in dticks]
    aticklbls = ["$%s$"%atick.replace('h', '^h').replace('m', '^m').replace('s','^s') for atick in aticks]

    ax.set_yticks(dt)
    ax.set_yticklabels(dticklbls)
    ax.set_xticks(at)
    ax.set_xticklabels(aticklbls)





for ax, c0, c1, tit in [(axS, SkyCoord("3h33m10s -27d58m00s"), SkyCoord("3h31m50s -27d38m00s"), 'GOODS-S'), (axN, SkyCoord("12h38m20s +62d04m45s"), SkyCoord("12h35m30s  +62d25m00s"), 'GOODS-N')]:


    yS0 = c0.dec.value
    yS1 = c1.dec.value

    if tit == "GOODS-N":
        xS0 = c0.ra.value - 360.
        xS1 = c1.ra.value - 360.

    if tit == "GOODS-S":
        xS0 = c0.ra.value
        xS1 = c1.ra.value


    ax.set_xlim(xS0, xS1)
    ax.set_ylim(yS1, yS0)
    ax.set_title(tit)










cs = [("3h32m20s -27d57m00s"), ("3h33m05s -27d54m00s"), ("3h32m40s -27d42m00s"), ("3h32m00s -27d44m00s")]

goods_coords = [(SkyCoord(c).ra.value, SkyCoord(c).dec.value) for c in cs]

gds_p = Polygon(goods_coords)
patch = PolygonPatch(gds_p, fc='grey', ec='black', alpha = 0.2, zorder=1)
axS.add_patch(patch)


cs = [("12h38m10s +62d15m00s"), ("12h37m00s +62d22m00s"), ("12h35m40s +62d10m00s"), (" 12h36m30s +62d05m00s")]
goodn_coords = [(SkyCoord(c).ra.value-360., SkyCoord(c).dec.value) for c in cs]
gdn_p = Polygon(goodn_coords)
patch = PolygonPatch(gdn_p, fc='grey', ec='black', alpha = 0.2, zorder=1)
axN.add_patch(patch)




for f in fits_fls:
    b = fits.open(f)
    footprints = b[1].data['footprint']
    BLUE = '#6699cc'
    for i, fp in enumerate(footprints):
        flt = b[1].data['filter'][i]
        if flt.startswith('G'):

            exptime = b[1].data['exptime'][i]/60./60.
            fp_s = fp.strip('POLYGON')
            l = []
            fp_split = fp_s.split(' ')
            try:
                fp_split= fp_split[1:len(fp_split)]
                coords = [(float(fp_split[2*i]), float(fp_split[2*i+1])) for i in arange(len(fp_split)/2)]
            except:
                fp_split= fp_split[2:len(fp_split)]
                coords = [(float(fp_split[2*i]), float(fp_split[2*i+1])) for i in arange(len(fp_split)/2)]
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


for ax, fig in [(axN, figN), (axS, figS)]:


    xr, yr = ax.get_xlim(), ax.get_ylim() 
    dx = (xr[1]-xr[0])*np.cos(yr[0]/180*np.pi)*60
    dy = (yr[1]-yr[0])*60
    ax.set_xlim(ax.get_xlim()[::-1])
    fig.set_size_inches(5,5*dy/dx)

    fig.subplots_adjust(left = 0.20, right = 0.95)



for ax in [axN, axS]:
    fs = 20
    ax.annotate('CANDELS WFC3', (0.06, 0.03), color = 'grey', xycoords = 'axes fraction', fontsize = fs)
    ax.annotate('G102', (0.06, 0.15), color = 'blue', fontweight = 'bold', xycoords = 'axes fraction', fontsize = fs)
    ax.annotate('G141',(0.06, 0.09),  color = 'red',  fontweight = 'bold',xycoords = 'axes fraction', fontsize = fs)

figN.savefig('/Users/rsimons/Dropbox/rcs_clear/figures/coverage_goodsN.png', dpi = 600)
figS.savefig('/Users/rsimons/Dropbox/rcs_clear/figures/coverage_goodsS.png', dpi = 600)




'''
'''