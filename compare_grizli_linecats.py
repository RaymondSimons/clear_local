from astropy.io import fits
import os
fields = ['GS1',
 'GS2',
 'GS3',
 'GS4',
 'GS5',
 'GN1',
 'GN2',
 'GN3',
 'GN4',
 'GN5',
 'GN7',
 'ERSPRIME']



O3 = []
eO3 = []
di = []
flds = []
for f, fld in enumerate(fields):
    #old_cat = fits.open('/Users/rsimons/Desktop/grizli_v2.1_cats_old/%s_lines_grizli.fits'%fld)
    fld_cat = fits.open('/Users/rsimons/Desktop/grizli_v2.1_cats/%s_lines_grizli.fits'%fld)
    O3_fld = fld_cat[1].data['OIII-4363_FLUX']
    eO3_fld = fld_cat[1].data['OIII-4363_FLUX_ERR']
    id_fld = fld_cat[1].data['ID']

    gd = where(O3_fld > 0)[0]

    for g in gd:
        O3.append(O3_fld[g])
        eO3.append(eO3_fld[g])
        di.append(id_fld[g])
        flds.append(fld)
O3  = array(O3 )
eO3 = array(eO3)
di  = array(di )
flds = array(flds)

gd = where(O3/eO3 > 5)[0]

for g in gd:
    print flds[g], di[g], '%.2f %.2f'%(O3[g], eO3[g])
    os.system('cp /Volumes/pegasus/clear/team_release_v2.1/grizli_extractions_v2.1/all_png/%s_%.5i* /Users/rsimons/Desktop/clear/figures/OIII_4363/'%(flds[g], di[g]))

