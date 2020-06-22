import astropy
from astropy.io import ascii
import numpy as np


fields = [('s', 'v4.3'), ('n', 'v4.4')]

cat_dir = '/Users/rsimons/Downloads'
for (f, v) in fields:
    cat_name = '%s/goods%s-F105W-astrodrizzle-%s_drz_sub_plus.cat'%(cat_dir, f, v)
    cat      = ascii.read(cat_name)
    print (max(cat['NUMBER']))
    cat_unmatched = ascii.read('%s/z67_in_CLEAR_G%s_unmatched_Finkelstein.txt'%(cat_dir,f))
    cat_matched   = ascii.read('%s/z67_in_CLEAR_G%s_matched_Finkelstein.txt'%(cat_dir,f))


    for i in np.arange(len(cat_unmatched)):
        di = cat_unmatched['ID_Finkelstein'][i]
        gd = np.where(di == cat['NUMBER'])[0]

        if len(gd) == 0:
            #no match found in astrodrizzle catalog, as expected
            cat.add_row(None) #add empty row to end of catalog
            cat[-1]['NUMBER'] = di
            cat[-1]['X_IMAGE'] = cat_unmatched['x-coord'][i]
            cat[-1]['Y_IMAGE'] = cat_unmatched['y-coord'][i]
            cat[-1]['X_WORLD'] = cat_unmatched['RA_Fin'][i]
            cat[-1]['Y_WORLD'] = cat_unmatched['DEC_Fin'][i]

        if len(gd) > 0:
            #this isn't supposed to happen, we found a match between the unmatched catalog and our existing catalog
            x = cat_unmatched['x-coord'][i]
            y = cat_unmatched['y-coord'][i]
            x_ad = cat['X_IMAGE'][gd[0]]
            y_ad = cat['Y_IMAGE'][gd[0]]

            ra   = cat_unmatched['RA_Fin'][i]
            dec  = cat_unmatched['DEC_Fin'][i]

            ra_ad   = cat['X_WORLD'][gd[0]]
            dec_ad  = cat['Y_WORLD'][gd[0]]

            print (f, di)
            print ('Unmatched')
            print ('Finkelstein Catalog      Astrodrizzle catalog')
            print ('X: ',x, x_ad)
            print ('Y: ',y, y_ad)

            print ('RA: ', ra,  ra_ad)
            print ('DEC: ', dec, dec_ad)

            print ('Diff: ', 3600 * np.sqrt((ra - ra_ad)**2. + (dec - dec_ad)**2))

            print ('\n\n')


    for i in np.arange(len(cat_matched)):
        di = cat_matched['ID_3DHST'][i]
        gd = np.where(di == cat['NUMBER'])[0]
        if len(gd) == 0:
            #this isn't supposed to happen, we don't have this object in our catalog
            print (f, di)

    print ('Writing output catalog for GOODS-%s...'%f)
    ascii.write(table = cat, output = cat_name.replace('plus.cat', 'plus_unmatched.cat'), overwrite=True,
                format = 'commented_header')















