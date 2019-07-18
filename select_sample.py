import astropy
from astropy.io import fits
import glob
from glob import glob
from numpy import *

#cat_fls = glob('/Users/rsimons/Desktop/clear/Catalogs/grizli_v2.1_cats/*_lines_grizli.fits')


if False: cat_dir = '/Users/rsimons/Dropbox/rcs_clear/catalogs'
else: cat_dir = '/user/rsimons/grizli_extractions/Catalogs'



cat_fls = glob(cat_dir + '/grizli_v2.1_cats/*_lines_grizli.fits')

#all of the strong line metallicity diagnostics
#not including those with N2 (which is blended in our data)
















f_R23   = open(cat_dir + '/sample_cats/R23_sample.cat', 'w+')
f_R2    = open(cat_dir + '/sample_cats/R2_sample.cat', 'w+')
f_R3    = open(cat_dir + '/sample_cats/R3_sample.cat', 'w+')
f_S2    = open(cat_dir + '/sample_cats/S2_sample.cat', 'w+')
f_O32   = open(cat_dir + '/sample_cats/O32_sample.cat', 'w+')
f_O3    = open(cat_dir + '/sample_cats/O3_sample.cat', 'w+')
f_O2    = open(cat_dir + '/sample_cats/O2_sample.cat', 'w+')
f_Ne3O2 = open(cat_dir + '/sample_cats/Ne3O2_sample.cat', 'w+')
f_O3S2  = open(cat_dir + '/sample_cats/O3S2_sample.cat', 'w+')
f_HaHb  = open(cat_dir + '/sample_cats/HaHb_sample.cat', 'w+')
f_any   = open(cat_dir + '/sample_cats/any_sample.cat', 'w+')

f_any.write('#these objects have two of the three at 5sigma: [OII], [OIII], Hb\n')

tot_R23    = 0
tot_R2     = 0
tot_R3     = 0
tot_S2     = 0
tot_O32    = 0
tot_O3    = 0
tot_O2    = 0
tot_Ne3O2  = 0
tot_O3S2   = 0
tot_HaHb   = 0


tot_any    = 0
tot_2      = 0
tot_3      = 0
tot_4      = 0
tot_5      = 0
tot_6      = 0
tot_7      = 0





for c, cat_fl in enumerate(cat_fls):
    cat = fits.open(cat_fl)
    #cat_data = cat[1].data
    cat_data = cat[1].data[cat[1].data['T_G102'] > 0.]
    sn = 5

    O2_f   = cat_data['OII_FLUX']
    O2_ef  = cat_data['OII_FLUX_ERR']
    O3_f   = cat_data['OIII_FLUX']
    O3_ef  = cat_data['OIII_FLUX_ERR']
    Ha_f   = cat_data['Ha_FLUX']
    Ha_ef  = cat_data['Ha_FLUX_ERR']
    Hb_f   = cat_data['Hb_FLUX']
    Hb_ef  = cat_data['Hb_FLUX_ERR']
    S2_f   = cat_data['SII_FLUX']
    S2_ef  = cat_data['SII_FLUX_ERR']
    Ne3_f  = cat_data['NeIII_FLUX']
    Ne3_ef = cat_data['NeIII_FLUX_ERR']

    O2_sn  = O2_f/O2_ef
    O3_sn  = O3_f/O3_ef  
    Ha_sn  = Ha_f/Ha_ef
    Hb_sn  = Hb_f/Hb_ef  
    S2_sn  = S2_f/S2_ef 
    Ne3_sn = Ne3_f/Ne3_ef



    bw_good_R23   = (O2_sn > sn)  & (O3_sn > sn) & (Hb_sn > sn)
    bw_good_R2    = (O2_sn > sn)  & (Hb_sn > sn)
    bw_good_R3    = (O3_sn > sn)  & (Hb_sn > sn)
    bw_good_S2    = (S2_sn > sn)  & (Ha_sn > sn)
    bw_good_O32   = (O2_sn > sn)  & (O3_sn > sn)
    bw_good_O3    = (O3_sn > sn)  & (Hb_sn > sn)
    bw_good_O2    = (O2_sn > sn)  & (Hb_sn > sn)
    bw_good_Ne3O2 = (Ne3_sn > sn) & (O2_sn > sn)
    bw_good_O3S2  = (O3_sn > sn) & (Hb_sn > sn) & (S2_sn > sn)  & (Ha_sn > sn)
    bw_good_HaHb  = (Ha_sn > sn) & (Hb_sn > sn)



    good_R23   = where(bw_good_R23  )[0]
    good_R2    = where(bw_good_R2   )[0]
    good_R3    = where(bw_good_R3   )[0]
    good_S2    = where(bw_good_S2   )[0]
    good_O32   = where(bw_good_O32  )[0]
    good_O3    = where(bw_good_O3  )[0]
    good_O2    = where(bw_good_O2  )[0]
    good_Ne3O2 = where(bw_good_Ne3O2)[0]
    good_O3S2  = where(bw_good_O3S2 )[0]
    good_HaHb  = where(bw_good_HaHb )[0]

    good_any = unique(concatenate((good_R2,   
                                  good_R3,   
                                  good_O32)))







    #gd_any = where(bw_good_R23 | bw_good_R2 | bw_good_R3 | bw_good_S2 | bw_good_O32 | bw_good_Ne3O2 | bw_good_O3S2)[0]
    gd_num = bw_good_R23.astype('int') + bw_good_R2.astype('int') + bw_good_R3.astype('int') + bw_good_S2.astype('int') + bw_good_O32.astype('int') + bw_good_O3.astype('int') + bw_good_O2.astype('int') + bw_good_Ne3O2.astype('int') + bw_good_O3S2.astype('int') + bw_good_HaHb.astype('int')
    gd_any = where(gd_num > 0)[0]
    gd_2 = where(gd_num > 1)[0]
    gd_3 = where(gd_num > 2)[0]
    gd_4 = where(gd_num > 3)[0]
    gd_5 = where(gd_num > 4)[0]
    gd_6 = where(gd_num > 5)[0]

    #print field, 'all 6 lines',  cat_data['ID'][gd_6]

    tot_any  += len(gd_any)
    tot_2 += len(gd_2)
    tot_3 += len(gd_3)
    tot_4 += len(gd_4)
    tot_5 += len(gd_5)
    tot_6 += len(gd_6)





    tot_R23   += len(good_R23  )
    tot_R2    += len(good_R2   )
    tot_R3    += len(good_R3   )
    tot_S2    += len(good_S2   )
    tot_O32   += len(good_O32  )
    tot_O3   += len(good_O3  )
    tot_O2   += len(good_O2  )
    tot_Ne3O2 += len(good_Ne3O2)
    tot_O3S2  += len(good_O3S2 )
    tot_HaHb  += len(good_HaHb )



    field = cat_fl.split('/')[-1].strip('lines_grizli.fits')
    print (field)

    print ('\tR23  :'), len(good_R23  )
    print ('\tR2   :'), len(good_R2   )
    print ('\tR3   :'), len(good_R3   )
    print ('\tS2   :'), len(good_S2   )
    print ('\tO32  :'), len(good_O32  )
    print ('\tO3  :'), len(good_O3  )
    print ('\tO2  :'), len(good_O2  )
    print ('\tNe3O2:'), len(good_Ne3O2)
    print ('\tO3S2 :'), len(good_O3S2 )
    print ('\tHaHb :'), len(good_HaHb )
    print ('\t>=1 :'), len(gd_any)
    print ('\t>=2 :'), len(gd_2)
    print ('\t>=3 :'), len(gd_3)
    print ('\t>=4 :'), len(gd_4)
    print ('\t>=5 :'), len(gd_5)
    print ('\t>=6 :'), len(gd_6)


    for g, gd in enumerate(good_R23  ): f_R23  .write('%s\t%.5i\n'%(field, cat_data['ID'][gd]))
    for g, gd in enumerate(good_R2   ): f_R2   .write('%s\t%.5i\n'%(field, cat_data['ID'][gd]))
    for g, gd in enumerate(good_R3   ): f_R3   .write('%s\t%.5i\n'%(field, cat_data['ID'][gd]))
    for g, gd in enumerate(good_S2   ): f_S2   .write('%s\t%.5i\n'%(field, cat_data['ID'][gd]))
    for g, gd in enumerate(good_O32  ): f_O32  .write('%s\t%.5i\n'%(field, cat_data['ID'][gd]))
    for g, gd in enumerate(good_O3  ): f_O3  .write('%s\t%.5i\n'%(field, cat_data['ID'][gd]))
    for g, gd in enumerate(good_O2  ): f_O2  .write('%s\t%.5i\n'%(field, cat_data['ID'][gd]))
    for g, gd in enumerate(good_Ne3O2): f_Ne3O2.write('%s\t%.5i\n'%(field, cat_data['ID'][gd]))
    for g, gd in enumerate(good_O3S2 ): f_O3S2 .write('%s\t%.5i\n'%(field, cat_data['ID'][gd]))
    for g, gd in enumerate(good_HaHb ): f_HaHb .write('%s\t%.5i\n'%(field, cat_data['ID'][gd]))
    for g, gd in enumerate(good_any ): f_any .write('%s\t%.5i\n'%(field, cat_data['ID'][gd]))



print ('total')
print ('R23  :'), tot_R23  
print ('R2   :'), tot_R2   
print ('R3   :'), tot_R3   
print ('S2   :'), tot_S2   
print ('O32  :'), tot_O32  
print ('O3  :'), tot_O3  
print ('O2  :'), tot_O2  
print ('Ne3O2:'), tot_Ne3O2
print ('O3S2 :'), tot_O3S2 
print ('HaHb  :'), tot_HaHb

print ('\t>=1:'), tot_any
print ('\t>=2 :'), tot_2
print ('\t>=3 :'), tot_3
print ('\t>=4 :'), tot_4
print ('\t>=5 :'), tot_5
print ('\t>=6 :'), tot_6



f_R23.close()   
f_R2.close()    
f_R3.close()    
f_S2.close()    
f_O32.close()   
f_O3.close()   
f_O2.close()   
f_Ne3O2.close()
f_O3S2.close()  
f_HaHb.close()  
f_any.close()













