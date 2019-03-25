import os
import glob



line_sets = ['R23',
             'R2',
             'R3',
             'S3',
             'O32',
             'O3S2',
             'HaHb']



sample_fulls_dir = '/user/rsimons/grizli_extractions/sample_fulls'



cat_dir = '/user/rsimons/grizli_extractions/Catalogs'
for line in line_sets:
    cat = np.loadtxt(cat_dir + '/sample_cats/%s_sample.cat'%line, dtype = 'str')
    for c in cat:
        fl = glob('/user/rsimons/grizli_extractions/%s/*/Prep/%s_%s.full.fits'%(c[0], c[0], c[1]))
        print fl
