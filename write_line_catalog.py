import astropy
from astropy.io import fits
import glob
from glob import glob
import numpy as np
from numpy import *
import joblib
from joblib import Parallel, delayed




fields = ['GS1','GS2', 'GS3', 'GS4', 'GS5', 'GN1', 'GN2', 'GN3', 'GN4', 'GN5', 'GN7', 'ERSPRIME']


def write_catalog(field):
    print (field)
    #fls = glob('/Volumes/pegasus/clear/grizli_extractions/%s/*/Prep/*%s*full.fits'%(field, field))
    #fits_name = '/Volumes/pegasus/clear/grizli_extractions/Catalogs/grizli_v2.1_cats/%s_lines_grizli.fits'%field
    fls = glob('/Users/rsimons/Desktop/clear/grizli_extractions_v3/%s/Prep/*%s*full.fits'%(field, field))
    fits_name = '/Users/rsimons/Desktop/clear/grizli_v3.0_cats/%s_lines_grizli.fits'%field


    #lines = ['Lya', 'CIV', 'MgII', 'OII',
    #         'Hd', 'Hg', 'OIIIx', 'HeII',
    #         'Hb', 'OIII', 'Ha', 'SII',
    #         'SIII', 'HeI', 'HeIb', 'NeIII',
    #         'NeV', 'NeVI', 'OI']

    lines =['ArIII-7138', 'CIII-1908', 'CIV-1549', 'H8', 'H9', 'Ha', 'Hb',
           'Hd', 'HeI-1083', 'HeI-5877', 'HeII-1640', 'Hg', 'Lya', 'MgII',
           'NIII-1750', 'NIV-1487', 'NV-1240', 'NeIII-3867', 'NeV-3346',
           'NeVI-3426', 'OI-6302', 'OII', 'OII-7325', 'OIII', 'OIII-1663',
           'OIII-4363', 'PaB', 'SII', 'SIII']


    fluxs = zeros((len(lines),2, len(fls))) - 99.
    exptime = zeros((2, len(fls)))

    zs = zeros((7, len(fls))) - 99.
    IDs = []
    ras = []
    decs = []
    nlines = []
    chimins = []
    dofs = []
    bics = []
    for f, fl in enumerate(fls):
        a = fits.open(fl)
        IDs.append(int(a[0].header['ID']))
        ras.append(a[0].header['ra'])
        decs.append(a[0].header['dec'])
        nlines.append(a[0].header['NUMLINES'])
        chimins.append(a[1].header['CHIMIN'])
        dofs.append(a[1].header['DOF'])
        bics.append(a[1].header['BIC_TEMP'])
        


        try:
            exptime[0,f] = float(a[0].header['T_G102'])
        except:
            exptime[0,f] = 0.0

        try:
            exptime[1,f] = float(a[0].header['T_G141'])
        except:
            exptime[1,f] = 0.0


        lines_f = a[0].header['HASLINES'].split(' ')
        if (len(lines_f) > 0) & (lines_f[0] !=''):
            zs[0, f] = a[1].header['Z50']
            zs[1, f] = a[1].header['Z02']
            zs[2, f] = a[1].header['Z16']
            zs[3, f] = a[1].header['Z84']
            zs[4, f] = a[1].header['Z97']
            zs[5, f] = a[1].header['Z_MAP']
            zs[6, f] = a[1].header['Z_RISK']


            for l, ln in enumerate(lines_f):
                flux_ln =  a[0].header['FLUX%.3i'%(l+1)]
                eflux_ln =  a[0].header['ERR%.3i'%(l+1)]
                ln_name =  a[0].header['LINE%.3i'%(l+1)]
                for ll, line in enumerate(lines):	
                    if ln_name == line:
                        fluxs[ll,0,f] = flux_ln * 1.e17
                        fluxs[ll,1,f] = eflux_ln  * 1.e17

    master_hdulist = []
    prihdr = fits.Header()

    prihdu = fits.PrimaryHDU(header=prihdr)    
    master_hdulist.append(prihdu)

    colhdr = fits.Header()

    col_list = [fits.Column(name='ID', format = 'J', array=array(IDs)),
    fits.Column(name='RA', format = 'D', array=array(ras)),
    fits.Column(name='DEC',format = 'D', array=array(decs)),
    fits.Column(name='nlines',format = 'J', array=array(nlines).astype('int')),

    fits.Column(name='z_50',format = 'D', array=zs[0,:]),
    fits.Column(name='z_02',format = 'D', array=zs[1,:]),
    fits.Column(name='z_16',format = 'D', array=zs[2,:]),
    fits.Column(name='z_84',format = 'D', array=zs[3,:]),
    fits.Column(name='z_97',format = 'D', array=zs[4,:]),
    fits.Column(name='z_MAP',format = 'D', array=zs[5,:]),
    fits.Column(name='z_RISK',format = 'D', array=zs[6,:])]




    for ll, line in enumerate(lines):
        col_list.append(fits.Column(name='%s_FLUX'%line,format = 'D', array=fluxs[ll,0,:]))
        col_list.append(fits.Column(name='%s_FLUX_ERR'%line,format = 'D', array=fluxs[ll,1,:]))
    col_list.append(fits.Column(name='T_G102',format = 'D', array=exptime[0,:]))
    col_list.append(fits.Column(name='T_G141',format = 'D', array=exptime[1,:]))
    col_list.append(fits.Column(name='BIC_TEMP',format = 'D', array=array(bics)))
    col_list.append(fits.Column(name='CHIMIN',format = 'D', array=array(chimins)))
    col_list.append(fits.Column(name='DOF',format = 'D', array=array(dofs)))        

    coldefs = fits.ColDefs(col_list)
    table_hdu = fits.BinTableHDU.from_columns(coldefs)
    master_hdulist.append(table_hdu)
    thdulist = fits.HDUList(master_hdulist)
    thdulist.writeto(fits_name, overwrite = True)



def make_master_cat():
    for f in ['S', 'N']:

        cat_dir = '/Users/rsimons/Desktop/clear/grizli_v3.0_cats'
        cat_names = glob(cat_dir + '/*%s*_lines_grizli.fits'%f)


        for c, cat_name in enumerate(cat_names):

            field = cat_name.split('/')[-1].split('_')[0]

            cat = fits.open(cat_name)

            if c == 0:
                print ('Creating empty dictionary for GD%s'%f)
                cat_all = {}

                formats = cat[1].columns.formats
                names = cat[1].columns.names
                for n in arange(len(formats)):
                    name = names[n]
                    form = formats[n]

                    cat_all[name] = {}
                    cat_all[name]['data'] = []
                    cat_all[name]['format'] = form



            for i in arange(len(cat[1].data)):
                if cat[1].data['ID'][i] not in cat_all['ID']['data']:
                    # this is the first recording of this object                
                    for name in names:
                        if (name == 'ID') | (name == 'DEC') | (name == 'RA'):
                            cat_all[name]['data'].append(cat[1].data[name][i])
                        else:
                            cat_all[name]['data'].append([cat[1].data[name][i]])

                else:
                    gd = where(cat[1].data['ID'][i]  == cat_all['ID']['data'])[0][0]
                    for name in names:

                        if (name == 'ID') | (name == 'DEC') | (name == 'RA'): continue
                        else:  cat_all[name]['data'][gd].append(cat[1].data[name][i])

        for i in arange(len(cat_all['ID']['data'])):
            t_g102 = array(cat_all['T_G102']['data'][i])
            t_g141 = array(cat_all['T_G141']['data'][i])

            if len(t_g102) > 1:
                where_both = where((t_g102 > 0) & (t_g141 > 0))[0]
                where_g102 = where(t_g102 > 0)[0]
                where_g141 = where(t_g141 > 0)[0]


                if len(where_both) > 0:
                    index_use = where_both[np.argmax(t_g102[where_both] + t_g141[where_both])]
                    #print ('\t multiple fits, at least one with both G102 and G141, using index %i'%index_use)
                elif len(where_g102) > 0:
                    index_use = where_g102[np.argmax(t_g102[where_g102])]
                    #print ('\t multiple fits, none with joint G102+G141, using max T G102, using index %i'%index_use)

                elif len(where_g141) > 0:
                    index_use = where_g141[np.argmax(t_g141[where_g141])]
                    #print ('\t multiple fits, none with G102, using max T G141, using index %i'%index_use)

            else: 
                index_use = 0
                #print ('\t single fit, using index %i'%index_use)

            for name in names:    
                if (name == 'ID') | (name == 'DEC') | (name == 'RA'): continue
                cat_all[name]['data'][i] = cat_all[name]['data'][i][index_use]
                
                
        print ('\t writing')

        
        master_hdulist = []
        prihdr = fits.Header()

        prihdu = fits.PrimaryHDU(header=prihdr)    
        master_hdulist.append(prihdu)

        colhdr = fits.Header()

        col_list = []
        for name in names:
            col_list.append(fits.Column(name=name, format = cat_all[name]['format'], array=array(cat_all[name]['data'])))
        coldefs = fits.ColDefs(col_list)
        table_hdu = fits.BinTableHDU.from_columns(coldefs)
        master_hdulist.append(table_hdu)
        thdulist = fits.HDUList(master_hdulist)
        thdulist.writeto(cat_dir + '/GD%s_lines_grizli_master.fits'%f, overwrite = True)
        
        print ('\t done')

if __name__ == '__main__':
    Parallel(n_jobs = -1, backend = 'threading')(delayed(write_catalog) (field = field) for field in fields)
    make_master_cat()



