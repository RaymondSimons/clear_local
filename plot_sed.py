import astropy
from astropy.io import fits
import grizli
from grizli.pipeline import photoz
from grizli import utils
import eazy
import os


cat_dir = '/Users/rsimons/Desktop/clear/Catalogs'
translate_file = '/Users/rsimons/Desktop/clear/Catalogs/goodsn_3dhst.v4.1.translate'

field = 'n'

params = {}
params['CATALOG_FILE'] = cat_dir + '/goods%s_3dhst.v4.1.cats/Catalog/goods%s_3dhst.v4.1.cat'%(field, field)
params['Z_STEP'] = 0.002
params['Z_MAX'] = 4
params['MAIN_OUTPUT_FILE'] = '{0}_3dhst.{1}.eazypy'.format('goodsn', 'v4.1')
params['PRIOR_FILTER'] = 205
params['MW_EBV'] = {'aegis':0.0066, 'cosmos':0.0148, 'goodss':0.0069, 
                    'uds':0.0195, 'goodsn':0.0103}['goodsn']
params['TEMPLATES_FILE'] = 'templates/fsps_full/tweak_fsps_QSF_12_v3.param'


eazy.symlink_eazy_inputs(path=os.path.dirname(eazy.__file__)+'/data', path_is_env=False)
templ0 = grizli.utils.load_templates(fwhm=1200, line_complexes=True, stars=False, 
                                     full_line_list=None,  continuum_list=None, 
                                     fsps_templates=True)


ez = eazy.photoz.PhotoZ(param_file=None, translate_file=translate_file, 
                        zeropoint_file=None, params=params, 
                        load_prior=True, load_products=False)


ep = photoz.EazyPhot(ez, grizli_templates=templ0, zgrid=ez.zgrid)

phot, ii, dd = ep.get_phot_dict(189.165441, 62.257319)




'''
figure(1)
clf()



cat_dir = '/Users/rsimons/Desktop/clear/Catalogs'
field = 'n'
di = 23186
if False:
    v41    = np.loadtxt(cat_dir + '/goods%s_3dhst.v4.1.cats/Catalog/goods%s_3dhst.v4.1.cat'%(field, field))
    v43    = np.loadtxt(cat_dir + '/goods%s_3dhst.v4.3.cat'%(field))
    v43nzp = np.loadtxt(cat_dir + '/goods%s_3dhst.v4.3.nzpcat'%(field))
    v41f   = fits.open(cat_dir + '/goods%s_3dhst.v4.1.cats/Catalog/goods%s_3dhst.v4.1.cat.FITS'%(field, field))




gd = where(v41f[1].data['id'] == di)[0][0]


filts = array([('U', 0.3593), ('B', 0.4448),  ('F435W', 0.4318), ('G', 0.4751), 
               ('V', 0.5470), ('F606W', 0.5919), ('Rs', 0.6819), ('R', 0.6276), 
               ('I', 0.7671), ('F775W', 0.7693), ('Z    ', 0.9028), ('F850LP', 0.9036), 
               ('F125W', 1.2471), ('J', 1.2517), ('H', 1.6347), ('F140W', 1.3924), 
               ('F160W', 1.5396), ('Ks', 2.1577), ('IRAC1', 3.5569), ('IRAC2', 4.5020), 
               ('IRAC3', 5.7450), ('IRAC4', 7.9158)])



for f, filt in enumerate(filts):
    print v41f[1].data['f_' + filt[0].lower()][gd]
    plot(float(filt[1]), v41f[1].data['f_' + filt[0].lower()][gd],'k.')
'''