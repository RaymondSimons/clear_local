import matplotlib as mpl
mpl.use('TkAgg')
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import glob
from glob import glob
from astropy.io import fits
from astropy.cosmology import Planck15 as cosmo
from numpy import *
import metal_calibs as calib
mpl.rcParams['text.usetex'] = True
plt.ioff()
plt.close('all')




'''
import argparse

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--ratio', default = 'R23')
args = vars(parser.parse_args())

rat = args['ratio']
'''


gds_cat = fits.open('/Volumes/pegasus/clear/catalogs/goods%s_3dhst.v4.4.cats/Eazy/goods%s_3dhst.v4.4.zout.fits'%('s', 's'))
gdn_cat = fits.open('/Volumes/pegasus/clear/catalogs/goods%s_3dhst.v4.4.cats/Eazy/goods%s_3dhst.v4.4.zout.fits'%('n', 'n'))

grizli_cat = fits.open('/Volumes/pegasus/clear/catalogs/grizli_v2.1_cats/master_cat.fits')

zmn = 10
zmx = 0
zs = []
with PdfPages('/Volumes/pegasus/clear/line_ratio_profiles/lineratio_profiles.pdf') as pdf:
    fig, axes = plt.subplots(2, 2, figsize = (12, 8))
    for r, rat in enumerate(['O32', 'R2', 'R3', 'R23']):

        if rat == 'R23':   cal =  calib.OH_R23
        if rat == 'R2':    cal =  calib.OH_R2
        if rat == 'R3':    cal =  calib.OH_R3
        if rat == 'O32':   cal =  calib.OH_O32




        fls = glob('/Volumes/pegasus/clear/line_ratio_profiles/%s/*_%s_profile.npy'%(rat, rat))
        ax = axes.ravel()[r]

        ms = []
        m_rs = []
        em_rl = []
        em_ru = []
        for f, fl in enumerate(fls):
            field = fl.split('/')[-1].split('_')[0]
            di = fl.split('/')[-1].split('_')[1]


            if (field == 'ERSPRIME') | ('GS' in field): cat = gds_cat
            if ('GN' in field): cat = gdn_cat

            gd = np.where(cat[1].data['id'] == int(di))[0][0]

            mass = cat[1].data['mass'][gd]
            gd_grizli = where(grizli_cat[field].data['id'] == int(di))[0][0]
            z = grizli_cat[field].data['z_50'][gd_grizli]

            zs.append(z)
            a = np.load(fl, allow_pickle = True)[()]
            if (len(a['r']) < 50) & (len(a['r']) > 10):# & (max(array(a['r']).max()) < 1.):
                zmn = min([z, zmn])
                zmx = max([z, zmx])


                ms.append(mass)
                arcsec_to_kpc = cosmo.arcsec_per_kpc_proper(z).value
                m_rs.append(a['m_mcmc'][0] * arcsec_to_kpc)
                em_rl.append(a['m_mcmc'][2] * arcsec_to_kpc)
                em_ru.append(a['m_mcmc'][1] * arcsec_to_kpc)

        ms  = array(ms )
        m_rs  = array(m_rs )
        em_rl = array(em_rl)
        em_ru = array(em_ru)

        gd = where(em_rl < 0.2)
        ax.errorbar(log10(ms[gd]), m_rs[gd], yerr = [em_rl[gd], em_ru[gd]], fmt = 'o', linewidth = 0.3, alpha = 0.5, markersize = 3, color = 'black')

        ax.set_ylabel(r'$\Delta$ %s/$\Delta$ r (dex kpc$^{-1}$)'%rat, fontsize = 15)
        ax.set_xlabel(r'$\log$ M$_*$ (M$_{\odot}$)', fontsize = 15)
        ax.axhline(y = 0.0, alpha = 0.2, zorder = 20, color = 'black', linestyle = '--')
        ax.set_xlim(8, 12)
        '''
        ax2 = ax.twinx()
        dOH_dr_ticks = arange(-3., 5., 0.5)
        dOH_dr_ticks_str = array(['%.1f'%oh for oh in dOH_dr_ticks])
        
        dR_dr_ticks = array([cal(oh) for oh in dOH_dr_ticks])
        ax2.set_yticks(dR_dr_ticks)
        ax2.set_yticklabels(dOH_dr_ticks_str)
        ax2.set_ylabel(r'$\Delta$ $\log$ OH/$\Delta$ r (dex kpc$^{-1}$)', fontsize = 15)
        '''
        ax.set_ylim(-0.50, 0.4)



    fig.tight_layout()
    pdf.savefig(fig)



