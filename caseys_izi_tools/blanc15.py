import os
import numpy as np
import pidly
import matplotlib.pyplot as plt

from scipy.integrate import trapz
from scipy.interpolate import interp1d

from calzetti import k as calk
from hri import hri

def izi(fluxes, errors, lines, idl=pidly.IDL(), dosave=False, savfile='res.sav', 
            grid=os.path.join(os.environ['IZI_DIR'],'grids','l09_high_csf_n1e2_6.0Myr.fits')) :

            #idl = pidly.IDL()
            idl('fluxes = {0}'.format(np.array2string(fluxes, separator=',',max_line_width=1000)))
            idl('errors = {0}'.format(np.array2string(errors, separator=',',max_line_width=1000)))
            idl('lines = {0}'.format(np.array2string(lines, separator=',',max_line_width=1000)))
            idl('forprint, fluxes, errors, lines')
            #print(grid, os.path.isfile(grid))
            #print('gridfile={0})'.format(grid))
            idl('res=izi(fluxes, errors, lines, NZ=100, gridfile="{0}")'.format(grid))
            if dosave :
                idl('save, file="{0}", res'.format(savfile))
            res = idl.ev('res', use_cache=True)
            return(res)
        
def cumstat(x, Px) :
    # given x, and Px, return a function such that f(0.50) = median, f(0.16) is 16th-tile, f(0.84) is 84th-tile, etc.
    # check to make sure Px is normalized so the integral is 1.0
    norm = trapz(Px, x)
    y = Px  / norm

    cum = np.zeros(len(y))
    for i in range(len(y)) : 
        if i==0 : 
            cum[i] = 0
        else : 
            ok = (x[i] > x)
            cum[i] = trapz(y[ok], x[ok])
        #print(i, x[i], cum[i])
    f = interp1d(cum, x)
    return(f)
        
def fit_izi(nZ, doplot=True, plotfile='nZ_IZI.png', plotformat='png', outdir='izi_out',
                doExtinctionCorrection=False) :
    # take nZ and calculate OH and errors using IZI.  Store those in nZ and return it.
    Zmode_izi = np.zeros(len(nZ))
    Zmode = np.zeros(len(nZ))
    Zmedian = np.zeros(len(nZ))
    Zlo68 = np.zeros(len(nZ))
    Zhi68 = np.zeros(len(nZ))
    Zlo68_izi = np.zeros(len(nZ))
    Zhi68_izi = np.zeros(len(nZ))
    npeaks = np.zeros(len(nZ))

    ny = np.ceil( len(nZ) / 4 ).astype(int)
    if doplot :
        fig, ax = plt.subplots(ny,4, figsize=(360,720))
        ax = ax.ravel()

    if not os.path.exists(outdir):
        os.mkdir(outdir)
        
    idl = pidly.IDL()

    # the next three things are for the extinction correction if called upon: 
    Vlam = 5470. # from Johnson Cousins_V 
    calkV = calk(Vlam)
    tlam = np.array([3727., 5007.,  4863.])
    
    for i in range(len(nZ)) : 
        fluxes = nZ[['OII_FLUX','OIII_FLUX','Hb_FLUX']].iloc[i].values
        errors = nZ[['OII_FLUX_ERR','OIII_FLUX_ERR','Hb_FLUX_ERR']].iloc[i].values
        
        lines = np.array(['oii3726;oii3729','oiii4959;oiii5007',  'hbeta'])

        if doExtinctionCorrection:
            # take AV from nZ:
            Av = nZ['Av'].iloc[i]
            
            for l in range(len(tlam)) :
                Alam = 1*Av * calk(tlam[l]) / calkV
                fluxes[l] = fluxes[l] * np.power(10, 0.4*Alam)
                errors[l] = errors[l] * np.power(10, 0.4*Alam)
                #print(Av, Alam, np.power(10,0.4*Alam))

        errors = errors / fluxes[2]
        fluxes = fluxes / fluxes[2]

        savfile = os.path.join(outdir,'izi_id{0}.sav'.format(nZ['ID'].iloc[i]))
        res = izi(fluxes, errors, lines, idl=idl, dosave=True, savfile=savfile,
                      grid=os.environ['IZI_DIR']+'/grids/d13_kappa20.fits')
#                      grid=os.environ['IZI_DIR']+'/grids/k01_SB99_K_csf_n1e1_8.0Myr.fits') # - k01 plots
#                      grid=os.environ['IZI_DIR']+'/grids/k01_SB99_L_csf_n1e1_8.0Myr.fits')
#                      grid=os.environ['IZI_DIR']+'/grids/k13_PN_csf_n3.5e2_4.0Myr.fits')

        if doplot :
            ax[i].plot(res['zarr'][0],res['zpdfmar'][0])
            ax[i].set_xlabel(r'$\log$(O/H) + 12')
            ax[i].set_title('ID={0}'.format(nZ['ID'].iloc[i]))
            ax[i].axvline(res['zgridmarmod'][0], linestyle=':', color='r',label='50%-tile')
            ax[i].axvline(res['zgridmarmod'][0]+res['eupzgridmarmod'][0],
                              linestyle=':', color='r',label='84%-tile')
            ax[i].axvline(res['zgridmarmod'][0]-res['edownzgridmarmod'][0],
                              linestyle='--', color='r',label='16%-tile')
        Zmode[i] = res['zgridmarmod'][0]

        
        Zmedian[i] = res['zgridmarmod'][0]

        #### Old way:
        ####
        #if False : 
        ftile = cumstat( res['zarr'][0], res['zpdfmar'][0])
        Zmode_izi[i] = ftile(0.50)
        Zlo68_izi[i] = ftile(0.16)
        Zhi68_izi[i] = ftile(0.84)

        (tZmod, tZlo, tZhi, tnpeaks) = hri( res['zarr'][0], res['zpdfmar'][0])
        Zmode[i] = tZmod
        Zlo68[i] = tZlo
        Zhi68[i] = tZhi
        npeaks[i] = tnpeaks
            
        #Zlo68[i] = res['zgridmarmod'][0]-res['edownzgridmarmod'][0]
        #Zhi68[i] = res['zgridmarmod'][0]+res['eupzgridmarmod'][0]
        print("IZI finished i={5}, ID={0} with Zmode, Z50, Z16, Z84 = {1:6.3}, {2:6.3}, {3:6.3}, {4:6.3}".format(nZ['ID'].iloc[i], Zmode[i], Zmedian[i], Zlo68[i], Zhi68[i], i))
                #ax[i].axvline(res['zgridmarmean'][0], linestyle='--', color='g', label='mean')
               # ax[i].legend()


    nZ.loc[:,'Zmode'] = Zmode # HRI
    nZ.loc[:,'Zmode_IZI'] = Zmode_izi
    nZ.loc[:,'npeaks'] = npeaks
    nZ['Z16_IZI'] = Zlo68_izi 
    nZ['Z84_IZI'] = Zhi68_izi
    nZ['Z16'] = Zlo68  # HRI
    nZ['Z84'] = Zhi68 # HRI
               
    if doplot :
        plt.savefig(plotfile, format=plotformat)
    idl.close()
    return(nZ)
