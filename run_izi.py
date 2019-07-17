import pidly
from calzetti import k as calk
import numpy as np

idl_path = '/grp/software/Linux/itt/idl/idl84/idl/bin/idl'
idl = pidly.IDL(idl_path)
Vlam = 5470. # from Johnson Cousins_V 
calkV = calk(Vlam)
tlam = np.array([3727., 5007.,  4863.])

lines = np.array(['oii3726;oii3729','oiii4959;oiii5007',  'hbeta'])
fluxes = np.array([0.4e-17, 0.8e-17, 1.e-17])
errors = np.array([0.1e-17, 0.2e-17, 1.e-18])

savfile = 'test.sav'

res = izi(fluxes, errors, lines, idl=idl, dosave=True, savfile=savfile,
              grid=os.environ['IZI_DIR']+'/grids/d13_kappa20.fits')

# take AV from nZ:
'''
Av = nZ['Av'].iloc[i]

for l in range(len(tlam)) :
    Alam = 1*Av * calk(tlam[l]) / calkV
    fluxes[l] = fluxes[l] * np.power(10, 0.4*Alam)
    errors[l] = errors[l] * np.power(10, 0.4*Alam)
    #print(Av, Alam, np.power(10,0.4*Alam))


errors = errors / fluxes[2]
fluxes = fluxes / fluxes[2]
'''



