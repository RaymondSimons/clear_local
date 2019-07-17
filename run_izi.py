import pidly
from calzetti import k as calk



idl = pidly.IDL()
Vlam = 5470. # from Johnson Cousins_V 
calkV = calk(Vlam)
tlam = np.array([3727., 5007.,  4863.])

# take AV from nZ:
Av = nZ['Av'].iloc[i]

for l in range(len(tlam)) :
    Alam = 1*Av * calk(tlam[l]) / calkV
    fluxes[l] = fluxes[l] * np.power(10, 0.4*Alam)
    errors[l] = errors[l] * np.power(10, 0.4*Alam)
    #print(Av, Alam, np.power(10,0.4*Alam))


errors = errors / fluxes[2]
fluxes = fluxes / fluxes[2]




