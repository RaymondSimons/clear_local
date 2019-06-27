
import numpy as np

def k(lam, gas=False, RV=4.045) :
    # returns the calzetti k(l) for a given lam in angstroms 
    # ref: Calzetti et al 2000
    # 
    # returns k(lambda) which gives the extinction (in mag) given an
    # E(B-V), i.e., A(lambda) = k(lambda) * E(B-V).

    # THIS IS THE EXTINCTION ON THE STELLAR CONTINUUM.
    # TO GET NEBULAR EXTINCTION, MULTIPLY BY 1.0/0.44

    # USAGE:
    try :
        lam
    except Exception :
        USAGE = 'from calzetti import k; X = calzetti.k(WAVELENGTH, gas=False)'
    #

    LL = lam / 1e4

    if np.size(np.shape(LL))==0: # then single value
        if (LL > 0.63) :
            cal=( 2.659*((-1.857) + 1.040/LL ) + RV)
        else :
            cal=(2.659*( -2.156 + 1.509/LL - 0.198/LL**2 + 0.011/LL**3 ) + RV)
    else :
        ok  = (LL > 0.63)
        cok = (LL<=0.63)
        cal = np.zeros(len(LL))
        cal[ok] = ( 2.659*((-1.857) + 1.040/LL[ok] ) + RV)
        cal[cok] = 2.659*( -2.156 + 1.509/LL[cok]- 0.198/LL[cok]**2 + 0.011/LL[cok]**3 ) + RV

    if gas :
        cal = cal * 0.44

    return(cal)

    
