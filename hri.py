import numpy as np
from scipy.integrate import trapz
from scipy.io import readsav
import matplotlib.pyplot as plt

def findpeaks(d1z, d2z, value=0.0, d2threshold=1e-2):
    
    # take the 1st and 2nd derivative and return the peaks... 
    
    diffseq = d1z - value
    signseq = np.sign(diffseq)
    zero_crossings = signseq[0:-2] != signseq[1:-1]
    indices = np.where((zero_crossings) & ((d2z[0:-2] > d2threshold) |(d2z[1:-1] > d2threshold)))[0]
    for i, v in enumerate(indices):
        #print(i,v,d1z[v], d2z[v])
        if ((abs(d1z[v + 1] - value) < abs(d1z[v] - value)) ):
            indices[i] = v + 1
    return indices

# now write a routine to measure the HDI
def hri_arrays(x, opx, pcut=0.68) : 
    ## make sure it is normalized: 
    norm = trapz(opx, x)
    px = opx / norm
    
    s = np.argsort( px)[::-1]
    sx = x[s]
    spx = px[s]
    
    intP = np.zeros(len(px))
    for i in range(len(px)) : 
        if i > 0 : 
            t = (px > spx[i])
            tpx = px
            tpx[t] = 0.0
            intP[i] = 1-  trapz(tpx, x)
    return(sx, spx, intP)

# now write a routine to measure the HDI
def hri(x, opx, pcut=0.6827, checkPeaks=True, doplot=False) : 

    #
    # Returns Zmode, Zlo, Zhi, npeaks
    #
    # checkPeaks means the code will test if Zlo and Zhi contain multiple peaks.  
    # If so, it will subtract the number of peaks (minus 1)  from npeaks
    #
    ## make sure it is normalized: 
    norm = trapz(opx, x)
    px = opx / norm
    
    sx, spx, intP = hri_arrays(x, opx, pcut=pcut)
    
    ## now, take derivatives and find peaks: 
    t = sx.argsort()
    intP = intP[t]
    d1z = (np.gradient(intP,x))
    d2z = np.gradient(np.gradient(intP,x),x)
    d1z = d1z / d1z.max()
    d2z = d2z / d2z.max()

    ans = findpeaks(d1z,d2z)
    #print("ans: ", ans)
 
    # number of peaks is here
    npeaks= len(ans)
    
    # now, take the highest peak and move down until you get to pcut: 
    
    # do lower bound: 
    print ('npeaks', npeaks)
    print ('intp', intP)
    print ('ans', ans)
    mni = np.argmin(intP[ans])
    i = ans[mni]
    Zmode=x[i]
    #print(i,z[i])
    done=False
    #print("Zlo")
    while not done :
        Zlo=x[i]
        i = i-1
        #print(i,z[i], Zlo, intP[i], len(intP), pcut)

        if ((i < 0) | (intP[i] > pcut)) : 
            done=True

    # do upper bound
    i = ans[mni]
    done=False
    while not done : 
        Zhi = x[i]
        i = i+1
        #print(i,z[i], Zhi, intP[i], len(intP), pcut)
        if ((i >= len(intP)) | (intP[i] > pcut)) : 
            done=True
        
    # check if the Zlo Zhi contains multiple peaks
    if checkPeaks : 
        count=0
        #print("starting: npeaks, count= ", npeaks, count)
        for a in ans : 
            if ((Zlo <= x[a]) & (Zhi >= x[a])) : 
                count=count+1
        npeaks = npeaks - (count-1) # count should be 1 if only 1 peak falls between Zlo and Zhi
        #print("ending: npeaks, count= ", npeaks, count)


    if doplot : 
        plt.plot(x,intP)
        plt.plot(x,d1z/d1z.max())
        plt.plot(x, d2z/d2z.max())
        plt.axvline(Zlo,ls=':')
        plt.axvline(Zhi,ls=':')

    return(Zmode, Zlo, Zhi, npeaks)
