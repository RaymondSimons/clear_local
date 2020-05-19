;########################################################################
;
; Copyright (C) 2014, Guillermo A. Blanc
; E-mail: gblancm@obs.carnegiescience.edu
;
; Updated versions of the software are available from my webpage
; http://users.obs.carnegiescience.edu/gblancm/izi
;
; Please reference Blanc et al. 2014 (in preparation) when publishing
; results that make use of this software.
;
; This software is provided as is without any warranty whatsoever.
; Permission to use, for non-commercial purposes is granted.
; Permission to modify for personal or internal use is granted,
; provided this copyright and disclaimer are included unchanged
; at the beginning of the file. All other rights are reserved.
;
; #######################################################################
;
;+
; NAME:
;
;   IZI 
;
; PURPOSE:
;
;   Compute the posterior PDF of the gas-phase metallicity and the
;   ionization parameter given a set of observed emission line fluxes
;   and errors, and a photo-ionization model grid. The code
;   interpolates the models to a finer grid and evaluates the
;   posterior PDF of the parameters given the data and the models.
;   
; CALLING SEQUENCE:
;
;  result = IZI(fluxin, errorin, idin,                                $
;  GRIDFILE=gridfile, PLOT=plot,                                      $
;  LOGOHSUN=logOHsun, EPSILON=epsilon, NZ=nz, NQ=nq,                  $
;  INTERGRIDFILE=intergridfile, OUTGRIDFILE=outgridfile,              $
;  LOGZLIMITS=logzlimits, LOGQLIMITS=logqlimits, LOGZPRIOR=logzprior, $
;  LOGQPRIOR=logqprior, NONORM=nonorm)
;
;
; INPUT PARAMETERS:
;
;    FLUXIN: float/double array of emission line fluxes in arbitrary units. Fluxes
;    are nornalized internally in the code to the flux of H-beta, or
;    in the case H-beta is not provided to the flux of the brightest
;    input emission line. Fluxes must be corrected for dust extinction.
;
;    ERRORIN: float/double array of emission line flux errors. Upper limits can be
;    provided by setting the a value of -666 in the FLUXIN array and
;    providing the 1 sigma upper limit in the ERRORIN parameter.
;
;    IDIN: array of strings containing the emission line IDs. The
;    names of the lines must comply with the IDs stored in the FITS
;    table storing the photo-ionization model grid. For all the grids that
;    are provided with IZI the following are some (not all) of the
;    adopted IDs:
;
;          'oii3726'    for [OII]-3726
;          'oii3729'    for [OII]-3729
;          'hgamma'     for H_gamma-4341
;          'oiii4959'   for [OIII]-4959
;          'hbeta'      for H_beta-4861
;          'oiii5007'   for [OIII]-5007
;          'oi6300'     for [OI]-6300
;          'nii6548'    for [NII]-6548
;          'halpha'     for H_alpha-6563
;          'nii6584'    for [NII]-6584
;          'sii6717'    for [SII]-6717
;          'sii6731'    for [SII]-6731 
;          'siii9068'   for [SII]-9068
;          'siii9532'   for [SII]-9532 
;
;
;    For a complete list of lines included in a particular grid read
;    the photo-ionization model grid into an IDL structure and check
;    the ID tag. 
;
;    Summed doublets (e.g. [OII]3726,9 at low resolution) can be
;    considered as a single line by providing the IDs of the
;    components separated by a semicolon in IDIN (e.g. 'oii3726;oii3729').
;
;        
; KEYWORD PARAMETERS:
;
;    GRIDFILE: string containing the filename for the photo-ionization
;    model grid to be used. If not provided IZI defaults to the
;    Levesque et al. 2010 models with high mass loss, constant star
;    formation history, n=1e2 cm^-3 and age of 6.0Myr. The format of
;    the grid files is a FITS table read as an IDL structure with an
;    entry for each (Z,q) element in the grid and the following tags for
;    each entry:
;      
;          NAME            STRING    name of the photo-ionization model
;          LOGOHSUN        FLOAT     assumed 12+log(O/H) solar abundance in the model 
;          LOGZ            FLOAT     log of metallicity (Z) in solar units
;          LOGQ            FLOAT     log of ion. parameter in log(cm/s) units
;          ID              STRING[Nlines]    emission line IDs
;          FLUX            FLOAT[Nlines]     emission line fluxes
;
;    /PLOT: set this keyword to produce plots (ps files) of the
;    parameters PDFs and of diagnostic line ratios versus log(O/H) and q.
;   
;    LOGOHSUN: set this keyword to a user suplied value which is used
;    instead of the LOGOHSUN value in the model grid file.
;
;    EPSILON: systematic uncertainty in dex for the emission line
;    fluxes in the model (see equation 3 in Blanc et al. 2014). If not
;    provided the code assumes a default value is 0.15 dex.
;
;    NZ: number of log(Z) elements in the interpolated grid, if not
;    provided the default is NZ=50
;   
;    NQ: number of log(q) elements in the interpolated grid, if not
;    provided the default is NZ=50
;
;    INTERGRIDFILE: string containing the filename for an already
;    interpolated photo-ionization model grid. If provided this file
;    is used intead of GRIDFILE and the interpolation step is
;    skiped. This interpolated grid file should be created using the
;    OUTGRIDFILE keyword. Useful for speeding computation time when
;    running IZI for a large sample of objects.
;
;    OUTGRIDFILE: string containing a filename to save the
;    interpolated grid file for latter use with INTERGRIDFILE
;
;    LOGZLIMITS: 2 element array conatining the lower and upper
;    log(Z) limits for the interpolated grid (section 3 of Blanc et
;    al. 2014) in units of 12+log(O/H)
;
;    LOGQLIMITS: 2 element array conatining the lower and upper
;    log(Z) limits for the interpolated grid (section 3 of Blanc et
;    al. 2014) in units of log(cm/s)
;
;    LOGZPRIOR: (Np x 2) element array especifying a prior for log(Z)
;    in 12+log(O/H) units. The first array [Np,0] contains Np
;    values for the metallicity and the second array [Np,1] contains
;    the probability for each value. This array is interpolated to the
;    NZ grid.
;
;    LOGQPRIOR: (Np x 2) element array especifying a prior for log(q)
;    in log(cm/s) units. The first array [Np,0] contains Np
;    values for the ionization parameter and the second array [Np,1] contains
;    the probability for each value. This array is interpolated to the
;    NQ grid.
;
;    /NONORM: set this keyword to avoid the normalization of the line
;    fluxes. This is useful when using model grids of line ratios
;    instead of line fluxes.
;
; 
; OUTPUT:
;
;    RESULT: The function outputs an IDL structure containing the following tags:
;
;      ID: names of the emission lines in the photo-ionization model
;      FLUX: nomalized input line fluxes 
;      ERROR: error in normalized line fluxes
;      CHI2: chi^2 between observed fluxes and model fluxes at mode of the joint PDF
;      ZGRID: joint mode best-fit metallicity in units of 12+log(O/H)
;      EUPZGRID: upper 1 sigma error in ZGRID
;      EDOWNZGRID: lower 1 sigma error in ZGRID
;      QGRID: joint mode best-fit ionization parameter in units of log(cm/s)
;      EUPQGRID: upper 1 sigma error in QGRID
;      EDOWNQGRID: lower 1 sigma error in QGRID
;      ZGRIDMARMOD: marginalized mode best-fit metallicity in units of 12+log(O/H)
;      EUPZGRIDMARMOD: upper 1 sigma error in ZGRIDMARMOD
;      EDOWNZGRIDMARMOD: lower 1 sigma error in ZGRIDMARMOD
;      QGRIDMARMOD: marginalized mode best-fit ionization parameter in units of log(cm/s)
;      EUPQGRIDMARMOD: upper 1 sigma error in QGRIDMARMOD
;      EDOWNQGRIDMARMOD: lower 1 sigma error in QGRIDMARMOD
;      ZGRIDMARMEAN: marginalized mean best-fit metallicity in units of 12+log(O/H)
;      EUPZGRIDMARMEAN: upper 1 sigma error in ZGRIDMARMEAN
;      EDOWNZGRIDMARMEAN: lower 1 sigma error in ZGRIDMARMEAN
;      QGRIDMARMEAN: marginalized mean best-fit ionization parameter in units of log(cm/s)
;      EUPQGRIDMARMEAN: upper 1 sigma error in QGRIDMARMEAN
;      EDOWNQGRIDMARMEAN: lower 1 sigma error in QGRIDMARMEAN
;      ZARR: array of interpolated metallicity values in units of 12+log(O/H)
;      ZPDFMAR: marginalized metallicity PDF as a function of ZARR 
;      QARR: array of interpolated ionization parameter values in units of log(cm/s)
;      QPDFMAR: marginalized metallicity PDF as a function of QARR 
;      FLAGS: 4 element array with flags [znpeaks, qnpeaks, zlimit, qlimit] with:
;             (znpeaks, qnpeaks): number of peaks in log(Z) and log(q) marginalized PDFs
;             (zlimit, qlimit): flags stating if the marginalized PDF is
;             bound (1), or is either an upper (2) or lowe (3) limit.
;
; MODIFICATION HISTORY: 
;
;    V1.0 - Created by G.A. Blanc
;
; =====================================================================================


FUNCTION uprior, xrange
; calculates a uniform prior for x
; recieves as input the range of allowed x values (2 element array)
  return, 1d/(xrange[1]-xrange[0])
END

; FUNCTION jprior, x, xrange
;; Jefferys Prior for x  - OBSOLETE BECAUSE ALL CALCULATIONS DONE IN LOG(Z) AND LOG(Q)
;; recieves as parameters the value of x and the allowed range of x (2 element array)
;  return, 1d/(x*alog(xrange[1]/xrange[0]))
;END

FUNCTION userprior, x, xarr, yarr
; interpolates a user provided prior (xarr, yarr) to x
; returns 0 if x is outside range of xarr
  if (x lt min(xarr) or x gt max(xarr)) then return, 0
  return, interpol(yarr, xarr, x, /NAN)
END

FUNCTION plotratioq, ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, range=range, title=title
; plot a line ratio vs the ionization parameter

  if not (keyword_set(range)) then range=[-2.0, 2.0]
  if not (keyword_set(title)) then title='log(ratio)'
  
  auxzarr=grid[uniq(grid.logz, sort(grid.logz))].logz
  auxzarr0=grid0[uniq(grid0.logz, sort(grid0.logz))].logz

  sel1=where(abs(grid.logz+logOHsun-d.zgrid) eq min(abs(grid.logz+logOHsun-d.zgrid))) 
  sel2=where(abs(grid.logz+logOHsun-d.zgrid-d.eupzgrid) eq min(abs(grid.logz+logOHsun-d.zgrid-d.eupzgrid)))
  sel3=where(abs(grid.logz+logOHsun-d.zgrid+d.edownzgrid) eq min(abs(grid.logz+logOHsun-d.zgrid+d.edownzgrid)))
  sel1=sel1[sort(grid[sel1].logq)]
  sel2=sel2[sort(grid[sel2].logq)]
  sel3=sel3[sort(grid[sel3].logq)]
                               ; plot
  gridratio=ga/gb
  gridratio0=ga0/gb0
 
  plot, grid[sel1].logq, alog10(gridratio[sel1]), yrange=range, xrange=[min(grid.logq)-0.1, max(grid.logq)+0.1], xthick=3, ythick=3, charthick=3, charsize=3, xtitle=textoidl('log(q)'), ytitle=title, thick=3, linestyle=2, xstyle=1, ystyle=1
  plotsym, 0, 1.2, /fill
;  for i=0, n_elements(auxzarr)-1 do begin
;     sel=where(grid.logz eq auxzarr[i]) 
;     sel=sel[sort(grid[sel].logq)]
;     oplot, grid[sel].logq, alog10(gridratio[sel]), color=cgcolor('gray');, psym=3
;  endfor
  for i=0, n_elements(auxzarr0)-1 do begin
     sel=where(grid0.logz eq auxzarr0[i]) 
     sel=sel[sort(grid0[sel].logq)]
     oplot, grid0[sel].logq, alog10(gridratio0[sel]), color=cgcolor('gray'), thick=3;, psym=3
  endfor
  oplot, grid[sel1].logq, alog10(gridratio[sel1]), thick=4, color=cgcolor('red')
  oplot, grid[sel2].logq, alog10(gridratio[sel2]), thick=4, color=cgcolor('red'), linestyle=2
  oplot, grid[sel3].logq, alog10(gridratio[sel3]), thick=4, color=cgcolor('red'), linestyle=2
  if (flaga eq 1 and flagb eq 1) then begin
     dratio=da/db
     edratio=sqrt((eda/db)^2+(da/db^2*edb)^2)/dratio/alog(10.)
     oploterror, [d.qgrid], [alog10(dratio)], [d.eupqgrid], [edratio], psym=3, thick=4, /hibar , color=cgcolor('black')
     oploterror, [d.qgrid], [alog10(dratio)], [d.edownqgrid], [edratio], psym=3, thick=4, /lobar , color=cgcolor('black')
  endif else if (flaga eq 2 and flagb eq 1) then begin
     dratio=eda/db
     edratio=eda/(db-edb)-eda/db
     oploterror, [d.qgrid], [alog10(dratio)], [d.eupqgrid], [edratio], psym=8, thick=4, /hibar , color=cgcolor('black')
     oploterror, [d.qgrid], [alog10(dratio)], [d.edownqgrid], [0], psym=8, thick=4, /lobar , color=cgcolor('black')
     plotsym, 1, 10, thick=4
     oplot, [d.qgrid], [alog10(dratio)], psym=8, thick=4 , color=cgcolor('black')
   endif else if (flaga eq 1 and flagb eq 2) then begin
     dratio=da/edb
     edratio=da/edb-(da-eda)/edb
     oploterror, [d.qgrid], [alog10(dratio)], [d.eupqgrid], [0], psym=8, thick=4, /hibar , color=cgcolor('black')
     oploterror, [d.qgrid], [alog10(dratio)], [d.edownqgrid], [edratio], psym=8, thick=4, /lobar , color=cgcolor('black')
     plotsym, 2, 10, thick=4
     oplot, [d.qgrid], [alog10(dratio)], psym=8, thick=4 , color=cgcolor('black')
  endif
END

FUNCTION plotratioz, ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, range=range, title=title
; plot line ratio vs metallicity

  if not (keyword_set(range)) then range=[-2.0, 2.0]
  if not (keyword_set(title)) then title='log(ratio)'

  auxqarr=grid[uniq(grid.logq, sort(grid.logq))].logq
  auxqarr0=grid0[uniq(grid0.logq, sort(grid0.logq))].logq
 
  sel1=where(abs(grid.logq-d.qgrid) eq min(abs(grid.logq-d.qgrid))) 
  sel2=where(abs(grid.logq-d.qgrid-d.eupqgrid) eq min(abs(grid.logq-d.qgrid-d.eupqgrid)))
  sel3=where(abs(grid.logq-d.qgrid+d.edownqgrid) eq min(abs(grid.logq-d.qgrid+d.edownqgrid)))

  gridratio=ga/gb
  gridratio0=ga0/gb0
  dratio=da/db
  edratio=sqrt((eda/db)^2+(da/db^2*edb)^2)/dratio/alog(10.)
  
  sel=where(abs(grid.logq-d.qgrid) eq min(abs(grid.logq-d.qgrid))) 

  plot, grid[sel].logz+logOHsun, alog10(gridratio[sel]), yrange=range, xrange=[min(grid.logz)+logOHsun-0.1, max(grid.logz)+logOHsun+0.1], xthick=3, ythick=3, charthick=3, charsize=3, xtitle=textoidl('12+log(O/H)'), ytitle=title, thick=3, linestyle=2, xstyle=1, ystyle=1
  plotsym, 0, 1.2, /fill
;  for i=0, n_elements(auxqarr)-1 do begin
;     sel=where(grid.logq eq auxqarr[i]) 
;     oplot, grid[sel].logz+logOHsun, alog10(gridratio[sel]), color=cgcolor('gray') ;, psym=3
;  endfor
  for i=0, n_elements(auxqarr0)-1 do begin
     sel=where(grid0.logq eq auxqarr0[i]) 
     oplot, grid0[sel].logz+logOHsun, alog10(gridratio0[sel]), color=cgcolor('gray'), thick=3 ;, psym=3
  endfor
  oplot, grid[sel1].logz+logOHsun, alog10(gridratio[sel1]), thick=4, color=cgcolor('red')
  oplot, grid[sel2].logz+logOHsun, alog10(gridratio[sel2]), thick=4, color=cgcolor('red'), linestyle=2
  oplot, grid[sel3].logz+logOHsun, alog10(gridratio[sel3]), thick=4, color=cgcolor('red'), linestyle=2
  if (flaga eq 1 and flagb eq 1) then begin
     dratio=da/db
     edratio=sqrt((eda/db)^2+(da/db^2*edb)^2)/dratio/alog(10.)
     oploterror, [d.Zgrid], [alog10(dratio)], [d.eupZgrid], [edratio], psym=3, thick=4, /hibar , color=cgcolor('black')
     oploterror, [d.Zgrid], [alog10(dratio)], [d.edownZgrid], [edratio], psym=3, thick=4, /lobar , color=cgcolor('black')
  endif else if (flaga eq 2 and flagb eq 1) then begin
     dratio=eda/db
     edratio=eda/(db-edb)-eda/db
     oploterror, [d.Zgrid], [alog10(dratio)], [d.eupZgrid], [edratio], psym=8, thick=4, /hibar , color=cgcolor('black')
     oploterror, [d.Zgrid], [alog10(dratio)], [d.edownZgrid], [0], psym=8, thick=4, /lobar , color=cgcolor('black')
     plotsym, 1, 10, thick=4
     oplot, [d.Zgrid], [alog10(dratio)], psym=8, thick=4 , color=cgcolor('black')
  endif else if (flaga eq 1 and flagb eq 2) then begin
     dratio=da/edb
     edratio=da/edb-(da-eda)/edb
     oploterror, [d.Zgrid], [alog10(dratio)], [d.eupZgrid], [0], psym=8, thick=4, /hibar , color=cgcolor('black')
     oploterror, [d.Zgrid], [alog10(dratio)], [d.edownZgrid], [edratio], psym=8, thick=4, /lobar , color=cgcolor('black')
     plotsym, 2, 10, thick=4
     oplot, [d.Zgrid], [alog10(dratio)], psym=8, thick=4 , color=cgcolor('black')
  endif

END

; ===============================================================
; =================== IZI =======================================
; ===============================================================

FUNCTION izi, fluxin, errorin, idin, gridfile=gridfile, plot=plot, logOHsun=logOHsun, epsilon=epsilon, nz=nz, nq=nq, intergridfile=intergridfile, outgridfile=outgridfile, logzlimits=logzlimits, logqlimits=logqlimits, logzprior=logzprior, logqprior=logqprior, nonorm=nonorm, linear=linear
  
; RENAME INPUT ARRAYS
  flux=fluxin
  error=errorin
  id=idin
  
; CHECK INPUT FOR CONSISTENCY
  nlines=n_elements(flux)
  if (n_elements(error) ne nlines or n_elements(id) ne nlines) then message, 'ERROR: Flux, Error, and ID arrays do not have the same number of elements'
  
; READ STRUCTURE CONTAINING PHOTO-IONIZATION MODEL GRID
  
; IF NOT SPECIFIED BY USER USE DEFAULT: 
; Levesque 2010, HIGH MASS LOSS, CSF 6Myr, n=100 cm^-3
  if not (keyword_set(gridfile)) then gridfile= getenv('IZI_DIR')+'/grids/l09_high_csf_n1e2_6.0Myr.fits'

; IF INTERPOLATED GRID IS PROVIDED USE IT
  if (keyword_set(intergridfile)) then gridfile=intergridfile

; READ GRID
  grid0=mrdfits(gridfile, 1)
  grid0.id=strtrim(grid0.id, 2)
  ngrid=n_elements(grid0.logz)

; GET LINE IDs IN THE GRID AND INDEX OF EACH LINE
  id0=grid0[0].id
  nlines0=n_elements(id0)
  for i=0, nlines0-1 do void=EXECUTE('in'+id0[i]+"=where(id0 eq '"+id0[i]+"')")

; TAKE SOLAR OXYGEN ABUNDANCE FROM MODEL GRID IF NOT PROVIDED
  if not (keyword_set(logOHsun)) then logOHsun=grid0[0].logOHsun      

; CUT GRID TO LOGZLIMITS AND LOGQLIMITS 
  if not (keyword_set(logzlimits)) then logzlimits=[min(grid0.logz+logOHsun), max(grid0.logz+logOHsun)]
  if not (keyword_set(logqlimits)) then logqlimits=[min(grid0.logq), max(grid0.logq)]

  grid0=grid0[where(grid0.logz+logOHsun ge logzlimits[0] and grid0.logz+logOHsun le logzlimits[1] and grid0.logq ge logqlimits[0] and grid0.logq le logqlimits[1])]

; CHANGE LOGZPRIOR TO SOLAR UNITS
  if (keyword_set(logZprior)) then logZprior[*,0]=logZprior-logOHsun

; INTERPOLATE STRUCTURE
  if not (keyword_set(intergridfile)) then begin
  ; MAKE NEW INTERPOLATED GRID WITH FINE SPACING

     nz1=long(50)
     nq1=long(50)
     if (keyword_set(nz)) then nz1=long(nz)
     if (keyword_set(nq)) then nq1=long(nq)
  
     zarr=range(min(grid0.logz), max(grid0.logz), nz1)
     qarr=range(min(grid0.logq), max(grid0.logq), nq1)

     dlogz=(zarr[nz1-1]-zarr[0])/(nz1-1)
     dlogq=(qarr[nq1-1]-qarr[0])/(nq1-1)

     fluxarr=dblarr(nlines0, nz1, nq1)

     nz0=n_elements(zarr0)
     nq0=n_elements(qarr0)


     if not (keyword_set(LINEAR)) then for i=0, nlines0-1 do fluxarr[i,*,*]=TRI_SURF(grid0.flux[i], grid0.logz, grid0.logq, bounds=[zarr[0], qarr[0], zarr[nz1-1], qarr[nq1-1]], nx=nz1, ny=nq1) $
else for i=0, nlines0-1 do fluxarr[i,*,*]=TRI_SURF(grid0.flux[i], grid0.logz, grid0.logq, bounds=[zarr[0], qarr[0], zarr[nz1-1], qarr[nq1-1]], nx=nz1, ny=nq1, /LINEAR) 

;; other interpolation methods I've played with
;     TRIANGULATE, grid0.logz, grid0.logq, tr, b 
;     for i=0, nlines0-1 do fluxarr[i,*,*]=TRIGRID(grid0.logz, grid0.logq, grid0.flux[i], tr, xout=zarr, yout=qarr, /QUINTIC) 
;     TRIANGULATE, grid0.logz, grid0.logq, tr, b 
;     for i=0, nlines0-1 do fluxarr[i,*,*]=TRIGRID(grid0.logz, grid0.logq, grid0.flux[i], tr, xout=zarr, yout=qarr) 
;     for i=0, nlines0-1 do fluxarr[i,*,*] = GRID_TPS(grid0.logz, grid0.logq, grid0.flux[i], NGRID=[nz1, nq1], START=[zarr[0], qarr[0]], DELTA=[dlogz, dlogq]) 
 ;    for i=0, nlines0-1 do fluxarr[i,*,*] = KRIG2D(grid0.flux[i], grid0.logz, grid0.logq, EXPONENTIAL=[10.0, 1.0,0.5], BOUNDS=[zarr[0], qarr[0], zarr[nz1-1], qarr[nq1-1]] , NX=nz1 , NY=nq1) 
;     for i=0, nlines0-1 do fluxarr[i,*,*] = MIN_CURVE_SURF(grid0.flux[i], grid0.logz, grid0.logq, BOUNDS=[zarr[0], qarr[0], zarr[nz1-1], qarr[nq1-1]] , NX=nz1 , NY=nq1, /TPS) 


     grid=replicate({name:grid0[0].name, logz:0., logq:0., id:strtrim(id0,2), flux:fltarr(nlines0)}, nz1*nq1)
  
     for i=0, nz1-1 do begin
        for j=0, nq1-1 do begin
           ngrid=i*nq1+j
           grid[ngrid].logz=zarr[i]
           grid[ngrid].logq=qarr[j]
           for k=0, nlines0-1 do grid[ngrid].flux[k]=fluxarr[k,i,j]
        endfor
     endfor
     ngrid=ngrid+1
  endif else begin     
     grid=grid0
     zarr=grid0[uniq(grid0.logz)].logz
     qarr=grid0[where(grid0.logz eq zarr[0])].logq
     nz1=n_elements(zarr)
     nq1=n_elements(qarr)
     dlogz=(zarr[nz1-1]-zarr[0])/(nz1-1)
     dlogq=(qarr[nq1-1]-qarr[0])/(nq1-1)
  endelse


  ; WRITE INTERPOLATED GRID IF USER WANTS TO
  if (keyword_set(outgridfile)) then mwrfits, grid, outgridfile, /create


  ; Check for summed sets of lines in input ID array and sum fluxes in grid
  ; All fluxes are summed to the first line and ID is set to that line 
  
  for i=0, n_elements(id)-1 do begin
     idsum=strsplit(id[i],';',/extract, count=nlinessum)
     if (nlinessum gt 1) then begin
        for j=1, nlinessum-1 do begin
            void=EXECUTE('grid.flux[in'+idsum[0]+',*]=grid.flux[in'+idsum[0]+',*]+grid.flux[in'+idsum[j]+',*]')
            void=EXECUTE('grid0.flux[in'+idsum[0]+',*]=grid0.flux[in'+idsum[0]+',*]+grid0.flux[in'+idsum[j]+',*]')
        endfor
        id[i]=idsum[0]
     endif
  endfor


  ; INCLUDE SYSTEMATIC UNCERTAINTY IN THE PHOTO-IONIZATION MODELS
  if not (keyword_set(epsilon)) then epsilon=0.15 ; default is 0.15 dex systematic uncertainty
  epsilon2=epsilon*alog(10) ; convert to scaling factor

  
  ; CREATE DATA STRUCTURE CONTAINING LINE FLUXES AND ESTIMATED PARAMETERS
  d={id:id0, $                       ; line id
     flux:dblarr(nlines0)-666, $     ; line flux      
     error:dblarr(nlines0)-666, $    ; line error
     chi2:0., $                      ; chi2 between Observed and Joint Mode fluxes.
     Zgrid:0., $                     ; Metallicity (Grid Interpolation)
     eupZgrid:0., $                  ; Error in Metallicity (Grid Interpolation)
     edownZgrid:0., $                ; Error in Metallicity (Grid Interpolation)
     qgrid:0., $                     ; Ionization Parameter (Grid Interpolation)
     eupqgrid:0., $                  ; Error in Ionization Parameter (Grid Interpolation)
     edownqgrid:0., $                ; Error in Ionization Parameter (Grid Interpolation)
     Zgridmarmod:0., $               ; Metallicity (Grid Interpolation)
     eupZgridmarmod:0., $            ; Error in Metallicity (Grid Interpolation)
     edownZgridmarmod:0., $          ; Error in Metallicity (Grid Interpolation)
     qgridmarmod:0., $               ; Ionization Parameter (Grid Interpolation)
     eupqgridmarmod:0., $            ; Error in Ionization Parameter (Grid Interpolation)
     edownqgridmarmod:0., $          ; Error in Ionization Parameter (Grid Interpolation)
     Zgridmarmean:0., $              ; Metallicity (Grid Interpolation)
     eupZgridmarmean:0., $           ; Error in Metallicity (Grid Interpolation)
     edownZgridmarmean:0., $         ; Error in Metallicity (Grid Interpolation)
     qgridmarmean:0., $              ; Ionization Parameter (Grid Interpolation)
     eupqgridmarmean:0., $           ; Error in Ionization Parameter (Grid Interpolation)
     edownqgridmarmean:0., $         ; Error in Ionization Parameter (Grid Interpolation)
     zarr:fltarr(nz1), $             ; Array of interpolated log(Z)
     zpdfmar:fltarr(nz1), $          ; Marginalized PDF for log(Z)
     qarr:fltarr(nq1), $             ; Array of interpolated log(q)
     qpdfmar:fltarr(nq1), $          ; Marginalized PDF for log(q)
     flags:[0,0,1,1] $              ; Flags: [znpeaks, qnpeaks, zlimit (2:upper 3:lower), qlimit (2:upper 3:lower)]
;     pdfjoint:fltarr(nz1, nq1), $    ; Joint PDF for log(Z) and log(q)
  } 
  
  
 ; FILL STRUCTURE WITH LINE FLUXES
  for i=0, nlines-1 do begin
     auxind=where(d.id eq id[i], nmatch)
     if (nmatch ne 1) then message, 'ERROR: ===== Line ID '+id[i]+'not recognized ====='
     d.flux[auxind]=flux[i]
     d.error[auxind]=error[i]
  endfor
  
  
  ; INDEX OF LINES WITH MEASUREMENTS
  good=where(d.error ne -666, ngood)
  measured=where(d.flux ne -666, nmeasured)
  upperlim=where(d.error ne -666 and d.flux eq -666, nupper)
  flag0=fltarr(nlines0)
  if (measured ne [-1]) then flag0[measured]=1      ; measured flux
  if (upperlim ne [-1]) then flag0[upperlim]=2      ; upper limit on flux
  flag=flag0[good]
  
 ; NORMALIZE LINE FLUXES TO H-BETA OR
 ; IF ABSENT NORMALIZE TO BRIGHTEST LINE

  if not (keyword_set(nonorm)) then begin ; use nonorm for line ratio fitting
  
     print, 'Normalizing Fluxes'

     idnorm='hbeta'
     if (d.flux[inhbeta] eq -666) then idnorm=(reverse((d.id[measured])[sort(d.flux[measured])]))[0]
                                ; normlaize data
     norm=d.flux[where(d.id eq idnorm)]
     d.flux[measured]=d.flux[measured]/norm[0]
     d.error[good]=d.error[good]/norm[0]
                                ; normalize grid
     for i=0, ngrid-1 do begin
        norm=grid[i].flux[where(grid[i].id eq idnorm)]
     grid[i].flux=grid[i].flux/norm[0]
  endfor
  
endif
  
  ; CALCULATE LIKELIHOOD AND POSTERIOR

  like=dblarr(ngrid)+1.0
  post=dblarr(ngrid)+1.0
  zrange=[min(grid.logz), max(grid.logz)]
  qrange=[min(grid.logq), max(grid.logq)]
  
  for i=0, ngrid-1 do begin
     for j=0, ngood-1 do begin
                                ; CALCULATE LIKELIHOOD
        if (flag[j] eq 1) then like[i]=like[i]*1.0/sqrt(2.0*!pi)/sqrt((d.error[good])[j]^2+(epsilon2*(grid[i].flux[good])[j])^2)*exp(-1.0*((d.flux[good])[j]-(grid[i].flux[good])[j])^2/2.0/((d.error[good])[j]^2+(epsilon2*(grid[i].flux[good])[j])^2)) ; if measured
        if (flag[j] eq 2) then like[i]=like[i]*0.5*(1+erf(((d.error[good])[j]-(grid[i].flux[good])[j])/(sqrt((d.error[good])[j]^2+(epsilon2*(grid[i].flux[good])[j])^2)*sqrt(2)))) ; if upper limit
     endfor   
                                ; CALCULATE POSTERIOR BY INCLUDING PRIORS AND NORMALIZING
     if (keyword_set(logzprior) eq 0 and keyword_set(logqprior) eq 0)  then  post[i]=uprior(zrange)*uprior(qrange)*like[i]
     if (keyword_set(logzprior) eq 1 and keyword_set(logqprior) eq 0)  then  post[i]=userprior(grid[i].logz, logzprior[*,0], logzprior[*,1])*uprior(qrange)*like[i]
     if (keyword_set(logzprior) eq 0 and keyword_set(logqprior) eq 1)  then  post[i]=uprior(zrange)*userprior(grid[i].logq, logqprior[*,0], logqprior[*,1])*like[i]
     if (keyword_set(logzprior) eq 1 and keyword_set(logqprior) eq 1)  then  post[i]=userprior(grid[i].logz, logzprior[*,0], logzprior[*,1])*userprior(grid[i].logq, logqprior[*,0], logqprior[*,1])*like[i]
   
  endfor
 
  like[where(finite(like) eq 0)]=0
  post[where(finite(post) eq 0)]=0
  
; SORT LIKELIHOOD AND POSTERIOR FOR GETTING BEST-FIT VALUES AND CONFIDENCE INTERVALS

  goodlike=where(finite(like))
  sortlike=(like[goodlike])[reverse(sort(like[goodlike]))]
  sortz=(grid[goodlike].logz)[reverse(sort(like[goodlike]))]
  sortq=(grid[goodlike].logq)[reverse(sort(like[goodlike]))]
  sumlike=dblarr(n_elements(sortlike))
  for i=0, n_elements(sortlike)-1 do sumlike[i]=total(sortlike[0:i])/total(sortlike)                                 
  goodpost=where(finite(post))  
  sortpost=(post[goodpost])[reverse(sort(post[goodpost]))]
  sortz=(grid[goodpost].logz)[reverse(sort(post[goodpost]))]
  sortq=(grid[goodpost].logq)[reverse(sort(post[goodpost]))]
  sumpost=dblarr(n_elements(sortpost))
  for i=0, n_elements(sortpost)-1 do sumpost[i]=total(sortpost[0:i])/total(sortpost) 

  
; CALCULATE BEST FIT METALLICITY, IONIZATION PARAMETER AND ERRORS

  post1sig=(sortpost[where(sumpost ge 0.683)])[0]
  post2sig=(sortpost[where(sumpost ge 0.955)])[0]
  post3sig=(sortpost[where(sumpost ge 0.997)])[0]

  like1sig=(sortlike[where(sumlike ge 0.683)])[0]
  like2sig=(sortlike[where(sumlike ge 0.955)])[0]
  like3sig=(sortlike[where(sumlike ge 0.997)])[0]


  d.Zgrid=sortz[0]+logOHsun
  d.edownZgrid=sortz[0]-min(sortz[where(sumpost le 0.683)])
  d.eupZgrid=max(sortz[where(sumpost le 0.683)])-sortz[0]

  d.qgrid=sortq[0]
  d.edownqgrid=sortq[0]-min(sortq[where(sumpost le 0.683)])
  d.eupqgrid=max(sortq[where(sumpost le 0.683)])-sortq[0]


  ; COMPUTE chi2
 
  bestgrid=where(grid.logz eq sortz[0] and grid.logq eq sortq[0])
  fobs=d.flux[where(d.flux ne -666)]  
  eobs=d.error[where(d.flux ne -666)]  
  fmod=grid[bestgrid].flux[where(d.flux ne -666)]
  emod=epsilon2*fmod
  d.chi2=total((fobs-fmod)^2/(eobs^2+emod^2))/n_elements(fobs)

; posterior for Z, marginalizing over q
  postz=dblarr(nz1)
  for j=0, nz1-1 do postz[j]=int_tabulated(sortq[where(sortz eq zarr[j])], sortpost[where(sortz eq zarr[j])], /SORT)
  postz=postz/total(postz)
  sumpostz=dblarr(n_elements(postz))
  for i=0, nz1-1 do sumpostz[i]=total(postz[0:i])
  d.Zgridmarmod=zarr[where(postz eq max(postz))]+logOHsun ; max of PDF
  d.Zgridmarmean=total(zarr*postz)/total(postz)+logOHsun ; first monet of PDF
  d.edownZgridmarmod=d.Zgrid-logOHsun-(zarr)[(where(sumpostz ge (1.0-0.683)/2.0))[0]]
  d.eupZgridmarmod=(zarr)[(where(sumpostz ge 1.0-(1.0-0.683)/2.0))[0]]-d.Zgrid+logOHsun
  d.edownZgridmarmean=d.Zgrid-logOHsun-(zarr)[(where(sumpostz ge (1.0-0.683)/2.0))[0]]
  d.eupZgridmarmean=(zarr)[(where(sumpostz ge 1.0-(1.0-0.683)/2.0))[0]]-d.Zgrid+logOHsun

  
; posterior for q, marginalizing over Z
  postq=dblarr(nq1)
  for j=0, nq1-1 do postq[j]=int_tabulated(sortz[where(sortq eq qarr[j])], sortpost[where(sortq eq qarr[j])], /SORT)
  postq=postq/total(postq)
  sumpostq=dblarr(n_elements(postq))
  for i=0, nq1-1 do sumpostq[i]=total(postq[0:i])
  d.qgridmarmod=qarr[where(postq eq max(postq))] ; max of PDF
  d.qgridmarmean=total(qarr*postq)/total(postq) ; first moment of PDF
  d.edownqgridmarmod=d.qgrid-(qarr)[(where(sumpostq ge (1.0-0.683)/2.0))[0]]
  d.eupqgridmarmod=(qarr)[(where(sumpostq ge 1.0-(1.0-0.683)/2.0))[0]]-d.qgrid
  d.edownqgridmarmean=d.qgrid-(qarr)[(where(sumpostq ge (1.0-0.683)/2.0))[0]]
  d.eupqgridmarmean=(qarr)[(where(sumpostq ge 1.0-(1.0-0.683)/2.0))[0]]-d.qgrid

; WRITE MARGINALIZED PDFS
  d.zarr=zarr+logOHsun
  d.zpdfmar=postz
  d.qarr=qarr
  d.qpdfmar=postq

; set FLAGS to warn bout multiple peaks and lower/upper limits
  dzpdf=deriv(postz)
  ddzpdf=deriv(dzpdf)
  auxpeak=intarr(nz1)
  for i=0, nz1-2 do if (dzpdf[i] gt 0 and dzpdf[i+1] le 0 and ddzpdf[i] lt 0) then auxpeak[i]=1
  zpeaks=where(auxpeak eq 1, auxct)
  d.flags[0]=auxct
  
  dqpdf=deriv(postq)
  ddqpdf=deriv(dqpdf)
  auxpeak=intarr(nq1)
  for i=0, nq1-2 do if (dqpdf[i] gt 0 and dqpdf[i+1] le 0 and ddqpdf[i] lt 0) then auxpeak[i]=1
  qpeaks=where(auxpeak eq 1, auxct)
  d.flags[1]=auxct
  
  if (max(postz[0:1]) ge 0.5*max(postz)) then d.flags[2]=2
  if (max(postz[nz1-2:nz1-1]) ge 0.5*max(postz)) then d.flags[2]=3
  if (max(postz[0:1]) ge 0.5*max(postz) and max(postz[nz1-2:nz1-1]) ge 0.5*max(postz)) then d.flags[2]=0
  if (max(postq[0:1]) ge 0.5*max(postq)) then d.flags[3]=2
  if (max(postq[nq1-2:nq1-1]) ge 0.5*max(postq)) then d.flags[3]=3
  if (max(postq[0:1]) ge 0.5*max(postq) and max(postq[nq1-2:nq1-1]) ge 0.5*max(postq)) then d.flags[3]=0


  print, '===== BEST FIT FROM JOINT PDF MODE ====='
  print, '===== Z =====', d.Zgrid, d.edownZgrid, d.eupZgrid, d.Zgrid-d.edownZgrid, d.Zgrid+d.eupZgrid 
  print, '===== q =====', d.qgrid, d.edownqgrid, d.eupqgrid, d.qgrid-d.edownqgrid, d.qgrid+d.eupqgrid
  c=2.99792458e10 ; cm/s 
  print, '=== U=q/c ===', d.qgrid-alog10(c), d.edownqgrid, d.eupqgrid, d.qgrid-alog10(c)-d.edownqgrid, d.qgrid-alog10(c)+d.eupqgrid
  print, '================ FLAGS ====================='
  print, d.flags
  print, '============================================'



; PLOT RESULTS
  if (keyword_set(plot)) then begin

      loadct, 39
    
      set_plot , 'PS'
      device, filename='izipdf.ps', /color, xsize=35, ysize=12, BITS_PER_PIXEL=8 
;      window, 1, xsize=1200, ysize=400
      !p.multi=[0,3,1,0]
      cgloadct, 8
      range=[0, 1e0]*max(post[where(finite(alog10(post)))])
      ncolors=256
      levels = range[0] + (range[1]-range[0])/(ncolors-1.0)*findgen(ncolors)
      contour, post, grid.logz+logOHsun, grid.logq, /fill, /irregular, levels=levels, charsize=3, xtitle=textoidl('12+log(O/H)'), ytitle=textoidl('log(q)'), xrange=[min(grid.logz)+logOHsun, max(grid.logz)+logOHsun], yrange=[min(grid.logq), max(grid.logq)], xstyle=1, ystyle=1, title=textoidl('P(Z,q)d(log Z)d(log q)'), xthick=3, ythick=3, charthick=3
      contour, post[where(finite(alog10(post)))], grid[where(finite(alog10(post)))].logz+logOHsun, grid[where(finite(alog10(post)))].logq, /irregular, levels=[post1sig], /overplot, c_thick=6, color=cgcolor('cyan')
      oplot, [logOHsun, logOHsun], [0,1e2], linestyle=2, thick=10, color=cgcolor('orange')
      oplot, [logOHsun, logOHsun]-1.0, [0,1e2], linestyle=2, thick=10, color=cgcolor('orange')
      xyouts, logOHsun+0.1, min(grid.logq)+0.1, textoidl('Z_{\odot}'), charsize=1.5, charthick=3, color=cgcolor('orange')
      xyouts, logOHsun-1+0.1, min(grid.logq)+0.1, textoidl('0.1Z_{\odot}'), charsize=1.5, charthick=3, color=cgcolor('orange')
      oplot, grid.logz+logOHsun, grid.logq, psym=3, color=cgcolor('black')
      plotsym, 0, 1.5, /fill
      oplot, [d.Zgrid], [d.qgrid], psym=8, color=cgcolor('orange')
      plotsym, 0, 1.5, thick=3
      oplot, [d.Zgrid], [d.qgrid], psym=8, color=cgcolor('black')
      
      
      loadct, 39
      sel=where(finite(postz))
      plot, zarr[sel]+logOHsun, postz[sel], xtitle='12+log(O/H)', ytitle='P(Z)d(log Z)', xthick=3, ythick=3, thick=6, charthick=3, charsize=3, xrange=[min(zarr)+logOHsun, max(zarr)+logOHsun], xstyle=1, ystyle=1, yrange=[0, 1.1*max(postz)]
      oplot, [d.Zgrid, d.Zgrid], [0,1e6], linestyle=0, thick=6, color=cgcolor('red')
      oplot, [d.Zgridmarmod, d.Zgridmarmod], [0,1e6], linestyle=0, thick=4, color=cgcolor('blue')
      oplot, [d.Zgridmarmean, d.Zgridmarmean], [0,1e6], linestyle=0, thick=4, color=cgcolor('green')
      oplot, [d.Zgrid-d.edownZgrid, d.Zgrid-d.edownZgrid], [0,1e6], linestyle=1, thick=2, color=cgcolor('red')
      oplot, [d.Zgrid+d.eupZgrid, d.Zgrid+d.eupZgrid], [0,1e6], linestyle=1, thick=2, color=cgcolor('red')
      oplot, [logOHsun, logOHsun], [0,1e2], linestyle=2, thick=4
      oplot, [logOHsun, logOHsun]-1.0, [0,1e2], linestyle=2, thick=4
      
      
      sel=where(finite(postq))
      plot, qarr[sel], postq[sel], xtitle='log(q)', ytitle='P(q)d(log q)', xthick=3, ythick=3, thick=6, charthick=3, charsize=3, xrange=[min(qarr), max(qarr)], xstyle=1, ystyle=1, yrange=[0, 1.1*max(postq)]
      oplot, [d.qgrid, d.qgrid], [0,1e6], linestyle=0, thick=6, color=cgcolor('red')
      oplot, [d.qgridmarmod, d.qgridmarmod], [0,1e6], linestyle=0, thick=4, color=cgcolor('blue')
      oplot, [d.qgridmarmean, d.qgridmarmean], [0,1e6], linestyle=0, thick=4, color=cgcolor('green')
      oplot, [d.qgrid-d.edownqgrid, d.qgrid-d.edownqgrid], [0,1e6], linestyle=1, thick=4, color=cgcolor('red')
      oplot, [d.qgrid+d.eupqgrid, d.qgrid+d.eupqgrid], [0,1e6], linestyle=1, thick=4, color=cgcolor('red')

       ;; xyouts, min(grid.logq)+1.1, 1.0*max(postq), 'Joint Mode', color=cgcolor('red'), charthick=3, charsize=1.5
       ;; xyouts, min(grid.logq)+1.1, 0.925*max(postq), 'Marg Mode', color=cgcolor('blue'), charthick=3, charsize=1.5
       ;; xyouts, min(grid.logq)+1.1, 0.85*max(postq), 'Marg. Mean', color=cgcolor('green'), charthick=3, charsize=1.5

      ;; xyouts, min(grid.logq)+1.1, 0.2275*max(postq), 'Joint Mode', color=cgcolor('red'), charthick=3, charsize=1.5
      ;; xyouts, min(grid.logq)+1.1, 0.15*max(postq), 'Marg Mode', color=cgcolor('blue'), charthick=3, charsize=1.5
      ;; xyouts, min(grid.logq)+1.1, 0.0725*max(postq), 'Marg. Mean', color=cgcolor('green'), charthick=3, charsize=1.5

      device, /close_file
      set_plot, 'X'
      

; ======== PLOTS VS METALLICITY ============

      !x.margin=[9,1]
      !y.margin=[4,1]

      set_plot , 'PS'
      device, filename='izi_zratios.ps', /color, xsize=40, ysize=20, BITS_PER_PIXEL=8 

;     window, 2, xsize=1000, ysize=600
      !p.multi=[0,4,2,0]
      
                                ; ======= R23 ========
                                ; define line ratio
      ga=grid.flux[inoiii5007]+grid.flux[inoii3726]
      gb=grid.flux[inhbeta]
      ga0=grid0.flux[inoiii5007]+grid0.flux[inoii3726]
      gb0=grid0.flux[inhbeta]
      da=d.flux[inoiii5007]+d.flux[inoii3726]
      db=d.flux[inhbeta]
      eda=sqrt(d.error[inoiii5007]^2+d.error[inoii3726]^2)
      edb=d.error[inhbeta]     
      flaga=1
      if (flag0[inoiii5007] eq 2 or flag0[inoii3726] eq 2) then begin
         flaga=2
         da=-666
         eda=total(([d.flux[inoiii5007], d.flux[inoii3726]])[where([flag0[inoiii5007],flag0[inoii3726]] eq 1)])+total(([d.error[inoiii5007],d.error[inoii3726]])[where([flag0[inoiii5007],flag0[inoii3726]] eq 2)])
     endif
      if (flag0[inoiii5007] eq 0 or flag0[inoii3726] eq 0) then flaga=0
      flagb=flag0[inhbeta]
      void=plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, title=textoidl('log(R23)'), range=[-1.5,1.5])

                                ; ======= N2O2 ========
                                ; define line ratio
     ga=grid.flux[innii6584]
     gb=grid.flux[inoii3726]
     ga0=grid0.flux[innii6584]
     gb0=grid0.flux[inoii3726]
     da=d.flux[innii6584]
     db=d.flux[inoii3726]
     eda=d.error[innii6584]
     edb=d.error[inoii3726]
     flaga=flag0[innii6584]
     flagb=flag0[inoii3726]
     void=plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, title='log(N2O2)', range=[-2,1])


                                ; ======= N2 ========
                                ; define line ratio
     ga=grid.flux[innii6584]
     gb=grid.flux[inhalpha]
     ga0=grid0.flux[innii6584]
     gb0=grid0.flux[inhalpha]
     da=d.flux[innii6584]
     db=d.flux[inhalpha]
     eda=d.error[innii6584]
     edb=d.error[inhalpha]
     flaga=flag0[innii6584]
     flagb=flag0[inhalpha]
     void=plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, title='log(N2)', range=[-3.0, 0.0])
        


                                ; ======= O3N2 ========
                                ; define line ratio
     ga=grid.flux[inoiii5007]
     gb=grid.flux[innii6584]
     ga0=grid0.flux[inoiii5007]
     gb0=grid0.flux[innii6584]
     da=d.flux[inoiii5007]
     db=d.flux[innii6584]
     eda=d.error[inoiii5007]
     edb=d.error[innii6584]
     flaga=flag0[inoiii5007]
     flagb=flag0[innii6584]
     void=plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, title='log(O3N2)', range=[-0.5,2.5])
        

                                ; ======= O3O2 ========
                                ; define line ratio
     ga=grid.flux[inoiii5007]
     gb=grid.flux[inoii3726]
     ga0=grid0.flux[inoiii5007]
     gb0=grid0.flux[inoii3726]
     da=d.flux[inoiii5007]
     db=d.flux[inoii3726]
     eda=d.error[inoiii5007]
     edb=d.error[inoii3726]
     flaga=flag0[inoiii5007]
     flagb=flag0[inoii3726]
     void=plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, title='log(O3O2)', range=[-1.5,1.5])
     

                                ; ======= R3 ========
                                ; define line ratio
     ga=grid.flux[inoiii5007]
     gb=grid.flux[inhbeta]
     ga0=grid0.flux[inoiii5007]
     gb0=grid0.flux[inhbeta]
     da=d.flux[inoiii5007]
     db=d.flux[inhbeta]
     eda=d.error[inoiii5007]
     edb=d.error[inhbeta]
     flaga=flag0[inoiii5007]
     flagb=flag0[inhbeta]
     void=plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb,  grid, grid0, logOHsun, d, title='log(R3)', range=[-2,1])

                                ; ======= N2S2 ========
                                ; define line ratio
     ga=grid.flux[innii6584]
     gb=grid.flux[insii6717]+grid.flux[insii6731]
     ga0=grid0.flux[innii6584]
     gb0=grid0.flux[insii6717]+grid0.flux[insii6731]
     da=d.flux[innii6584]
     db=d.flux[insii6717]
     eda=d.error[innii6584]
     edb=d.error[insii6717]
     flaga=flag0[innii6584]
     flagb=flag0[insii6717]
     void=plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, title='log(N2S2)', range=[-1.5,1.5])

                                ; ======= S2 ========
                                ; define line ratio
     ga=grid.flux[insii6717]+grid.flux[insii6731]
     gb=grid.flux[inhalpha]
     ga0=grid0.flux[insii6717]+grid0.flux[insii6731]
     gb0=grid0.flux[inhalpha]
     da=d.flux[insii6717]
     db=d.flux[inhalpha]
     eda=d.error[insii6717]
     edb=d.error[inhalpha]
     flaga=flag0[insii6717]
     flagb=flag0[inhalpha]
     void=plotratioz(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, title='log(S2)', range=[-3.0,0.0])


      device, /close_file
      set_plot, 'X'

; ============= PLOTS VERSUS IONIZATION PARAMETER ===========


;     window, 3, xsize=1000, ysize=600
      set_plot , 'PS'
      device, filename='izi_qratios.ps', /color, xsize=40, ysize=20, BITS_PER_PIXEL=8 
      
      !p.multi=[0,4,2,0]
     

                                ; ======= R23 ========
                                ; define line ratio
     ga=grid.flux[inoiii5007]+grid.flux[inoii3726]
     gb=grid.flux[inhbeta]
     ga0=grid0.flux[inoiii5007]+grid0.flux[inoii3726]
     gb0=grid0.flux[inhbeta]
     da=d.flux[inoiii5007]+d.flux[inoii3726]
     db=d.flux[inhbeta]
     eda=sqrt(d.error[inoiii5007]^2+d.error[inoii3726]^2)
     edb=d.error[inhbeta]
     flaga=1
     if (flag0[inoiii5007] eq 2 or flag0[inoii3726] eq 2) then begin
        flaga=2
        da=-666
        eda=total(([d.flux[inoiii5007],d.flux[inoii3726]])[where([flag0[inoiii5007],flag0[inoii3726]] eq 1)])+total(([d.error[inoiii5007], d.error[inoii3726]])[where([flag0[inoiii5007],flag0[inoii3726]] eq 2)])
     endif
     if (flag0[inoiii5007] eq 0 or flag0[inoii3726] eq 0) then flaga=0
     flagb=flag0[inhbeta]
     void=plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, title='log(R23)')
     
                                ; ======= N2O2 ========
                                ; define line ratio
     ga=grid.flux[innii6584]
     gb=grid.flux[inoii3726]
     ga0=grid0.flux[innii6584]
     gb0=grid0.flux[inoii3726]
     da=d.flux[innii6584]
     db=d.flux[inoii3726]
     eda=d.error[innii6584]
     edb=d.error[inoii3726]
     flaga=flag0[innii6584]
     flagb=flag0[inoii3726]
     void=plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, title='log(N2O2)', range=[-3,1])

                                ; ======= N2 ========
                                ; define line ratio
     ga=grid.flux[innii6584]
     gb=grid.flux[inhalpha]
     ga0=grid0.flux[innii6584]
     gb0=grid0.flux[inhalpha]
     da=d.flux[innii6584]
     db=d.flux[inhalpha]
     eda=d.error[innii6584]
     edb=d.error[inhalpha]
     flaga=flag0[innii6584]
     flagb=flag0[inhalpha]
     void=plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, title='log(N2)', range=[-3.5, 0.5])

                                ; ======= O3N2 ========

     ga=grid.flux[inoiii5007]
     gb=grid.flux[innii6584]
     ga0=grid0.flux[inoiii5007]
     gb0=grid0.flux[innii6584]
     da=d.flux[inoiii5007]
     db=d.flux[innii6584]
     eda=d.error[inoiii5007]
     edb=d.error[innii6584]
     flaga=flag0[inoiii5007]
     flagb=flag0[innii6584]
     void=plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, title='log(O3N2)', range=[-2.0,2.0])

                                ; ======= O3O2 ========
                                ; define line ratio
     ga=grid.flux[inoiii5007]
     gb=grid.flux[inoii3726]
     ga0=grid0.flux[inoiii5007]
     gb0=grid0.flux[inoii3726]
     da=d.flux[inoiii5007]
     db=d.flux[inoii3726]
     eda=d.error[inoiii5007]
     edb=d.error[inoii3726]
     flaga=flag0[inoiii5007]
     flagb=flag0[inoii3726]
     void=plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, title='log(O3O2)', range=[-2.5,1.5])
     
                                ; ======= O3 ========
                                ; define line ratio
     ga=grid.flux[inoiii5007]
     gb=grid.flux[inhbeta]
     ga0=grid0.flux[inoiii5007]
     gb0=grid0.flux[inhbeta]
     da=d.flux[inoiii5007]
     db=d.flux[inhbeta]
     eda=d.error[inoiii5007]
     edb=d.error[inhbeta]
     flaga=flag0[inoiii5007]
     flagb=flag0[inhbeta]
     void=plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, title='log(R3)', range=[-2.5,1.5])


                                ; ======= N2S2 ========
                                ; define line ratio
     ga=grid.flux[innii6584]
     gb=grid.flux[insii6717]+grid.flux[insii6731]
     ga0=grid0.flux[innii6584]
     gb0=grid0.flux[insii6717]+grid0.flux[insii6731]
     da=d.flux[innii6584]
     db=d.flux[insii6717]
     eda=d.error[innii6584]
     edb=d.error[insii6717]
     flaga=flag0[innii6584]
     flagb=flag0[insii6717]
     void=plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, title='log(N2S2)', range=[-2.5,1.5])


                                ; ======= S2 ========
                                ; define line ratio
     ga=grid.flux[insii6717]+grid.flux[insii6731]
     gb=grid.flux[inhalpha]
     ga0=grid0.flux[insii6717]+grid0.flux[insii6731]
     gb0=grid0.flux[inhalpha]
     da=d.flux[insii6717]
     db=d.flux[inhalpha]
     eda=d.error[insii6717]
     edb=d.error[inhalpha]
     flaga=flag0[insii6717]
     flagb=flag0[inhalpha]
     void=plotratioq(ga, gb, ga0, gb0, da, db, eda, edb, flaga, flagb, grid, grid0, logOHsun, d, title='log(S2)', range=[-3.0,1.0])



      device, /close_file
      set_plot, 'X'


  endif

  ; RETURN RESULTS
  return, d
  
  
END


