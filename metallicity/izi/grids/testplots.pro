pro testplots

dinf=mrdfits('d13_kappaINF.fits',1)
d20=mrdfits('d13_kappa20.fits',1)
l10=mrdfits('l09_high_csf_n1e2_6.0Myr.fits',1)
k01=mrdfits('d00_SB99_ins_n3.5e2_0.0Myr.fits',1)




zarr=d20[uniq(d20.logz, sort(d20.logz))].logz
qarr=d20[uniq(d20.logq, sort(d20.logq))].logq


; Figure 20
window, 1, xsize=800, ysize=800
plot, [1e6], [1e6], xrange=[-1.1, 0.8], yrange=[-2.0, 1.5], charsize=2, /xs, /ys, xtitle='NII/SII', ytitle='NII/OII'

; Dopita: kappa=INF (white) and kappa=20 (green)

inoii3726=where(strtrim(d20[0].id,2) eq 'oii3726')
inoii3729=where(strtrim(d20[0].id,2) eq 'oii3729')
inoiii5007=where(strtrim(d20[0].id,2) eq 'oiii5007')
innii6584=where(strtrim(d20[0].id,2) eq 'nii6584')
insii6717=where(strtrim(d20[0].id,2) eq 'sii6717')
insii6731=where(strtrim(d20[0].id,2) eq 'sii6731')

n2s2=alog10(dinf.flux[innii6584]/(dinf.flux[insii6717]+dinf.flux[insii6731]))
n2o2=alog10(dinf.flux[innii6584]/(dinf.flux[inoii3726]+dinf.flux[inoii3729]))
for i=0, n_elements(zarr)-1 do oplot, n2s2[where(dinf.logz eq zarr[i])], n2o2[where(dinf.logz eq zarr[i])]
for i=0, n_elements(qarr)-1 do oplot, n2s2[where(dinf.logq eq qarr[i])], n2o2[where(dinf.logq eq qarr[i])]

n2s2=alog10(d20.flux[innii6584]/(d20.flux[insii6717]+d20.flux[insii6731]))
n2o2=alog10(d20.flux[innii6584]/(d20.flux[inoii3726]+d20.flux[inoii3729]))
for i=0, n_elements(zarr)-1 do oplot, n2s2[where(d20.logz eq zarr[i])], n2o2[where(d20.logz eq zarr[i])], color=cgcolor('green')
for i=0, n_elements(qarr)-1 do oplot, n2s2[where(d20.logq eq qarr[i])], n2o2[where(d20.logq eq qarr[i])], color=cgcolor('green')



; Figure 17
window, 0, xsize=800, ysize=800
plot, [1e6], [1e6], xrange=[-2.5, 2], yrange=[-2.0, 1.0], charsize=2, /xs, /ys, xtitle='NII/OII', ytitle='OIII/OII'

; Dopita: kappa=INF (white) and kappa=20 (green)
inoii3726=where(strtrim(d20[0].id,2) eq 'oii3726')
inoii3729=where(strtrim(d20[0].id,2) eq 'oii3729')
inoiii5007=where(strtrim(d20[0].id,2) eq 'oiii5007')
innii6584=where(strtrim(d20[0].id,2) eq 'nii6584')
insii6717=where(strtrim(d20[0].id,2) eq 'sii6717')
insii6731=where(strtrim(d20[0].id,2) eq 'sii6731')


o3o2=alog10(dinf.flux[inoiii5007]/(dinf.flux[inoii3726]+dinf.flux[inoii3729]))
n2o2=alog10(dinf.flux[innii6584]/(dinf.flux[inoii3726]+dinf.flux[inoii3729]))
for i=0, n_elements(zarr)-1 do oplot, n2o2[where(dinf.logz eq zarr[i])], o3o2[where(dinf.logz eq zarr[i])]
for i=0, n_elements(qarr)-1 do oplot, n2o2[where(dinf.logq eq qarr[i])], o3o2[where(dinf.logq eq qarr[i])]

o3o2=alog10(d20.flux[inoiii5007]/(d20.flux[inoii3726]+d20.flux[inoii3729]))
n2o2=alog10(d20.flux[innii6584]/(d20.flux[inoii3726]+d20.flux[inoii3729]))
for i=0, n_elements(zarr)-1 do oplot, n2o2[where(d20.logz eq zarr[i])], o3o2[where(d20.logz eq zarr[i])], color=cgcolor('green')
for i=0, n_elements(qarr)-1 do oplot, n2o2[where(d20.logq eq qarr[i])], o3o2[where(d20.logq eq qarr[i])], color=cgcolor('green')

; Levesque (blue)

inoii3726=where(strtrim(l10[0].id,2) eq 'oii3726')
inoii3729=where(strtrim(l10[0].id,2) eq 'oii3729')
inoiii5007=where(strtrim(l10[0].id,2) eq 'oiii5007')
innii6584=where(strtrim(l10[0].id,2) eq 'nii6584')
insii6717=where(strtrim(l10[0].id,2) eq 'sii6717')
insii6731=where(strtrim(l10[0].id,2) eq 'sii6731')

zarr=l10[uniq(l10.logz, sort(l10.logz))].logz
qarr=l10[uniq(l10.logq, sort(l10.logq))].logq

o3o2=alog10(l10.flux[inoiii5007]/(l10.flux[inoii3726]+l10.flux[inoii3729]))
n2o2=alog10(l10.flux[innii6584]/(l10.flux[inoii3726]+l10.flux[inoii3729]))
for i=0, n_elements(zarr)-1 do oplot, n2o2[where(l10.logz eq zarr[i])], o3o2[where(l10.logz eq zarr[i])], color=cgcolor('blue')
for i=0, n_elements(qarr)-1 do oplot, n2o2[where(l10.logq eq qarr[i])], o3o2[where(l10.logq eq qarr[i])], color=cgcolor('blue')

; Kewley (red)

inoii3726=where(strtrim(k01[0].id,2) eq 'oii3726')
inoii3729=where(strtrim(k01[0].id,2) eq 'oii3729')
inoiii5007=where(strtrim(k01[0].id,2) eq 'oiii5007')
innii6584=where(strtrim(k01[0].id,2) eq 'nii6584')
insii6717=where(strtrim(k01[0].id,2) eq 'sii6717')
insii6731=where(strtrim(k01[0].id,2) eq 'sii6731')

zarr=k01[uniq(k01.logz, sort(k01.logz))].logz
qarr=k01[uniq(k01.logq, sort(k01.logq))].logq

print, qarr

o3o2=alog10(k01.flux[inoiii5007]/(k01.flux[inoii3726]+k01.flux[inoii3729]))
n2o2=alog10(k01.flux[innii6584]/(k01.flux[inoii3726]+k01.flux[inoii3729]))
for i=0, n_elements(zarr)-1 do oplot, n2o2[where(k01.logz eq zarr[i])], o3o2[where(k01.logz eq zarr[i])], color=cgcolor('red')
for i=0, n_elements(qarr)-1 do oplot, n2o2[where(k01.logq eq qarr[i])], o3o2[where(k01.logq eq qarr[i])], color=cgcolor('red')


stop
end
