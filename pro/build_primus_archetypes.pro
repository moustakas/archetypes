code is fine (although incomplete right at the end), but made obsolete by build_primus_lp_archetypes 
    
;+
; NAME:
;   BUILD_PRIMUS_ARCHETYPES
;
; PURPOSE:
;   Build the PRIMUS archetypes (template set).
;
; OPTIONAL INPUTS: 
;   del_chisq - delta-chi^2 value, used to assess whether a new
;     archetype is needed
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Nov, UCSD - rewritten, based largely on
;     K. Wong's PRIMUS_BUILD_TEMPSET
;
; Copyright (C) 2009, John Moustakas
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

pro build_primus_archetypes, del_chisq=del_chisq

    common test, info, ages, fopt_primus

    ftweak = 0.1 ; refluxing prior
    if (n_elements(del_chisq) eq 0) then del_chisq = 3.0

; fiducial redshift; the fiducial mask and rerun is needed to get a
; "typical" PRIMUS inverse variance array 
    fiducial_rerun = '0021'
    fiducial_mask = 'vvds0240'
    zbase = 0.25     
    
; read the parent sample of AGES spectra
    outpath = getenv('PRIMUS_DATA')+'/template_basis/'
    if (n_elements(info) eq 0) or (n_elements(ages) eq 0) then $
      ages = primus_readin_parenttemplate(version='2.0',$
      templatepar=info)
    nobj = n_elements(info)
    
; sort by decreasing *dereddened* u-g color, putting the red galaxies
; first 
    srt = reverse(sort(info.ug_dered))
    info = info[srt]
    ages.flux = ages.flux[*,srt]

; generate a fiducial wavelength and inverse variance array using a
; "typical" mask; we fit everything between 4500 and 9500 A
    splog, 'Generating fiducial inverse variance array'
    wave1 = primus_default_wavevec(wrange=[4500.0,9500.0])
    npix = n_elements(wave1) 

    primus_read, fiducial_mask, fiducial_rerun, extract=extract1
    good = where(extract1.badextract eq 0)
    extract = extract1[good]

;interpolate object ivars to match wave1    
    fivar_interp = fltarr(npix,n_elements(extract)) 
    for ii = 0, n_elements(extract)-1 do fivar_interp[*,ii] = $
      interpol(extract[ii].fivar1,extract[ii].wave1,wave1)
    ivar_primus = djs_median(fivar_interp,2)

; normalize and convolve each model spectrum to PRIMUS resolution at a
; fiducial redshift ZBASE
    splog, 'Normalizing and convolving to PRIMUS resolution'
    ww = where((ages.wave gt 6800.0) and (ages.wave lt 7200.0))
    for ii = 0, nobj-1 do ages.flux[*,ii] = ages.flux[*,ii]*$
      median(ages.flux[ww,0])/median(ages.flux[ww,ii])

    if (n_elements(fopt_primus) eq 0) then begin
    fopt_primus = fltarr(npix,nobj)
    for ii = 0L, nobj-1 do begin 
       fopt_primus[*,ii] = template2primus(ages.wave*(1.0+zbase),ages.flux[*,ii],$
         wave1,ccd=5,slitwidth=8,pair=1,/default_throughput,airmass=1.1);,$
;        objpsf=fltarr(npix)+2.25,instpsf=fltarr(npix)+0.9)
    endfor
    endif

; assign each spectrum a constant S/N~20 between 6800 and 7200 A
    waverange = where((wave1 gt 6800.0) and (wave1 lt 7200.0))
    fopt1 = fopt_primus*0.0
    ivar1 = fopt_primus*0.0
    for ii = 0, nobj-1 do begin 
       fopt1[*,ii] = fopt_primus[*,ii]/median(fopt_primus[waverange,ii])
       snr = mean(fopt1[waverange,ii] * sqrt(ivar_primus[waverange]))
       ivar1[*,ii] = ivar_primus * (20.0 / snr)^2
    endfor

; generate the refluxing vectors
    nreflux = 2L
    reflux = fltarr(npix,nreflux)
    for ii = 0, nreflux-1 do reflux[*,ii] = (findgen(npix)/float(npix-1)*2.0-1.0)^(ii+1)

; generate the basis set
    soname = filepath('libprimus.'+idlutils_so_ext(),root_dir=$
      getenv("PRIMUS_DIR"),subdirectory='lib')

    basis1 = {$
      wave:         wave1, $
      flux:  fltarr(npix), $
      ivar:  fltarr(npix), $
      index:           -1, $ ; relative to the parent sample
      nfit:             1}
    basis = basis1
    basis.flux = fopt1[*,0] ; initialize
    basis.ivar = ivar1[*,0]
    basis.index = 0
    
;   for ii = 0, 100-1 do begin 
    for ii = 0, nobj-1 do begin 

       nbasis = n_elements(basis)
       chisq = fltarr(nbasis)
       coeffs = fltarr(nreflux+1,nbasis)
       index = 0L
       s = call_external(soname,'primus_iterfit_call',$
         float(fopt1[*,ii]),float(ivar1[*,ii]),float(basis.flux),$
         float(reflux),float(ftweak),long(nbasis),long(nreflux),$
         long(npix),long(index),float(chisq),float(coeffs))
       
; generate the model PRIMUS spectrum for testing purposes
;      model[*,ii] = coeffs[0,index]*reflux[*,0]*fopt1[*,ii] + $
;        coeffs[1,index]*reflux[*,1]*fopt1[*,ii] + $
;        coeffs[2,index]*tempset[*,index]

; if min(chisq)>del_chisq based on the current basis set, then add
; this template to the basis set
       if (min(chisq) gt del_chisq) then begin 
          splog, 'Adding template '+string(ii,format='(I0)')+$
            ' to the basis set'
          newbasis = basis1
          newbasis.flux = fopt1[*,ii]
          newbasis.ivar = ivar1[*,ii]
          newbasis.index = ii
          basis = [temporary(basis),newbasis]
       endif else begin
          these = where(chisq lt del_chisq)
          basis[these].nfit = basis[these].nfit + 1
;         splog, 'Template '+string(ii,format='(I0)')+$
;           ' is fitted by basis template '+string(ii,format='(I0)')+$
;           ' with chi^2='+min(
;         basis[index].nfit = basis[index].nfit + total(chisq lt del_chisq)
       endelse

;      djs_plot, ages.wave, ages.flux[*,ii]/median(ages.flux[*,ii]), xsty=3, ysty=3, xrange=[3500,9900]
;      djs_oplot, ages.wave, ages.flux[*,index]/median(ages.flux[*,index]), color='cyan'
;      cc = get_kbrd(1)
       
; if the lowest chi2 between the current model in the parent sample
; and the best-fit template in the existing set basis is > DEL_CHISQ,
; then add the current model to the basis set
    endfor

    plothist, nfit, bin=0.5, xx, yy

stop
    
    if keyword_set(model) then $
      model = model[*,temp_index]
    
;create output files
    nbasis = n_elements(temp_index)
    basis = replicate({ages_id : long(0),$
      highwave : fltarr(n_elements(wavedata[0].restwave)),$
      highflux : fltarr(n_elements(wavedata[0].restflux))},nbasis)
    basis[0].ages_id = 0.0
    basis[0].highwave = red_temp.highwave
    basis[0].highflux = red_temp.highflux_mean
    for i = 1, nbasis - 1 do begin $
      basis[i].ages_id = data[temp_index[i]].ages_id
       basis[i].highwave = wavedata[temp_index[i]].restwave
       basis[i].highflux = wavedata[temp_index[i]].restflux
    endfor

stop    
    
;create plots
    !p.multi=[0,1,2]
set_plot,'ps'

outps='ages_d4000_halpha_'+strtrim(string(del_chisq,format='(f4.1)'),1)+'.ps'
device,filename=outps,xsize=7,ysize=10,yoffset=0.5,/inches
plotsym,8,0.4
plot,data.d4000_narrow_model[0],data.h_alpha_ew[0],xr=[0.9,2.9],yr=[0.01,1E3],$
	title='AGES Data after cuts, del_chi^2='+strtrim(string(del_chisq,format='(f4.1)'),1),xtitle='D4000',ytitle='Log H-alpha EW (A)',$
	charsize=1.5,xthick=4.0,ythick=4.0,xs=1,ys=1,/ylog,psym=8
plot,data[temp_index].d4000_narrow_model[0],data[temp_index].h_alpha_ew[0],xr=[0.9,2.9],yr=[0.01,1E3],title='Basis set',$
	xtitle='D4000',ytitle='H-alpha EW (A)',psym=8,charsize=1.5,xthick=4.0,ythick=4.0,$
	xs=1,ys=1,/ylog
device,/close
spawn,'gzip '+outps

u = fltarr(n_elements(wavedata))
g = fltarr(n_elements(wavedata))
r = fltarr(n_elements(wavedata))
u_basis = fltarr(n_elements(basis))
g_basis = fltarr(n_elements(basis))
r_basis = fltarr(n_elements(basis))
for i = 0, n_elements(wavedata) - 1 do begin $
	maggies = k_project_filters(wavedata[i].restwave,wavedata[i].restflux,/silent)
	u[i] = -2.5 * alog10(maggies[0,0,0])
	g[i] = -2.5 * alog10(maggies[0,0,1])
	r[i] = -2.5 * alog10(maggies[0,0,2])
endfor
for i = 0, n_elements(basis) - 1 do begin $
	basis_maggies = k_project_filters(basis[i].highwave,basis[i].highflux,/silent)
	u_basis[i] = -2.5 * alog10(basis_maggies[0,0,0])
	g_basis[i] = -2.5 * alog10(basis_maggies[0,0,1])
	r_basis[i] = -2.5 * alog10(basis_maggies[0,0,2])
endfor
outps='ages_ug_gr_'+strtrim(string(del_chisq,format='(f4.1)'),1)+'.ps'
device,filename=outps,xsize=7,ysize=10,yoffset=0.5,/inches
plot,u-g,g-r,title='AGES Data after cuts, del_chi^2='+strtrim(string(del_chisq,format='(f4.1)'),1),xtitle='u-g',ytitle='g-r',$
	psym=8,xs=1,ys=1,xr=[-1,3],yr=[-1,2],charsize=1.5,xthick=4.0,ythick=4.0
plot,u_basis-g_basis,g_basis-r_basis,title='Basis set',xtitle='u-g',ytitle='g-r',psym=8,xs=1,ys=1,$
	xr=[-1,3],yr=[-1,2],charsize=1.5,xthick=4.0,ythick=4.0
device,/close
spawn,'gzip '+outps

!p.multi = [0,1,5]
outps = 'ages_basis_set_'+strtrim(string(del_chisq,format='(f4.1)'),1)+'.ps'
device,filename=outps,xsize=7,ysize=10,yoffset=0.5,/inches
for i = 0, nbasis - 1 do begin $
	plot,basis[i].highwave,basis[i].highflux,xr=[3000,8000],title=string(i)+'    AGES ID#: '+strtrim(string(basis[i].ages_id,format='(i5)'),1),$
		xtitle='Rest Wavelength (A)',ytitle='Flux',charsize=1.2
endfor
device,/close
spawn,'gzip '+outps

outps = 'ages_basis_set_'+strtrim(string(del_chisq,format='(f4.1)'),1)+'_detail.ps'
device,filename=outps,xsize=7,ysize=10,yoffset=0.5,/inches
!p.multi = [0,1,4]
for i = 0, nbasis - 1 do begin $
	m = median(basis[i].highflux[where(basis[i].highwave ge 3000 and basis[i].highwave le 8000)])
	s = stddev(basis[i].highflux[where(basis[i].highwave ge 3000 and basis[i].highwave le 8000)])
	plot,basis[i].highwave,basis[i].highflux,xs=1,xr=[3000,8000],charsize=1.2,$
		title=string(i)+'    AGES ID#: '+strtrim(string(basis[i].ages_id,format='(i5)'),1),xtitle='Rest Wavelength (A)',ytitle='Flux'
	plot,basis[i].highwave,basis[i].highflux,xs=1,ys=1,xr=[3000,8000],yr=[0,m+2*s],$
		title='Continuum detail view',xtitle='Rest Wavelength (A)',ytitle='Flux',charsize=1.2
endfor
device,/close
spawn,'gzip '+outps
set_plot,'x'
!p.multi = [0,1,1]

;create fits file
outfits = 'ages_basis_set_'+strtrim(string(del_chisq,format='(f4.1)'),1)+'.fits'
mwrfits,basis,outfits,/create
spawn,'gzip '+outfits

return
end
