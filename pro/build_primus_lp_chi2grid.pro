;+
; NAME:
;   BUILD_PRIMUS_LP_CHI2GRID
;
; PURPOSE:
;   Using the parent sample of AGES galaxies, build the chi^2 grid
;   that BUILD_PRIMUS_LP_ARCHETYPES expects.
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;   clobber - overwrite the existing output file
;
; OUTPUTS: 
;   A chi^2 grid in binary FITS format.
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Nov, UCSD - largely based on existing code (see
;     BUILD_PRIMUS_AGES_ARCHETYPES, which was based on code by
;     K. Wong)
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

; build_primus_lp_chi2grid, version='v2.0', ftweak=0.1
; build_primus_lp_chi2grid, version='v3.0', ftweak=0.0

pro build_primus_lp_chi2grid, version=version, ftweak=ftweak, clobber=clobber

    common lp_chi2grid, info, ages, fopt_primus

; output file names
    if (n_elements(version) eq 0) then version = $
      primus_parenttemplate_defaultversion()
    grid_outfile = getenv('PRIMUS_DATA')+'/ages/'+$
      'ages_chi2grid_'+version+'.fits'

    if file_test(grid_outfile+'.gz',/reg) and $
      (keyword_set(clobber) eq 0) then begin
       splog, 'File '+grid_outfile+'.gz exists; use /CLOBBER'
       return
    endif
       
; fiducial redshift; the fiducial mask and rerun is needed to get a
; "typical" PRIMUS inverse variance array 
    if (n_elements(ftweak) eq 0) then ftweak = 0.1

    fiducial_rerun = '0021'
    fiducial_mask = 'vvds0240'
    zbase = 0.25     

; read the parent sample of AGES spectra (see
; BUILD_PRIMUS_AGES_PARENT) 
    if (n_elements(info) eq 0) or (n_elements(ages) eq 0) then $
      ages = primus_readin_parenttemplate(version=version,$
      templatepar=info)
    nobj = n_elements(info)

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
    ww = where((ages.wave gt 6800.0) and (ages.wave lt 7200.0))
    for ii = 0, nobj-1 do ages.flux[*,ii] = ages.flux[*,ii]*$
      median(ages.flux[ww,0])/median(ages.flux[ww,ii])

    if (n_elements(fopt_primus) eq 0) then begin       
       splog, 'Convolving to PRIMUS resolution'
       t0 = systime(1)
       fopt_primus = fltarr(npix,nobj)
       for ii = 0L, nobj-1 do begin 
          fopt_primus[*,ii] = template2primus(ages.wave*(1.0+zbase),ages.flux[*,ii],$
            wave1,ccd=5,slitwidth=8,pair=1,/default_throughput,airmass=1.1)
       endfor
       splog, 'Total time = ', (systime(1)-t0)/60.0
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

;; pack the convolved spectra into a structure and write out    
;    splog, 'Writing '+convolved_outfile+'.gz'
;    mwrfits, chi2grid, convolved_outfile, /create
;    spawn, 'gzip -f '+convolved_outfile, /sh
    
; generate the refluxing vectors
    nreflux = 2L
    reflux = fltarr(npix,nreflux)
    for ii = 0, nreflux-1 do reflux[*,ii] = $
      (findgen(npix)/float(npix-1)*2.0-1.0)^(ii+1)

; build the chi^2 grid; for every object in the parent set, ii, compute
; the chi^2 between that object and every other spectrum in the parent
; set, jj; obviously chi^2=0.0 for ii=jj
    soname = filepath('libprimus.'+idlutils_so_ext(),root_dir=$
      getenv("PRIMUS_DIR"),subdirectory='lib')

    out = {$
      rerun:    fiducial_rerun,$
      mask:     fiducial_mask,$
      zbase:    zbase,$
      ftweak:   ftweak,$
      nobj:     fix(nobj),$
      nreflux:  fix(nreflux),$
      reflux:   reflux,$
      wave:     wave1,$
      flux:     fopt1,$
      ivar:     ivar1,$
      coeffs:   fltarr(nreflux+1,nobj),$
      chi2grid: fltarr(nobj,nobj)}

    t0 = systime(1)
    for ii = 0, nobj-1 do begin 
       print, format='("Object ",I0,"/",I0,A10,$)', ii+1, $
         nobj-1, string(13b)
       chisq = fltarr(nobj)
       coeffs = fltarr(nreflux+1,nobj)
       index = 0L
       
       s = call_external(soname,'primus_iterfit_call',$
         float(fopt1[*,ii]),float(ivar1[*,ii]),float(fopt1),$
         float(reflux),float(ftweak),long(nobj),long(nreflux),$
         long(npix),long(index),float(chisq),float(coeffs))
; in BUILD_PRIMUS_LP_ARCHETYPES we need the first index to correspond
; to the data and the second index to correspond to each model
       out.chi2grid[ii,*] = chisq
       out.coeffs = coeffs
    endfor
    splog, 'Total time = ', (systime(1)-t0)/60.0
    
; write out
    splog, 'Writing '+grid_outfile+'.gz'
    mwrfits, out, grid_outfile, /create
    spawn, 'gzip -f '+grid_outfile, /sh

return
end
