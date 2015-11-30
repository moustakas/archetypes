;+
; NAME:
;   BUILD_PRIMUS_LP_ARCHETYPES
;
; PURPOSE:
;   Build the PRIMUS archetypes (template set) using LP CPLEX.
;
; OPTIONAL INPUTS: 
;   del_chisq - delta-chi^2 value, used to assess whether a new
;     archetype is needed (can be a vector)
;
; KEYWORD PARAMETERS: 
;   lp_solve - solve the LP problem (takes a long time!)
;   lp_parse - parse the LP output and write out the final set of
;     archetypes 
;   final - once the final DEL_CHISQ value has been chosen, use
;   /LP_PARSE and /FINAL to write out the final archetypes 
;   clobber - overwrite existing files of the same name
;   qaplot - build a QAplot of the archetypes
;
; OUTPUTS: 
;   Various.
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Nov, UCSD - largely based on existing code (see
;     BUILD_PRIMUS_AGES_ARCHETYPES, which was based on code by
;     K. Wong); essentially *all* the LP code was written by D. Hogg
;     (see the astrometry code repository)
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

pro build_primus_lp_archetypes, del_chisq, version=version, lp_format=lp_format, $
  lp_solve=lp_solve, lp_parse=lp_parse, final=final, clobber=clobber, $
  qaplot=qaplot, hogg=hogg, basis=basis, addmgii=addmgii

    common test, info, ages

    light = 2.99792458D5        ; speed of light [km/s]
    fwhm2sig = 2D0*sqrt(2.0*alog(2.0))

    if (n_elements(version) eq 0) then version = $
      primus_parenttemplate_defaultversion()
    if (n_elements(del_chisq) eq 0) then del_chisq = 2.0

;   del_chisq = [0.5,0.7,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.0,10.0]
;   build_primus_lp_archetypes, del_chisq, /lp_format, /clobber
    
; call this routine recursively    
    ndel_chisq = n_elements(del_chisq)
    if (ndel_chisq gt 1) then begin
       if keyword_set(final) then begin
          splog, 'DEL_CHISQ must be a scalar when /FINAL'
          return
       endif
       for ii = 0, ndel_chisq-1 do begin
          build_primus_lp_archetypes, del_chisq[ii], version=version, lp_format=lp_format, $
            lp_solve=lp_solve, lp_parse=lp_parse, final=final, clobber=clobber, $
            qaplot=qaplot, hogg=hogg, basis=basis, addmgii=addmgii
       endfor
; make a tarball for Sam
       if keyword_set(lp_format) then begin
          spawn, 'tar czvf '+getenv('PRIMUS_DATA')+'/ages/'+$
            'ages.lp.tar.gz ages_??.??_'+version+'.lp', /sh
       endif
       return
    endif
    
    if (del_chisq lt 10.0) then chistr = $
      '0'+string(del_chisq,format='(F4.2)') else $
        chistr = string(del_chisq,format='(F5.2)') 

    gridfile = getenv('PRIMUS_DATA')+'/ages/'+$
      'ages_chi2grid_'+version+'.fits.gz'
    basis_outfile = getenv('PRIMUS_DATA')+'/ages/'+$
      'ages_basis_'+chistr+'_'+version+'.fits'
    final_basis_outfile = getenv('PRIMUS_DATA')+'/template_basis/'+$
      'ages_basis_'+chistr+'_'+version+'_'+$
      repstr(strmid(hogg_iso_date(),2,8),'-','')+'.fits.gz'
    if file_test(basis_outfile,/reg) and $
      (keyword_set(clobber) eq 0) then begin
       splog, 'File '+basis_outfile+' exists; use /CLOBBER'
       return
    endif

    lp_outfile = getenv('PRIMUS_DATA')+'/ages/'+$
      'ages_'+chistr+'_'+version+'.lp'
    if keyword_set(hogg) then begin
       cplex_outfile = lp_outfile+'.output'
    endif else begin
       cplex_outfile = getenv('PRIMUS_DATA')+'/ages/'+$
         'ages_'+chistr+'_'+version+'.cplex.out'
;        'ages_'+chistr+'_'+version+'.binary.cplex.out'
    endelse

; read the parent sample of AGES spectra (see
; BUILD_PRIMUS_AGES_PARENT) 
    if (n_elements(info) eq 0) or (n_elements(ages) eq 0) then $
      ages = primus_readin_parenttemplate(version=version,$
      templatepar=info)
    nobj = n_elements(info)
    
; read the output from BUILD_PRIMUS_LP_CHI2GRID
    splog, 'Reading '+gridfile
    gridinfo = mrdfits(gridfile,1)
    grid = gridinfo.chi2grid
    
; the first index in GRID corresponds to the data and the second index
; corresponds to the models    
    yesnogrid = grid lt del_chisq
    foo = size(yesnogrid,/dimens)
    nspectra = foo[0]
    nmodel = foo[1]

; ##################################################
; part 1 - write out the input .LP file
    if keyword_set(lp_format) then begin
       if file_test(lp_outfile,/reg) and $
         (keyword_set(clobber) eq 0) then begin
          splog, 'File '+lp_outfile+' exists; use /CLOBBER'
          return
       endif
       lp_format, yesnogrid, lp_outfile, /binary
    endif

; ##################################################
; part 2 - solve the LP problem using Hogg's method, or send the
; LP file to Sam! 
    if keyword_set(lp_solve) and keyword_set(hogg) then begin
       cmd = 'glpsol --cpxlp '+lp_outfile+' -o '+$
         cplex_outfile+' --mipgap 0.01'
       splog, cmd
       spawn, cmd, /sh
    endif

; ##################################################
; part 3 - parse the LP output and write out the final archetypes;
; need to treat the output from CPLEX (via Sam) and GLPSOL (via Hogg)
; different 
    if keyword_set(lp_parse) then begin
       if keyword_set(hogg) then begin ; output from GLPSOL/Hogg
          grepfilename = cplex_outfile+'output.grep'
          cmd = 'grep "[0-9].a[0-9][0-9][0-9][0-9][0-9][0-9]" '+$
            cplex_outfile+' > '+grepfilename
          splog, cmd
          spawn, cmd, /sh
          readcol, grepfilename, foo, name, st, activity, format='I,A,A,F'
          indx = long(strmid(name,1))
          amplitude = fltarr(max(indx)+1)
          amplitude[indx] = activity
          rmfile, grepfilename
       endif else begin ; output from CPLEX/Sam
          splog, 'Reading '+cplex_outfile
          file = djs_readlines(cplex_outfile)
          file = file[where(strmatch(file,'*variable name*'))]
          file = repstr(repstr(repstr(repstr(repstr(file,'"',' '),$
            '<variable name=',' '),'index=',' '),'value=',' '),'/>','')
          file = strtrim(file,2)
          amplitude = fltarr(nmodel)
          for ii = 0, nmodel-1 do amplitude[ii] = $
            (strsplit(file[ii],' ',/extract))[2]
          amplitude = float(fix(amplitude)) ; note sure if this is right!
       endelse

; sort and pick through templates based on LP amplitude
       splog, 'Sorting by LP amplitude...'
       sindx = reverse(sort(amplitude))
       jj = 0L
       cover = 0L
       use = bytarr(nmodel)
       repeat begin
          subgrid = yesnogrid[*,sindx[0:jj]]
          if (jj GT 0) then subgrid = total(subgrid,2)
          oldcover = cover
          cover = round(total(subgrid GT 0))
          dcover = cover-oldcover
          if (dcover GT 0) then use[sindx[jj]] = 1
          jj = jj+1
       endrep until ((cover EQ nspectra) or (jj GE nmodel))
       useindx = where(use,nuse)
       
; sort and pick through templates based on responsibility
       splog, 'Sorting by responsibility...'
       nonzero = where(amplitude GT 0.0)
       resp = round(total(yesnogrid[*,nonzero],1))
       sindx = nonzero[reverse(sort(resp))]
       jj = 0L
       cover = 0L
       ruse = bytarr(nmodel)
       repeat begin
          subgrid = yesnogrid[*,sindx[0:jj]]
          if (jj GT 0) then subgrid = total(subgrid,2)
          oldcover = cover
          cover = round(total(subgrid GT 0))
          dcover = cover-oldcover
          if (dcover GT 0) then ruse[sindx[jj]] = 1
          jj = jj+1
       endrep until ((cover EQ nspectra) or (jj GE nmodel))
       ruseindx = where(ruse,nruse)
       
; choose and save shorter list
       splog, 'nuse:',nuse,' nruse:',nruse
       if (nruse LT nuse) then begin
          splog, 'responsibility sorting wins'
          use = ruse
       endif else begin
          splog, 'amplitude sorting wins'
       endelse
       useindx = where(use,nbasis)

; final basis set; RESP is the "responsibility" of each template
; (i.e., the number of parent galaxies that this template represents),
; RESPINDX is the corresponding index array, and WEIGHT is the final
; summed weight
       basisinfo = info[useindx]
       basis = {$
         del_chisq: 0.0, $
         basisindx: fix(useindx), $
         nbasis:      0, $
         resp:       intarr(nbasis), $
         respindx:   intarr(nspectra,nbasis)-1, $
         chi2grid:   fltarr(nspectra,nbasis)+1E6,$
         multiplicity: intarr(nspectra),$ 
         ages_weight:  fltarr(nspectra,nbasis),$
         chi2_weight:  fltarr(nspectra,nbasis),$
         weight:     fltarr(nbasis), $
         unweighted_weight:   fltarr(nbasis), $
         highwave: ages.wave, $
         highflux: ages.flux[*,useindx]}
       basis.del_chisq = del_chisq
       basis.nbasis = nbasis
       basis.resp = round(total(yesnogrid[*,useindx],1))

       for ii = 0, nbasis-1 do begin
          respindx = where(yesnogrid[*,useindx[ii]] gt 0,nresp)
          basis.respindx[0:nresp-1,ii] = respindx
          basis.chi2grid[respindx,ii] = grid[respindx,useindx[ii]]
          basis.ages_weight[respindx,ii] = info[respindx].ages_weight*$
            info[respindx].field_weight ; =W_g
       endfor

; D. Eisenstein - 2009-Dec-20:
;       
; You have ~3000 galaxies, each of which has a weight W_g.  Some ~200 of
; them were selected as archetypes; I'll enumerate those with the index a.
; You have the chi^2 computed between every pair of galaxies, in particular
; between every galaxy and every archetype, chi2_{g,a}.
; 
; The original version was to assign every galaxy to exactly one archetype
; (choosing the minimum chi2).  Then the weight for every archetype is
; the sum of all of the weights of its galaxies.
; 
; What I'm suggesting is that each galaxy in fact could be matched by
; multiple archetypes, with the chi2's as probabilities.
;        
; p_ga = exp(-0.5*chi2_ga)
; P_ga = p_ga / sum_a(p_ga)    [ so sum_a P_ga = 1 ]
; w_a = sum_g (P_ga*W_g)       [ so sum_a w_a = sum_g W_g ]

; first compute the normalized probability that each galaxy is
; represented by each archetype based on chi^2 (note that many/most of
; these probabilities will be zero, except for the handful of
; archetypes that are most similar to a given object); also store the
; number of archetypes that have a non-zero probability of
; representing each spectrum (MULTIPLICITY)
       basis.chi2_weight = exp(-0.5*basis.chi2grid)
       for jj = 0, nspectra-1 do begin
          norm = total(basis.chi2_weight[jj,*],/double)
          basis.chi2_weight[jj,*] = basis.chi2_weight[jj,*]/$
            (norm+(norm eq 0.0))*(norm ne 0.0) ; sum_a P_ga = 1
          basis.multiplicity[jj] = fix(total(basis.chi2_weight[jj,*] gt 0.0))
       endfor
;      print, basis.chi2_weight[*,100], total(basis.chi2_weight[*,100])
;      print, total(basis.chi2_weight,2) ; sum_a P_ga = 1

       basis.weight = total(basis.chi2_weight*basis.ages_weight,1)
       basis.unweighted_weight = total(basis.ages_weight,1)
;      print, total(basis.weight), total(info.ages_weight*info.field_weight)
       
; write out       
       splog, 'Writing '+basis_outfile
       mwrfits, basis, basis_outfile, /create
       mwrfits, basisinfo, basis_outfile
       spawn, 'gzip -f '+basis_outfile, /sh

; write out the final basis set, if desired
       if keyword_set(final) then spawn, 'rsync -auv '+$
         basis_outfile+'.gz '+final_basis_outfile, /sh
    endif 

; ##################################################
; add Mg II - this should have been run after
; BUILD_PRIMUS_LP_ARCHETYPES, /LP_PARSE, /FINAL
    if keyword_set(addmgii) then begin
       mgii = mrdfits(getenv('PRIMUS_DATA')+'/ages/ages_mgii.fits.gz',1)
       splog, 'Fixing the template version!!!'
       vv = 'ages_basis_10.00_v3.0_100113'
;      vv = 'ages_basis_02.00_v2.0_100112'
;      vv = get_template_version()
       archfile = getenv('PRIMUS_DATA')+'/template_basis/'+vv+'.fits.gz'
       arch = mrdfits(archfile,1)
       archinfo = mrdfits(archfile,2)

       djs_plot, archinfo.ug, archinfo.ew_oii, psym=6, $
         xtitle='^{0.1}(u-g)', ytitle='EW([O II]) (\AA)'
       djs_oplot, mgii.ugriz_absmag[0]-mgii.ugriz_absmag[1], $
         mgii.oii_3727_ew[0], psym=6, color='red'
       djs_oplot, [0.0,1.0], [15,15], color='cyan', thick=2
       djs_oplot, [1.0,1.0], [-20,15], color='cyan', thick=2

       these = where((archinfo.ug lt 1.0) and $
         (archinfo.ew_oii lt 15),nthese)

       wave = arch.highwave
       flux = arch.highflux[*,these]
       mgflux = flux

       lineew = 25.0       ; EW [Angstrom]
       linewave = 2799.495
       linesigma = 1200.0  ; [km/s]

;      debug = 1
       for jj = 0L, nthese-1L do begin

          ww = where((wave gt linewave-100.0) and $
            (wave lt linewave+100.0))
          cflux = median(flux[ww,jj])
          lineflux = lineew*cflux
          peakflux = lineflux/alog(10.0)/linewave
          
          sigma = linesigma/light/alog(10.0) ; log-lambda
          term1 = exp(-0.5*(alog10(wave)-alog10(linewave))^2.0/sigma^2.0)
          linemodel = peakflux*term1/(sqrt(2.0*!pi)*sigma)
          
          mgflux[*,jj] = flux[*,jj] + linemodel
          
          if keyword_set(debug) then begin
             djs_plot, wave, mgflux[*,jj], xr=[2700,4500], color='cyan'
             djs_oplot, wave, flux[*,jj]
             djs_oplot, wave[ww], flux[ww,jj], color='red'
             djs_oplot, linewave+[-100,100], cflux*[1,1], line=0, color='green'
             cc = get_kbrd(1)
          endif
          
       endfor

;      arch.highflux[*,these] = mgflux

       newarchinfo = [archinfo,archinfo[these]]
       newarch = {$
         del_chisq:  arch.del_chisq,$
         basisindx:  [arch.basisindx,arch.basisindx[these]],$
         nbasis:     arch.nbasis + nthese,$
         resp:       [arch.resp,arch.resp[these]],$
         respindx:   [[arch.respindx],[arch.respindx[*,these]]],$
         chi2grid:   [[arch.chi2grid],[arch.chi2grid[*,these]]],$
         multiplicity: [arch.multiplicity,arch.multiplicity[these]],$
         ages_weight:  [[arch.ages_weight],[arch.ages_weight[*,these]]],$
         chi2_weight:  [[arch.chi2_weight],[arch.chi2_weight[*,these]]],$
         weight:       [arch.weight,arch.weight[these]],$
         unweighted_weight: [arch.unweighted_weight,arch.unweighted_weight[these]],$
         highwave:          arch.highwave,$
         highflux:          [[arch.highflux],[mgflux]]}

; write out!         
         outfile = getenv('PRIMUS_DATA')+'/template_basis/'+$
           strmid(vv,0,strlen(vv)-7)+'_mgii_'+$
           repstr(strmid(hogg_iso_date(),2,8),'-','')+'.fits'

         splog, 'Writing '+outfile
         mwrfits, newarch, outfile, /create
         mwrfits, newarchinfo, outfile
         spawn, 'gzip -f '+outfile, /sh
    endif

return
end


;;; now go back through and compute the total weight
;;       for ii = 0, nbasis-1 do begin
;;          good = where(basis.respindx[ii,*] ne -1,nresp)
;;          respindx = reform(basis.respindx[ii,good])
;;          total_weight = info[respindx].ages_weight*info[respindx].field_weight ; =W_g
;;
;;          chi2weight = fltarr(nresp)+1.0
;;          for jj = 0, nresp-1 do begin
;;             these = where(respindx[jj] eq basis.respindx,nthese)
;;             if (nthese ne 0) then begin
;;;               niceprint, basis.respindx[these], basis.chi2grid[these] & print
;;                allchi2 = basis.chi2grid[these]
;;                thischi2 = basis.chi2grid[ii,jj]
;;                chi2min = min(allchi2)
;;                chi2weight[jj] = exp(-0.5*(thischi2-chi2min))/$
;;                  total(exp(-0.5*(allchi2-chi2min)),/double)
;;             endif
;;          endfor
;;          basis.weight[ii] = total(chi2weight*total_weight)/total(chi2weight)
;;          basis.unweighted_weight[ii] = total(total_weight)
;;;         splog, 'Note!!!!'
;;          basis.weight[ii] = total(total_weight)
;;
;;;         print, basis.unweighted_weight[ii], basis.weight[ii]
;;;         basis.weight[ii] = total(total_weight)
;;;         chi2weight = exp(-0.5*grid[respindx,useindx[ii]]) ; wrong!!
;;       endfor       

