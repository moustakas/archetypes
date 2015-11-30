;+
; NAME:
;   QAPLOT_LP_ARCHETYPES
;
; PURPOSE:
;   Build some handy QAplots of the LP archetypes.
;
; INPUTS: 
;
; OPTIONAL INPUTS: 
;
; KEYWORD PARAMETERS: 
;
; OUTPUTS: 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Dec 16, UCSD - written
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

pro qaplot_lp_archetypes, del_chisq, version=version

    common qaplot_archetypes, parent, parentinfo, gandalf, grid

    if (n_elements(version) eq 0) then version = $
      primus_parenttemplate_defaultversion()
    if (n_elements(del_chisq) eq 0) then del_chisq = 2.0

;   del_chisq = [0.5,0.7,1.0,1.5,2.0,2.5,3.0,4.0,5.0,7.0,10.0]
; call this routine recursively    
    ndel_chisq = n_elements(del_chisq)
    if (ndel_chisq gt 1) then begin
       for ii = 0, ndel_chisq-1 do qaplot_lp_archetypes, del_chisq[ii]
       return
    endif

    if (del_chisq lt 10.0) then chistr = $
      '0'+string(del_chisq,format='(F4.2)') else $
        chistr = string(del_chisq,format='(F5.2)') 

    agespath = getenv('PRIMUS_DATA')+'/ages/'
    
; read in all the various structures
    if (n_elements(parent) eq 0) or (n_elements(parentinfo) eq 0) then $
      parent = primus_readin_parenttemplate(version=version,$
      templatepar=parentinfo)
    
    if (n_elements(grid) eq 0) then begin
       gridfile = agespath+'ages_chi2grid_'+version+'.fits.gz'
       splog, 'Reading '+gridfile
       grid = mrdfits(gridfile,1)
    endif

    if (n_elements(gandalf) eq 0) then begin
       gandalffile = agespath+'gandalf_parent_'+version+'.fits.gz'
;      ancillaryfile = agespath+'ancillary_parent_'+version+'.fits.gz'
       splog, 'Reading '+gandalffile
       gandalf = mrdfits(gandalffile,1)
;      splog, 'Reading '+ancillaryfile
;      ancillary = mrdfits(ancillaryfile,1)
    endif

; read the basis set    
    archfile = agespath+'ages_basis_'+$
      chistr+'_'+version+'.fits.gz'
    splog, 'Reading '+archfile
    arch = mrdfits(archfile,1)
    archinfo = mrdfits(archfile,2)

    archindx = arch.basisindx

; make the various plots    
    qafile = agespath+'qaplots/qaplot_lp_archetypes_'+chistr+'_'+version+'.ps'
    im_plotconfig, 0, pos, psfile=qafile, xmargin=[1.3,0.4], $
      width=6.8

; --------------------------------------------------
; apparent magnitude versus redshift
    djs_plot, parentinfo.z, parentinfo.i, $
      psym=6, symsize=0.3, yrange=[14,20.2], xrange=[0,0.4], $
      xtitle='Redshift', ytitle='I_{NDWFS} (mag)', $
      position=pos, xsty=1, ysty=1
    symsize1 = 6.0*sqrt(arch.weight/max(arch.weight))
    djs_oplot, archinfo.z, archinfo.i, psym=symcat(6,thick=5), $
      color='red', symsize=symsize1

; --------------------------------------------------
; absolute magnitude versus redshift
    djs_plot, parentinfo.z, parentinfo.mg, $
      psym=6, symsize=0.3, yrange=[-16,-23.0], xrange=[0,0.4], $
      xtitle='Redshift', ytitle='M_{0.1g}', $
      position=pos, xsty=1, ysty=1
    symsize1 = 6.0*sqrt(arch.weight/max(arch.weight))
    djs_oplot, archinfo.z, archinfo.mg, psym=symcat(6,thick=5), $
      color='red', symsize=symsize1

; --------------------------------------------------
; color-magnitude diagram    
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xrange=[-16,-23.0], yrange=[0.1,2.3], xtitle='M_{0.1g}', $
      ytitle='^{0.1}(u - g)', position=pos
    djs_oplot, parentinfo.mg, parentinfo.ug, psym=symcat(16), symsize=0.5
    symsize1 = 6.0*sqrt(arch.weight/max(arch.weight))
    djs_oplot, archinfo.mg, archinfo.ug, psym=symcat(6,thick=5), $
      color='red', symsize=symsize1

; --------------------------------------------------
; color-magnitude diagram - dust-corrected
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xrange=[-16,-23.0], yrange=[0.1,2.3], xtitle='M_{0.1g, cor}', $
      ytitle='^{0.1}(u - g)_{cor}', position=pos
    djs_oplot, parentinfo.mg, parentinfo.ug_dered, psym=symcat(16), symsize=0.5
    symsize1 = 6.0*sqrt(arch.weight/max(arch.weight))
    djs_oplot, archinfo.mg, archinfo.ug_dered, psym=symcat(6,thick=5), $
      color='red', symsize=symsize1

; --------------------------------------------------
; D(4000) vs EW([OII])
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xrange=[0.7,2.5], yrange=[0.5,200], xtitle='D_{n}(4000)', $
      ytitle='EW([O II]) (\AA)', position=pos, /ylog
    djs_oplot, parentinfo.d4000, parentinfo.ew_oii>0.8, psym=symcat(16), symsize=0.5
    symsize1 = 6.0*sqrt(arch.weight/max(arch.weight))
    djs_oplot, archinfo.d4000, archinfo.ew_oii>0.8, psym=symcat(6,thick=5), $
      color='red', symsize=symsize1

; --------------------------------------------------
; BPT diagram
;   arch_gandalf = gandalf[arch.basisindx]
    class = iclassification(gandalf,ratios=ratios,snrcut_class=3.0)
    
    good1 = where((ratios.nii_ha gt -900.0) and (ratios.oiii_hb gt -900.0))
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xtitle=textoidl('log ([N II] \lambda6584)'), $
      ytitle=textoidl('log ([O III] \lambda5007)'), $
      xrange=[-1.6,0.5], yrange=[-1.1,1.2], position=pos
    djs_oplot, ratios[good1].nii_ha, ratios[good1].oiii_hb, $
      psym=symcat(16), symsize=0.5
;   hogg_scatterplot, ratios[good1].nii_ha, ratios[good1].oiii_hb, $
;     xtitle=textoidl('log ([N II] \lambda6584)'), $
;     ytitle=textoidl('log ([O III] \lambda5007)'), $
;     xrange=[-2,1], yrange=[-1.5,1.5], position=pos
    good2 = where((ratios[archindx].nii_ha gt -900.0) and $
      (ratios[archindx].oiii_hb gt -900.0))
    symsize1 = 6*sqrt(arch.weight[good2]/max(arch.weight[good2]))
    djs_oplot, ratios[archindx[good2]].nii_ha, ratios[archindx[good2]].oiii_hb, $
      psym=symcat(6,thick=5), color='red', symsize=symsize1

; --------------------------------------------------
; luminosity vs Balmer decrement
    dust = iunred_linedust(gandalf,snrcut=3.0)
    good1 = where(dust.ebv_hahb_err gt 0)
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xrange=[-16,-23.0], yrange=[-0.05,1.5], xtitle='M_{0.1g}', $
      ytitle='E(B-V) from H\alpha/H\beta (mag)', position=pos
    djs_oplot, parentinfo[good1].mg, dust[good1].ebv_hahb, $
      psym=symcat(16), symsize=0.5
    good2 = where(dust[archindx].ebv_hahb_err gt 0)
    symsize1 = 6*sqrt(arch.weight[good2]/max(arch.weight[good2]))
    djs_oplot, archinfo[good2].mg, dust[archindx[good2]].ebv_hahb, $
      psym=symcat(6,thick=5), color='red', symsize=symsize1

; --------------------------------------------------
; luminosity vs continuum reddening
    djs_plot, [0], [0], /nodata, xsty=1, ysty=1, $
      xrange=[-16,-23.0], yrange=[-0.05,0.75], xtitle='M_{0.1g}', $
      ytitle='E(B-V) of the Continuum (mag)', position=pos
    djs_oplot, parentinfo.mg, parentinfo.ebv, $
      psym=symcat(16), symsize=0.5
    symsize1 = 6*sqrt(arch.weight/max(arch.weight))
    djs_oplot, archinfo.mg, archinfo.ebv, $
      psym=symcat(6,thick=5), color='red', symsize=symsize1

; close the plot    
    im_plotconfig, psfile=qafile, /psclose, /gzip
    spawn, 'rsync -auv '+qafile+'.gz ~', /sh

; plot the archetypes themselves as well as the objects that each
; archetype represents, sorted by total weight
    qafile = agespath+'qaplots/qaplot_lp_archetypes_'+chistr+'_'+version+'_sed.ps'
    im_plotconfig, 8, pos, psfile=qafile, charsize=1.6
;   pos = im_getposition(nx=1,ny=2,xmargin=[1.1,0.4],ymargin=[0.4,1.1],$
;     xspace=0.0,yspace=0.0,width=7.0,height=[4.0,3.0],xpage=8.5,ypage=8.5)

    srt = reverse(sort(arch.weight))
    normwave = (k_lambda_eff(filterlist='ndwfs_I.par'))[0]
    highwave = arch.highwave
    wave = parent.wave
    npix = n_elements(wave)
    inrange = where((wave gt normwave-150.0) and $
      (wave lt normwave+150.0))
    xrange = [2500.0,1.9E4]
    plotrange = where((wave gt xrange[0]) and $
      (wave lt xrange[1]))
;   xrange = minmax(wave)
;   for jj = 25, 35 do begin
    for jj = 0, arch.nbasis-1 do begin
; first deal with the full range of galaxies that this archetype
; represents so that we can get the yrange
       good = where(arch.respindx[*,srt[jj]] ne -1,nthese)
       these = reform(arch.respindx[good,srt[jj]])
       flux = -2.5*alog10(parent.flux[*,these])
       for kk = 0, nthese-1 do flux[*,kk] = flux[*,kk] - $
         djs_median(flux[inrange,kk])
       if (nthese gt 1L) then begin
;         ugmin = min(parentinfo[these].ug,minindx,$
;           subscript_max=maxindx)
;         fmin = flux[*,minindx] ; reddest SED
;         fmax = flux[*,maxindx] ; bluest SED
          fmin = min(flux,dim=2)
          fmax = max(flux,dim=2)
       endif       

       yrange = [max(flux[plotrange,*]),min(flux[plotrange,*])]
; plot the archetype itself in the top panel       
       highflux = -2.5*alog10(arch.highflux[*,srt[jj]])
       highflux = highflux - djs_median(highflux[inrange])

;       if jj eq 9 then begin
;          dfpsclose
;          djs_plot, highwave, highflux, xrange=xrange, yrange=yrange
;          djs_oplot, wave, fmin, color='red'
;          djs_oplot, wave, fmax, color='blue'
;          print, 
;stop
;       endif
; now make a grey shaded region of the full range
       djs_plot, [0], [0], /nodata, xsty=3, ysty=3, $
         position=pos, /xlog, xrange=xrange, yrange=yrange, $
         xtitle='Rest Wavelength (\AA)', ytitle='Relative Flux (AB mag)'
;      xyouts, pos[0,0]-0.1, pos[1,0], 'Relative Flux (AB mag)', align=0.5, $
;        /norm, orientation=90.0
       if (nthese gt 1L) then begin
;         polyfill, [wave[0],wave,wave[npix-1],wave], [fmin[0],fmax,fmin[npix-1],fmin], $
;           /data, /fill, color=fsc_color('medium grey',100), noclip=0
          polyfill, [wave[0],wave,wave[npix-1],wave], [fmin[0],fmax,fmin[npix-1],fmin], $
            /data, /line_fill, orientation=135, spacing=0.1, linestyle=0, $
            color=fsc_color('grey',100), noclip=0, thick=1
          polyfill, [wave[0],wave,wave[npix-1],wave], [fmin[0],fmax,fmin[npix-1],fmin], $
            /data, /line_fill, orientation=45, spacing=0.1, linestyle=0, $
            color=fsc_color('grey',100), noclip=0, thick=1
       endif
; overlay the archetype and make a legend
       djs_oplot, highwave, highflux
       label1 = [$
         '#'+string(srt[jj]+1,format='(I0)'),$
         'N_{resp}='+string(arch.resp[srt[jj]],format='(I0)'),$
         'Weight='+strtrim(string(arch.weight[srt[jj]],format='(F12.2)'),2)]
       label2 = [$
         '<Age>='+strtrim(string(archinfo[srt[jj]].age,format='(F12.2)'),2)+' Gyr',$
         'E(B-V)='+strtrim(string(archinfo[srt[jj]].ebv,format='(F12.3)'),2),$
         '^{0.1}(u-g)='+strtrim(string(archinfo[srt[jj]].ug,format='(F12.3)'),2),$
         '^{0.1}(u-g)_{cor}='+strtrim(string(archinfo[srt[jj]].ug_dered,format='(F12.3)'),2)]
       im_legend, label1, /right, /top, box=0, charsize=1.3
       im_legend, label2, /right, /bottom, box=0, charsize=1.3
;if jj eq 18 then stop
    endfor
    
; close the plot    
    im_plotconfig, psfile=qafile, /psclose, /gzip
    spawn, 'rsync -auv '+qafile+'.gz ~', /sh

stop    
    
return
end
    
