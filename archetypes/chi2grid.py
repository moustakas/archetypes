#!/usr/bin/env python

"""Build the binary chi2 grid based on the parent sample.

"""
from __future__ import division, print_function

import os
import sys
import argparse
import logging
import numpy as np

from archetypes import io

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Build the chi2 grid.')
                                         
    parser.add_argument('-o','--objtype', type=str, default=None, metavar='', 
                        help='object type (ELG, LRG, STAR)') 
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='toggle on verbose output')
    parser.add_argument('--snr', type=float, default=20.0, metavar='', 
                        help='signal-to-noise ratio')
    parser.add_argument('--outfile', type=str, default='OBJTYPE-chi2grid.txt', metavar='', 
                        help='output ASCII file name')

    args = parser.parse_args()
    if args.objtype is None:
        parser.print_help()
        sys.exit(1)

    # Set the debugging level
    if args.verbose:
        lvl = logging.DEBUG
    else:
        lvl = logging.INFO
    logging.basicConfig(format='%(message)s',level=lvl,stream=sys.stdout)
    log = logging.getLogger('__name__')

    log.info('Building chi2 grid for {}s.'.format(objtype))

    # Set default output file name.
    ltype = objtype.lower()
    if args.outfile:
        outfile = args.outfile
        if outfile == 'OBJTYPE-chi2grid.txt':
            outfile = ltype+'-chi2grid.txt'
    else: 
        outfile = ltype+'-chi2grid.txt'
    
    # Read the parent sample of spectra.
    key = 'DESI_BASIS_TEMPLATES'
    if key not in os.environ:
        log.fatal('Required ${} environment variable not set'.format(key))
        raise EnvironmentError

    objpath = os.getenv(key)

    objfile_wild = glob(os.path.join(objpath,ltype+'_templates_*.fits'))
    nfile = len(objfile)

    if nfile>0:
        objfile_latest = objfile[nfile-1] # latest version
        if os.path.isfile(objfile_latest):
            log.info('Reading {}'.format(objfile_latest))
        else: 
            log.error('Parent spectra file {} not found'.format(objfile_latest))
            raise IOError()
    else:
        log.error('Parent spectra file {} not found'.format(objfile_wild))
        raise IOError()

    flux, hdr = fits.getdata(objfile_latest, 0, header=True)
    meta = Table(fits.getdata(objfile_latest, 1))
    #wave = 10**(hdr['CRVAL1'] + np.arange(hdr['NAXIS1'])*hdr['CDELT1'])
    wave = fits.getdata(objfile_latest, 2)

    npix = flux.shape[1]
    ntemp = flux.shape[0]

    # Assign each spectrum a constant S/N~20 between 6800 and 7200 A.
    ivar = np.zeros_like(flux)

    waverange = np.where(((wave>6800) & (wave<7200))*1)
    for ib in range(ntemp):
        scale = np.median(flux[ib,waverange])
        scale = flux[ib,waverange]/args.snr
        ivar[ib,:] = (args.snr/flux[ib,:])^2


    waverange = where((wave1 gt 6800.0) and (wave1 lt 7200.0))
    fopt1 = fopt_primus*0.0
    ivar1 = fopt_primus*0.0
    for ii = 0, nobj-1 do begin 
       fopt1[*,ii] = fopt_primus[*,ii]/median(fopt_primus[waverange,ii])
       snr = mean(fopt1[waverange,ii] * sqrt(ivar_primus[waverange]))
       ivar1[*,ii] = ivar_primus * (20.0 / snr)^2
    endfor


    # Build the chi2 grid
    chi2grid = np.zeros_like(flux)

if __name__ == '__main__':
    main()
