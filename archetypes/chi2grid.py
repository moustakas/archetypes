#!/usr/bin/env python

"""Build the chi2 grid based on the parent sample.

"""
from __future__ import division, print_function

import os
import sys
import logging
import argparse
import numpy as np
from glob import glob
from astropy.io import fits
from astropy.table import Table

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

    objtype = args.objtype
    ltype = objtype.lower()
    log.info('Building chi2 grid for {}s.'.format(objtype))

    # Set default output file name.
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

    objfile_wild = os.path.join(objpath,ltype+'_templates_*.fits')
    objfile = glob(objfile_wild)
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

    flux1, hdr = fits.getdata(objfile_latest, 0, header=True)
    meta = Table(fits.getdata(objfile_latest, 1))
    #wave = 10**(hdr['CRVAL1'] + np.arange(hdr['NAXIS1'])*hdr['CDELT1'])
    wave1 = fits.getdata(objfile_latest, 2)

    # Restrict the wavelength range to something reasonable.
    wavecut = np.where(((wave1>3500) & (wave1<1E4))*1)[0]
    #wavecut = np.where(((wave1>1200) & (wave1<5E4))*1)[0]
    flux = flux1[:,wavecut]
    wave = wave1[wavecut]

    flux = flux[:40,:] # testing!
    npix = flux.shape[1]
    ntemp = flux.shape[0]

    ivar = np.ones_like(flux)
    #ivar = (20.0/flux)**2

    # Normalize to the median flux around 6800-7200 A.
    waverange = np.where(((wave>6800) & (wave<7200))*1)
    for ib in range(ntemp):
        flux[ib,:] /= np.median(flux[ib,waverange])
        ivar[ib,:] = (10/flux[ib,:])**2

    # Testing!
    if args.verbose:
        import matplotlib.pyplot as plt
        for ib in range(ntemp):
            plt.plot(wave,flux[ib,:])
        plt.show()

    # Build the chi2 grid
    chi2grid = np.zeros([ntemp,ntemp])
    for ii in range(ntemp):
        for jj in range(ntemp):
            chi2grid[jj,ii] = np.sum(ivar[ii,:]*(flux[ii,:]-flux[jj,:])**2)
    chi2grid[chi2grid>0] -= chi2grid[chi2grid>0].min()

    print(chi2grid)
    print((chi2grid<2)*1)
    np.savetxt(outfile,chi2grid)

if __name__ == '__main__':
    main()
