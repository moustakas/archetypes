#!/usr/bin/env python

"""Build the chi2 grid from the parent sample for a given object type.

"""
from __future__ import division, print_function

import os
import sys
import logging
import argparse
import numpy as np

from ioarchetypes import read_parent

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Build the chi2 matrix.')
    parser.add_argument('-o','--objtype', type=str, default=None, metavar='', 
                        help='object type (ELG, LRG, STAR)') 
    parser.add_argument('--minwave', type=float, default=1500.0, 
                        help='minimum wavelength [Angstrom]')
    parser.add_argument('--maxwave', type=float, default=50000.0,
                        help='maximum wavelength [Angstrom]')
    parser.add_argument('-c','--chunksize', type=long, default=500,
                        help='number of spectra to include in each chunk')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='toggle on verbose output')
    parser.add_argument('-t', '--test', action='store_true', 
                        help='test the code by running on a mini dataset')
    parser.add_argument('--outfile', type=str, default='OBJTYPE-chi2.txt', metavar='', 
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
    log.info('Building chi2 for {}s.'.format(objtype.upper()))

    # Read the parent sample and restrict the wavelength range.
    flux, wave, meta = read_parent('ELG')
    nspec, npix = flux.shape
    
    log.info('Restricting the wavelength range to {}, {}'.format(args.minwave,args.maxwave))
    trim = (wave>args.minwave)*(wave<args.maxwave)
    flux = flux[:,trim]
    wave = wave[trim]
    nspec, npix = flux.shape

    log.info('Found {} spectra, each with {} pixels.'.format(nspec,npix))

    # Reduce the sample for testing purposes.
    if args.test:
        #keep = np.arange(0,5)
        keep = (meta['D4000']>1.4)*(meta['D4000']<2)
        log.info('Keeping {} spectra in the test sample.'.format(np.sum(keep*1)))
        flux = flux[keep,:]
    nspec, npix = flux.shape

    # Build chi2 in chunks because the array operations can be very large.
    chi2 = np.zeros((nspec,nspec))

    chunksize = np.min((args.chunksize,nspec))
    nchunk = np.ceil(nspec/chunksize).astype(long)
    log.info('Computing chi2 in {} chunk(s).'.format(nchunk))

    for ichunk in range(nchunk):
        log.info('\r  Working on chunk {:d}/{:d}'.format(ichunk+1,nchunk))
        
        i1 = ichunk*chunksize
        i2 = np.min(((ichunk+1)*chunksize,nspec))
        fdata = flux[i1:i2,None,:]
        fmodel = flux[None,i1:i2,:]

        # See equations 7-9 in Benitez+00
        fdd = np.sum(fdata*fdata, axis=2)
        fmm = np.sum(fmodel*fmodel, axis=2)
        fdm = np.sum(fdata*fmodel, axis=2)
        amp = np.divide(fdm,fmm)
        chi2[i1:i2, i1:i2] = fdd - amp*fdm

    ## This brute-force loop demonstrates that the we're computing chi^2
    ## correctly.
    chi2test = np.zeros((nspec,nspec))
    for ii in range(nspec):
        #for jj in range(nspec):
        for jj in range(ii+1,nspec):
            amp = np.sum(flux[ii,:]*flux[jj,:])/np.sum(flux[jj,:]**2)
            chi2test[ii,jj] = np.sum((flux[ii,:]-amp*flux[jj,:])**2)
    print(chi2test/chi2)

    # Write out
    outfile = os.path.join(os.getenv('ARCHETYPES_DIR'),objtype.lower()+'-chi2.txt')
    log.info('Writing {}'.format(outfile))
    np.savetxt(outfile,chi2)

    # Plot some spectra, for fun.
    if args.verbose:
        import matplotlib.pyplot as plt
        for ib in range(nspec):
            plt.plot(wave,flux[ib,:])
        plt.show()

if __name__ == '__main__':
    main()
