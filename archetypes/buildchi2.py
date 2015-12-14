#!/usr/bin/env python

"""Build the chi2 grid from the parent sample for a given object type.

"""
from __future__ import division, print_function

import os
import sys
import logging
import argparse
import numpy as np

from archetypes.io import read_parent

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Build the chi2 matrix.')
    parser.add_argument('-o','--objtype', type=str, default=None, metavar='', 
                        help='object type (ELG, LRG, STAR)') 
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='toggle on verbose output')
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
    log.info('Building chi2 grid for {}s.'.format(objtype))

    # Set default output file name.
    if args.outfile:
        outfile = args.outfile
        if outfile == 'OBJTYPE-chi2grid.txt':
            outfile = objtype.lower()+'-chi2grid.txt'
    else: 
        outfile = objtype.lower()+'-chi2grid.txt'

    # Read the parent sample.
    ff, ww, mm = read_parent('ELG')

    flux1, hdr = fits.getdata(objfile_latest, 0, header=True)
    meta = Table(fits.getdata(objfile_latest, 1))
    #wave = 10**(hdr['CRVAL1'] + np.arange(hdr['NAXIS1'])*hdr['CDELT1'])
    wave1 = fits.getdata(objfile_latest, 2)

    # Restrict the wavelength range to something reasonable.
    wavecut = np.where(((wave1>3500) & (wave1<1E4))*1)[0]
    #wavecut = np.where(((wave1>1200) & (wave1<5E4))*1)[0]
    flux = flux1[:,wavecut]
    wave = wave1[wavecut]

    keep = np.where((meta['D4000']>1.4)*1 & (meta['D4000']<2)*1)[0]
    flux = flux[keep,:] # testing!
    #flux = flux[:40,:] # testing!
    flux = flux[np.argsort(meta['D4000'][keep]),:]
    
    npix = flux.shape[1]
    ntemp = flux.shape[0]
    print(ntemp, npix)

    #ivar = np.ones_like(flux)
    #ivar = (1.0/flux)**2

    prec = 0.1 # desired precision 

    # Normalize to the median flux around 6800-7200 A.
    waverange = np.where(((wave>6800) & (wave<7200))*1)
    for ib in range(ntemp):
        #print(ib)
        flux[ib,:] /= np.median(flux[ib,waverange])
        #ivar[ib,:] = (10/flux[ib,:])**2

    # Testing!
    if args.verbose:
        import matplotlib.pyplot as plt
        for ib in range(ntemp):
            plt.plot(wave,flux[ib,:])
        plt.show()

    # Build the chi2 grid
    # The next line assumes that all amp=1 and all ivar=1!! Hack
    print('Computing chi^2')
    chi2grid = np.sum((flux[:,np.newaxis,:] - flux[np.newaxis,:,:])**2, axis=2)

    ## Build the chi2 grid
    #chi2grid = np.zeros([ntemp,ntemp])
    #for ii in range(ntemp):
    #    print(ii)
    #    for jj in range(ii+1,ntemp):
    #        #amp = np.sum(ivar[ii,:]*flux[ii,:]*flux[jj,:])/np.sum(ivar[ii,:]*flux[jj,:]**2)
    #        #chi2grid[jj,ii] = np.sum(ivar[ii,:]*(flux[ii,:]-amp*flux[jj,:])**2)
    #        chi2grid[jj,ii] = np.sum((flux[ii,:]-flux[jj,:])**2)
    #        chi2grid[ii,jj] = chi2grid[jj,ii]

    #chi2grid[chi2grid>0] -= chi2grid[chi2grid>0].min()

    matrix = (chi2grid<(npix*prec**2))*1
    #print(chi2grid)
    print(matrix)
    print(np.sum(matrix,axis=0), np.sum(matrix,axis=1))
    np.savetxt(outfile,chi2grid)

    #outfile = '{}-{:.1f}.lp'.format(ltype,args.chi2cut)
    outfile = '{}.lp'.format(ltype)
    np.savetxt(outfile,matrix,fmt='%g')

if __name__ == '__main__':
    main()
