#!/usr/bin/env python

"""Stack the SEQUELS LRG spectra.

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

from speclite import redshift, 

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Build the chi2 grid.')
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='toggle on verbose output')
    parser.add_argument('--zmin', type=float, default=0.5, metavar='', 
                        help='minimum redshift')
    parser.add_argument('--zmax', type=float, default=1.0, metavar='', 
                        help='maximum redshift')

    args = parser.parse_args()

    # Set the debugging level
    if args.verbose:
        lvl = logging.DEBUG
    else:
        lvl = logging.INFO
    logging.basicConfig(format='%(message)s',level=lvl,stream=sys.stdout)
    log = logging.getLogger('__name__')

    topdir = os.getenv('SEQUELS_DIR')
    outdir = os.path.join(topdir,'stacks')

    # Loop on all the plates
    plates = ['7280']
    for pl in plates:
        platefile = os.path.join(topdir,pl,'spPlate-'+pl+'-56709.fits')
        zbestfile = os.path.join(topdir,pl,'spZbest-'+pl+'-56709.fits')

        zbest = fits.getdata(zbestfile,1)
        plug = fits.getdata(platefile,5)
        flux, hdr = fits.getdata(platefile,0,header=True)
        wave = 10**(hdr['CRVAL1'] + np.arange(hdr['NAXIS1'])*hdr['CD1_1'])
        ivar = fits.getdata(platefile,1)

        # Pick out the "good" LRGs
        iz_wise = (plug['EBOSS_TARGET0'] & 2)>0
        ri_wise = (plug['EBOSS_TARGET0'] & 4)>0
        zcut = (zbest['Z']>args.zmin)*(zbest['Z']<args.zmax)*(zbest['RCHI2DIFF_NOQSO']>0.005)
        lrg = np.where(np.logical_and(np.logical_or(iz_wise,ri_wise),zcut))[0]

        for ib in range(len([0,1])):
            zflux = flux[ib,:]
            rules = [dict(name='wlen', exponent=+1
            
        



if __name__ == '__main__':
    main()
