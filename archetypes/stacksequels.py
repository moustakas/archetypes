#!/usr/bin/env python

"""Stack the SEQUELS LRG spectra.

From Abhi: I am attaching the file containing visual inspection of some of these
LRGs. CLASS of 3 or 4 are good spectra. I wouldn't believe CLASS 1 at all. This
classification is same as DEEP2. Some of the text files will have
classifications from 0 to 3. Only god knows why eBOSS people chose to set up the
classification this way.

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

from speclite import redshift, accumulate

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Build the chi2 grid.')
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='toggle on verbose output')
    parser.add_argument('--zmin', type=float, default=0.5, metavar='', 
                        help='minimum redshift')
    parser.add_argument('--zmax', type=float, default=1.0, metavar='', 
                        help='maximum redshift')
    parser.add_argument('--zbin', type=float, default=0.25, metavar='', 
                        help='redshift bin size')

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

    # Build the redshift bin (centers). Only works if ZBIN is an integer
    # multiple of [ZMIN,ZMAX].
    zbins = np.concatenate(([zmin],np.arange(zmin,zmax,zbin)+zbin/2,[zmax]))

    # Loop on all the plates
    plates = ['7280']
    for pl in plates:
        platefile = glob(os.path.join(topdir,pl,'spPlate-'+pl+'-?????.fits'))[0]
        zbestfile = glob(os.path.join(topdir,pl,'spZbest-'+pl+'-?????.fits'))[0]

        zall = fits.getdata(zbestfile,1)
        plug = fits.getdata(platefile,5)
        fluxall, hdr = fits.getdata(platefile,0,header=True)
        wave = 10**(hdr['CRVAL1'] + np.arange(hdr['NAXIS1'])*hdr['CD1_1'])
        ivarall = fits.getdata(platefile,1)

        # Pick out the "good" LRGs
        iz_wise = (plug['EBOSS_TARGET0'] & 2)>0
        ri_wise = (plug['EBOSS_TARGET0'] & 4)>0
        zcut = (zall['Z']>args.zmin)*(zall['Z']<args.zmax)*(zall['RCHI2DIFF_NOQSO']>0.005)
        lrg = np.where(np.logical_and(np.logical_or(iz_wise,ri_wise),zcut))[0]

        if len(lrg)>0:
            flux = fluxall[lrg,:]
            ivar = ivarall[lrg,:]
            zlrg = zall[lrg]
        
            # Assign each object to the right redshift bin.
            idx  = np.digitize(x,bins)
            #print(nbin, bins, xmin, xmax)

        
    
    




        for ib in range(len([0,1])):
            zflux = flux[lrg[ib],:]
            zivar = ivar[lrg[ib],:]
            rules = [dict(name='wave', exponent=+1, array_in=wave),
                     dict(name='flux', exponent=-1, array_in=zflux),
                     dict(name='ivar', exponent=+2, array_in=zivar)]
            out = redshift(zbest['z'][lrg[ib]],0.2,rules=rules)
            print(out.shape)
        

if __name__ == '__main__':
    main()
