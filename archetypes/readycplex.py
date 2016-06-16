#!/usr/bin/env python

"""Get ready to run CPLEX.  Code operates on the output of chi2grid.py. 

"""
from __future__ import division, print_function

import os
import sys
import logging
import argparse
import numpy as np

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Build the files needed by CPLEX.')
    parser.add_argument('-o','--objtype', type=str, default=None, metavar='', 
                        help='object type (ELG, LRG, STAR)') 
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='toggle on verbose output')
    parser.add_argument('--chi2cut', type=float, default=1000.0, metavar='', 
                        help='chi^2 threshold')
    parser.add_argument('--chi2file', type=str, default='OBJTYPE-chi2grid.txt', metavar='', 
                        help='OUTFILE from chi2grid.py')

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
    log.info('Building CPLEX files for {}s.'.format(objtype))

    # Set default output file name.
    outfile = '{}-{:.1f}.lp'.format(ltype,args.chi2cut)

    # Read the output of chi2grid.py
    if args.chi2file:
        chi2file = args.chi2file
        if chi2file == 'OBJTYPE-chi2grid.txt':
            chi2file = ltype+'-chi2grid.txt'
    else: 
        chi2file = ltype+'-chi2grid.txt'
    chi2grid = np.loadtxt(chi2file)

    matrix = (chi2grid<args.chi2cut)*1
    print(chi2grid, matrix)
    print(np.sum(matrix))

    np.savetxt(outfile,matrix,fmt='%g')

    # Write out the input .LP file
    '''
    mm = bytarr(40,40)
    openr, lun, 'elg.lp', /get_lun
    readf, lun, mm
    free_lun, lun
    lp_format, mm, 'elg.output', /binary
    '''

if __name__ == '__main__':
    main()
