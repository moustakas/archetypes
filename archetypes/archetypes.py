#!/usr/bin/env python

"""Get ready to run CPLEX.  Code operates on the output of chi2grid.py. 

"""
from __future__ import division, print_function

import os
import sys
import logging
import argparse
import numpy as np

from glob import glob

import cplex
import matplotlib.pyplot as plt

log = logging.getLogger('__name__')

def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description='Build the files needed by CPLEX.')
    parser.add_argument('-o','--objtype', type=str, default='ELG', metavar='', 
                        help='object type (ELG, LRG, STAR)') 
    parser.add_argument('-v', '--verbose', action='store_true', 
                        help='toggle on verbose output')
    parser.add_argument('--chi2cut', type=float, default=1000.0, metavar='', 
                        help='chi^2 threshold')
    parser.add_argument('--chi2file', type=str, default='OBJTYPE-chi2grid.txt', metavar='', 
                        help='OUTFILE from chi2grid.py')

    args = parser.parse_args()
    #if args.objtype is None:
    #    parser.print_help()
    #    sys.exit(1)

    # Set the debugging level
    if args.verbose:
        lvl = logging.DEBUG
    else:
        lvl = logging.INFO
    logging.basicConfig(format='%(message)s',level=lvl,stream=sys.stdout)
    log = logging.getLogger('__name__')

    objtype = args.objtype
    ltype = objtype.lower()
    log.info('Solving using CPLEX for {}s.'.format(objtype))

    # Set default output file name.
    outfile = '{}-{:.1f}.lp'.format(ltype,args.chi2cut)

    # Read the output of chi2grid.py
    if args.chi2file:
        chi2file = args.chi2file
        if chi2file == 'OBJTYPE-chi2grid.txt':
            chi2file = ltype+'-chi2grid.txt'
    else: 
        chi2file = ltype+'-chi2grid.txt'
    #chi2grid = np.loadtxt(chi2file)
    chi2grid = np.loadtxt(chi2file)

    fix me --------

    prec = 0.1 # desired precision 
    matrix = (chi2grid<(npix*prec**2))*1
    #print(chi2grid)
    print(matrix)
    print(np.sum(matrix,axis=0), np.sum(matrix,axis=1))
    np.savetxt(outfile,chi2grid)

    #outfile = '{}-{:.1f}.lp'.format(ltype,args.chi2cut)
    outfile = '{}.lp'.format(ltype)
    np.savetxt(outfile,matrix,fmt='%g')

    -------------

    ntemp_objective = 30

    ndata, ntemp = chi2grid.shape
    logchi2cuts = np.arange(0.5,1.0,0.01)
    n_templates = np.zeros_like(logchi2cuts).astype("int")
    binaries = np.zeros((len(logchi2cuts), ntemp)).astype("int")
    for ii,logchi2cut in enumerate(logchi2cuts):
        binaries[ii] = gimme_template_binaries(chi2grid, 10.**logchi2cut)
        n_templates[ii] = np.sum(binaries[ii])

    this = np.argmin(np.abs(n_templates-ntemp_objective))
    indices = np.where(binaries[this,:])[0]
    thismatrix = chi2grid < 10.**logchi2cuts[this]
    responsibilities = np.sum(thismatrix[:,indices], axis=0)
    print(indices.shape, responsibilities.shape)
    for ii, rr in zip(indices, responsibilities):
        print(ii, rr)

    plt.scatter(logchi2cuts, n_templates)
    plt.show()

    sys.exit(1)

    wave, flux, meta = read_parent()
    plt.plot(meta['D4000'],meta['OII_3727_EW'],'r.',alpha=0.25)
    plt.plot(meta['D4000'][indices],meta['OII_3727_EW'][indices],'bo')
    plt.show()

    for ii in indices:
        plt.plot(wave,flux[ii])
        plt.show()

    plt.imshow(binaries.T, cmap="gray", interpolation="nearest", origin="lower")
    plt.show()

def gimme_template_binaries(chi2grid, chi2cut):

    matrix = (chi2grid < chi2cut) * 1
    # note that in principle the rows and columns of the matrix are different!
    # you might need a ".T" in here...
    ndata, ntemp = matrix.shape

    # each template contributes an amplitude to the objective function
    # each data point contributes a constraint

    my_obj      = np.ones(ntemp) # coefficients for the objective function
    my_rhs      = np.ones(ndata) # lower bounds for the constraints
    my_sense    = "G" * ndata    # need to make sure this is >=

    # in cplex terminology, rows are constraints, and columns are amplitudes
    rows, cols = np.where(matrix) # matrix is just ones and zeros!
    vals = np.ones_like(rows)

    prob = cplex.Cplex()
    prob.objective.set_sense(prob.objective.sense.minimize)
    prob.linear_constraints.add(rhs = my_rhs, senses = my_sense)
    prob.variables.add(obj = my_obj, types = "B" * ntemp)
    prob.linear_constraints.set_coefficients(zip(rows, cols, vals))
    prob.solve()
    return prob.solution.get_values()

if False:    
    print("Solution status = " , prob.solution.get_status(), ":", end=' ')
    print(prob.solution.status[prob.solution.get_status()])
    print("Solution value  = ", prob.solution.get_objective_value())
    slack = prob.solution.get_linear_slacks()
    pi    = prob.solution.get_dual_values()
    x     = prob.solution.get_values()
    dj    = prob.solution.get_reduced_costs()
    for i in range(numrows):
        print("Row %d:  Slack = %10f  Pi = %10f" % (i, slack[i], pi[i]))
    for j in range(numcols):
        print("Column %d:  Value = %10f Reduced cost = %10f" % (j, x[j], dj[j]))

    prob.write('elg.output')

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
