#!/usr/bin/env python

"""I/O routines for archetypes.

TODO (@moustakas): Remove dependence on desispec.interpolation.resample_flux 

"""
from __future__ import division, print_function

import os
import sys
import logging
import numpy as np

from desispec.interpolation import resample_flux

logging.basicConfig(format='%(message)s',level=logging.INFO,stream=sys.stdout)
log = logging.getLogger('__name__')

#- Utility function to wrap resample_flux for multiprocessing map
def _resample_flux(args):
    return resample_flux(*args)

def read_parent(objtype, outwave=None, nspec=None, infile=None):
    """Read and return the parent spectral templates for a given object type.
       Optionally returns a randomly selected subset of nspec spectra sampled at
       wavelengths outwave.

    Args:
      objtype (str): object type to read (e.g., ELG, LRG, QSO, STAR, FSTD, WD).
      outwave (numpy.array, optional): array of wavelength at which to sample 
        the spectra.
      nspec (int, optional): number of templates to return
      infile (str, optional): full path to input template file to read,
        over-riding the contents of the $ARCHETYPES_DIR environment
        variable.

    Returns:
      outflux (numpy.ndarray): Array [ntemplate,npix] of flux values [erg/s/cm2/A]. 
      outwave (numpy.ndarray): Array [npix] of wavelengths for FLUX [Angstrom].
      meta (astropy.Table): Meta-data table for each object.  The contents of this
        table varies depending on what OBJTYPE has been read.

    Raises:
      EnvironmentError: If the required $ARCHETYPES_DIR environment
        variable is not set.
      IOError: If the requisite FITS file is not found.

    """
    from glob import glob
    from astropy.io import fits
    from astropy.table import Table

    key = 'ARCHETYPES_DIR'
    if key not in os.environ:
        log.fatal('Required ${} environment variable not set'.format(key))
        raise EnvironmentError
    objpath = os.getenv(key)

    ltype = objtype.lower()
    if infile is None:
        objfile_wild = os.path.join(objpath,ltype+'_templates_*.fits')
    else:
        objfile_wild = infile
        
    objfile = glob(objfile_wild)
    nfile = len(objfile)

    if nfile>0:
        objfile_latest = objfile[nfile-1] # latest version
        if os.path.isfile(objfile_latest):
            log.info('Reading {}'.format(objfile_latest))
        else: 
            log.error('Parent spectra {} not found'.format(objfile_latest))
            raise IOError()
    else:
        log.error('Parent spectra {} not found'.format(objfile_wild))
        raise IOError()

    flux, hdr = fits.getdata(objfile_latest, 0, header=True)
    meta = Table(fits.getdata(objfile_latest, 1))
    #wave = 10**(hdr['CRVAL1'] + np.arange(hdr['NAXIS1'])*hdr['CDELT1'])
    wave = fits.getdata(objfile_latest, 2)

    # Optionally choose a random subset of spectra. There must be a fast way to
    # do this using fitsio.
    ntemplates = flux.shape[0]
    if nspec is not None:
        these = np.random.choice(np.arange(ntemplates),nspec)
        flux = flux[these,:]
        meta = meta[these]

    # Optionally resample the templates at specific wavelengths.  Use
    # multiprocessing to speed this up.
    if outwave is None:
        outflux = flux # Do I really need to copy these variables!
        outwave = wave
    else:
        args = list()
        for jj in range(nspec):
            args.append((outwave, wave, flux[jj,:]))

        ncpu = multiprocessing.cpu_count() // 2   #- avoid hyperthreading
        pool = multiprocessing.Pool(ncpu)
        outflux = pool.map(_resample_flux, args)
        outflux = np.array(outflux)    

    return outflux, outwave, meta
