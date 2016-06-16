def read_parent():

    from astropy.io import fits
    from astropy.table import Table

# Read the parent sample of spectra.
    key = 'DESI_BASIS_TEMPLATES'
#    if key not in os.environ:
#        log.fatal('Required ${} environment variable not set'.format(key))
#        raise EnvironmentError

    objpath = os.getenv(key)
    objtype = 'ELG'
    ltype = objtype.lower()

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

    return wave1, flux1, meta

