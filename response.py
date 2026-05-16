import re
import numpy as np
import os, errno
import glob
try:
    import astropy.io.fits
except:
    pass
import sys

import util

# Taken from https://stackoverflow.com/questions/2545532/python-analog-of-natsort-function-sort-a-list-using-a-natural-order-algorithm
def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

def rmf_directories():
    try:
        return os.environ['RMFPATH'].split(':')
    except:
        return ['/data/legs/rpete/flight/rmfs']

def arf_directories():
    return os.environ['ARFPATH'].split(':')


def rmf_files(detector, grating, maxorder=None, archive=False):
    if archive:
        files = archive_rmf_files(detector, grating, maxorder)
    else:
        files= generic_files(rmf_directories(), '{}-'.format(detector), '.rmf', grating, maxorder)
    return files

def zeroth_arf_file(obsid, archive=False):
    if archive:
        obsid = f'{int(obsid):05d}'
        globstr = f'{util.archivedir}/[is]/{obsid}/analysis/tg/hrcf{obsid}_0th.arf'
        #sys.stderr.write(globstr+"\n")
        return glob.glob(globstr)[0]

    files = []
    for d in arf_directories():
        files.extend(glob.glob('{}/{}_0th.arf'.format(d, obsid)))
    if len(files) != 1:
        raise ValueError('{} files found for obsid {}'.format(len(files), obsid))
    thefile = files[0]
    return thefile

def archive_rmf_files(detector, grating, maxorder):
    files = { 'neg' : None, 'pos' : None }
    prefixes = { 'neg' : 'm', 'pos' : 'p' }
    detectors = { 'HRC-S' : 'hrcs', 'HRC-I' : 'hrci' }
    grating = grating.lower()

    for orders in files:
        globstr = f'{util.archivedir}/../hrc/rmfs/{detectors[detector]}_{grating}_{prefixes[orders]}[1-9]*.rmf'
        #sys.stderr.write(globstr+"\n")
        files[orders] = sorted(glob.glob(globstr), key=natural_key)[:maxorder]

    return files

def archive_garf_files(obsid, grating, maxorder):
    files = { 'neg' : None, 'pos' : None }
    prefixes = { 'neg' : 'm', 'pos' : 'p' }
    grating = grating.lower()

    obsid = f'{int(obsid):05d}'

    for orders in files:
        globstr = f'{util.archivedir}/[is]/{obsid}/analysis/tg/hrcf{obsid}_{grating}_{prefixes[orders]}[1-9]*.arf'
        #sys.stderr.write(globstr+"\n")
        files[orders] = sorted(glob.glob(globstr), key=natural_key)[:maxorder]

    return files

def garf_files(obsid, grating, maxorder=None, archive=False):
    if archive:
        return archive_garf_files(obsid, grating, maxorder)

    return generic_files(arf_directories(), '{}_'.format(obsid), '_garf.fits', grating, maxorder)

def generic_files(directories, prefix, postfix, grating, maxorder):
    files = { 'neg' : [], 'pos' : [] }
    signs = { 'neg' : '-', 'pos' : '' }

    for d in directories:
        for arm in files:
            tmp_files = glob.glob('{}/{}{}_{}[1-9]*{}'.format(d, prefix, grating, signs[arm], postfix))
            if len(tmp_files):
                if len(files[arm]):
                    raise ValueError('Found multiple directories for files with prefix={}, postfix={}'.format(prefix, postfix))
                else:
                    files[arm] = sorted(tmp_files, key=natural_key)[:maxorder]

    if not len(files['pos']) or not len(files['neg']) or len(files['neg']) is not len(files['pos']):
        raise ValueError('response.generic_files() - found the following number of files for prefix={}, postfix={}: N_neg={}, N_pos={}'.format(prefix, postfix, len(files['neg']), len(files['pos'])))

    return files

def read_rmf(filename):
    if 'rmfs' not in read_rmf.__dict__:
        read_rmf.rmfs = {}
    hdulist = astropy.io.fits.open(filename)
    hdr = hdulist['matrix'].header
    datasum = hdr['datasum']
    if datasum not in read_rmf.rmfs:
        sys.stderr.write("response.read_rmf() - reading {}\n".format(filename))
        data = hdulist['matrix'].data
        assert 0==np.size(np.where(data.field('n_grp')!=1))
        energ_lo = data.field('energ_lo')
        energ_hi = data.field('energ_hi')
        f_chan = data.field('f_chan')
        n_chan = data.field('n_chan')
        matrix = data.field('matrix')
        n = f_chan.shape[0]
        resp = np.zeros(n)
        for i in range(n):
            resp[f_chan[i][0]-1:f_chan[i][0]+n_chan[i][0]-1] += matrix[i]

        read_rmf.rmfs[datasum] = { 'bin_lo':e2w(energ_hi), 'bin_hi':e2w(energ_lo), 'response':resp }
    hdulist.close()
    return read_rmf.rmfs[datasum]['bin_lo'], read_rmf.rmfs[datasum]['bin_hi'], read_rmf.rmfs[datasum]['response']

def read_arf(filename):

    sys.stderr.write("response.read_arf() - reading {}\n".format(filename))

    hdulist = astropy.io.fits.open(filename)
    header = hdulist['specresp'].header
    data = hdulist['specresp'].data
    hdulist.close()

    energ_lo = data.field('energ_lo')
    energ_hi = data.field('energ_hi')
    specresp = data.field('specresp')

    return e2w(energ_hi), e2w(energ_lo), specresp

def get_rmfs(detector, grating, maxorder=None, archive=False):
    if 'rmfs' not in get_rmfs.__dict__ or True:
        files = rmf_files(detector, grating, maxorder, archive=archive)
        rmfs = {}
        bins_set = False

        for arm in files:
            bin_lo_, bin_hi_, response = read_rmf(files[arm][0])

            rmfs[arm] = np.zeros((len(files[arm]), response.size))
            rmfs[arm][0] = response

            if not bins_set:
                bin_lo = np.zeros(rmfs[arm].shape)
                bin_hi = bin_lo.copy()
                bin_lo[0] = bin_lo_
                bin_hi[0] = bin_hi_

            for i in range(1, len(files[arm])):
                bin_lo_, bin_hi_, rmfs[arm][i] = read_rmf(files[arm][i])
                if not bins_set:
                    bin_lo[i] = bin_lo_
                    bin_hi[i] = bin_hi_

            bins_set = True

        get_rmfs.bin_lo = bin_lo
        get_rmfs.bin_hi = bin_hi
        get_rmfs.rmfs = rmfs

    return get_rmfs.bin_lo, get_rmfs.bin_hi, get_rmfs.rmfs

def get_garfs(obsid, grating, maxorder=None, archive=False):
    files = garf_files(obsid, grating, maxorder, archive=archive)
    garfs = {}
    bins_set = False
    
    for order in files:
        bin_lo_, bin_hi_, specresp = read_arf(files[order][0])
        garfs[order] = np.zeros((len(files[order]), specresp.size))
        garfs[order][0] = specresp

        if not bins_set:
            bin_lo = np.zeros(garfs[order].shape)
            bin_hi = bin_lo.copy()
            bin_lo[0] = bin_lo_
            bin_hi[0] = bin_hi_

        for i in range(1, len(files[order])):
            bin_lo_, bin_hi_, garfs[order][i] = read_arf(files[order][i])
            if not bins_set:
                bin_lo[i] = bin_lo_
                bin_hi[i] = bin_hi_

        bins_set=True

    return bin_lo, bin_hi, garfs

def read_header(filename):
    hdulist = astropy.io.fits.open(filename)
    header = hdulist[1].header
    hdulist.close()
    return header

def get_response(obsid, arm, detnam=None, maxorder=None, tg_reprocess='tg_reprocess', archive=False):
    if detnam is None:
        detnam = util.detnam(obsid, tg_reprocess=tg_reprocess, archive=archive)
    if detnam[:4] == 'ACIS': detnam = 'ACIS-S'
    bin_lo, bin_hi, rmfs = get_rmfs(detnam, arm, maxorder, archive=archive)
    bin_lo, bin_hi, garfs = get_garfs(obsid, arm, maxorder, archive=archive)
    response = {}
    for order in rmfs:
        response[order] = rmfs[order] * garfs[order]
    return bin_lo, bin_hi, response

def e2w(energy):
    return 12.39854/energy

def w2e(wavelength):
    return e2w(wavelength)

def bins2energy(lo, hi):
    return w2e(0.5 * (lo + hi))

