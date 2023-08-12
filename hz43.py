import glob
import numpy as np
import astropy.io.fits

import flux
import util
import response
 
basedir='/data/legs/rpete/flight/hz43'
datadir=basedir+'/data'

def model():
    return np.loadtxt('/data/legs/rpete/flight/hz43/model/non_LTE/15.model', unpack=True)

def obsids(exclude=None):
    global basedir
    obsids = np.loadtxt(basedir+'/obsids', usecols=(0,), unpack=True).astype(int)
    if exclude:
        for o in exclude:
            obsids = obsids[obsids!=o]
    return obsids

def predicted_rates(obsids):
    rates = []
    model_flux = None
    for obsid in obsids:
        w1, w2, specresp = response.read_arf(response.zeroth_arf_file(obsid))
        wav = 0.5*(w1+w2)

        if model_flux is None:
            model_wav, model_flux = model()
            model_flux = np.interp(wav, model_wav, model_flux) * (w2-w1)

        rates.append( np.sum(flux.flux2rate(wav, model_flux, specresp)) )
    return np.array(rates)

# return numpy arrays of HZ 43
def obsids_years(detector=None, offaxis=False, exclude=None):

    if 'data' not in obsids_years.__dict__:

        obsids_, years, dets = [], [], []

        for obsid in obsids(exclude=exclude):
            pha2 = util.pha2_file(obsid)
            hdr = util.read_header(pha2)
            obsids_.append(obsid)
            years.append(hdr['year'])
            dets.append(hdr['detnam'])

        obsids_ = np.array(obsids_)
        years = np.array(years)
        dets = np.array(dets)

        si = np.argsort(years)

        obsids_ = obsids_[si]
        years = years[si]
        dets = dets[si]

        onaxis_mask = (obsids_!=1516) & (obsids_!=1004) & (obsids_!=1005) & (obsids_!=2601) & (obsids_!=2603) & (obsids_!=18415)

        hrci_mask = dets == 'HRC-I'
        hrcs_mask = dets == 'HRC-S'

        obsids_years.data = {
            False : { # do not include offaxis
                None : { 'obsids' : obsids_[onaxis_mask], 'years' : years[onaxis_mask] },
                'HRC-I' : { 'obsids' : obsids_[hrci_mask&onaxis_mask], 'years' : years[hrci_mask&onaxis_mask] },
                'HRC-S' : { 'obsids' : obsids_[hrcs_mask&onaxis_mask], 'years' : years[hrcs_mask&onaxis_mask] }
            },
            True : { # include all observations, including offaxis
                None : { 'obsids' : obsids_, 'years' : years },
                'HRC-I' : { 'obsids' : obsids_[hrci_mask], 'years' : years[hrci_mask] },
                'HRC-S' : { 'obsids' : obsids_[hrcs_mask], 'years' : years[hrcs_mask] }
            },
        }

    return obsids_years.data[offaxis][detector]['obsids'], obsids_years.data[offaxis][detector]['years']
