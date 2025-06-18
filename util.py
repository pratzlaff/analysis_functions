import glob
try:
    import astropy.io.fits
except:
    pass
import numpy as np
import errno, os, sys

import flux
import response

basedir='/data/legs/rpete/flight/hz43'
archivedir='/data/loss/rpete/hrc'
datadir=basedir+'/data'

tg_parts = { 'HEG' : 1, 'MEG' : 2, 'LEG' : 3 }

def calc_rate(src, bg, hdr):
    mult =  bg_mult(hdr)

    net = src - bg * mult
    net_var = src + bg * mult**2

    rate = net / hdr['exposure']
    rate_var = net_var / hdr['exposure']**2

    return rate, np.sqrt(rate_var)

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def e2w(energy):
    return 12.39854/energy

def w2e(wavelength):
    return e2w(wavelength)

def bins2energy(lo, hi):
    return w2e(0.5 * (lo + hi))

def evt2_file(obsid, tg_reprocess='tg_reprocess', archive=False):

    if archive:
        global archivedir
        obsid = f'{int(obsid):05d}'
        globstr = f'{archivedir}/[is]/{obsid}/analysis/hrcf{obsid}_evt2.fits'
        #sys.stderr.write(globstr+"\n")
        return glob.glob(globstr)[0]

    global datadir
    globstr = '{}/{}/{}/*_evt2.fits'.format(datadir, obsid, tg_reprocess)
    #sys.stderr.write(globstr+"\n")
    return glob.glob(globstr)[0]

def pha2_file(obsid, tg_reprocess='tg_reprocess', archive=False):

    if archive:
        global archivedir
        obsid = f'{int(obsid):05d}'
        globstr = f'{archivedir}/[is]/{obsid}/analysis/hrcf{obsid}_pha2.fits'
        #sys.stderr.write(globstr+"\n")
        return glob.glob(globstr)[0]

    global datadir
    globstr = '{}/{}/{}/*_pha2.fits'.format(datadir, obsid, tg_reprocess)
    #sys.stderr.write(globstr+"\n")
    return glob.glob(globstr)[0]

def read_pha2(filename):
        hdulist = astropy.io.fits.open(filename)
        hdr = hdulist['spectrum'].header
        data = hdulist['spectrum'].data
        hdulist.close()

        try:
            hdr['year'] = 1998 + (hdr['mjd_obs'] - 50814)/365.2425
        except:
            hdr['year'] = 1998 + (hdr['mjd-obs'] - 50814)/365.2425

        return data, hdr

def read_header(filename):
    hdulist = astropy.io.fits.open(filename)
    header = hdulist[1].header
    hdulist.close()

    # evt2 files used to have mjd_obs, now they have mjd-obs???
    try:
        header['year'] = 1998 + (header['mjd_obs'] - 50814)/365.2425
    except:
        header['year'] = 1998 + (header['mjd-obs'] - 50814)/365.2425
    return header

def zeroth_rates(obsids, tg_reprocess='tg_reprocess', archive=False):

    rates = np.zeros(obsids.size)
    rate_errs = rates.copy()
    exposures = rates.copy()

    for i in range(obsids.size):
        obsid = obsids[i]
        evt2 = evt2_file(obsid, tg_reprocess=tg_reprocess, archive=archive)
        sys.stderr.write('Using '+evt2+"\n")
        hdr = read_header(evt2)

        x0, y0, r0, x, y, jnk, jnk = read_evt2(evt2)
        r0 = r0*2.5
        #sys.stderr.write("r0={}".format(r0)+"\n")

        # background annulus
        bgratio = 2.0
        bg_r1 = r0
        bg_r2 = np.sqrt(r0**2 * bgratio + bg_r1**2)

        dist = np.sqrt((x-x0)**2 + (y-y0)**2)

        src = np.where(dist < r0)[0].size
        bg = np.where((dist > bg_r1) & (dist < bg_r2))[0].size

        src_err = np.sqrt(src)
        bg_err = np.sqrt(bg)

        bg /= bgratio
        bg_err /= bgratio

        net = src - bg
        net_err = np.sqrt(src_err**2 + bg_err**2)

        exposure = hdr['exposure']

        rates[i] = net/exposure
        rate_errs[i] = net_err/exposure
        exposures[i] = exposure

    return rates, rate_errs, exposures

def read_evt2(filename):
    hdulist = astropy.io.fits.open(filename)
    header = hdulist['events'].header
    data = hdulist['events'].data

    x = data.field('x')
    y = data.field('y')

    chipx = data.field('chipx')
    chipy = data.field('chipy')

    data_region = hdulist['region'].data
    x0 = data_region['x'][0]
    y0 = data_region['y'][0]
    r0 = data_region['r'][0][0]

    hdulist.close()

    return x0, y0, r0, x, y, chipx, chipy

def bg_mult(hdr):
    return hdr['backscal'] / (hdr['backscup']+hdr['backscdn'])

def grating(obsid, tg_reprocess='tg_reprocess', archive=False):
    pha2 = pha2_file(obsid, tg_reprocess, archive=archive)
    hdr = read_header(pha2)
    return hdr['detnam']

def detnam(obsid, tg_reprocess='tg_reprocess', archive=False):
    pha2 = pha2_file(obsid, tg_reprocess, archive=archive)
    hdr = read_header(pha2)
    return hdr['detnam']

def detnam(obsid, tg_reprocess='tg_reprocess', archive=False):
    pha2 = pha2_file(obsid, tg_reprocess, archive=archive)
    hdr = read_header(pha2)
    return hdr['detnam']

def get_spectrum(obsid, grating, combine=1, tg_reprocess='tg_reprocess', archive=False):
    spectrum = get_spectra((obsid,), grating, combine, tg_reprocess, archive=archive)
    return spectrum

def get_spectra2(obsids, arm, combine=1, tg_reprocess='tg_reprocess', archive=False):

    global tg_parts

    spectra = {'neg':{}, 'pos':{}}
    wav = None

    orders = {'neg':-1, 'pos':+1}

    rate = None
    bin_lo, bin_hi = {}, {}
    exposure = np.zeros(len(obsids))

    for i in range(len(obsids)):
 
        obsid = obsids[i]

        pha2 = pha2_file(obsid, tg_reprocess=tg_reprocess, archive=archive)
        jnk, jnk, response_ = get_response(obsid, arm, maxorder=1, archive=archive)
        data, hdr = read_pha2(pha2)
        exposure[i] = hdr['exposure']

        if rate is None:
            rate, rate_var, flux_, flux_var, response = ({} for j in range(5))
            for order in orders:
                rate[order] = np.zeros((len(obsids), response_[order][0].size))
                rate_var[order] = rate[order].copy()
                flux_[order] = rate[order].copy()
                flux_var[order] = rate[order].copy()
                response[order] = rate[order].copy()

        for order in spectra:
                
            row = np.where((data['tg_m']==orders[order]) & (data['tg_part']==tg_parts[arm]))[0][0]
            rate_, rate_err_ = calc_rate(data['counts'][row], data['background_up'][row] + data['background_down'][row], hdr)
            # reverse arrays, to be ordered by increasing wavelength
            rate_=rate_[::-1]
            rate_err_=rate_err_[::-1]
            response[order][i] = response_[order][0][::-1]

            if order not in bin_lo:
                bin_lo[order] = data['bin_lo'][row][::-1]
                bin_hi[order] = data['bin_hi'][row][::-1]

            flux__, flux_err_ = flux.rate2flux(rate_, rate_err_, response[order][i], bin_lo[order], bin_hi[order])

            rate[order][i] = rate_
            rate_var[order][i] = rate_err_**2
            flux_[order][i] = flux__
            flux_var[order][i] = flux_err_**2

    # end of for obsids

    # rebin everything
    if combine != 1:

        for order in rate:
            mod = np.mod(rate[order].shape[1], combine)
            if mod:
                # sys.stderr.write("mod={}".format(mod)+"\n")
                rate[order]     = rate[order][:,:-mod]
                rate_var[order] = rate_var[order][:,:-mod]
                flux_[order]     = flux_[order][:,:-mod]
                flux_var[order] = flux_var[order][:,:-mod]
                bin_lo[order]   = bin_lo[order][:-mod]
                bin_hi[order]   = bin_hi[order][:-mod]
                response[order] = response[order][:-mod]

            newshape = (rate[order].shape[0], rate[order].shape[1]/combine, combine)

            rate[order]     = np.reshape(rate[order], newshape).sum(axis=-1)
            rate_var[order] = np.reshape(rate_var[order], newshape).sum(axis=-1)
            flux_[order]     = np.reshape(flux_[order], newshape).sum(axis=-1)
            flux_var[order] = np.reshape(flux_var[order], newshape).sum(axis=-1)
            response[order] = np.reshape(response[order], newshape).sum(axis=-1)

            bin_lo[order] = bin_lo[order][::combine]
            bin_hi[order] = bin_hi[order][combine-1::combine]

            # delta_lambda = np.abs(
            #     np.median(s[order]['lambda'][:-1] - s[order]['lambda'][1:]) )

            # for key in 'flux', 'flux_err':
            #     s[order][key] /= delta_lambda

    retval = { }
    weights = (exposure / exposure.sum())[:, np.newaxis]
    for order in spectra:

        retval[order] = {}
        retval[order]['lambda'] = 0.5*(bin_lo[order]+bin_hi[order])

        retval[order]['flux'] = (weights * flux[order]).sum(axis=0)
        retval[order]['flux_err'] = np.sqrt((weights**2 * flux_var[order]).sum(axis=0))

        retval[order]['rate'] = (weights * rate[order]).sum(axis=0)
        retval[order]['rate_err'] = np.sqrt((weights**2 * rate_var[order]).sum(axis=0))

        retval[order]['response'] = (weights * response[order]).sum(axis=0)

    return retval

def get_spectra(obsids, arm, combine=1, tg_reprocess='tg_reprocess', archive=False):

    global tg_parts

    spectra = {'neg':{}, 'pos':{}}
    wav = None

    orders = {'neg':-1, 'pos':+1}

    rates = None

    for i in range(len(obsids)):

        obsid = obsids[i]
        pha2 = pha2_file(obsid, tg_reprocess=tg_reprocess, archive=archive)
        jnk, jnk, response_ = response.get_response(obsid, arm, maxorder=1, archive=archive)

        pha2_d, pha2_h = read_pha2(pha2)

        s = { 'neg':{}, 'pos':{} }

        for order in spectra:
            row = np.where((pha2_d['tg_m']==orders[order]) & (pha2_d['tg_part']==tg_parts[arm]))[0][0]
            src = np.copy(pha2_d['counts'][row])
            bg = pha2_d['background_up'][row] + pha2_d['background_down'][row]
            bin_lo = np.copy(pha2_d['bin_lo'][row])
            bin_hi = np.copy(pha2_d['bin_hi'][row])
            resp = np.ma.masked_invalid(np.copy(response_[order][0]))
            resp = np.copy(response_[order][0])

            # reverse the arrays to be ordered by increasing wavelength
            src=src[::-1]
            bg=bg[::-1]
            bin_lo=bin_lo[::-1]
            bin_hi=bin_hi[::-1]
            resp=resp[::-1]

            # rebin everything
            if combine != 1:

                mod = np.mod(src.shape[0], combine)

                if mod:
                    sys.stderr.write("mod={}".format(mod)+"\n")
                    src = src[:-mod]
                    bg = bg[:-mod]
                    bin_lo = bin_lo[:-mod]
                    bin_hi = bin_hi[:-mod]
                    resp = resp[:-mod]

                newshape = (int(src.shape[0]/combine), combine)

                src = np.reshape(src, newshape).sum(axis=1)
                bg = np.reshape(bg, newshape).sum(axis=1)
                resp = np.reshape(resp, newshape).sum(axis=1)

                bin_lo = bin_lo[::combine]
                bin_hi = bin_hi[combine-1::combine]

            wav = 0.5*(bin_lo+bin_hi)

            rate, rate_err = calc_rate(src, bg, pha2_h)
            flux_, flux_err = flux.rate2flux(rate, rate_err, resp, bin_lo, bin_hi)

            s[order]['counts'] = src
            # s[order]['bg'] = bg
            # s[order]['bg_err'] = bg_err
            s[order]['rate'] = rate
            s[order]['rate_err'] = rate_err
            s[order]['flux'] = flux_
            s[order]['flux_err'] = flux_err
            s[order]['response'] = resp
            s[order]['bin_lo'] = bin_lo
            s[order]['bin_hi'] = bin_hi
            s[order]['lambda'] = wav

            # delta_lambda = np.abs(
            #     np.median(s[order]['lambda'][:-1] - s[order]['lambda'][1:]) )

            # for key in 'flux', 'flux_err':
            #     s[order][key] /= delta_lambda

            # create the data arrays if they don't yet exist
            if not spectra[order].keys():
                wav = s[order]['lambda']
                zeros = np.zeros((len(obsids), s[order]['flux'].size))
                for key in 'rate', 'rate_err', 'flux', 'flux_err', 'response':
                    spectra[order][key] = np.copy(zeros)
                spectra[order]['exposure'] = np.zeros(len(obsids))

            spectra[order]['exposure'][i] = pha2_h['exposure']
            for key in 'rate', 'rate_err', 'flux', 'flux_err', 'response':
                spectra[order][key][i] = s[order][key]

    retval = { }
    for order in spectra:
        weights = spectra[order]['exposure'] / spectra[order]['exposure'].sum()
        weights = weights[:,np.newaxis]

        retval[order] = { }
        retval[order]['lambda'] = wav

        retval[order]['flux'] = (weights * spectra[order]['flux']).sum(axis=0)
        retval[order]['flux_err'] = np.sqrt((weights**2 * spectra[order]['flux_err']**2).sum(axis=0))

        retval[order]['rate'] = (weights * spectra[order]['rate']).sum(axis=0)
        retval[order]['rate_err'] = np.sqrt((weights**2 * spectra[order]['rate_err']**2).sum(axis=0))

        retval[order]['response'] = (weights * spectra[order]['response']).sum(axis=0)

    return retval
