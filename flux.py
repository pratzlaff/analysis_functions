import numpy as np
import util

def rate2flux(rate, rate_err, response, bin_lo, bin_hi):
    energy = util.bins2energy(bin_lo, bin_hi)
    mult = energy * 1.602e-9 / response
    flux_ = rate * mult
    flux_err = rate_err * mult
    return flux_, flux_err

def flux2rate(wav, flux_, response):
    rate = flux_ / 1.602e-9 / util.w2e(wav) * response
    return rate

def rate_summed(src, bg, bin_lo, bin_hi, hdr, factor=1):

    bg_mult = util.bg_mult(hdr)

    net = src - bg * bg_mult
    net_var = src + bg * bg_mult**2

    rate = net / hdr['exposure']
    rate_err = np.sqrt(net_var) / hdr['exposure']

    rate *= factor
    rate_err *= factor

    rate = rate.sum()
    rate_err = np.sqrt((rate_err**2).sum())

    return rate, rate_err

def flux_summed(src, bg, resp, bin_lo, bin_hi, hdr, factor=1):

    bg_mult = util.bg_mult(hdr)

    net = src - bg * bg_mult
    net_var = src + bg * bg_mult**2

    rate = net / hdr['exposure']
    rate_err = np.sqrt(net_var) / hdr['exposure']

    rate *= factor
    rate_err *= factor

    flux_, flux_err = rate2flux(rate, rate_err, resp, bin_lo, bin_hi)

    rate = rate.sum()
    rate_err = np.sqrt((rate_err**2).sum())

    flux_ = flux_.sum()
    flux_err = np.sqrt((flux_err**2).sum())

    return rate, rate_err, flux_, flux_err

def flux_summed2(src, bg, resp, bin_lo, bin_hi, hdr):

    bg_mult = util.bg_mult(hdr)

    src_sum = src.sum()
    bg_sum = bg.sum()
    resp_sum = resp.mean()

    net = src_sum - bg_sum * bg_mult
    net_var = src_sum + bg_sum * bg_mult**2

    rate = net / hdr['exposure']
    rate_err = np.sqrt(net_var) / hdr['exposure']

    flux_, flux_err = rate2flux(rate, rate_err, resp_sum, bin_lo.min(), bin_hi.max())
    return rate, rate_err, flux_, flux_err
