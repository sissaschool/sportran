# -*- coding: utf-8 -*-

import numpy as np
from scipy.signal import lfilter


def filter_and_sample(y_big, W, DT, window='rectangular', even_NSTEPS=True, detrend=False, drop_first=True):
    """Filter signal with moving average window of width W and then sample it
    with time step DT."""

    if (W > 1):
        if (window == 'rectangular'):
            y_f = lfilter((1. / W) * np.ones(W), 1., y_big, axis=0)
        else:
            raise NotImplementedError('Not implemented window type.')

        # drop first W steps (initial conditions)
        if drop_first:
            y = y_f[(W - 1)::DT]
        else:
            y = y_f[::DT]
    else:
        y = y_big[::DT]

    # remove the mean
    if detrend:
        y = y - np.mean(y, axis=0)

    # keeps an even number of points
    if even_NSTEPS:
        if (y.shape[0] % 2 == 1):
            return y[:-1]
    return y


def resample_psd(freqs, psd, cutfrequency):
    if (cutfrequency >= freqs[-1]):
        return freqs, psd
    NFREQS = freqs.size - 1
    cutidx = (np.abs(freqs - cutfrequency)).argmin()
    if (NFREQS % cutidx == 0):   # cut frequency is sub-multiple of max freq
        DT = NFREQS / cutidx
        if (DT > 2):
            raise Warning('DT Not implemented.')
        newpsd = psd.copy()[:cutidx + 1]
        newpsd = newpsd + psd[:-cutidx - 2:-1]
        newfreqs = freqs[:cutidx + 1]
        #log.write_log(cutidx, DT, freqs[cutidx], newpsd.size)
    else:
        raise NotImplementedError('Not implemented.')
    return newfreqs, newpsd
