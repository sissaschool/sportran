# -*- coding: utf-8 -*-

import numpy as np
from .tools.resample import filter_and_sample
from sportran.utils import log


def resample_timeseries(x, TSKIP=None, fstar_THz=None, FILTER_W=None, plot=False, PSD_FILTER_W=None, freq_units='THz',
                        FIGSIZE=None, verbose=True):   # yapf: disable
    """
    Simulate the resampling of a time series x.

    Parameters
    ----------
    x            = the time series object, a (subclass of) MDSample
    TSKIP        = sampling time [steps]
    fstar_THz    = target cutoff frequency [THz]
    TSKIP and fstar_THZ are mutually exclusive.

    FILTER_W     = pre-sampling filter window width [steps]
    plot         = plot the PSD [False]
    PSD_FILTER_W = PSD filtering window width [chosen frequency units]
    freq_units   = 'thz'  [THz]
                   'red'  [omega*DT/(2*pi)]
    FIGSIZE      = plot figure size
    verbose      = print log [True]

    Returns
    -------
    xf : a filtered & resampled time series object
    ax : an array of plot axes, optional (if plot=True)
    """
    from .mdsample import MDSample

    if not isinstance(x, MDSample):
        raise ValueError('x must be a MDSample object or a subclass.')
    if (TSKIP is not None) and (fstar_THz is not None):
        raise ValueError('Please specify either TSKIP or fstar_THz.')
    elif (TSKIP is None) and (fstar_THz is None):
        raise ValueError('Please specify either TSKIP or fstar_THz.')
    if x.psd is None or PSD_FILTER_W:
        x.compute_psd(PSD_FILTER_W, freq_units)   # this ensures that Nyquist_f_THz is defined
    if fstar_THz:
        # find TSKIP that gives the closest fstar
        TSKIP = int(round(x.Nyquist_f_THz / fstar_THz))
    if FILTER_W is None:
        FILTER_W = TSKIP
    if plot:
        raise NotImplementedError()

    fstar_THz = x.Nyquist_f_THz / TSKIP
    fstar_idx = np.argmin(x.freqs_THz < fstar_THz)

    # get the builder of the original time series object
    # builder is a dictionary containing the parameters to define the new time series
    TimeSeries, builder = x._get_builder()
    trajectory = builder['traj']

    # filter & resample time series -- for 3D time series, apply to each row
    if (len(trajectory.shape) == 3):
        new_trajectory = np.array([filter_and_sample(traj, FILTER_W, TSKIP, 'rectangular') for traj in trajectory])
    elif (len(trajectory.shape) < 3):
        new_trajectory = filter_and_sample(trajectory, FILTER_W, TSKIP, 'rectangular')
    else:
        raise ValueError('Trajectory shape not valid.')

    # define new time series by updating the builder
    builder.update(traj=new_trajectory, DT_FS=(x.DT_FS * TSKIP))
    xf = TimeSeries(**builder)

    # define new filtering window width to be equal to the original one
    new_PSD_FILTER_W = x.PSD_FILTER_W_THZ if freq_units in ('THz', 'thz') else x.PSD_FILTER_W * TSKIP
    if xf.psd is None:
        xf.compute_psd(new_PSD_FILTER_W, freq_units)   # this ensures that Nyquist_f_THz is defined

    # write log
    xf.resample_log = \
        '-----------------------------------------------------\n' +\
        '  RESAMPLE TIME SERIES\n' +\
        '-----------------------------------------------------\n' +\
        ' Original Nyquist freq  f_Ny =  {:12.5f} THz\n'.format(x.Nyquist_f_THz) +\
        ' Resampling freq          f* =  {:12.5f} THz\n'.format(fstar_THz) +\
        ' Sampling time         TSKIP =  {:12d} steps\n'.format(TSKIP) +\
        '                             =  {:12.3f} fs\n'.format(TSKIP * x.DT_FS) +\
        ' Original  n. of frequencies =  {:12d}\n'.format(x.NFREQS) +\
        ' Resampled n. of frequencies =  {:12d}\n'.format(xf.NFREQS)
    if x.fpsd is not None and xf.fpsd is not None:   # TODO: maybe substitute with if x.PSD_FILTER_W != 0
        xf.resample_log += \
            ' PSD      @cutoff  (pre-filter&sample) ~ {:12.5f}\n'.format(x.fpsd[fstar_idx]) +\
            '                  (post-filter&sample) ~ {:12.5f}\n'.format(xf.fpsd[-1]) +\
            ' log(PSD) @cutoff  (pre-filter&sample) ~ {:12.5f}\n'.format(x.flogpsd[fstar_idx]) +\
            '                  (post-filter&sample) ~ {:12.5f}\n'.format(xf.flogpsd[-1])
    xf.resample_log += \
        ' min(PSD)          (pre-filter&sample) = {:12.5f}\n'.format(x.psd_min) +\
        ' min(PSD)         (post-filter&sample) = {:12.5f}\n'.format(xf.psd_min) +\
        ' % of original PSD Power f<f* (pre-filter&sample)  = {:6.3f} %\n'.format(np.trapz(x.psd[:fstar_idx+1]) / x.psd_power * 100.)
    if x.fpsd is None or xf.fpsd is None:
        xf.resample_log += ' fPSD not calculated before resampling\n '
    xf.resample_log += '-----------------------------------------------------\n'

    if verbose:
        log.write_log(xf.resample_log)

    if plot:
        axes = x.plot_resample(xf, PSD_FILTER_W, freq_units=freq_units, axes=axes, FIGSIZE=FIGSIZE, mode='log')
        return xf, axes
    else:
        return xf
