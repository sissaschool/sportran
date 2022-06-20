# -*- coding: utf-8 -*-

from sportran.current import Current
from sportran.utils import log
from sportran.plotter.plotter import plot_fstar_analysis

__all__ = ['fstar_analysis']

def fstar_analysis(x, TSKIP_LIST, aic_type='aic', aic_Kmin_corrfactor=1.0, manual_cutoffK=None, plot=True, axes=None,
                   FIGSIZE=None, verbose=False, **plot_kwargs):   # yapf: disable
    """
    Perform cepstral analysis on a set of resampled time series, to study the effect of f*.
    For each TSKIP in TSKIP_LIST, the HeatCurrent x is filtered & resampled, and then cesptral-analysed.

    Parameters
    ----------
    TSKIP_LIST    = list of sampling times [steps]
    aic_type      = the Akaike Information Criterion function used to choose the cutoff ('aic', 'aicc')
    aic_Kmin_corrfactor = correction factor multiplied by the AIC cutoff (cutoffK = aic_Kmin * aic_Kmin_corrfactor)
    manual_cutoffK = (P*-1) = manual cutoff. If set, the AIC cutoff will be ignored.

    plot          = plot the PSD (default: True)
    axes          = matplotlib.axes.Axes object (if None, create one)
    FIGSIZE       = plot figure size
    verbose       = verbose output (default: False)
    **plot_kwargs = other parameters passed to plot function

    Returns
    -------
    xf : array_like
        an array of HeatCurrents, corresponding the values of TSKIP_LIST
    ax : array_like, optional (if plot=True)
        array of plot axes
    fig : figure, optional (if axes=None)
        plot figure
    """

    if not isinstance(x, Current):
        raise ValueError('x must be a Current object or a subclass.')

    xf = []
    for TSKIP in TSKIP_LIST:
        log.write_log('TSKIP = {:4d} - FSTAR = {:8g} THz'.format(TSKIP, x.Nyquist_f_THz / TSKIP))
        xff = x.resample(TSKIP=TSKIP, plot=False, verbose=verbose)
        xff.cepstral_analysis(aic_type=aic_type, aic_Kmin_corrfactor=aic_Kmin_corrfactor, manual_cutoffK=manual_cutoffK)
        xf.append(xff)
    FSTAR_THZ_LIST = [xff.Nyquist_f_THz for xff in xf]

    if plot:
        try:
            return plot_fstar_analysis(xf, FSTAR_THZ_LIST, original_current=x, axes=axes, FIGSIZE=FIGSIZE,
                                       **plot_kwargs)
        except AttributeError:
            print('Plotter does not support the plot_resample method')
    else:
        return xf
