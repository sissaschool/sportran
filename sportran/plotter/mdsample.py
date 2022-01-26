# -*- coding: utf-8 -*-
"""
Defines the MDSamplePlotter class.
"""

__all__ = ['MDSamplePlotter']

from . import plotter


class MDSamplePlotter(plotter.Plotter):
    """
    A Plotter subclass containing the plot functions used by an object of type Current.
    """
    from .plotter import (plot_trajectory, plot_resample)
    _plot_style = 'api_style.mplstyle'

    def plot_periodogram(current, PSD_FILTER_W=None, *, freq_units='THz', freq_scale=1.0, axes=None, FIGSIZE=None,
                         mode='log', **plot_kwargs):
        """
        Plots the current's periodogram (psd)
        :param current:         current object to plot periodogram
        :param PSD_FILTER_W:    width of the filtering window
        :param freq_units:      'thz'  [THz]
                                'red'  [omega*DT/(2*pi)]
        :param freq_scale:      rescale red frequencies by this factor (e.g. 2 --> freq = [0, 0.25])
        :param axes:            plot periodograms in units of kappa (default: True) - NB: log-psd not converted
        :param FIGSIZE:         size of the plot

        :return: a matplotlib.axes.Axes object
        """
        # kappa_units is not supported by MDSample
        from .plotter import plot_periodogram
        plot_kwargs.pop('kappa_units')
        return plot_periodogram(current, PSD_FILTER_W=PSD_FILTER_W, freq_units=freq_units, freq_scale=freq_scale,
                                axes=axes, kappa_units=False, FIGSIZE=FIGSIZE, mode=mode, **plot_kwargs)

    pass
