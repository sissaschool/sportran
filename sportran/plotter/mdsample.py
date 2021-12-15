# -*- coding: utf-8 -*-
"""
Defines the MDSamplePlotter class.
"""

__all__ = ('MDSamplePlotter')

from . import plotter


class MDSamplePlotter(plotter.Plotter):
    """
    A Plotter subclass containing the plot functions used by an object of type Current.
    """
    from .plotter import (plot_trajectory, plot_periodogram, plot_resample)
    _plot_style = 'api_style.mplstyle'
    pass
