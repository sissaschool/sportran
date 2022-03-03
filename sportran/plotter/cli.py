# -*- coding: utf-8 -*-
"""
Defines the CLIPlotter class.
"""

__all__ = ['CLIPlotter']

from . import plotter


class CLIPlotter(plotter.Plotter):
    """
    A Plotter subclass containing the plot functions used by the command-line interface.
    """
    from .plotter import (plot_periodogram, plot_ck, plot_L0_Pstar, plot_kappa_Pstar, plot_cepstral_spectrum,
                          plot_fstar_analysis, plot_resample, plot_psd, plot_cospectrum_component)
    _plot_style = 'cli_style.mplstyle'   # TODO define style for CLI
    pass


# probably we should decorate all these functions with addPlotToPdf and other decorators that allow any changes of
# style needed
