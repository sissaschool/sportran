# -*- coding: utf-8 -*-
"""
Utilities to pretty-plot the results
"""

__all__ = ['plt', 'use_plot_style', 'Plotter', 'MDSamplePlotter', 'CurrentPlotter', 'CLIPlotter']

import matplotlib
# matplotlib.use('Agg')  # if needed use force=True, warn=False
import matplotlib.pyplot as plt

from .plotter import Plotter, addPlotToPdf, PdfPages
from .style import use_plot_style
from .mdsample import MDSamplePlotter
from .current import CurrentPlotter
from .cli import CLIPlotter
