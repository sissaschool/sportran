# -*- coding: utf-8 -*-

from os.path import isfile
from . import plt
from warnings import warn

DEFAULT_PLOT_STYLE = 'api_style.mplstyle'


def use_plot_style(plot_style_filename=None):
    """
    Use a matplotlib plot style file.
    """
    if plot_style_filename is None:
        plot_style_filename = DEFAULT_PLOT_STYLE
    #print('Using {} plot style.'.format(plot_style_filename))

    # try to import matplotlib style settings
    if isfile(plot_style_filename):
        pltstyle_file = plot_style_filename
    else:
        try:
            from importlib.resources import files
            pltstyle_file = str(files('sportran.plotter.styles').joinpath(plot_style_filename))
        except Exception:
            try:
                import pkg_resources
                pltstyle_file = pkg_resources.resource_filename('sportran.plotter.styles', plot_style_filename)
            except Exception:
                import os
                tc_path = os.path.dirname(os.path.abspath(__file__))
                pltstyle_file = os.path.join(tc_path, 'styles', plot_style_filename)

    try:
        # print('using style ', plot_style_filename)
        plt.style.use(pltstyle_file)
    except:
        warn('The plot style {} could not be loaded.'.format(pltstyle_file))
