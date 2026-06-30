# -*- coding: utf-8 -*-

import os
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
        pltstyle_file = None
        try:
            from importlib import resources

            with resources.as_file(
                resources.files('sportran.plotter.styles').joinpath(plot_style_filename)
            ) as style_path:
                pltstyle_file = str(style_path)
        except Exception:
            pass

        if pltstyle_file is None:
            # fallback (if sportran is not installed...)
            pltstyle_file = os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                'styles',
                plot_style_filename,
            )

    try:
        # print('using style ', plot_style_filename)
        plt.style.use(pltstyle_file)
    except:
        warn('The plot style {} could not be loaded.'.format(pltstyle_file))
