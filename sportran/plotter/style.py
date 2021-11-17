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
            import pkg_resources
            pltstyle_file = pkg_resources.resource_filename('sportran.plotter.styles', plot_style_filename)
        except:
            # fallback (if sportran is not installed...)
            try:
                abs_path = os.path.abspath(__file__)
                tc_path = abs_path[:abs_path.rfind('/')]
                os.path.append(tc_path[:tc_path.rfind('/')])
            except:
                abs_path = '.'
            pltstyle_file = tc_path + '/plotter/styles/' + plot_style_filename

    try:
        # print('using style ', plot_style_filename)
        plt.style.use(pltstyle_file)
    except:
        warn('The plot style {} could not be loaded.'.format(plot_style_filename))
