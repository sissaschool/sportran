# -*- coding: utf-8 -*-
"""
==========================
Sportran library
==========================

This is the core library containing all the code necessary to perform a nice
and fast cepstral analysis on your time series.
"""

from . import md
from . import current
from .current import *
from . import utils
from . import plotter
from . import i_o

__all__ = [current.__all__ + md.__all__]

__license__ = 'GPL-3.0 license, see LICENSE.txt file.'
__version__ = '1.0.0rc1'
__authors__ = 'Loris Ercole, Riccardo Bertossa, Sebastiano Bisacchi'
__paper__ = (
    'L. Ercole, R. Bertossa, S. Bisacchi, S.Baroni, "SporTran: a code to estimate transport coefficients from the '
    'cepstral analysis of (multivariate) current time series", arXiV:2202.11571 (2022), submitted to Computer'
    'Physics Communications, https://doi.org/10.48550/arXiv.2202.11571')
__paper_short__ = 'L. Ercole et al., arXiv:2202.117571 (2022)'
