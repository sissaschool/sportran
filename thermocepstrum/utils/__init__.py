# -*- coding: utf-8 -*-

from .loadAfterPlt import Plt
from .logger import PrintMethod

__all__ = ('loadAfterPlt', 'Plt')

try:
    log = PrintMethod()
except:
    raise RuntimeError('PrintMethod not defined.')

try:
    plt = Plt()
except:
    raise Warning('Plotter undefined')
