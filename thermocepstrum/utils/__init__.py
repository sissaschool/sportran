# -*- coding: utf-8 -*-

from .logger import PrintMethod

try:
    log = PrintMethod()
except:
    raise RuntimeError('PrintMethod not defined.')

__all__ = ('log')
