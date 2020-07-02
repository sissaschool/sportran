# -*- coding: utf-8 -*-

#from .current import Current
#from .heat import HeatCurrent

#__all__ = ['Current', 'HeatCurrent']

from .current import *
from .heat import *
from .electric import *

__all__ = ('Current', 'HeatCurrent', 'ElectricCurrent',)

#import os, pkgutil, importlib
#pkg_dir = os.path.dirname(__file__)
#print(pkg_dir)
#for (module_loader, name, ispkg) in pkgutil.iter_modules([pkg_dir]):
#    print('.' + name)
#    importlib.import_module('.' + name, __package__)

#    # Add the class to this package's variables
#    globals()[attribute_name] = attribute
#all_modules = [x[1] for x in pkgutil.iter_modules(path=[pkg_dir])]

#__all__ = [cls.__name__ for cls in current.Current.__subclasses__()]
#print(all_modules)
#print(__file__, __package__)

