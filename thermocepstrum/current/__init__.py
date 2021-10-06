# -*- coding: utf-8 -*-

from .current import *
from .heat import *
from .electric import *
from .stress import *

__all__ = ('Current', 'HeatCurrent', 'ElectricCurrent', 'StressCurrent',)

# define list of all classes with units defined
import inspect
_all = dir()


def _get_currents_with_units():
    """Inspect all the classes accessible from this module, and detect the ones that contains the member `_current_type`. Then calls the methods to inspect the units implemented for each discovered class
    :return: ( {'_current_type': (CurrentClass, ['unit_list'], ['parameter_list])}  )
    """
    currents_with_units = {}
    all_units = []
    all_parameters = []
    for k in _all:
        v = globals()[k]
        att = getattr(v, '_current_type', None)
        if att is not None:
            parameters = []
            units = []
            for unit, funct in v._get_units().items():
                param = list(inspect.signature(funct).parameters.keys())
                units.append(unit)
                parameters += param
            parameters = list(set(parameters))
            currents_with_units[att] = (v, units, parameters)
            all_units += units
            all_parameters += parameters
    all_units = list(set(all_units))
    all_parameters = list(set(all_parameters))
    return currents_with_units, all_units, all_parameters


def _list_of_currents_and_units(verbose=False):
    s = ''
    for k, v_ in all_currents.items():
        v = v_[0]
        s += f"'{k}': {v._input_parameters}\n"
        for u in v.get_units_list():
            s += f"   - '{u}'\n"
            if verbose:
                s += v._get_units()[u].__doc__
    return s


#: list of currents classes, units implemented and parameters that are found dinamically when the module is imported
all_currents, all_units, all_parameters = _get_currents_with_units()


def build_currents_units_table(col=9):
    """print a table with the Current classes and the units implemented for each class"""
    table = ''
    table += ' ' * col
    for u in all_units:
        table += u[:col].ljust(col)
    table += '\n'
    for k, v in all_currents.items():
        table += k[:col].ljust(col)
        for unit in all_units:
            if unit in v[1]:
                table += 'X'.ljust(col)
            else:
                table += ' ' * col
        table += '\n'
    return table
