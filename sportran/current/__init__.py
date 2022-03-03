# -*- coding: utf-8 -*-
"""
This module handles all the various type of currents and units of measure that can be handled natively by this library.

The currently supported units and currents can be viewed by calling the function :py:func:`build_currents_units_table`

To add a new current type and/or unit you have to do the following.
1. define in a new file a class derived from :py:class:`.current.Current` with the following class variables:
  - `_current_type` that you can set to the name of the current that will be used in all the user interfaces
  - `_input_parameters` is the set of the parameters that your current needs. Usually `{'DT_FS', 'UNITS', 'TEMPERATURE', 'VOLUME'}`
  - `_KAPPA_SI_UNITS` string that describes the units of the final result to be displayed where appropriate
2. define `__init__` like
  ``
  def __init__(self, traj, **params):
     super().__init__(traj, **params)
  ``
3. define `_get_builder` that must return the current class type and a list of parameters that applied to the class constructor generate an instance identical to the current one, so that it can be used like
   ``
   CurrentType, builder = self._get_builder()
   new_ts = CurrentType(**builder)
   ``
4. define the units for your class. The units are automatically retrieved from the module :py:mod:`.current.units`. You must do the following:
   - add a file named `_current_type`.py (a file with the same name of your class attribute value)
   - write inside the file some functions `scale_kappa_*` the part of the name after `scale_kappa_` will be the name of the unit used in the user interfaces.
     The parameters must be consistent with `_input_parameters`. The parameters (with the same name) will be provided by the user in the class constructor or in the various interfaces

"""
from .current import *
from .heat import *
from .electric import *
from .stress import *

__all__ = ['Current', 'HeatCurrent', 'ElectricCurrent', 'StressCurrent']

# define list of all classes with units defined
import inspect

_all = dir()


def _get_currents_with_units():
    """
    Inspect all the classes accessible from this module, and detect the ones that contains the attribute `_current_type`.
    Then call the `get_units` method to inspect the units implemented for each discovered class.
    :return: ( {'_current_type': (CurrentClass, ['unit_list'], ['parameter_list'])} )
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


# list of currents classes, units implemented and parameters that are found dynamically when the module is imported
all_currents, all_units, all_parameters = _get_currents_with_units()


def build_currents_units_table(col=9):
    """Print a table with the Current classes and the units implemented for each class"""
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
