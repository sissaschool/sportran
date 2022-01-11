# -*- coding: utf-8 -*-

import numpy as np
from . import Current
from . import units

__all__ = ['StressCurrent']


class StressCurrent(Current):
    """
    StressCurrent API for thermo-cepstral analysis.
    Defines a StressCurrent object with useful tools to perform analysis.

    INPUT parameters:
     - traj          the current time series (N, N_EQUIV_COMPONENTS) array
       For a multi-component fluid use a (N_CURRENTS, N, N_EQUIV_COMPONENTS) array
     - DT_FS         MD time step [fs]
     - UNITS         the units of current ('metal', 'real', ...) - use the method `get_units_list()` to get a list of supported units
     - TEMPERATURE   average temperature [K]
     - VOLUME        simulation cell volume [A^3]

    OPTIONAL parameters:
     - PSD_FILTER_W  PSD filter window [freq_units] (optional)
     - FREQ_UNITS    frequency units   [THz or red] (optional)
     - MAIN_CURRENT_INDEX for a multi-current time series, the index of the "main" current (e.g. energy) [0]
     - MAIN_CURRENT_FACTOR factor to be multiplied by the main current [1.0]
    """
    _current_type = 'stress'
    KAPPA_SI_UNITS = 'Pa*s'
    _input_parameters = {'DT_FS', 'UNITS', 'TEMPERATURE', 'VOLUME'}

    # _optional_parameters = {'PSD_FILTER_W', 'FREQ_UNITS', 'MAIN_CURRENT_INDEX', 'MAIN_CURRENT_FACTOR'}

    def __init__(self, traj, **params):
        # params: (DT_FS, UNITS, TEMPERATURE, VOLUME, PSD_FILTER_W=None, FREQ_UNITS='THz')
        super().__init__(traj, **params)

    @property
    def _builder(self):
        """
        Returns a dictionary of all keyworded parameters needed to rebuild an identical object of the same class.
        The trajectory is excluded. Used by self._get_builder().
        """
        return dict(DT_FS=self.DT_FS, UNITS=self.UNITS, TEMPERATURE=self.TEMPERATURE, VOLUME=self.VOLUME,
                    PSD_FILTER_W=self.PSD_FILTER_W_THZ, FREQ_UNITS='THz')
