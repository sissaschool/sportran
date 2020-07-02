# -*- coding: utf-8 -*-

import numpy as np
from . import Current
from . import units

from thermocepstrum.utils.loadAfterPlt import plt
from thermocepstrum.utils.utils import PrintMethod
log = PrintMethod()

try:
    plt
except:
    log.write_log('Warning: plt undefined')

__all__ = ('ElectricCurrent',)

class ElectricCurrent(Current):
    """
    HeatCurrent API for thermo-cepstral analysis.
    Defines a HeatCurrent object with useful tools to perform analysis.

    INPUT:
     - traj          the heat current time series (N * N_COMPONENTS array)
       For a multi-component fluid use a (N_FLUID_COMPONENTS * N * N_COMPONENTS array)
     - UNITS         the units of current ('metal', 'real')
     - DT_FS         MD time step [fs]
     - TEMPERATURE   average temperature [K]
     - VOLUME        simulation cell volume [A^3]
     - PSD_FILTER_W  PSD filter window [freq_units] (optional)
     - FREQ_UNITS    frequency units   [THz or red] (optional)
    """
    _current_type = 'electric'
    _input_parameters = {'DT_FS', 'UNITS', 'TEMPERATURE', 'VOLUME'}

    # _optional_parameters = {'PSD_FILTER_W', 'FREQ_UNITS', 'MAIN_CURRENT_INDEX', 'MAIN_CURRENT_FACTOR'}

    def __init__(self, traj, **params):
        # params: (DT_FS, UNITS, TEMPERATURE, VOLUME, PSD_FILTER_W=None, FREQ_UNITS='THz')
        super().__init__(traj, **params)

    def _get_builder(self):
        """
        Get a tuple (class, builder) that can be used to build a new object with same parameters:
          TimeSeries, builder = self._get_builder()
          new_ts = TimeSeries(**builder)
        """
        if self.MANY_CURRENTS:
            traj_array = np.row_stack(([self.traj], [j.traj for j in self.otherMD]))
        else:
            traj_array = self.traj
        builder = dict(traj=traj_array, DT_FS=self.DT_FS, UNITS=self.UNITS, TEMPERATURE=self.TEMPERATURE,
                       VOLUME=self.VOLUME, PSD_FILTER_W=self.PSD_FILTER_W_THZ, FREQ_UNITS='THz')
        return type(self), builder

