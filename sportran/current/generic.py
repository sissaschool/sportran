# -*- coding: utf-8 -*-

from . import Current

__all__ = ['GenericCurrent']


class GenericCurrent(Current):
    """
    GenericCurrent API for thermo-cepstral analysis.
    Defines a HeatCurrent object with useful tools to perform analysis.

    INPUT parameters:
     - traj          the current time series (N, N_EQUIV_COMPONENTS) array
       For a multi-component fluid use a (N_CURRENTS, N, N_EQUIV_COMPONENTS) array
     - DT_FS         MD time step [fs]
     - KAPPA_SCALE   the GK conversion factor, multiplies the GK integral

    OPTIONAL parameters:
     - PSD_FILTER_W  PSD filter window [freq_units] (optional)
     - FREQ_UNITS    frequency units   [THz or red] (optional)
     - MAIN_CURRENT_INDEX for a multi-current time series, the index of the "main" current (e.g. energy) [0]
     - MAIN_CURRENT_FACTOR factor to be multiplied by the main current [1.0]
    """
    _current_type = None
    _input_parameters = {'DT_FS', 'KAPPA_SCALE'}
    _KAPPA_SI_UNITS = ''
    # _optional_parameters = {'PSD_FILTER_W', 'FREQ_UNITS', 'MAIN_CURRENT_INDEX', 'MAIN_CURRENT_FACTOR'}

    @property
    def _builder(self):
        """
        Returns a dictionary of all keyworded parameters needed to rebuild an identical object of the same class.
        The trajectory is excluded. Used by self._get_builder().
        """
        return dict(DT_FS=self.DT_FS, KAPPA_SCALE=self.KAPPA_SCALE, PSD_FILTER_W=self.PSD_FILTER_W_THZ,
                    FREQ_UNITS='THz')
