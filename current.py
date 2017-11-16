################################################################################
###   "current" API
################################################################################

import numpy as np
import md


class HeatCurrent(MDSample):
   """
   HeatCurrent API for thermo-cepstral analysis.
   """

   def __init__(self, j, units, DT_FS, TEMPERATURE, VOLUME):
      MDSample.__init__(self, traj=j, DT_FS=DT_FS)
      self.initialize_units(units, TEMPERATURE, VOLUME, DT_FS)
      return

   def initialize_units(self, units, TEMPERATURE, VOLUME, DT_FS):
      self.units = units 
      if (self.units == 'metal'):
         self.kappa_scale = md.units.scale_kappa_METALtoSI(TEMPERATURE, VOLUME, DT_FS)
      elif (self.units == 'real'):
         self.kappa_scale = md.units.scale_kappa_REALtoSI(TEMPERATURE, VOLUME, DT_FS)
      else:
         raise ValueError('Units not supported.')
      return

   def 
