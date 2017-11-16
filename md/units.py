

def scale_kappa_REALtoSI ( temp, volume, timestep ):
  """Conversion factor for the thermal conductivity from REAL LAMMPS units to SI units.
  INPUT:    temp      =  temperature [ K ]
            volume    =  cell volume [ A^3 ]
            timestep  =  integration time step [ fs ]"""

  kB = 1.3806504
  NA = 6.02214
  massunit = 1.660538921
  charge = 1.6021765;
  return (4184./NA/temp)**2/kB/volume*timestep*100.


def scale_kappa_METALtoSI ( temp, volume, timestep ):
  """Conversion factor for the thermal conductivity from METAL LAMMPS units to SI units.
  INPUT:    temp      =  temperature [ K ]
            volume    =  cell volume [ A^3 ]
            timestep  =  integration time step [ fs ]"""

  kB = 1.3806504
  NA = 6.02214
  massunit = 1.660538921
  charge = 1.6021765;
  return (charge/temp)**2/kB/volume*timestep*10000.

