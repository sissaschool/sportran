

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

def scale_kappa_DLPOLYtoSI ( temp, volume, timestep ):
  """Conversion factor for the thermal conductivity from DL_POLY units to SI units.
  INPUT:    temp      =  temperature [ K ]
            volume    =  cell volume [ A^3 ]
            timestep  =  integration time step [ fs ]"""

  kB = 1.3806504
  NA = 6.022140857;
  return (1.0/NA/temp)**2/kB/volume*timestep*1e10

def scale_kappa_CHARGEtoSI ( temp, volume, timestep ):
  """Conversion factor for the thermal conductivity from CHARGE DL_POLY units to SI units.
  INPUT:    temp      =  temperature [ K ]
            volume    =  cell volume [ A^3 ]
            timestep  =  integration time step [ fs ]"""

  kB = 1.3806504;
  charge = 1.6021765;
  return charge**2/temp/kB/volume*timestep*100

def scale_kappa_VELtoSI ( temp, volume, timestep ):
  """Conversion factor for the thermal conductivity from VEL DL_POLY units to SI units.
  INPUT:    temp      =  temperature [ K ]
            volume    =  cell volume [ A^3 ]
            timestep  =  integration time step [ fs ]"""

  kB = 1.3806504;
  return timestep*100
