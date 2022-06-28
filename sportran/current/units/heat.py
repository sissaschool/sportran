# -*- coding: utf-8 -*-

from . import constants


def scale_kappa_real(TEMPERATURE, VOLUME):
    """
    Conversion factor for the thermal conductivity from REAL LAMMPS units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    The input current is in units of kcal/mole * Angstrom/femtosecond. Current is EXTENSIVE
    """
    return (constants.kcal / constants.NA / TEMPERATURE)**2 / constants.kB / VOLUME * 100.


def scale_kappa_metal(TEMPERATURE, VOLUME):
    """
    Conversion factor for the thermal conductivity from METAL LAMMPS units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    The input current is in units of eV * Angstrom/picosecond. Current is EXTENSIVE
    """
    return (constants.charge / TEMPERATURE)**2 / constants.kB / VOLUME * 10000.


def scale_kappa_qepw(TEMPERATURE, VOLUME):
    """
    Conversion factor for the thermal conductivity from Quantum Espresso PW/HARTREE-HEATCURRENT units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    The input current is in units of Rydberg atomic units. Current is EXTENSIVE
    Rydberg atomic units for the energy current are: Ry * a_0 / tau_{a.u.} where
      Ry = 13.606eV = 2.1799 10^{-18} J
      a0 = 5.2918 10^{-11}m
      tau_{a.u.} = 4.8378 10^{-17} s
    """
    return (constants.charge / TEMPERATURE)**2 / constants.kB / VOLUME * 10000. * (constants.Ry *
                                                                                   constants.J_PWtoMETAL)**2


def scale_kappa_gpumd(TEMPERATURE, VOLUME):
    """
    Conversion factor for the thermal conductivity from GPUMD units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    Note that units for time are derived from energy [eV], mass [amu] and position [A].
    Therefore, units for square velocity are [eV/amu].
    The input current is in units of eV ^ {3/2} / atomic_mass_unit ^ {1/2}. Current is EXTENSIVE
    """
    return (constants.charge)**3 / (TEMPERATURE)**2 / constants.massunit / constants.kB / VOLUME * 1.0e8
