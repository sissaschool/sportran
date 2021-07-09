# -*- coding: utf-8 -*-

from . import constants

# Note: the factor 10000 appearing in metal and qepw units is due to the following:
# 1e(2*(-19) + 2*(-10) -2*(-12) -(-23) -3*(-10) -15) = 1e4
#       e_SI      ang       ps    kB       ang   fs
# since the multiplicative factor would be
# [(charge*1e-19)**2] [velocity**2] / (kB*1e-23) / [TEMPERATURE] / [VOLUME] * [integration DT_FS]
#
# The additional "1.0e6" factor appearing in real units is due to the velocity units,
# which are ang/fs in real units, rather than ang/ps as in metal units.


def scale_kappa_real(TEMPERATURE, VOLUME):
    """
    Conversion factor for the electrical conductivity from REAL LAMMPS units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    """
    return constants.charge**2 / TEMPERATURE / constants.kB / VOLUME * 10000. * 1.0e6


def scale_kappa_metal(TEMPERATURE, VOLUME):
    """
    Conversion factor for the thermal conductivity from METAL LAMMPS units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    """
    return constants.charge**2 / TEMPERATURE / constants.kB / VOLUME * 10000.


def scale_kappa_qepw(TEMPERATURE, VOLUME):
    """
    Conversion factor for the electrical conductivity from Quantum Espresso PW units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    """
    return constants.charge**2 / TEMPERATURE / constants.kB / VOLUME * 10000. * constants.J_PWtoMETAL**2


def scale_kappa_gpumd(TEMPERATURE, VOLUME):
    """
    Conversion factor for the thermal conductivity from GPUMD units to SI units.
    Note that units for time are derived from energy [eV], mass [amu] and position [A].
    Therefore, units for square velocity are [eV/amu].
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    """
    return constants.charge**3 / TEMPERATURE / constants.massunit / constants.kB / VOLUME * 1.0e8


def scale_kappa_dlpoly(TEMPERATURE, VOLUME):
    """
    Conversion factor for the thermal conductivity from DL_POLY units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    """
    return constants.charge**2 / TEMPERATURE / constants.kB / VOLUME * 10000.
