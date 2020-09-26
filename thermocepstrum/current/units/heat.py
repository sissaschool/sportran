# -*- coding: utf-8 -*-

kB = 1.3806504
NA = 6.02214
massunit = 1.660538921
charge = 1.6021765


def scale_kappa_real(TEMPERATURE, VOLUME):
    """
    Conversion factor for the thermal conductivity from REAL LAMMPS units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    """
    return (4184. / NA / TEMPERATURE)**2 / kB / VOLUME * 100.


def scale_kappa_metal(TEMPERATURE, VOLUME):
    """
    Conversion factor for the thermal conductivity from METAL LAMMPS units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    """
    return (charge / TEMPERATURE)**2 / kB / VOLUME * 10000.


def scale_kappa_qepw(TEMPERATURE, VOLUME):
    """
    Conversion factor for the thermal conductivity from Quantum Espresso PW/HARTREE-HEATCURRENT units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    """
    J_PWtoMETAL = 1.0 / 6.719329152e-6   # [Ry*Bohr/tau_PW] --> [ev*A/ps]  (1tau_PW = 4.8378e-5 ps)
    return (charge / TEMPERATURE)**2 / kB / VOLUME * 10000. * J_PWtoMETAL**2


def scale_kappa_gpumd(TEMPERATURE, VOLUME):
    """
    Conversion factor for the thermal conductivity from GPUMD units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    """
    return (charge)**3 / (TEMPERATURE)**2 / massunit / kB / VOLUME * 1.0e8


def scale_kappa_dlpoly(TEMPERATURE, VOLUME):
    """
    Conversion factor for the thermal conductivity from DL_POLY units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    """
    return (1.0 / NA / TEMPERATURE)**2 / kB / VOLUME * 1e10
