# -*- coding: utf-8 -*-

from . import constants


def scale_kappa_GPa(TEMPERATURE, VOLUME):
    """
    Conversion factor for the viscosity to SI units when stress is in GPa.
    Notice that GPa are the units used for stress in cp.x *.str file.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    Input stress in GPa.
    """
    return 1.0e-4 * VOLUME / TEMPERATURE / constants.kB


def scale_kappa_real(TEMPERATURE, VOLUME):
    """
    Conversion factor for the viscosity from REAL LAMMPS units to SI units.
    Stress is computed with https://docs.lammps.org/compute_pressure.html
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    Input stress in atmospheres.
    """
    return constants.atm**2 * 1.0e-12 * VOLUME / TEMPERATURE / constants.kB


def scale_kappa_metal(TEMPERATURE, VOLUME):
    """
    Conversion factor for the viscosity from METAL LAMMPS units to SI units.
    Stress is computed with https://docs.lammps.org/compute_pressure.html
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    Input stress in bars.
    """
    return 1.0e-12 * VOLUME / TEMPERATURE / constants.kB


def scale_kappa_qepw(TEMPERATURE, VOLUME):
    """
    Conversion factor for the viscosity from Quantum Espresso PW units to SI units.
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    Input stress is in units of Ry / bohr_radius^3 .
    """
    return (constants.charge * constants.Ry_per_bohr3)**2 * VOLUME / TEMPERATURE / constants.kB


def scale_kappa_gpumd(TEMPERATURE, VOLUME):
    """
    Conversion factor for the viscosity from GPUMD units to SI units.
    Note that units for stress and pressure, derived from energy [eV] and volume [A^3],
    are [eV/A^3].
    INPUT:
    TEMPERATURE [K]
    VOLUME      cell VOLUME [A^3]
    Input stress is in units of eV / Angstrom^3 .
    """
    return constants.charge**2 * VOLUME / TEMPERATURE / constants.kB
