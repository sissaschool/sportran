# -*- coding: utf-8 -*-

kB = 1.3806504
NA = 6.02214
massunit = 1.660538921
charge = 1.6021765
kcal = 4184.   # Joule
Ha = 27.211386   # eV
Ry = Ha / 2.0   # eV
bohr = 0.529177   # A
tau_PW = 4.8378e-5   # ps
atm = 1.01325   # bars

J_PWtoMETAL = bohr / tau_PW   # [Bohr/tau_PW] --> [A/ps]
Ry_per_bohr3 = Ry / bohr**3   # [Ry/bohr^3] --> [eV/A^3]
dlpoly_press = massunit * 1.0e2   # bars

A_to_m = 1e-10   # [A] --> [m]
A_per_fs_to_A_per_ps = 1.0e-3   # [A/fs] --> [A/ps]
A_per_ps_to_m_per_s = 1.0e2   # [A/ps] --> [m/s]
