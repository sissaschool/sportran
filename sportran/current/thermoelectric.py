# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import CubicSpline

from . import Current
from .results import ThermoelectricResult
from .units import electric as electric_units
from .units import heat as heat_units
from sportran.md.maxlike import scale_matrix

__all__ = ['ThermoElectricCurrent']


class ThermoElectricCurrent(Current):
    """
    Thermoelectric current API for coupled heat/charge transport.

    The class expects already prepared heat and charge currents.

    INPUT parameters:
     - traj          coupled current time series as a (2, N, N_EQUIV_COMPONENTS) array,
                     where index 0 is heat current and index 1 is charge current.
                     Alternatively, pass heat_current and charge_current named arguments.
     - DT_FS         MD time step [fs]
     - UNITS         units of current ('metal', 'real', ...)
     - TEMPERATURE   average temperature [K]
     - VOLUME        simulation cell volume [A^3]
    """

    _current_type = 'thermoelectric'
    _input_parameters = {'DT_FS', 'UNITS', 'TEMPERATURE', 'VOLUME'}
    _KAPPA_SI_UNITS = 'W/m/K'

    def __init__(
        self,
        traj=None,
        heat_current=None,
        charge_current=None,
        **params,
    ):
        if traj is None:
            if heat_current is None or charge_current is None:
                raise ValueError(
                    'Provide `traj` or both `heat_current` and `charge_current`.'
                )
            traj = np.array([heat_current, charge_current])
        elif (heat_current is not None) or (charge_current is not None):
            raise ValueError(
                'Use either `traj` or (`heat_current`, `charge_current`), not both.'
            )

        super().__init__(traj, **params)

        if self.N_CURRENTS != 2:
            raise ValueError(
                'ThermoElectricCurrent requires exactly two currents: heat and charge.'
            )

        self.thermoelectric_result = None

    @property
    def _builder(self):
        return dict(
            DT_FS=self.DT_FS,
            UNITS=self.UNITS,
            TEMPERATURE=self.TEMPERATURE,
            VOLUME=self.VOLUME,
            PSD_FILTER_W=self.PSD_FILTER_W_THZ,
            FREQ_UNITS='THz',
        )

    @staticmethod
    def _default_wishart_model(x, y):
        xx = np.concatenate([-x[::-1], x[1:]])
        yy = np.asarray(y).reshape(x.size, 3)
        yy = np.concatenate([yy[::-1], yy[1:]])
        return CubicSpline(xx, yy)

    def analyze_wishart(
        self,
        n_parameters='AIC',
        model=None,
        mask=None,
        solver='BFGS',
        guess_runave_window=50,
        minimize_kwargs=None,
        limits=None,
        omega_fixed=None,
        mc_samples=1000,
    ):
        """
        Run one-shot Wishart analysis and return thermoelectric coefficients.
        """
        if model is None:
            model = self._default_wishart_model

        self.maxlike_estimate(
            model=model,
            n_parameters=n_parameters,
            mask=mask,
            likelihood='wishart',
            solver=solver,
            guess_runave_window=guess_runave_window,
            minimize_kwargs=minimize_kwargs,
            ext_guess=None,
            limits=limits,
            omega_fixed=omega_fixed,
        )

        self.thermoelectric_result = self._build_thermoelectric_result(mc_samples)
        return self.thermoelectric_result

    def coefficients(self):
        if self.thermoelectric_result is None:
            raise RuntimeError('Run `analyze_wishart` first.')
        return {
            'sigma': self.thermoelectric_result.sigma,
            'sigma_std': self.thermoelectric_result.sigma_std,
            'kappa': self.thermoelectric_result.kappa,
            'kappa_std': self.thermoelectric_result.kappa_std,
            'seebeck': self.thermoelectric_result.seebeck,
            'seebeck_std': self.thermoelectric_result.seebeck_std,
        }

    def spectra(self):
        if self.thermoelectric_result is None:
            raise RuntimeError('Run `analyze_wishart` first.')
        return {
            'omega': self.thermoelectric_result.omega,
            'spectral_matrix_mean': self.thermoelectric_result.spectral_matrix_mean,
            'spectral_matrix_std': self.thermoelectric_result.spectral_matrix_std,
        }

    def coefficients_vs_frequency(
        self,
        source='wishart',
        units='si',
        with_uq=True,
        uq='mc',
        mc_samples=256,
        random_seed=None,
    ):
        """
        Compute sigma, kappa and Seebeck as functions of frequency.

        Parameters
        ----------
        source : {'raw', 'filtered', 'wishart'}
            Source spectrum used to build coefficients.
        units : {'si', 'plot'}
            Output units. 'si' gives (S/m, W/m/K, V/K), while 'plot' gives
            (S/cm, W/m/K, mV/K).
        with_uq : bool
            If True, include uncertainty if available.
        uq : {'mc'}
            UQ backend. Currently only Monte Carlo from parameter covariance.
        mc_samples : int
            Number of Monte Carlo samples for uncertainty estimation.
        random_seed : int or None
            Seed for Monte Carlo sampler.
        """
        source = source.lower()
        units = units.lower()
        uq = uq.lower()

        S, omega = self._spectral_tensor_by_source(source)
        coeff = self._coefficients_from_spectral_tensor(S, units=units)

        out = {
            'omega': omega,
            'sigma': coeff['sigma'],
            'kappa': coeff['kappa'],
            'seebeck': coeff['seebeck'],
        }

        if not with_uq:
            return out

        if uq != 'mc':
            raise NotImplementedError('Only MC uncertainty is currently supported.')

        if source == 'wishart':
            sampled = self._sample_coefficients_vs_frequency_mc(
                mc_samples=mc_samples,
                units=units,
                random_seed=random_seed,
            )
            if sampled is not None:
                out['sigma_std'] = sampled['sigma_std']
                out['kappa_std'] = sampled['kappa_std']
                out['seebeck_std'] = sampled['seebeck_std']

        return out

    def _build_thermoelectric_result(self, mc_samples):
        scales = self._thermoelectric_scales()
        omega = np.copy(self.maxlike.omega)
        spectral_mean = np.copy(self.maxlike.NLL_mean)

        spectral_std = None
        try:
            spectral_std = np.copy(self.maxlike.NLL_std)
        except AttributeError:
            pass

        S0 = spectral_mean[0]
        S0_std = None if spectral_std is None else spectral_std[0]

        coeff = self._coefficients_from_spectral_matrix(S0, scales)
        sigma = coeff['sigma']
        kappa = coeff['kappa']
        seebeck = coeff['seebeck']
        onsager_zero_mean = coeff['onsager_zero']

        sigma_std = None
        kappa_std = None
        seebeck_std = None
        coeff_cov = None
        onsager_zero_std = S0_std

        sampled = self._sample_coefficients_mc(scales, mc_samples)
        if sampled is not None:
            sigma_std = sampled['sigma_std']
            kappa_std = sampled['kappa_std']
            seebeck_std = sampled['seebeck_std']
            coeff_cov = sampled['coeff_cov']

        coeff_omega = self._coefficients_from_spectral_tensor(spectral_mean, units='si')

        sigma_omega_std = None
        kappa_omega_std = None
        seebeck_omega_std = None
        sampled_omega = self._sample_coefficients_vs_frequency_mc(
            mc_samples=mc_samples,
            units='si',
            random_seed=None,
        )
        if sampled_omega is not None:
            sigma_omega_std = sampled_omega['sigma_std']
            kappa_omega_std = sampled_omega['kappa_std']
            seebeck_omega_std = sampled_omega['seebeck_std']

        return ThermoelectricResult(
            omega=omega,
            spectral_matrix_mean=spectral_mean,
            spectral_matrix_std=spectral_std,
            onsager_zero_mean=onsager_zero_mean,
            onsager_zero_std=onsager_zero_std,
            sigma=float(sigma),
            sigma_std=None if sigma_std is None else float(sigma_std),
            kappa=float(kappa),
            kappa_std=None if kappa_std is None else float(kappa_std),
            seebeck=float(seebeck),
            seebeck_std=None if seebeck_std is None else float(seebeck_std),
            coefficients_cov=coeff_cov,
            sigma_omega=coeff_omega['sigma'],
            kappa_omega=coeff_omega['kappa'],
            seebeck_omega=coeff_omega['seebeck'],
            sigma_omega_std=sigma_omega_std,
            kappa_omega_std=kappa_omega_std,
            seebeck_omega_std=seebeck_omega_std,
            scales=scales,
            metadata={
                'temperature': self.TEMPERATURE,
                'volume': self.VOLUME,
                'units': self.UNITS,
                'n_currents': self.N_CURRENTS,
                'n_components': self.N_EQUIV_COMPONENTS,
                'n_parameters': self.maxlike.n_parameters,
            },
        )

    def _sample_coefficients_mc(self, scales, mc_samples):
        params_mean = getattr(self.maxlike, 'parameters_mean', None)
        params_cov = getattr(self.maxlike, 'parameters_cov', None)
        if params_mean is None or params_cov is None:
            return None

        try:
            samples = params_mean + np.random.multivariate_normal(
                mean=np.zeros_like(params_mean),
                cov=params_cov,
                size=mc_samples,
            )
        except (ValueError, np.linalg.LinAlgError):
            return None

        coeff_list = []
        for sample in samples:
            Sm = (
                scale_matrix(
                    self.maxlike.model,
                    sample,
                    self.maxlike.omega,
                    self.maxlike.omega_fixed,
                    self.N_CURRENTS,
                )
                / self.N_CURRENTS
            )
            coeff = self._coefficients_from_spectral_matrix(Sm[0], scales)
            coeff_list.append([coeff['sigma'], coeff['kappa'], coeff['seebeck']])

        coeff_arr = np.asarray(coeff_list)
        return {
            'sigma_std': coeff_arr[:, 0].std(axis=0),
            'kappa_std': coeff_arr[:, 1].std(axis=0),
            'seebeck_std': coeff_arr[:, 2].std(axis=0),
            'coeff_cov': np.cov(coeff_arr.T),
        }

    def _sample_coefficients_vs_frequency_mc(self, mc_samples, units='si', random_seed=None):
        params_mean = getattr(self.maxlike, 'parameters_mean', None)
        params_cov = getattr(self.maxlike, 'parameters_cov', None)
        if params_mean is None or params_cov is None:
            return None

        rng = np.random.default_rng(random_seed)
        try:
            samples = params_mean + rng.multivariate_normal(
                mean=np.zeros_like(params_mean),
                cov=params_cov,
                size=mc_samples,
            )
        except (ValueError, np.linalg.LinAlgError):
            return None

        sigma_list = []
        kappa_list = []
        seebeck_list = []
        for sample in samples:
            Sm = (
                scale_matrix(
                    self.maxlike.model,
                    sample,
                    self.maxlike.omega,
                    self.maxlike.omega_fixed,
                    self.N_CURRENTS,
                )
                / self.N_CURRENTS
            )
            coeff = self._coefficients_from_spectral_tensor(Sm, units=units)
            sigma_list.append(coeff['sigma'])
            kappa_list.append(coeff['kappa'])
            seebeck_list.append(coeff['seebeck'])

        sigma_arr = np.asarray(sigma_list)
        kappa_arr = np.asarray(kappa_list)
        seebeck_arr = np.asarray(seebeck_list)

        return {
            'sigma_std': sigma_arr.std(axis=0),
            'kappa_std': kappa_arr.std(axis=0),
            'seebeck_std': seebeck_arr.std(axis=0),
        }

    def _spectral_tensor_by_source(self, source):
        if source == 'raw':
            S = self.cospectrum.real.transpose((2, 0, 1)) / self.N_EQUIV_COMPONENTS
            omega = self.freqs_THz
            return S, omega
        if source == 'filtered':
            if self.fcospectrum is None:
                raise RuntimeError(
                    'Filtered cospectrum is not available. Set PSD_FILTER_W and recompute PSD.'
                )
            S = self.fcospectrum.real.transpose((2, 0, 1))
            omega = self.freqs_THz
            return S, omega
        if source == 'wishart':
            if self.thermoelectric_result is None:
                raise RuntimeError('Run `analyze_wishart` first.')
            S = self.thermoelectric_result.spectral_matrix_mean
            omega = self.freqs_THz[self.maxlike.mask[2]]
            return S, omega

        raise ValueError('`source` must be one of: raw, filtered, wishart')

    def _coefficients_from_spectral_tensor(self, spectral_tensor, units='si'):
        units = units.lower()
        if units not in ('si', 'plot'):
            raise ValueError("`units` must be either 'si' or 'plot'.")

        scales = self._thermoelectric_scales()
        kappa_scale = scales['kappa_scale']
        sigma_scale = scales['sigma_scale']
        mixed_scale = scales['mixed_scale']

        S = np.asarray(spectral_tensor)
        Lqq = S[..., 0, 0] * kappa_scale * self.TEMPERATURE
        Lqc = S[..., 0, 1] * mixed_scale
        Lcc = S[..., 1, 1] * sigma_scale

        eps = np.finfo(float).eps
        Lcc_safe = np.where(np.abs(Lcc) > eps, Lcc, np.nan)

        sigma = Lcc
        kappa = (Lqq - Lqc**2 / Lcc_safe) / self.TEMPERATURE
        seebeck = Lqc / (self.TEMPERATURE * Lcc_safe)

        if units == 'plot':
            sigma = sigma / 100.0
            seebeck = seebeck * 1.0e3

        return {
            'sigma': np.asarray(sigma),
            'kappa': np.asarray(kappa),
            'seebeck': np.asarray(seebeck),
        }

    def _thermoelectric_scales(self):
        unit = self.UNITS.lower()
        hscale = getattr(heat_units, 'scale_kappa_{}'.format(unit))
        escale = getattr(electric_units, 'scale_kappa_{}'.format(unit))

        heat_scale = 0.5 * hscale(self.TEMPERATURE, self.VOLUME)
        sigma_scale = 0.5 * escale(self.TEMPERATURE, self.VOLUME)
        mixed_scale = 0.5 * np.sqrt(
            hscale(self.TEMPERATURE, self.VOLUME)
            * escale(self.TEMPERATURE, self.VOLUME)
            * self.TEMPERATURE
        )
        return {
            'kappa_scale': float(heat_scale),
            'sigma_scale': float(sigma_scale),
            'mixed_scale': float(mixed_scale),
        }

    def _coefficients_from_spectral_matrix(self, spectral_matrix, scales):
        kappa_scale = scales['kappa_scale']
        sigma_scale = scales['sigma_scale']
        mixed_scale = scales['mixed_scale']

        Lqq = spectral_matrix[0, 0] * kappa_scale * self.TEMPERATURE
        Lqc = spectral_matrix[0, 1] * mixed_scale
        Lcc = spectral_matrix[1, 1] * sigma_scale

        onsager_zero = np.array([[Lqq, Lqc], [Lqc, Lcc]], dtype=float)

        kappa = (Lqq - Lqc**2 / Lcc) / self.TEMPERATURE
        sigma = Lcc
        seebeck = Lqc / (self.TEMPERATURE * Lcc)

        return {
            'onsager_zero': onsager_zero,
            'sigma': float(sigma),
            'kappa': float(kappa),
            'seebeck': float(seebeck),
        }
