# -*- coding: utf-8 -*-

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np


@dataclass
class ThermoelectricResult:
    """
    Container for thermoelectric Wishart analysis outputs.
    """

    omega: np.ndarray
    spectral_matrix_mean: np.ndarray
    spectral_matrix_std: Optional[np.ndarray]
    onsager_zero_mean: np.ndarray
    onsager_zero_std: Optional[np.ndarray]
    sigma: float
    sigma_std: Optional[float]
    kappa: float
    kappa_std: Optional[float]
    seebeck: float
    seebeck_std: Optional[float]
    coefficients_cov: Optional[np.ndarray]
    sigma_omega: np.ndarray
    kappa_omega: np.ndarray
    seebeck_omega: np.ndarray
    sigma_omega_std: Optional[np.ndarray]
    kappa_omega_std: Optional[np.ndarray]
    seebeck_omega_std: Optional[np.ndarray]
    scales: Dict[str, float]
    metadata: Dict[str, object]

    def to_dict(self):
        return {
            'omega': self.omega,
            'spectral_matrix_mean': self.spectral_matrix_mean,
            'spectral_matrix_std': self.spectral_matrix_std,
            'onsager_zero_mean': self.onsager_zero_mean,
            'onsager_zero_std': self.onsager_zero_std,
            'sigma': self.sigma,
            'sigma_std': self.sigma_std,
            'kappa': self.kappa,
            'kappa_std': self.kappa_std,
            'seebeck': self.seebeck,
            'seebeck_std': self.seebeck_std,
            'coefficients_cov': self.coefficients_cov,
            'sigma_omega': self.sigma_omega,
            'kappa_omega': self.kappa_omega,
            'seebeck_omega': self.seebeck_omega,
            'sigma_omega_std': self.sigma_omega_std,
            'kappa_omega_std': self.kappa_omega_std,
            'seebeck_omega_std': self.seebeck_omega_std,
            'scales': self.scales,
            'metadata': self.metadata,
        }
