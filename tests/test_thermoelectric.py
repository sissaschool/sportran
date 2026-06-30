# -*- coding: utf-8 -*-

import numpy as np


def _toy_currents(n_steps=1024, n_components=3, seed=7):
    rng = np.random.default_rng(seed)
    base = rng.normal(size=(n_steps, n_components))
    heat = 2.0 * base + 0.3 * rng.normal(size=(n_steps, n_components))
    charge = 0.5 * base + 0.5 * rng.normal(size=(n_steps, n_components))
    return heat, charge


def test_thermoelectric_current_init_named_inputs():
    import sportran as st

    j_q, j_c = _toy_currents()
    current = st.ThermoElectricCurrent(
        heat_current=j_q,
        charge_current=j_c,
        DT_FS=1.0,
        UNITS='metal',
        TEMPERATURE=900.0,
        VOLUME=1000.0,
    )
    assert current.N_CURRENTS == 2
    assert current.MANY_CURRENTS is True


def test_thermoelectric_wishart_analysis_outputs():
    import sportran as st

    j_q, j_c = _toy_currents(n_steps=2048, seed=11)
    current = st.ThermoElectricCurrent(
        heat_current=j_q,
        charge_current=j_c,
        DT_FS=1.0,
        UNITS='metal',
        TEMPERATURE=900.0,
        VOLUME=1000.0,
    )

    result = current.analyze_wishart(
        n_parameters=4,
        mask=(slice(None), slice(None), slice(None, None, 32)),
        minimize_kwargs={
            'options': {'maxiter': 80, 'disp': False},
        },
        mc_samples=64,
    )

    assert result is current.thermoelectric_result
    assert np.isfinite(result.sigma)
    assert np.isfinite(result.kappa)
    assert np.isfinite(result.seebeck)
    assert result.onsager_zero_mean.shape == (2, 2)
    assert result.spectral_matrix_mean.shape[1:] == (2, 2)


def test_thermoelectric_coefficients_vs_frequency_api():
    import sportran as st

    j_q, j_c = _toy_currents(n_steps=2048, seed=17)
    current = st.ThermoElectricCurrent(
        heat_current=j_q,
        charge_current=j_c,
        DT_FS=1.0,
        UNITS='metal',
        TEMPERATURE=900.0,
        VOLUME=1000.0,
        PSD_FILTER_W=0.1,
        FREQ_UNITS='THz',
    )

    current = current.resample(fstar_THz=40.0, plot=False, freq_units='thz')
    current.analyze_wishart(
        n_parameters=4,
        mask=(slice(None), slice(None), slice(None, None, 4)),
        minimize_kwargs={'options': {'maxiter': 80, 'disp': False}},
        mc_samples=32,
    )

    raw_plot = current.coefficients_vs_frequency(source='raw', units='plot', with_uq=False)
    filt_plot = current.coefficients_vs_frequency(
        source='filtered', units='plot', with_uq=False
    )
    wishart_si = current.coefficients_vs_frequency(
        source='wishart',
        units='si',
        with_uq=True,
        uq='mc',
        mc_samples=32,
        random_seed=123,
    )
    wishart_plot = current.coefficients_vs_frequency(
        source='wishart',
        units='plot',
        with_uq=True,
        uq='mc',
        mc_samples=32,
        random_seed=123,
    )

    assert raw_plot['omega'].shape[0] == raw_plot['sigma'].shape[0]
    assert filt_plot['omega'].shape[0] == filt_plot['kappa'].shape[0]
    assert wishart_si['omega'].shape[0] == wishart_si['seebeck'].shape[0]

    assert 'sigma_std' in wishart_si
    assert wishart_si['sigma_std'].shape == wishart_si['sigma'].shape
    assert wishart_si['kappa_std'].shape == wishart_si['kappa'].shape
    assert wishart_si['seebeck_std'].shape == wishart_si['seebeck'].shape

    # unit consistency: plot units are S/cm and mV/K
    assert np.allclose(wishart_plot['sigma'] * 100.0, wishart_si['sigma'])
    assert np.allclose(wishart_plot['seebeck'] / 1000.0, wishart_si['seebeck'])
