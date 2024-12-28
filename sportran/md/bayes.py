# -*- coding: utf-8 -*-
# Methods to perform a bayesian estimation of the transport coefficients

import numpy as np
import emcee
import scipy.special as sp
from . import aic
from .cepstral import (dct_coefficients, dct_filter_psd, dct_filter_tau, CepstralFilter, multicomp_cepstral_parameters,)
from .tools.filter import runavefilter
from sportran.utils import log
from multiprocessing import Pool
import time

__all__ = ['BayesFilter']
EULER_GAMMA = (
    0.57721566490153286060651209008240243104215933593992   # Euler-Mascheroni constant
)
LOG2 = np.log(2)


class BayesFilter(object):
    """
    BAYESIAN ANALYSIS based filtering.

    ** INPUT VARIABLES:
    spectrum        = the original periodogram (if single-component) of spectral matrix (if multi-component)
    is_offdiag      = If True, estimate the off-diagonal matrix element of the spectral matrix (default = True)
    is_diag         = If True, estimate the diagonal matrix elements of the spectral matrix (default = False)
    model           = the function that models the data (for now only spline)
    n_parameters    = the number of parameters to be used for the fit

    ** INTERNAL VARIABLES:
    samplelogpsd  = the original sample log-PSD - logpsd_THEORY_mean

    logpsdK  = the cepstrum of the data, \\hat{C}_n (i.e. the DCT of samplelogpsd)
    aic_min  = minimum value of the AIC
    aic_Kmin = cutoffK that minimizes the AIC
    aic_Kmin_corrfactor = aic_Kmin cutoff correction factor (default: 1.0)
    cutoffK  = (P*-1) = cutoff used to compute logtau and logpsd (by default = aic_Kmin * aic_Kmin_corrfactor)
    manual_cutoffK_flag = True if cutoffK was manually specified, False if aic_Kmin is being used

    logtau          = filtered log(tau) as a function of cutoffK, L_0(P*-1)
    logtau_cutoffK  = filtered log(tau) at cutoffK, L*_0
    logtau_var_cutoffK = theoretical L*_0 variance
    logtau_std_cutoffK = theoretical L*_0 standard deviation
    logpsd          = filtered log-PSD at cutoffK

    tau          = filtered tau as a function of cutoffK, S_0(P*-1)
    tau_cutoffK  = filtered tau at cutoffK, S*_0
    tau_var_cutoffK = theoretical S*_0 variance
    tau_std_cutoffK = theoretical S*_0 standard deviation
    psd          = filtered PSD at the specified cutoffK

    p_aic... = Bayesian AIC weighting stuff
    """

    def __init__(self, spectrum, model, n_parameters, n_components, is_restart=False, n_steps=2000000,
                 backend='chain.h5', burn_in=None, thin=None, mask=None,
                ):

        if not isinstance(spectrum, np.ndarray):
            raise TypeError('spectrum should be an object of type numpy.ndarray')
        if spectrum.shape[0] != 2 or spectrum.shape[1] != 2:
            raise TypeError('spectrum should be a 2x2xN numpy.ndarray')

        self.spectrum = spectrum / n_components
        self.model = model
        self.mask = mask
        self.n_components = n_components
        self.n_parameters = n_parameters
        self.is_restart = is_restart
        self.n_steps = n_steps
        self.backend = backend
        self.burn_in = burn_in
        self.thin = thin

    def __repr__(self):
        msg = 'BayesFilter:\n'   # + \
        #   '  AIC type  = {:}\n'.format(self.aic_type) + \
        #   '  AIC min   = {:f}\n'.format(self.aic_min) + \
        #   '  AIC_Kmin  = {:d}\n'.format(self.aic_Kmin)
        # if self.cutoffK is not None:
        #     msg += \
        #         '  AIC_Kmin_corrfactor = {:f}\n'.format(self.aic_Kmin_corrfactor) + \
        #         '  cutoffK = (P*-1) = {:d} {:}\n'.format(self.cutoffK, '(manual)' if self.manual_cutoffK_flag else '(auto)') + \
        #         '  L_0*   = {:15f} +/- {:10f}\n'.format(self.logtau_cutoffK, self.logtau_std_cutoffK) + \
        #         '  S_0*   = {:15f} +/- {:10f}\n'.format(self.tau_cutoffK, self.tau_std_cutoffK)
        return msg

    ################################################

    def run_mcmc(self, n_parameters=None, n_steps=None, is_restart=None, mask=None, filename=None, n_walkers=None,
                 log_like='off',
                ):

        # Initialize the parameters if undefined
        if n_parameters is None:
            n_parameters = self.n_parameters
        if n_steps is None:
            n_steps = self.n_steps
        if is_restart is None:
            is_restart = self.is_restart
        if mask is None:
            mask = self.mask
        if filename is None:
            filename = self.backend

        # Initialize the parameters for cepstral analysis
        ck_THEORY_var, psd_THEORY_mean = multicomp_cepstral_parameters(self.spectrum.shape[2], self.n_components)
        # Cepstral analysis for the diagonal elements
        cepf1 = CepstralFilter(np.log(self.spectrum[0, 0].real), ck_theory_var=ck_THEORY_var,
                               psd_theory_mean=psd_THEORY_mean, aic_type='aic',
                              )
        cepf1.scan_filter_tau(cutoffK=None, aic_Kmin_corrfactor=1.0)
        cepf2 = CepstralFilter(np.log(self.spectrum[1, 1].real), ck_theory_var=ck_THEORY_var,
                               psd_theory_mean=psd_THEORY_mean, aic_type='aic',
                              )
        cepf2.scan_filter_tau(cutoffK=None, aic_Kmin_corrfactor=1.0)
        self.sigma1 = cepf1.psd[0]
        self.sigma2 = cepf2.psd[0]

        nu = 2
        ell = self.n_components
        # Define noisy data
        if mask is not None:
            noisy_data = (self.spectrum.real[0, 1] / np.sqrt(cepf1.psd * cepf2.psd))[mask]
        else:
            noisy_data = self.spectrum.real[0, 1] / np.sqrt(cepf1.psd * cepf2.psd)

        self.factor = np.sqrt(cepf1.psd[0] / cepf2.psd[0])

        # Define initial points for the MCMC
        try:
            guess_data = runavefilter(noisy_data, 100)
        except:
            guess_data = runavefilter(noisy_data, 10)

        args = np.int32(np.linspace(0, len(noisy_data) - 1, n_parameters, endpoint=True))

        # MCMC sampling
        # number of walkers must be larger than twice the number of parameters (and often a power of 2)
        if n_walkers is None:
            n_walkers = int(2**np.ceil(np.log2(2 * n_parameters)))

        log.write_log('MCMC with {} parameters and {} walkers'.format(n_parameters, n_walkers))
        log.write_log(f'Running up to {n_steps} steps')

        p0 = guess_data
        p0 = np.clip(p0[args][np.newaxis, :n_parameters] + np.random.normal(0, 0.1, (n_walkers, n_parameters)), -0.98,
                     0.98,
                    )

        omega = np.arange(noisy_data.size)
        omega_fixed = omega[args]

        self.omega = omega
        self.omega_fixed = omega_fixed

        # Set up the backend
        # Don't forget to clear it in case the file already exists
        backend = emcee.backends.HDFBackend(filename)
        if not is_restart:
            backend.reset(n_walkers, n_parameters)

        # Initialize the sampler
        if log_like == 'off':
            sampler = emcee.EnsembleSampler(n_walkers, n_parameters, self.log_posterior_offdiag,
                                            args=(omega, omega_fixed, noisy_data, nu, ell), backend=backend,
                                           )
        elif log_like == 'normal':
            sampler = emcee.EnsembleSampler(n_walkers, n_parameters, self.log_posterior_normal,
                                            args=(omega, omega_fixed, noisy_data, nu, ell), backend=backend,
                                           )

        # Run MCMC
        # We'll track how the average autocorrelation time estimate changes
        index = 0
        autocorr = np.empty(n_steps)

        # This will be useful to testing convergence
        old_tau = np.inf

        # Now we'll sample for up to max_n steps
        if is_restart:
            coord = backend.get_chain()[-1]
            good_idx = None
            tau = backend.get_autocorr_time(discard=500)
            burn_in = int(2 * np.max(tau))
            thin = np.max([1, int(0.5 * np.min(tau))])
            log.write_log('MCMC autocorrelation time = {}'.format(tau))
            if good_idx is None:
                samples = backend.get_chain(discard=burn_in, flat=True, thin=thin)
            else:
                samples = backend.get_chain(discard=burn_in, flat=True, thin=thin)[:, good_idx]
                n_parameters_ = samples.shape[1]
                self.n_parameters = n_parameters_

            # Compute marginalized errors
            rho = []
            rho_min = []
            rho_max = []
            print(sampler.iteration)
            self.acceprance_fraction = np.mean(sampler.acceptance_fraction)
            for i in range(self.n_parameters):
                mcmc = np.percentile(samples[:, i], [16, 50, 84])
                q = np.diff(mcmc)
                rho.append(mcmc[1])
                rho_min.append(mcmc[1] - q[0])
                rho_max.append(mcmc[1] + q[1])

            # The estimated parameters are the values of rho as a function of frequency
            self.parameters_mean = np.array(rho)
            self.parameters_args = args
            self.parameters_std = 0.5 * (np.array(rho_max) - np.array(rho_min))
            self.sampler = sampler
            self.noisy_data = noisy_data
            if log_like == 'normal':
                self.aic = (
                    2 * self.log_likelihood_normal(self.parameters_mean, omega, omega_fixed, noisy_data, nu, ell) -
                    2 * self.n_parameters)
            elif log_like == 'off':
                self.aic = (
                    2 * self.log_likelihood_offdiag(self.parameters_mean, omega, omega_fixed, noisy_data, nu, ell) -
                    2 * self.n_parameters)
                burn_in = int(2 * np.max(tau))
                thin = np.max([1, int(0.5 * np.min(tau))])
            return
        else:
            coord = np.copy(p0)

        todo = True
        if is_restart:
            todo = False
        disc = 0
        for sample in sampler.sample(coord, iterations=n_steps, progress=True, store=True):
            # Only check convergence every 100 steps
            if sampler.iteration % 250:
                continue

            # Compute the autocorrelation time so far
            # Using tol=0 means that we'll always get an estimate even
            # if it isn't trustworthy
            tau = sampler.get_autocorr_time(tol=0, discard=disc)
            autocorr[index] = np.mean(tau)
            index += 1
            if sampler.iteration % 2000 == 0:
                print(tau)

            if todo and sampler.iteration > 1000:
                s_old = np.ones(n_parameters)
                for i in range(100, int(sampler.iteration / 2) + 1, 100):

                    s = sampler.get_autocorr_time(tol=0, discard=i)

                    if np.all(abs((s - s_old) / s) * 100 < 5):
                        disc = i
                        todo = False
                        print(s)
                        break
                    s_old = s

            # Check convergence
            converged = np.all(tau * 100 < sampler.iteration)
            converged &= np.all(np.abs(old_tau - tau) / tau < 0.015)
            if converged:
                break
            old_tau = tau

        # Compute chains auto-correlation time to estimate convergence
        # If AutocorrError, probably the chain is too short. You can still use ~2*max(tau) as burn_in
        good_idx = None
        try:
            tau = sampler.get_autocorr_time(discard=disc)
            burn_in = int(6 * np.max(tau))
            thin = np.max([1, int(1.5 * np.min(tau))])
            log.write_log('MCMC autocorrelation time = {}'.format(tau))
        except emcee.autocorr.AutocorrError:
            log.write_log('The chain is probably too short')
            burn_in = int(sampler.iteration * 0.3)
            thin = int(np.max([int(0.05 * sampler.iteration), 10]))
        except ValueError:
            log.write_log(f'There is something wrong with tau: tau = {tau}')
            good_idx = ~np.isnan(tau)
            tau = tau[good_idx]
            log.write_log('Fixed MCMC autocorrelation time = {}'.format(tau))
            burn_in = int(6 * np.max(tau))
            thin = np.max([1, int(1.5 * np.min(tau))])

        log.write_log('MCMC burn in = {}; thin = {}'.format(burn_in, thin))

        if good_idx is None:
            samples = sampler.get_chain(discard=burn_in, flat=True, thin=thin)
        else:
            samples = sampler.get_chain(discard=burn_in, flat=True, thin=thin)[:, good_idx]
            n_parameters_ = samples.shape[1]
            self.n_parameters = n_parameters_

        # Compute marginalized errors
        rho = []
        rho_min = []
        rho_max = []
        for i in range(self.n_parameters):
            mcmc = np.percentile(samples[:, i], [16, 50, 84])
            q = np.diff(mcmc)
            rho.append(mcmc[1])
            rho_min.append(mcmc[1] - q[0])
            rho_max.append(mcmc[1] + q[1])

        # The estimated parameters are the values of rho as a function of frequency
        self.parameters_mean = np.array(rho)
        self.parameters_args = args
        self.parameters_std = 0.5 * (np.array(rho_max) - np.array(rho_min))
        self.sampler = sampler
        self.noisy_data = noisy_data
        if log_like == 'normal':
            self.aic = (2 * self.log_likelihood_normal(self.parameters_mean, omega, omega_fixed, noisy_data, nu, ell) -
                        2 * self.n_parameters)
        elif log_like == 'off':
            self.aic = (2 * self.log_likelihood_offdiag(self.parameters_mean, omega, omega_fixed, noisy_data, nu, ell) -
                        2 * self.n_parameters)
            self.dic = (
                -2 * self.log_likelihood_offdiag(self.parameters_mean, omega, omega_fixed, noisy_data, nu, ell) +
                4 * sampler.get_log_prob(discard=burn_in, flat=True, thin=thin).mean())
        with open('aic_{}'.format(self.n_parameters), 'w+') as g:
            g.write('{}\t{}\n'.format(self.n_parameters, self.aic))
        with open('dic_{}'.format(self.n_parameters), 'w+') as g:
            g.write('{}\t{}\n'.format(self.n_parameters, self.dic))

    ################################################
    # Helper functions

    # The log-likelihood function
    # def log_likelihood_offdiag(self, w, omega, omega_fixed, data_, nu, ell):
    #     spline = self.model(omega_fixed, w)
    #     rho = np.clip(spline(omega), -0.99, 0.99)

    #     one_frac_rho2 = 1/(1-rho**2)

    #     # Data is distributed according to a Variance-Gamma distribution with parameters:
    #     # mu = 0; alpha = 1/(1-rho**2); beta = rho/(1-rho**2); lambda = ell*nu/2
    #     # Its expectation value is ell*nu*rho
    #     data = ell*data_
    #     z = data - ell*nu*rho

    #     log_pdf = -np.log(sp.gamma(0.5*nu)) + 0.5*(nu-1)*np.log(np.abs(z)) - 0.5*np.log(2**(nu-1)*np.pi/one_frac_rho2) + rho*z*one_frac_rho2 +\
    #             np.log(sp.kv(0.5*(nu-1), np.abs(z)*one_frac_rho2))
    #     return np.sum(log_pdf)

    def run_mcmc_scratch(self, n_parameters=None, n_steps=None, is_restart=None, mask=None, filename=None,
                         n_walkers=None, log_like='off',
                        ):

        # Initialize the parameters if undefined
        if n_parameters is None:
            n_parameters = self.n_parameters
        else:
            self.n_parameters = n_parameters
        if n_steps is None:
            n_steps = self.n_steps
        if is_restart is None:
            is_restart = self.is_restart
        if mask is None:
            mask = self.mask
        if filename is None:
            filename = self.backend

        # Initialize the parameters for cepstral analysis

        nu = 2
        ell = self.n_components
        # Define noisy data
        if mask is not None:
            noisy_data = (self.spectrum[0, 1])[mask]
        else:
            noisy_data = self.spectrum[0, 1]

        # Define initial points for the MCMC
        try:
            guess_data = runavefilter(noisy_data, 200)
        except:
            guess_data = runavefilter(noisy_data, 100)

        args = np.int32(np.linspace(0, len(noisy_data) - 1, n_parameters, endpoint=True))

        # MCMC sampling
        # number of walkers must be larger than twice the number of parameters (and often a power of 2)
        if n_walkers is None:
            n_walkers = int(2**np.ceil(np.log2(2 * n_parameters)))

        log.write_log('MCMC with {} parameters and {} walkers'.format(n_parameters, n_walkers))
        log.write_log(f'Running up to {n_steps} steps')

        p0 = guess_data
        if log_like == 'normal' or log_like == 'off':
            p0 = np.clip(p0[args][np.newaxis, :n_parameters] + np.random.normal(0, 0.2, (n_walkers, n_parameters)), -40,
                         40,
                        )

        if log_like == 'diag':
            p0 = np.clip(p0[args][np.newaxis, :n_parameters] + np.random.normal(0, 10, (n_walkers, n_parameters)), 0,
                         1e6,
                        )

        omega = np.arange(noisy_data.size)
        omega_fixed = omega[args]

        self.omega = omega
        self.omega_fixed = omega_fixed

        # Set up the backend
        # Don't forget to clear it in case the file already exists
        backend = emcee.backends.HDFBackend(filename)
        if not is_restart:
            backend.reset(n_walkers, n_parameters)

        # Initialize the sampler
        if log_like == 'off':
            sampler = emcee.EnsembleSampler(n_walkers, n_parameters, self.log_posterior_offdiag,
                                            args=(omega, omega_fixed, noisy_data, nu, ell), backend=backend,
                                           )
        elif log_like == 'normal':
            sampler = emcee.EnsembleSampler(n_walkers, n_parameters, self.log_posterior_normal,
                                            args=(omega, omega_fixed, noisy_data, nu, ell), backend=backend,
                                           )

        elif log_like == 'diag':
            sampler = emcee.EnsembleSampler(n_walkers, n_parameters, self.log_posterior_diag,
                                            args=(omega, omega_fixed, noisy_data, ell), backend=backend,
                                           )

        # Run MCMC
        # We'll track how the average autocorrelation time estimate changes
        index = 0
        autocorr = np.empty(n_steps)

        # This will be useful to testing convergence
        old_tau = np.inf

        # Now we'll sample for up to max_n steps
        if is_restart:
            coord = backend.get_chain()[-1]
            good_idx = None
            tau = backend.get_autocorr_time(discard=1000)
            burn_in = int(6 * np.max(tau))
            thin = np.max([1, int(1.5 * np.min(tau))])
            log.write_log('MCMC autocorrelation time = {}'.format(tau))
            if good_idx is None:
                samples = backend.get_chain(discard=burn_in, flat=True, thin=thin)
            else:
                samples = backend.get_chain(discard=burn_in, flat=True, thin=thin)[:, good_idx]
                n_parameters_ = samples.shape[1]
                self.n_parameters = n_parameters_

            # Compute marginalized errors
            rho = []
            rho_min = []
            rho_max = []
            for i in range(self.n_parameters):
                mcmc = np.percentile(samples[:, i], [16, 50, 84])
                # mcmc = np.percentile(samples[:, i], [50-95/2, 50, 50+95/2])
                q = np.diff(mcmc)
                # print(mcmc[1], q[0], q[1])
                rho.append(mcmc[1])
                rho_min.append(mcmc[1] - q[0])
                rho_max.append(mcmc[1] + q[1])
            # The estimated parameters are the values of rho as a function of frequency
            self.parameters_mean = np.array(rho)
            self.parameters_args = args
            self.parameters_std = 0.5 * (np.array(rho_max) - np.array(rho_min))
            self.sampler = sampler
            self.noisy_data = noisy_data
            if log_like == 'normal':
                self.aic = (
                    2 * self.log_likelihood_normal(self.parameters_mean, omega, omega_fixed, noisy_data, nu, ell) -
                    2 * self.n_parameters)
            elif log_like == 'off':
                self.aic = (
                    2 * self.log_likelihood_offdiag(self.parameters_mean, omega, omega_fixed, noisy_data, nu, ell) -
                    2 * self.n_parameters)
                burn_in = int(2 * np.max(tau))
                thin = np.max([1, int(0.5 * np.min(tau))])
            elif log_like == 'diag':
                self.aic = (2 * self.log_likelihood_diag(self.parameters_mean, omega, omega_fixed, noisy_data, ell) -
                            2 * self.n_parameters)
                burn_in = int(2 * np.max(tau))
                thin = np.max([1, int(0.5 * np.min(tau))])
            return
        else:
            coord = np.copy(p0)

        todo = True
        disc = 0
        for sample in sampler.sample(coord, iterations=n_steps, progress=True, store=True):
            # Only check convergence every 100 steps
            if sampler.iteration % 250:
                continue

            # Compute the autocorrelation time so far
            if sampler.iteration % 1000 == 0:
                print(tau)

            # Using tol=0 means that we'll always get an estimate even
            # if it isn't trustworthy
            tau = sampler.get_autocorr_time(tol=0, discard=disc)
            autocorr[index] = np.mean(tau)
            index += 1

            # Check convergence
            converged = np.all(tau * 100 < sampler.iteration)
            converged &= np.all(np.abs(old_tau - tau) / tau < 0.015)

            if todo and sampler.iteration % 500 == 0 and sampler.iteration > 1000:
                s_old = np.ones(n_parameters)
                for i in range(100, int(sampler.iteration / 2) + 1, 100):

                    s = sampler.get_autocorr_time(tol=0, discard=i)

                    if np.all(abs((s - s_old) / s) * 100 < 2):
                        disc = i
                        todo = False
                        break
                    s_old = s
            if converged:
                break
            old_tau = tau

        # Compute chains auto-correlation time to estimate convergence
        # If AutocorrError, probably the chain is too short. You can still use ~2*max(tau) as burn_in
        good_idx = None
        try:
            print(disc, ' discard')
            tau = sampler.get_autocorr_time(discard=disc)
            burn_in = max(int(6 * np.max(tau)), disc)
            thin = np.max([1, int(1.5 * np.min(tau))])
            log.write_log('MCMC autocorrelation time = {}'.format(tau))
        except emcee.autocorr.AutocorrError:
            log.write_log('The chain is probably too short')
            burn_in = max(int(sampler.iteration * 0.3), disc)
            thin = int(np.max([int(1.5 * sampler.iteration), 10]))
        except ValueError:
            log.write_log(f'There is something wrong with tau: tau = {tau}')
            good_idx = ~np.isnan(tau)
            tau = tau[good_idx]
            log.write_log('Fixed MCMC autocorrelation time = {}'.format(tau))
            burn_in = max(int(2 * np.max(tau)), disc)
            thin = np.max([1, int(1.5 * np.min(tau))])
        # if self.burn_in is not None:
        #     burn_in = self.burn_in
        # else:
        #     self.burn_in = burn_in
        # if self.thin is not None:
        #     thin = self.thin
        # else:
        #     self.thin = thin
        log.write_log('MCMC burn in = {}; thin = {}'.format(burn_in, thin))

        if good_idx is None:
            samples = sampler.get_chain(discard=burn_in, flat=True, thin=thin)
        else:
            samples = sampler.get_chain(discard=burn_in, flat=True, thin=thin)[:, good_idx]
            n_parameters_ = samples.shape[1]
            self.n_parameters = n_parameters_

        # Compute marginalized errors
        rho = []
        rho_min = []
        rho_max = []
        for i in range(self.n_parameters):
            mcmc = np.percentile(samples[:, i], [16, 50, 84])
            # mcmc = np.percentile(samples[:, i], [50-95/2, 50, 50+95/2])
            q = np.diff(mcmc)
            # print(mcmc[1], q[0], q[1])
            rho.append(mcmc[1])
            rho_min.append(mcmc[1] - q[0])
            rho_max.append(mcmc[1] + q[1])

        # The estimated parameters are the values of rho as a function of frequency
        self.parameters_mean = np.array(rho)
        self.parameters_args = args
        self.parameters_std = 0.5 * (np.array(rho_max) - np.array(rho_min))
        self.sampler = sampler
        self.noisy_data = noisy_data
        if log_like == 'normal':
            self.aic = (2 * self.log_likelihood_normal(self.parameters_mean, omega, omega_fixed, noisy_data, nu, ell) -
                        2 * self.n_parameters)
        elif log_like == 'off':
            self.aic = (2 * self.log_likelihood_offdiag(self.parameters_mean, omega, omega_fixed, noisy_data, nu, ell) -
                        2 * self.n_parameters)
            self.dic = (
                -2 * self.log_likelihood_offdiag(self.parameters_mean, omega, omega_fixed, noisy_data, nu, ell) +
                4 * sampler.get_log_prob(discard=burn_in, flat=True, thin=thin).mean())
        elif log_like == 'diag':
            self.aic = (2 * self.log_likelihood_diag(self.parameters_mean, omega, omega_fixed, noisy_data, ell) -
                        2 * self.n_parameters)
            self.dic = (-2 * self.log_likelihood_diag(self.parameters_mean, omega, omega_fixed, noisy_data, ell) +
                        4 * sampler.get_log_prob(discard=burn_in, flat=True, thin=thin).mean())
        with open('aic_{}'.format(self.n_parameters), 'w+') as g:
            g.write('{}\t{}\n'.format(self.n_parameters, self.aic))
        with open('dic_{}'.format(self.n_parameters), 'w+') as g:
            g.write('{}\t{}\n'.format(self.n_parameters, self.dic))

    ################################################
    # Helper functions

    # The log-likelihood function
    # def log_likelihood_offdiag(self, w, omega, omega_fixed, data_, nu, ell):
    #     spline = self.model(omega_fixed, w)
    #     rho = np.clip(spline(omega), -0.99, 0.99)

    #     one_frac_rho2 = 1/(1-rho**2)

    #     # Data is distributed according to a Variance-Gamma distribution with parameters:
    #     # mu = 0; alpha = 1/(1-rho**2); beta = rho/(1-rho**2); lambda = ell*nu/2
    #     # Its expectation value is ell*nu*rho
    #     data = ell*data_
    #     z = data - ell*nu*rho

    #     log_pdf = -np.log(sp.gamma(0.5*nu)) + 0.5*(nu-1)*np.log(np.abs(z)) - 0.5*np.log(2**(nu-1)*np.pi/one_frac_rho2) + rho*z*one_frac_rho2 +\
    #             np.log(sp.kv(0.5*(nu-1), np.abs(z)*one_frac_rho2))
    #     return np.sum(log_pdf)

    def log_likelihood_wishart(self, w, omega, omega_fixed, data_, nu, ell):
        """
        Logarithm of the Wishart probability density function.
        """
        from scipy.special import multigammaln

        spline = self.model(omega_fixed, w)
        V = spline(omega)
        # TODO: convert the notation from wikipedia to Baroni
        X = data_ * ell * nu
        n = ell
        p = 2
        a, b, d = X[..., 0, 0], X[..., 0, 1], X[..., 1, 1]
        detX = a * d - b**2

        # detX = np.linalg.det(X)

        a, b, d = V[..., 0, 0], V[..., 0, 1], V[..., 1, 1]
        invV = (1 / (a * d - b**2) * np.array([[d, -b], [-b, a]])).transpose(2, 0, 1)
        trinvV_X = np.trace(invV @ X, axis1=1, axis2=2)
        detV = a * d - b**2
        # trinvV_X = np.trace(np.linalg.inv(V)@X)

        log_pdf = 0.5 * (
            (n - p - 1) * np.log(detX) - trinvV_X - n * p * LOG2 - n * np.log(detV)) - multigammaln(0.5 * n, 2)

        return np.sum(log_pdf)

    def log_likelihood_offdiag(self, w, omega, omega_fixed, data_, nu, ell):
        """
        Logarithm of the Variance-Gamma probability density function.
        """
        spline = self.model(omega_fixed, w)
        rho = np.clip(spline(omega), -0.98, 0.98)
        _alpha = 1 / (1 - rho**2)
        _beta = rho / (1 - rho**2)
        _lambda = 0.5 * ell * nu
        _gamma2 = _alpha**2 - _beta**2
        _lambda_minus_half = _lambda - 0.5

        # Data is distributed according to a Variance-Gamma distribution with parameters (notation as in Wikipedia):
        # mu = 0; alpha = 1/(1-rho**2); beta = rho/(1-rho**2); lambda = ell*nu/2
        # Its expectation value is ell*nu*rho
        z = data_ * ell * nu
        absz = np.abs(z)
        # z = data
        log_pdf = (_lambda * np.log(_gamma2) + _lambda_minus_half * np.log(absz) +
                   np.log(sp.kv(_lambda_minus_half, _alpha * absz)) + _beta * z - 0.5 * np.log(np.pi) -
                   np.log(sp.gamma(_lambda)) - _lambda_minus_half * np.log(2 * _alpha))

        res = np.sum(log_pdf)
        return res

    def log_likelihood_diag(self, w, omega, omega_fixed, data_, ell):
        spline = self.model(omega_fixed, w)
        rho = np.clip(spline(omega), 1e-6, 1e6)

        # Data is distributed according to a Chi-squared distribution with parameters (notation as in Wikipedia):
        # Its expectation value is ell*rho
        z = data_ * ell / rho
        absz = np.abs(z)
        # z = data
        log_pdf = (ell / 2 - 1) * np.log(absz) - absz / 2 - np.log(rho)

        res = np.sum(log_pdf)
        return res

    def log_likelihood_normal(self, w, omega, omega_fixed, data_, nu, ell):
        spline = self.model(omega_fixed, w)
        rho = np.clip(spline(omega), -0.98, 0.98)

        log_pdf = -((data_ - rho)**2)
        return np.sum(log_pdf)

    # The log-prior function
    def log_prior_offdiag(self, w):
        # Uniform prior
        if np.all((w >= -1) & (w <= 1)):
            return 1
        else:
            return -np.inf

    # The log-prior function
    def log_prior_diag(self, w):
        # Uniform prior
        if np.all((w >= 1e-6) & (w <= 1e6)):
            return 1
        else:
            return -np.inf

    # The log-posterior function
    def log_posterior_offdiag(self, w, omega, omega_fixed, data, nu=6, ell=3):
        return self.log_prior_offdiag(w) + self.log_likelihood_offdiag(w, omega, omega_fixed, data, nu, ell)

    # The log-posterior function
    def log_posterior_diag(self, w, omega, omega_fixed, data, ell=3):
        return self.log_prior_diag(w) + self.log_likelihood_diag(w, omega, omega_fixed, data, ell)

    # The log-posterior function
    def log_posterior_normal(self, w, omega, omega_fixed, data, nu=6, ell=3):
        return self.log_prior_offdiag(w) + self.log_likelihood_normal(w, omega, omega_fixed, data, nu, ell)

    def initialize_cepstral_distribution(self, ck_theory_var=None, psd_theory_mean=None):
        """
        Initialize the theoretical distribution of the cepstral coefficients.
        The samplelogpsd must has been already set.

        Input parameters:
            ck_theory_var   = the theoretical variance of cepstral coefficients, \\sigma*^2(P*,N)
            psd_theory_mean = the theoretical bias of log-PSD, \\lambda_l

        If ck_theory_var and/or psd_theory_mean are not specified, the default theoretical values will be used.
        """
        NF = self.samplelogpsd.size
        N = 2 * (NF - 1)

        if psd_theory_mean is None:
            # by default the THEORETICAL means are the one component ones:
            # ck THEORY mean:
            #    - EULER_GAMMA - log(2)   for k = {0, N/2}
            #    - EULER_GAMMA            otherwise
            self.logpsd_THEORY_mean = -EULER_GAMMA * np.ones(NF)
            self.logpsd_THEORY_mean[0] = -EULER_GAMMA - np.log(2)
            self.logpsd_THEORY_mean[-1] = -EULER_GAMMA - np.log(2)
        else:
            self.logpsd_THEORY_mean = psd_theory_mean

        # set theoretical errors
        if ck_theory_var is None:
            # by default the THEORETICAL variances are the one component ones:
            # ck THEORY variances:
            #    (pi^2)/3/N   for k = {0, N/2}
            #    (pi^2)/6/N   otherwise
            self.logpsdK_THEORY_var = (1.0 / N * np.concatenate(
                ([np.pi**2 / 3], [np.pi**2 / 6.0] * (NF - 2), [np.pi**2 / 3])))
            self.logpsdK_THEORY_std = np.sqrt(self.logpsdK_THEORY_var)
            # logtau THEORY variances:  (we assume to be summing ck up to K, included)
            #    (pi^2)/3/N*(2*K+1)   for K = {0, N/2-1}
            #    (pi^2)/3             for K = N/2
            self.logtau_THEORY_var = (1.0 / N * np.concatenate(
                (np.pi**2 / 3.0 * (2 * np.arange(NF - 1) + 1), [np.pi**2 / 3.0 * N])))
            self.logtau_THEORY_std = np.sqrt(self.logtau_THEORY_var)
        else:
            self.logpsdK_THEORY_var = ck_theory_var
            self.logpsdK_THEORY_std = np.sqrt(self.logpsdK_THEORY_var)
            self.logtau_THEORY_var = np.zeros(NF)
            self.logtau_THEORY_var[0] = self.logpsdK_THEORY_var[0]
            for K in range(1, NF - 1):
                self.logtau_THEORY_var[K] = (self.logtau_THEORY_var[K - 1] + 4.0 * self.logpsdK_THEORY_var[K])
            self.logtau_THEORY_var[-1] = (self.logtau_THEORY_var[-2] + self.logpsdK_THEORY_var[-1])
            self.logtau_THEORY_std = np.sqrt(self.logtau_THEORY_var)

    def scan_filter_tau(self, cutoffK=None, aic_Kmin_corrfactor=1.0, correct_mean=True):
        """
        Computes tau as a function of the cutoffK (= P*-1).
        Also computes psd and logpsd for the given cutoffK.
        If cutoffK is None, aic_Kmin * aic_Kmin_corrfactor will be used.

        Input parameters:
            cutoffK = (P*-1) = cutoff used to compute logtau and logpsd (by default = aic_Kmin * aic_Kmin_corrfactor)
            aic_Kmin_corrfactor = aic_Kmin cutoff correction factor (default: 1.0)
            correct_mean = fix the bias introduced by the log-distribution (default: True)

        self.tau_cutoffK will contain the value of tau for the specified cutoff cutoffK

        If cutoffK is out of range, the maximum K will be used.
        """
        if cutoffK is not None:
            if not isinstance(cutoffK, int) or (cutoffK < 0):
                raise ValueError('cutoffK must be a positive integer.')
            if aic_Kmin_corrfactor != 1.0:
                raise ValueError(
                    'If you specify cutoffK manually, the AIC will not be used, hence aic_Kmin_corrfactor will be ignored.'
                )
        self.aic_Kmin_corrfactor = aic_Kmin_corrfactor

        if cutoffK is None:
            self.cutoffK = int(round(self.aic_Kmin * self.aic_Kmin_corrfactor))
            self.manual_cutoffK_flag = False
        else:
            self.cutoffK = cutoffK
            self.manual_cutoffK_flag = True

        if self.cutoffK >= self.samplelogpsd.size:
            log.write_log('! Warning:  cutoffK ({:}) is out of range.'.format(self.cutoffK))
            # log.write_log('! Warning:  cutoffK ({:}) is out of range. The maximum frequency ({:}) will be used.'.format(self.cutoffK, self.samplelogpsd.size - 1))
            # self.cutoffK = self.samplelogpsd.size - 1

        # COS-filter analysis with frequency cutoff K
        self.logtau = dct_filter_tau(self.samplelogpsd)
        self.logpsd = dct_filter_psd(self.samplelogpsd, self.cutoffK)   # that is log(psd) for the chosen cutoffK
        self.psd = np.exp(self.logpsd)
        self.tau = np.exp(self.logtau)
        self.tau_THEORY_std = self.tau * self.logtau_THEORY_std

        if self.cutoffK < self.samplelogpsd.size:
            self.logtau_cutoffK = self.logtau[self.cutoffK]
            self.logtau_var_cutoffK = self.logtau_THEORY_var[self.cutoffK]
            self.logtau_std_cutoffK = self.logtau_THEORY_std[self.cutoffK]
            self.tau_cutoffK = self.tau[self.cutoffK]
            self.tau_std_cutoffK = self.tau_THEORY_std[self.cutoffK]
            self.tau_var_cutoffK = self.tau_std_cutoffK**2
        else:
            self.logtau_cutoffK = np.NaN
            self.logtau_var_cutoffK = np.NaN
            self.logtau_std_cutoffK = np.NaN
            self.tau_cutoffK = np.NaN
            self.tau_var_cutoffK = np.NaN
            self.tau_std_cutoffK = np.NaN

        if correct_mean:
            self.logpsd = self.logpsd + self.logpsd_THEORY_mean
            self.logtau = self.logtau + self.logpsd_THEORY_mean[0]
            self.logtau_cutoffK = self.logtau_cutoffK + self.logpsd_THEORY_mean[0]

    def scan_filter_psd(self, cutoffK_LIST, correct_mean=True):
        """Computes the psd and tau as a function of the cutoff K.
        Repeats the procedure for all the cutoffs in cutoffK_LIST."""
        self.cutoffK_LIST = cutoffK_LIST
        self.logpsd_K_LIST = np.zeros((self.samplelogpsd.size, len(self.cutoffK_LIST)))
        self.psd_K_LIST = np.zeros((self.samplelogpsd.size, len(self.cutoffK_LIST)))
        self.logtau_K_LIST = np.zeros(len(self.cutoffK_LIST))   # DEFINED AS log(PSD[0]), no factor 0.5 or 0.25
        self.tau_K_LIST = np.zeros(len(self.cutoffK_LIST))

        for k, K in enumerate(self.cutoffK_LIST):
            # COS-filter analysis with frequency cutoff K
            self.logpsd_K_LIST[:, k] = dct_filter_psd(self.samplelogpsd, K)
            self.logtau_K_LIST[k] = self.logpsd_K_LIST[0, k]
            self.psd_K_LIST[:, k] = np.exp(self.logpsd_K_LIST[:, k])
            self.tau_K_LIST[k] = np.exp(self.logtau_K_LIST[k])

            if correct_mean:
                self.logpsd_K_LIST[:, k] = (self.logpsd_K_LIST[:, k] + self.logpsd_THEORY_mean)
                self.logtau_K_LIST[k] = (self.logtau_K_LIST[k] + self.logpsd_THEORY_mean[0])
