# -*- coding: utf-8 -*-
import numpy as np
import scipy.special as sp

# from scipy.special import multigammaln
from scipy.optimize import minimize
from scipy.linalg import cholesky
import opt_einsum
from sportran.utils import log
from .tools.filter import runavefilter

EULER_GAMMA = (
    0.57721566490153286060651209008243104215933593992   # Euler-Mascheroni constant
)
LOG2 = np.log(2)


class MaxLikeFilter:
    """
    Maximum-likelihood estimate of the Onsager or transport coefficient.

    Parameters:
    - data: The noisy data (spectral matrix or one of its components).
    - model: Function that models the data (e.g., spline function).
    - n_parameters: Number of parameters for the fit or 'AIC' for automatic selection.
    - n_components: Number of independent samples the data is generated from.
    - n_currents: Number of independent flux types.
    - likelihood: Type of likelihood function to use (
                                                        'wishart',
                                                        'chisquare',
                                                        'variancegamma'
                                                      ).
    - solver: Optimization solver (e.g., 'BFGS').
    - omega_fixed: Fixed frequencies for the model nodes.
    """

    def __init__(self, data=None, model=None, n_parameters=None, n_components=None, n_currents=None, likelihood=None,
                 solver=None, omega_fixed=None, ext_guess=None, alpha=10**(np.linspace(-10, -3, 10000)),
                ):
        """
        Initialize the MaxLikeFilter class with the provided parameters.
        """
        log.write_log('MaxLikeFilter Initialization')

        self.data = data
        self.model = model
        self.alpha = alpha
        self.n_parameters = n_parameters
        self.n_components = n_components
        self.n_currents = n_currents
        self.solver = solver
        self.omega_fixed = omega_fixed
        self.omega = np.arange(data.shape[-1]) if data is not None else None
        self._data_prepared = False

        # Set likelihood function
        self.log_like = self._get_likelihood_function(likelihood)

        # Store optimization results
        self.parameters_mean = None
        self.parameters_std = None
        self.parameters_cov = None
        self.optimizer_res = None
        self.log_likelihood_value = None
        self.aic_values = None
        self.optimal_nparameters = None

        # DEBUG
        self.ext_guess = ext_guess

    def _get_likelihood_function(self, likelihood):
        """
        Get the likelihood function based on the provided likelihood type.
        """
        if likelihood is None:
            return None
        likelihood = likelihood.lower()
        likelihoods = {
            'wishart': self.log_likelihood_wishart,
            'chisquare': self.log_likelihood_diag,
            'chisquared': self.log_likelihood_diag,
            'variancegamma': self.log_likelihood_offdiag,
            'variance-gamma': self.log_likelihood_offdiag,
        }
        if likelihood in likelihoods:
            return likelihoods[likelihood]
        else:
            raise NotImplementedError('Supported likelihoods: wishart, chisquare, variance-gamma')

    def _validate_parameters(self):
        """
        Ensure that all necessary parameters are set before running maxlike.
        """
        assert (self.n_parameters is not None), 'Number of parameters (n_parameters) must be provided'
        assert self.solver is not None, 'Solver must be provided'
        assert self.data is not None, '`data` must be provided'
        assert self.log_like is not None, 'Likelihood must be set'
        assert self.model is not None, 'Model must be provided'

    def maxlike(self, data=None, model=None, n_parameters=None, likelihood=None, solver=None, mask=None,
                n_components=None, n_currents=None, guess_runave_window=50, omega_fixed=None, write_log=True,
                minimize_kwargs=None, limits=None,
               ):
        """
        Perform the maximum-likelihood estimation.
        """
        minimize_kwargs = minimize_kwargs or {}

        # FIXME
        self.limits = limits

        # Update instance variables if provided
        self._update_parameters(data, model, n_parameters, likelihood, solver, n_components, n_currents, omega_fixed,)

        # Validate necessary variables
        self._validate_parameters()

        if write_log:
            log.write_log(f'Maximum-likelihood estimation with n_parameters = {self.n_parameters}')

        # Prepare data
        self._prepare_data(mask)

        # Prepare n_parameters
        n_parameters_list = self._prepare_n_parameters(self.n_parameters)

        if isinstance(n_parameters_list, int):
            # Run with fixed number of parameters
            self._run_maxlike_fixed_parameters(n_parameters_list, guess_runave_window, minimize_kwargs, write_log)
        else:
            # Run over a range of n_parameters and select optimal one via AIC
            self._run_maxlike_aic(n_parameters_list, guess_runave_window, minimize_kwargs, write_log)

    def _update_parameters(self, data, model, n_parameters, likelihood, solver, n_components, n_currents, omega_fixed,):
        """
        Update class parameters with new values if provided.
        """
        if data is not None:
            self.data = data
            self.omega = np.arange(data.shape[-1])
        if model is not None:
            self.model = model
        if n_parameters is not None:
            self.n_parameters = n_parameters
        if likelihood is not None:
            self.log_like = self._get_likelihood_function(likelihood)
        if solver is not None:
            self.solver = solver
        if n_components is not None:
            self.n_components = n_components
        if n_currents is not None:
            self.n_currents = n_currents
        if omega_fixed is not None:
            self.omega_fixed = omega_fixed

    def _prepare_data(self, mask):
        """
        Prepare data for processing, applying mask if provided.
        """
        if mask is not None:
            self.data = self.data[mask]

        # Validate data shape
        self._validate_data_shape()

        self.omega = np.arange(self.data.shape[-1])

        if not self._data_prepared:
            # Perform the axis transformation only once
            self._orig_data = np.copy(self.data)
            self.data = np.moveaxis(self.data, -1, 0)
            self._data_prepared = True

    def _validate_data_shape(self):
        """
        Validate the shape of the input data.
        """
        if len(self.data.shape) == 3:
            assert (self.log_like == self.log_likelihood_wishart), 'Misshaped `data` for likelihood'
            assert (self.data.shape[0] == self.data.shape[1]), 'Data for a Wishart estimate must be a (n,n,N) array.'
        elif len(self.data.shape) != 1:
            raise ValueError('`data` should be a 1D or 3D array')

    def _prepare_n_parameters(self, n_parameters):
        """
        Prepare the number of parameters array based on input.
        """
        if isinstance(n_parameters, int):
            return n_parameters
        elif isinstance(n_parameters, (list, np.ndarray)):
            n_parameters = np.asarray(n_parameters)
            assert np.issubdtype(n_parameters.dtype, np.integer), '`n_parameters` must be an integer array-like'
            log.write_log((f'Optimal number of parameters between {n_parameters.min()} '
                           f'and {n_parameters.max()} chosen by AIC'))
            return n_parameters
        elif isinstance(n_parameters, str) and n_parameters.lower() == 'aic':
            n_parameters = np.arange(3, 40)
            log.write_log('Optimal number of parameters between 3 and 40 chosen by AIC')
            return n_parameters
        else:
            raise ValueError('Invalid value for n_parameters')

    def _run_maxlike_fixed_parameters(self, n_parameters, guess_runave_window, minimize_kwargs, write_log):
        """
        Run maximum likelihood estimation with a fixed number of parameters.
        """
        self.n_parameters = n_parameters
        self._initialize_spline_nodes(write_log)
        if self.ext_guess is not None:
            guess_data = self.ext_guess
        else:
            guess_data = self.guess_data(self._orig_data, self.omega, self.omega_fixed, self.n_components,
                                         self.n_currents, window=guess_runave_window,
                                        )

        # Perform optimization
        self._optimize_parameters(guess_data, minimize_kwargs, write_log)

    def _run_maxlike_aic(self, n_parameters_list, guess_runave_window, minimize_kwargs, write_log):
        """
        Run maximum likelihood estimation over a range of parameters and choose the best
        one with AIC.
        """
        _aic = []
        _aic_max = -np.inf
        _steps_since_last_aic_update = 0

        for n_par in n_parameters_list:
            log.write_log(f'n_parameters = {n_par}')
            self.n_parameters = int(n_par)
            self.omega_fixed = None   # Reset omega_fixed to recompute spline nodes
            self._initialize_spline_nodes(write_log)
            guess_data = self.guess_data(self._orig_data, self.omega, self.omega_fixed, self.n_components,
                                         self.n_currents, window=guess_runave_window,
                                        )

            self._optimize_parameters(guess_data, minimize_kwargs, write_log=False)
            _new_aic = self.log_likelihood_value - n_par
            _aic.append(_new_aic)

            if _new_aic > _aic_max:
                _aic_max = _new_aic
                self.optimal_nparameters = n_par
                self.best_parameters_mean = self.parameters_mean.copy()
                self.best_parameters_std = self.parameters_std.copy()
                self.best_parameters_cov = self.parameters_cov.copy()
                self.best_omega_fixed = self.omega_fixed.copy()
                self.best_log_likelihood_value = self.log_likelihood_value
                _steps_since_last_aic_update = 0
            else:
                _steps_since_last_aic_update += 1

            if _steps_since_last_aic_update > 5:
                break

            if write_log:
                log.write_log((f'AIC: {_new_aic}; Steps since last AIC '
                               f'update: {_steps_since_last_aic_update}'))

        self.aic_values = np.array(_aic)
        # After the loop, set the parameters to the best found
        self.n_parameters = self.optimal_nparameters
        self.parameters_mean = self.best_parameters_mean
        self.parameters_std = self.best_parameters_std
        self.parameters_cov = self.best_parameters_cov
        self.omega_fixed = self.best_omega_fixed
        self.log_likelihood_value = self.best_log_likelihood_value

    def _initialize_spline_nodes(self, write_log):
        """
        Initialize spline nodes for the model.
        """
        if self.omega_fixed is None:
            if write_log:
                log.write_log('Spline nodes are equispaced from 0 to the Nyquist frequency.')
            args = np.int32(np.linspace(0, self.data.shape[0] - 1, self.n_parameters, endpoint=True))
            self.omega_fixed = self.omega[args]
        assert self.omega_fixed.shape[0] == self.n_parameters

    def guess_data(self, data, omega, omega_fixed, ell, nu, window=10):
        """
        Moving average of the input data as initial guess for the parameter estimation.
        """
        try:
            guess_data = self._compute_moving_average(data, window)
        except Exception as e:
            log.write_log(f'Guessing data failed with exception: {e}. Window changed to 10.')
            guess_data = self._compute_moving_average(data, 10)

        return self._process_guess_data(guess_data, omega, omega_fixed, ell, nu)

    def _compute_moving_average(self, data, window):
        """
        Compute the moving average of the data.
        """
        shape = data.shape
        return np.array([runavefilter(c, window) for c in data.reshape(-1, shape[-1])]).reshape(shape)

    def _process_guess_data(self, guess_data, omega, omega_fixed, ell, nu):
        """
        Process guess data based on the likelihood.
        """
        indices = [np.argmin(np.abs(omega - omega_fixed[i])) for i in range(len(omega_fixed))]
        guess_data = np.array([guess_data[..., idx] for idx in indices])

        if self.log_like == self.log_likelihood_wishart:
            guess_data = self._process_wishart_guess_data(guess_data, nu)
        elif self.log_like == self.log_likelihood_diag:
            guess_data = np.log(guess_data)

        return guess_data.flatten()

    def _process_wishart_guess_data(self, guess_data, nu):
        """
        Process the guess data for Wishart likelihood.
        """
        guess_data = np.array([cholesky(g, lower=False) for g in guess_data])
        upper_triangle_indices = np.triu_indices(nu)
        nw = self.omega_fixed.shape[0]
        return self._flatten_wishart_parameters(guess_data, upper_triangle_indices, nw, nu)

    def _flatten_wishart_parameters(self, guess_data, upper_triangle_indices, nw, nu):
        """
        Flatten the Wishart parameters for optimization.
        """
        guess_params = np.zeros((nw, nu * (nu + 1) // 2))
        for idx, (i, j) in enumerate(zip(*upper_triangle_indices)):
            guess_params[:, idx] = guess_data[:, i, j]
        return guess_params

    def _optimize_parameters(self, guess_data, minimize_kwargs, write_log):
        """
        Perform the optimization to find the parameters that maximize the likelihood.
        """

        self._used_guess_data = guess_data, self.data
        res = minimize(fun=self.log_like, x0=guess_data,
                       args=(self.model, self.omega, self.omega_fixed, self.data, self.n_currents, self.n_components,
                            ), bounds=self.limits, method=self.solver, **(minimize_kwargs or {}),
                      )
        self._store_optimization_results(res, write_log)

    def _optimize_alpha(self, res):
        """
        Assume a gaussian prior (alpha/pi)**(P/2) e**(-alpha*||w||**2) and maximize
        the marginal distribution of alpha. We sample the posterior distribution
        assuming it is Gaussian
        (see https://en.wikipedia.org/wiki/Bernstein–von_Mises_theorem)
        and compute p(D|alpha) reweighting the posterior at alpha=0: see
        reweight_alpha and reweight_logev_alpha_vec.
        """

        w = res.x

        cov = res.hess_inv
        if self.solver == 'L-BFGS-B':
            cov = cov.todense()

        samples = generate_samples_mc_alpha(w, cov)
        dic_alpha, self.alpha_plot = reweight_logev_alpha_vec(alpha=self.alpha, samples=samples)
        parameters_mean, parameters_cov = reweight_alpha(alpha=dic_alpha['alpha_s'], samples=samples)

        return dic_alpha, parameters_mean, parameters_cov

    def _store_optimization_results(self, res, write_log):
        """
        Store the results of the optimization.
        """

        if hasattr(res, 'hess_inv'):
            if write_log:
                log.write_log(f'The {self.solver} solver provides Hessian. '
                              'Covariance matrix estimated through Laplace approximation.')

            self.best_alpha, self.parameters_mean, self.parameters_cov = (self._optimize_alpha(res=res))

            try:
                cov = res.hess_inv.todense()
            except AttributeError:
                cov = res.hess_inv
            self.parameters_cov = cov

            self.parameters_mean = res.x
            self.parameters_std = np.sqrt(np.abs(self.parameters_cov.diagonal()))
        else:
            if write_log:
                log.write_log((f'The {self.solver} solver does not provide Hessian. '
                               'No covariance matrix output.'))
            self.parameters_mean = res.x
            self.parameters_std = None
            self.parameters_cov = None

        self.optimizer_res = res
        self.log_likelihood_value = -self.log_like(res.x, self.model, self.omega, self.omega_fixed, self.data,
                                                   self.n_currents, self.n_components,
                                                  )

    def log_likelihood_wishart(self, w, model, omega, omega_fixed, data_, nu, ell, eps=1e-3):
        """
        Logarithm of the Wishart probability density function.
        """
        n = ell
        p = nu

        # Compute scale matrix from the model
        V = scale_matrix(model, w, omega, omega_fixed, p)
        X = data_

        if n < p:
            S = np.linalg.svd(X, hermitian=True, compute_uv=False)
            detX = np.array([np.prod(s[abs(s) > eps]) for s in S])
        else:
            detX = np.linalg.det(X)

        invV = np.linalg.inv(V)
        detV = np.linalg.det(V)

        trinvV_X = opt_einsum.contract('wab,wba->w', invV, X)

        log_pdf = 0.5 * ((n - p - 1) * np.log(detX) - trinvV_X - n * np.log(detV))

        return -np.sum(log_pdf)

    def log_likelihood_diag(self, w, model, omega, omega_fixed, data, M, ell):
        """
        Negative of the logarithm of the Chi-squared probability density function.

        :param array like: data is the PSD
        :param int: nu is the number of thermodynamically independent fluxes (usally
                    useful for thermal transport)
        """

        # Number of degrees of freedom
        dof = int(ell - M + 1)

        # Model for the spectrum
        spline = model(omega_fixed, w)
        spectrum = np.exp(spline(omega))

        # Log-likelihood of data = spectrum * chi^2(dof) / dof
        log_pdf = 0.5 * (np.log(data**(dof - 2) / spectrum**dof) - dof * data / spectrum)

        # Return the negative log-likelihood
        return -np.sum(log_pdf)

    def log_likelihood_offdiag(self, w, model, omega, omega_fixed, data_, nu, ell):
        """
        Negative of the logarithm of the Variance-Gamma probability density function.
        """
        spline = model(omega_fixed, w)
        rho = np.clip(spline(omega), -0.98, 0.98)
        _alpha = 1 / (1 - rho**2)
        _beta = rho / (1 - rho**2)
        _lambda = 0.5 * ell * nu
        _gamma2 = _alpha**2 - _beta**2
        _lambda_minus_half = _lambda - 0.5

        z = data_ * nu * ell
        absz = np.abs(z)
        term1 = _lambda * np.log(_gamma2)
        term2 = _lambda_minus_half * np.log(absz)
        term3 = np.log(sp.kv(_lambda_minus_half, _alpha * absz))
        term4 = _beta * z
        term5 = -_lambda_minus_half * np.log(2 * _alpha)
        log_pdf = term1 + term2 + term3 + term4 + term5
        return -np.sum(log_pdf)

    def extract_and_scale_results(self):
        """
        Extract results and scale matrices according to the likelihood.
        """
        omega_fixed = self.omega_fixed
        params = self.parameters_mean
        params_cov = self.parameters_cov
        omega = self.omega

        if self.log_like == self.log_likelihood_wishart:
            self.NLL_mean = (scale_matrix(self.model, params, omega, omega_fixed, self.n_currents) / self.n_currents)

            self.NLL_std = (scale_matrix_std_mc(self.model, params, omega, omega_fixed, self.n_currents,
                                                self.parameters_cov, size=1000,
                                               ) / self.n_currents)
        else:
            _NLL_spline = self.model(omega_fixed, params)
            self.NLL_mean = np.exp(_NLL_spline(omega))

            try:
                err = np.random.multivariate_normal(np.zeros(params_cov.shape[0]), params_cov, size=1000)
                samples = [self.model(omega_fixed, params + e)(omega) for e in err]
                self.NLL_std = np.std([np.exp(s) for s in samples], axis=0)
            except TypeError:
                pass
            except AttributeError:
                pass

    def __repr__(self):
        return 'MaxLikeFilter:\n'


def scale_matrix(model, w, omega, omega_fixed, n):
    """
    Compute the scale matrix from the model.
    """
    elements = model(omega_fixed, w)(omega)
    is_complex = elements.dtype == np.complex128

    # Precompute the upper triangle indices
    triu_indices = np.triu_indices(n)

    # Preallocate L
    L = np.zeros((n, n, omega.shape[0]), dtype=elements.dtype)

    # Assign elements to L using vectorized operations
    for idx, (i, j) in enumerate(zip(*triu_indices)):
        L[i, j] = elements[:, idx]

    # Compute the scale matrix S using einsum
    S = np.einsum('jiw,jkw->wik', L.conj() if is_complex else L, L)

    return S


def scale_matrix_std_mc(model, w, omega, omega_fixed, n, cov_w, size=1000):
    """
    Compute the standard deviation of the scale matrix via Monte Carlo sampling.
    """
    sample = w + np.random.multivariate_normal(mean=np.zeros_like(w), cov=cov_w, size=size)
    sample_S = np.stack([scale_matrix(model, ww, omega, omega_fixed, n) for ww in sample])
    S_std = sample_S.std(axis=0)
    return S_std


def reweight_logev_alpha_vec(samples, alpha):
    """
    samples: shape is (N, P): N number of samples, P number of parameters
    array: array of alpha to test
    """
    M = samples.shape[1]
    means = np.mean(np.exp(-alpha[:, None] * np.linalg.norm(samples, axis=1)**2), axis=1)
    l = np.where(means > 1e-300)[0]
    truth_mean = np.log(means[l]) + M / 2 * np.log(alpha[l] * 2 / np.pi)
    dic_alpha = {}
    dic_alpha['lev_s'] = truth_mean
    dic_alpha['alpha_s'] = alpha[np.argmax(dic_alpha['lev_s'])]

    return dic_alpha, alpha[l]


def reweight_alpha(alpha, samples):
    """
    samples: shape is (N, P): N number of samples, P number of parameters
    alpha: scalar
    """
    # Compute the squared norms
    norm_samples = np.linalg.norm(samples, axis=1)**2

    # Use log-sum-exp trick to prevent underflow in the exponential terms
    max_exp_term = np.max(-alpha * norm_samples)

    exp_term = np.exp(-alpha * norm_samples - max_exp_term)

    # Weight denominator
    weight_denominator = np.mean(exp_term, axis=0)

    # Compute the weighted mean
    truth_mean = (np.mean(samples.T[:, :] * exp_term, axis=1,)) / weight_denominator

    # Compute the weighted covariance with log-sum-exp normalization
    weighted_samples = samples.T[:, None, :] * samples.T[None, :, :]

    # Adjust the covariance by scaling with the exponential
    truth_cov = (np.mean(weighted_samples * exp_term, axis=-1,)) / weight_denominator

    # Subtract the outer product of the means
    truth_cov = truth_cov - np.outer(truth_mean, truth_mean)

    return truth_mean, truth_cov


def generate_samples_mc_alpha(w, cov_w, size=1000):
    """
    samples shape is (N, P): N number of samples, P number of parameters
    w: parameters mean as estimated by self.maxlike
    cov_w: array PxP
    """

    sample = w + np.random.multivariate_normal(mean=np.zeros_like(w), cov=cov_w, size=size)

    return sample
