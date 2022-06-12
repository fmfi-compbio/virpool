import logging
import math
from statistics import stdev
import argh
import numba
import numpy as np
import yaml
from numba.experimental import jitclass
from numpy import average
from scipy.optimize import Bounds, minimize
from scipy.special import softmax
from scipy.stats import chi2

from helpers import load_posteriors


# silence NumbaPerformanceWarning
import warnings
from numba.core.errors import NumbaPerformanceWarning
warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)


@jitclass([
    ('posteriors', numba.float64[:, :]),
])
class Model:
    def __init__(self, posteriors):
        self.posteriors = posteriors

    def loglikelihood(self, weights):
        read_count = self.posteriors.shape[0]
        result = 0
        for r in np.arange(read_count):
            partial = np.log(max(
                np.float64(1e-300),
                np.float64(np.dot(weights, self.posteriors[r, :]))
            )
            )
            result += partial
        return result / read_count

    @staticmethod
    def softmax(x):
        expd = np.exp(x)
        result = expd / np.sum(expd)
        return result

    def target(self, x):
        weights = self.softmax(x)
        return -self.loglikelihood(weights)

    def jacobian(self, x):
        weights = self.softmax(x)
        jac_weights = self.jacobian_weights(weights)
        jac_softmax = self.jacobian_softmax(x)
        result = np.zeros(x.shape, dtype=np.float64)
        for i in np.arange(x.shape[0]):
            result[i] = np.dot(jac_weights, jac_softmax[:, i])
        return -result

    def jacobian_weights(self, weights):
        result = np.zeros(weights.shape, dtype=np.float64)
        read_count = self.posteriors.shape[0]
        for r in np.arange(read_count):
            common = 1 / max(np.float64(1e-300), np.float64(np.dot(weights, self.posteriors[r, :])))
            subresult = common * self.posteriors[r, :]
            result += subresult
        result /= read_count
        return result

    @staticmethod
    def jacobian_softmax(x):
        # result[i][j] = d w_i / d x_j
        k = x.shape[0]
        result = np.zeros((k, k), dtype=np.float64)
        g = max(np.sum(np.exp(x)), np.float64(1e-150))
        for i in np.arange(k):
            for j in np.arange(k):
                p1 = np.exp(x[i]) * g if i == j else 0
                p2 = np.exp(x[i]) * np.exp(x[j])
                result[i][j] = (p1 - p2) / (g ** 2)
        return result


def estimate_weights(posteriors):
    variant_count = posteriors.shape[1]
    rng = np.random.default_rng()
    x0 = [math.log(rng.uniform(1, 20)) for _ in np.arange(variant_count)]
    model = Model(posteriors)

    bounds = Bounds(
        lb=[0.0 for _ in np.arange(variant_count)],
        ub=[50.0 for _ in np.arange(variant_count)]
    )

    values_during_run = []

    def logger(x):
        nonlocal values_during_run
        weights = softmax(x)
        row = {"weights": list(map(float, weights))}
        values_during_run.append(row)
        # print(weights)

    minimizer_args = {
        "fun": model.target,
        "x0": np.array(x0),
        "jac": model.jacobian,
        "bounds": bounds,
        "options": {
            "disp": False,
            "gtol": 1e-08,
            "ftol": 1e-13,
            "maxcor": 50,
        },
        "callback": logger
    }

    optimisation_result = minimize(**minimizer_args)

    result_weights = softmax(optimisation_result['x'])
    result = {
        "estimated_weights": list(map(float, result_weights)),
        "fun": float(optimisation_result['fun']),
        "raw_output": str(optimisation_result),
        "values_during_run": values_during_run
    }

    return result


def estimate_weights_multiple_tries(posteriors, tries):
    subresults = []
    for t in range(tries):
        subresult = estimate_weights(posteriors)
        subresults.append(subresult)
    return subresults


def find_best_weights(results):
    best_weights, best_fun = None, math.inf
    for result in results:
        if result['fun'] < best_fun:
            best_weights, best_fun = result['estimated_weights'], result['fun']
    return best_weights, best_fun


def estimator_characteristics(results):
    variant_count = len(results[0]['estimated_weights'])
    means = [average([r['estimated_weights'][i] for r in results]) for i in range(variant_count)]
    std = [stdev([r['estimated_weights'][i] for r in results]) for i in range(variant_count)]
    return means, std


def pvalues_simple(posteriors, weights):
    """Pvalues (Bonferroni corrected) of a particular weight being equal to zero,
     but the P-value is just an upper bound"""
    read_count = posteriors.shape[0]
    variant_count = posteriors.shape[1]
    model = Model(posteriors)
    total_ll = model.loglikelihood(weights)
    pvalues = np.zeros(weights.shape)
    for v in range(variant_count):
        subweights = np.copy(weights)
        subweights[v] = 0
        subweights /= sum(subweights)
        ll = model.loglikelihood(subweights)
        s = -2 * read_count * (ll - total_ll)
        pvalue = chi2.sf(s, df=1)
        pvalue_adj = min(1.0, pvalue * variant_count)
        pvalues[v] = float(pvalue_adj)
    return pvalues


def pvalues_optimisation(posteriors, full_weights):
    """Pvalues (Bonferroni corrected) of a particular weight being equal to zero."""
    read_count = posteriors.shape[0]
    variant_count = posteriors.shape[1]
    model = Model(posteriors)
    total_ll = model.loglikelihood(full_weights)
    pvalues = np.zeros(full_weights.shape)
    for v in range(variant_count):
        subposteriors = np.delete(posteriors, v, 1)
        subresult = estimate_weights_multiple_tries(subposteriors, 10)
        subweights, _ = find_best_weights(subresult)
        ll = Model(subposteriors).loglikelihood(np.array(subweights))
        s = -2 * read_count * (ll - total_ll)
        pvalue = chi2.sf(s, df=1)
        pvalue_adj = min(1.0, pvalue * variant_count)
        pvalues[v] = float(pvalue_adj)
    return pvalues


@argh.arg("-n", "--tries", type=int)
@argh.arg("--pvalue-threshold", type=float)
@argh.arg("-l", "--logfile")
def main(posteriors_filename,
         output_filename,
         tries: int = 10,
         pvalue_threshold=None,
         tries_pvalues: int = 1,
         logfile=None):
    logging_level = logging.INFO
    if logfile is None:
        logging.basicConfig(level=logging_level)
    else:
        logging.basicConfig(filename=logfile, level=logging_level)

    with open(posteriors_filename, newline='') as f:
        variant_names, posteriors = load_posteriors(f)

    results = estimate_weights_multiple_tries(posteriors, tries)
    for result in results:
        logging.info(result['raw_output'])

    best_weights, best_fun = find_best_weights(results)
    # weights_mean, weights_std = estimator_characteristics(results)
    # print(weights_mean, weights_std)
    logging.info(f"best weights (full model): {list(best_weights)}")

    if pvalue_threshold is not None:
        pvalues = pvalues_optimisation(posteriors, np.array(best_weights))
        logging.info(f"P-values: {list(pvalues)}")
        pvalue_threshold = max(pvalue_threshold, np.amin(pvalues))
        columns_to_remove = np.array([i for i in range(pvalues.shape[0])
                                      if pvalues[i] > pvalue_threshold], dtype=int)
        columns_to_keep = [i for i in range(pvalues.shape[0])
                           if pvalues[i] <= pvalue_threshold]

        subposteriors = np.delete(posteriors, columns_to_remove, axis=1)
        subresults = estimate_weights_multiple_tries(subposteriors, tries_pvalues)
        best_subweights, best_subfun = find_best_weights(subresults)

        best_weights = np.zeros(pvalues.shape)
        for i in range(len(columns_to_keep)):
            best_weights[columns_to_keep[i]] = best_subweights[i]

    result = {variant_names[i]: float(best_weights[i]) for i in range(len(variant_names))}
    with open(output_filename, "w") as f:
        yaml.dump(result, f)


if __name__ == "__main__":
    argh.dispatch_command(main)
