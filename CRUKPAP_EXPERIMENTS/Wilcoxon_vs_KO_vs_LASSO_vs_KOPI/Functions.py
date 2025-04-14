import numpy as np
import sanssouci as sa
from scipy.stats import hmean
from joblib import Parallel, delayed
import warnings
from sklearn.utils.validation import check_memory
from hidimstat.knockoffs.knockoff_aggregation import _empirical_pval


# Function to perform variable selection with KOPI
# all the functions are based on (or exact copies of) function written in https://github.com/alexblnn/KOPI

def get_null_pis(B, p):

    pi0 = np.zeros((B, p))
    for b in range(B):
        signs = (np.random.binomial(1, 0.5, size=p) * 2) - 1
        Z = 0
        for j in range(p):
            if signs[j] < 0:
                pi0[b][j] = 1
                Z += 1
            else:
                pi0[b][j] = (1 + Z) / p
    pi0 = np.sort(pi0, axis=1)
    return pi0

def get_pi_template(B, p):
    """
    Build a template by sampling from the joint distribution 
    of pi statistics under the null.

    Parameters
    ----------

    B: int
        Number of samples
    p: int
        Number of variables

    Returns
    -------
    template : array-like of shape (B, p)
        Sorted set of candidate threshold families
    
    """
    pi0 = get_null_pis(B, p)
    template = np.sort(pi0, axis=0) # extract quantile curves
    return template




def aggregate_list_of_matrices(pi0_raw):
    """
    Provided "draws" null pi-statistics matrices pi0,
    aggregate them to retain a single pi-statistics matrix.
    """
    pi0_raw = np.array(pi0_raw)
    draws, B, p = pi0_raw.shape
    pi0_raw_ = np.reshape(pi0_raw, (B, draws, p))
    pi0 = np.vstack([hmean(pi0_raw_[i], axis=0) for i in range(B)])

    return pi0




def create_hmean_p_values_learned_tpl(p, draws = 100, B = 100, n_jobs = 4):
    parallel = Parallel(n_jobs)
    pi0_raw = np.array(parallel(delayed(get_null_pis)(B, p) for draw in range(draws)))
    
    pi0_hmean = aggregate_list_of_matrices(pi0_raw)

    learned_tpl_raw = np.array(parallel(delayed(get_pi_template)(B, p) for draw in range(draws)))
    learned_tpl_hmean_ = aggregate_list_of_matrices(learned_tpl_raw)

    learned_tpl_hmean = np.sort(learned_tpl_hmean_, axis=0)
    
    return pi0_hmean, learned_tpl_hmean




def report_selected(p_values, region_size):
    """
    For p-values or e-values and a given region size, 
    compute FDP and power.
    """
    p = np.shape(p_values)[0]
    
    if region_size == 0:
        return None, None

    else:
        selected = np.argsort(p_values)[:region_size]
        
    prediction = np.array([0] * p)
    prediction[selected] = 1


    return selected, prediction






def KOPI_given_all (ko_stats, draws, alpha, target_fdp, learned_tpl_hmean, pi0_hmean, k_max):
    """
    Compute results for a single run.
    """
    warnings.filterwarnings("error")

    try:
        calibrated_thr_hmean = sa.calibrate_jer(alpha, learned_tpl_hmean, pi0_hmean, k_max)
    except UserWarning:
        print("No threshold family controls the JER at level alpha ( all family control the JER at level alpha for extreme alpha values), increase B. Aborting simulation", alpha)
        return None

    warnings.filterwarnings("default")
    
    pvals = np.array([_empirical_pval(ko_stats[i], 1)for i in range(draws)])

    p_values_hmean = hmean(pvals, axis=0)

    size_hmean = sa.find_largest_region(p_values_hmean, calibrated_thr_hmean, 1 - target_fdp)
    
    selected_hmean, prediction_hmean = report_selected(p_values_hmean, size_hmean)

    return selected_hmean, prediction_hmean





#def _empirical_pval(test_score, offset=1):

#    pvals = []
#   n_features = test_score.size

#    if offset not in (0, 1):
#        raise ValueError("'offset' must be either 0 or 1")

#    test_score_inv = -test_score
#    for i in range(n_features):
#        if test_score[i] <= 0:
#            pvals.append(1)
#        else:
#            pvals.append(
#               (offset + np.sum(test_score_inv >= test_score[i])) /
#                n_features
#            )

#    return np.array(pvals)


# Quel est l'effet du paramètre k_max ?
# Quel est l'effet du paramètre alpha ?
