from __future__ import print_function
import numpy as np
import time, sys

def time_string () :
    return time.strftime("[%Y-%m-%d %H:%M:%S]", time.localtime ())

def error_message (program_name, message) :
    full_message = "%s %s: ERROR: %s"%(time_string (), program_name, message)
    print (full_message, file=sys.stderr)

def log_message (program_name, message) :
    full_message = "%s %s: INFO: %s"%(time_string (), program_name, message)
    print (full_message)

def warn_message (program_name, message) :
    full_message = "%s %s: WARNING: %s"%(time_string (), program_name, message)
    print (full_message)

def model_r2 (population, expr_binned, bias=0.01, zerorows=None) :
    """
    Returns the R2 of the population, compared with the binned reporter
    expression
    """
    notnan = np.isnan (expr_binned)
    if zerorows is not None :
        mask = np.logical_and (~notnan, ~zerorows)
    else :
        mask = ~notnan
    x = expr_binned [mask]
    y = np.log2 (bias + population [mask])
    return np.corrcoef (x,y)[0,1]**2

def ps (s,alpha) :
    """
    theoretical p(s) with exponent alpha given
    """
    if s==0. :
        return 0.
    else :
        return np.power (s,-alpha)

def fill_p_with_ps (P,i,alpha) :
    """
    Fills the i-th row and column of P with values given
    by the p(s) function
    """
    dim = P.shape[0]
    for k in range (dim) :
        s = abs(i-k)
        P[i,k] = ps(s,alpha)
        P[k,i] = ps(s,alpha)

def row_normalize_matrix (M) :
    n = np.sum (M,axis=1)
    N = M.shape [0]
    Mnorm = M.copy ()
    for i in range (N) :
        if n[i]!=0. :
            Mnorm[i] /= n[i]
    return Mnorm

def P_hop_plus_gene (P_hop, P_gene, p_firing) :
    N = P_hop.shape[0]
    P = np.zeros ([N,N])
    for i in range (N) :
        n = np.sum (P_gene [i,:])
        if n!=0. :
            P [i,:] = p_firing*P_gene[i,:]/n + (1.-p_firing)*P_hop[i,:]
        else :
            P [i,:] = P_hop [i,:]
    return P

def get_promoters_and_terminators (genes, N, hic_res=2000) :
    """
    Given the gene array "genes", extract positions of promoters
    and terminators, and both, into an array of N bins.
    """
    promoters = np.zeros (N)
    terminators = np.zeros (N)
    both = np.zeros (N)
    for gene in genes :
        startsite = gene ['start']
        endsite = gene ['end']
        strand = gene ['strand']
        # IMPORTANT NOTE: here we suppose that we are using Python 2! The division
        # between integers in Python 3 yields a float.
        if strand == '+' :
            i = startsite/hic_res
            j = endsite/hic_res
        elif strand == '-' :
            i = endsite/hic_res
            j = startsite/hic_res
        else :
            raise ValueError
        promoters [i] += 1
        terminators [j] += 1
        both [i] += 1
        both [j] += 1
    return promoters, terminators, both
