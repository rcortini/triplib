import numpy as np

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
    n = np.sum (M,axis=0)
    return (M/n).T
