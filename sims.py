import numpy as np
from scipy.sparse.linalg import eigs

def jump_to (p):
    """
    Select a site according to the probability vector p
    """
    value = np.random.random ()
    return np.argmax (value-p<0.)

def do_the_search (i0, nsteps, P) :
    """
    Perform a search of nsteps on the system described by the transition
    matrix P, with a starting point given by i0.
    """
    nsites = P.shape [0]
    visits = np.zeros (nsites).astype (int)
    # init
    i = i0
    # cycle on steps
    for t in range (nsteps) :
        visits [i] += 1
        i = jump_to (P [i])
    return visits

def get_equilibrium_distribution (P,
                                  from_eigs=False,
                                  nsteps=100000,
                                  ntrials=10) :
    """
    Calculate equilibrium distribution of a random walk on the row-normalized
    graph described by the matrix P. Return the vector that corresponds to the
    found solution
    """
    if from_eigs :
        hw, hv = eigs (P.T,k=1,which='LM')
        # select the index of the largest eigenvalue
        imax = hw.argmax ()
        # select the eigenvector corresponding to the largest eigenvalue
        population = hv [:,imax].real
    else :
        nsites = P.shape[0]
        i0 = np.random.randint (0,nsites,ntrials)
        Pn = np.cumsum (P,axis=1)
        visits = np.array ([do_the_search (i, nsteps, Pn)\
                            for i in range (ntrials)])
        # get the population from average visits
        population = np.mean (visits, axis=0)
    # final result
    mean = np.mean (population)
    return population/mean

def get_life_and_death (P, P_gene, tau, where='prom', ntrials=10, hic_res=2000) :
    """
    Get the population of the graph if including the following hypothesis: the
    particles may start only at the promoters, terminators, or both.
    """
    nsites = P.shape[0]
    # select starting sites depending on the "where" parameters, passed to the
    # function
    promoters, terminators = np.where (P_gene)
    if where=='prom' :
        i0 = np.random.choice (promoters,size=ntrials)
    elif where=='term' :
        i0 = np.random.choice (terminators,size=ntrials)
    elif where=='both' :
        i0 = np.random.choice (np.concatenate((promoters,terminators)),size=ntrials)
    # calculate the cumulative sum probability matrix
    Pn = np.cumsum (P,axis=1)
    visits = np.zeros (nsites)
    for i in range (ntrials) :
        t = int (np.random.exponential (tau))+1
        visits += do_the_search (i0[i],t,Pn)
    # get the population from average visits
    mean = np.mean (visits)
    return visits/mean

def model_r2 (population, expr_binned) :
    """
    Returns the R2 of the population, compared with the binned reporter
    expression
    """
    mask = [~np.isnan (expr_binned)]
    x = expr_binned [mask]
    y = population [mask]
    return np.corrcoef (x,y)[0,1]**2
