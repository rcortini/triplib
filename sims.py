import numpy as np
from scipy.sparse.linalg import eigs
from scipy.optimize      import minimize

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
                                  ntrials=10,
                                  divide_by_sum=True) :
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
    if divide_by_sum :
        return population/np.sum(population)
    else :
        return population

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

def propagate_dirac_comb (startsites, P, nsteps=100, with_identity=False) :
    """
    Propagates a solution of elements starting at sites described
    by the array of indices startsites, using the probability
    matrix P, with nsteps
    """
    nsites = P.shape [0]
    # initialize matrix
    X = np.zeros ((nsteps,nsites))
    if with_identity :
        X[0,:] = startsites
    else :
        X[0,:] = np.dot (startsites,P)
    for i in range (1,nsteps) :
        # get state at previous step
        x = X[i-1,:]
        # feed it to the evolution matrix
        X[i,:] = np.dot (x,P)
    return X

def weigh_with_exponential (X, tau, with_identity=False) :
    """
    Given a matrix of time steps describing the convergence to the
    equilibrium distribution X, this function returns the average population
    of the sites, supposing that there is an exponential process with
    half-life tau
    """
    nsteps = X.shape [0]
    nsites = X.shape [1]
    if with_identity :
        steps = np.arange (0,nsteps,1)
    else :
        steps = np.arange (1,nsteps+1,1)
    p = np.exp (-steps/tau)
    predicted = np.zeros (nsites)
    for i, x in enumerate (X) :
        predicted += p[i] * x
    return predicted

def propagate_dirac_comb_k (k, startsites, P, nsteps=100) :
    """
    A parallel-ready version of propagate_dirac_comb
    """
    sites = np.append (startsites, k)
    return propagate_dirac_comb (sites, P, nsteps=nsteps)

def obj_diffusion (x,A,expr,mask) :
    """
    Objective function to minimize when optimizing the start vector for a
    diffusion on a Hi-C matrix. Return the gradient along with the function for
    an efficient calculation of both.
    """
    xA = np.dot (x,A)
    mxA = np.mean (xA)
    h = np.log2 (xA/mxA)
    h_minus_expr = h-expr 
    F = np.sum ((h_minus_expr**2)[mask])
    mAl = np.mean (A,axis=1)
    dh = ((A/xA).T - mAl/mxA).T
    df = np.sum ((h_minus_expr*dh)[:,mask],axis=1)
    # 2/ln(2) = 2.8853900817779268
    return F, 2.8853900817779268*df

# the constraints
cons = ({'type': 'eq',
                 'fun' : lambda x: np.sum (x)-1,
                 'jac' : lambda x: np.ones (len (x))},
        {'type': 'ineq',
                 'fun' : lambda x: 1-x,
                 'jac' : lambda x: np.identity (len(x))},
        {'type': 'ineq',
                 'fun' : lambda x: x,
                 'jac' : lambda x: np.identity (len(x))})

def optimize_start_diffusion (xstart,A,expr,mask,disp=True) :
    """
    Optimize the start vector of the diffusion of factors in a Hi-C matrix. The
    matrix A must be precomputed to be the exp-weighed powers of P, which was
    previously row-normalized. Provide the full expr_binned vector along with a
    pre-calculated mask describing the valid values of the expr_binned vector
    """
    res = minimize (obj_diffusion,
                    xstart,
                    args=(A,expr,mask),
                    jac=True,
                    constraints=cons,
                    method='SLSQP',
                    options={'disp': disp,'maxiter' : 100})
    return res

# the function to minimize
def obj_contacts (x,A,expr,mask) :
    """
    Objective function to minimize when optimizing the start vector for a
    contact vector with Hi-C matrix.
    """
    xA = np.dot (A,x)
    mxA = np.mean (xA)
    h = np.log2 (xA/mxA)
    h_minus_expr = h-expr 
    F = np.sum ((h_minus_expr**2)[mask])
    return F

def optimize_start_contacts (xstart,P,expr,mask,disp=True) :
    """
    Optimize the start vector of the contacts of sites in a Hi-C matrix. P was
    previously row-normalized. Provide the full expr_binned vector along with a
    pre-calculated mask describing the valid values of the expr_binned vector
    """
    res = minimize (obj_contacts,
                    xstart,
                    args=(P,expr,mask),
                    jac=False,
                    constraints=cons,
                    method='SLSQP',
                    options={'disp': disp,'maxiter' : 100})
    return res

def get_xstart (N,in_xstart_file=None,out_xstart_file=None) :
    """
    Get the xstart vector. If in_xstart_file is given, load the file and return.
    Otherwise, generate a random start vector for the optimization, and if requested save
    it to a file
    """
    if in_xstart_file is not None :
        return np.loadtxt (in_xstart_file)
    else :
        xstart = np.random.random (N)
        xstart /= np.sum (xstart)
        if out_xstart_file is not None :
            np.savetxt (out_xstart_file,xstart)
        return xstart

def optimize_model (what,matrix,expr,mask,target_fval=None,
                    in_xstart_file=None,
                    out_xstart_file=None,
                    disp=True) :
    """
    This function allows to get to the correct optimizing procedure, and iterate
    until the target f value is lower than the supplied one
    """
    # what shall we optimize?
    if what == 'contacts' :
        target_f = optimize_start_contacts
    elif what == 'diffusion' :
        target_f = optimize_start_diffusion
    # start vector
    N = matrix.shape[0]
    xstart = get_xstart (N,in_xstart_file,out_xstart_file)
    # optimize the first time
    res = target_f (xstart,matrix,expr,mask,disp=disp)
    if target_fval is not None :
        while (res.fun > target_fval) :
            xstart = get_xstart (N,in_xstart_file=None,out_xstart_file=out_xstart_file)
            res = target_f (xstart,matrix,expr,mask,disp=disp)
    return res.x
