import numpy as np

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
