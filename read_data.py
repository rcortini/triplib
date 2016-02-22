import numpy as np

# module-wide variables
datadir = "/home/rcortini/work/data/"
hic_file = datadir + "chr2L.mtx"

# load hi-c data
def load_hic () :
    return np.loadtxt (hic_file)
