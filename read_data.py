import numpy as np

# module-wide variables
datadir = "/home/rcortini/work/data/"
hic_file = datadir + "chr2L.mtx"
gene_file = datadir + "Kc_exp_color.bed"

# load hi-c data
def load_hic () :
    return np.loadtxt (hic_file)

# load the gene data
def load_genes () :
    gene_dtype = {'names' : ['chr','start','end','strand','expr','color','gene','state9'],
                             'formats' : ['S10','i8','i8','S2','f','S10','S10','i4']}
    return np.genfromtxt (gene_file, dtype=np.dtype (gene_dtype), skip_header=1)
