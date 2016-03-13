import numpy as np
from .data_process import warn_message

# assign the promoters
def site_idx (n, hic_res=2000) :
    return int (n/hic_res)

def gene_expression_probability_matrix (nsites, genes, Pfunc,
                                        hic_res=2000, promoter_distance=500,
                                       n_exclude_start=0,
                                       n_exclude_end=0) :
    Pg = np.zeros ([nsites,nsites])
    pfunc = Pfunc[0]
    pfunc_args = Pfunc[1]
    for gene in genes :
        startsite = gene ['start']
        endsite = gene ['end']
        strand = gene ['strand']
        if strand == '+' :
            i = startsite - promoter_distance
            j = endsite
        else :
            i = endsite + promoter_distance
            j = startsite
        i = site_idx (i, hic_res)
        j = site_idx (j, hic_res)
        try :
            Pg [i,j] += pfunc (gene,pfunc_args)
        except IndexError :
            warn_message ("gene_exp", "Index error with %d/%d" % (startsite,endsite))
            continue
    return Pg [n_exclude_start:n_exclude_end,n_exclude_start:n_exclude_end]

# the firing function default: fire if gene present
def always_fire (gene,args) :
    return 1.

# fire if gene is active
def fire_if_active (gene,args) :
    if gene ['expr'] != 0. :
        return 1.
    else :
        return 0.

# fire proportionally to gene expression level
def linear_fire (gene,args) :
    expr = gene ['expr']
    min_expression = args [0]
    max_expression = args [1]
    return (expr-min_expression)/(max_expression-min_expression)
