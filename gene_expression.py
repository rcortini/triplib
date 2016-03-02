import numpy as np

# assign the promoters
def site_idx (n, hic_res=2000) :
    return int (n/hic_res)

def gene_expression_probability_matrix (nsites, genes, pfunc,
                                        hic_res=2000, promoter_distance=500,
                                       n_exclude=0) :
    Pg = np.zeros ([nsites,nsites])
    for gene in genes :
        startsite = gene ['start']
        endsite = gene ['end']
        strand = gene ['strand']
        if strand == '+' :
            startsite -= promoter_distance
        else :
            endsite = startsite
            startsite = endsite + promoter_distance
        startsite = site_idx (startsite, hic_res)
        endsite = site_idx (endsite, hic_res)
        try :
            Pg [startsite,endsite] = pfunc (gene)
        except IndexError :
            continue
    return Pg [n_exclude:,n_exclude:]

# the firing function default: fire if gene present
def always_fire (gene) :
    return 1.

# fire if gene is active
def fire_if_active (gene) :
    if gene ['expr'] != 0. :
        return 1.
    else :
        return 0.
