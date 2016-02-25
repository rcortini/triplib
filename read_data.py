import numpy as np

# module-wide variables
datadir = "/home/rcortini/work/data/"
hic_file = datadir + "chr2L.mtx"
hic_file_normalized = datadir + "chr2L_normalized.mtx"
gene_file = datadir + "Kc_exp_color.bed"
colors_file = datadir + "colors.bed"
reporter_file = datadir + "allprom.txt"

# load hi-c data
def load_hic (normalized=True) :
    if normalized :
        return np.loadtxt (hic_file_normalized)
    else :
        return np.loadtxt (hic_file)


# load the gene data
def load_genes () :
    gene_dtype = {'names' : ['chr','start','end','strand','expr','color','gene','state9'],
                             'formats' : ['S10','i8','i8','S2','f','S10','S10','i4']}
    return np.genfromtxt (gene_file, dtype=np.dtype (gene_dtype), skip_header=1)

# load the colors data
def load_colors () :
    colors_dtype = {'names' : ['chr','start','end','color','size'],
                             'formats' : ['S10','i8','i8','S8','i8']}
    return np.genfromtxt (colors_file, dtype=np.dtype (colors_dtype))

# load the reporter data
def load_reporters () :
    with open (reporter_file, 'r') as f :
        for line in f :
            if line [0] == '#' :
                continue
            break
    keys = line.strip ('\n').split ()
    types = ['S30','S2','S2','i12','f','S4','i2']
    for i in range (len (types), len (keys)) :
        types.append ('f')
    reporter_dtype = {'names' : keys, 'formats' : types}
    return np.genfromtxt (reporter_file, dtype=np.dtype (reporter_dtype), skip_header=6)

# get hic-binned average reporter expression
def load_expr_binned (reporters, nsites,
                      responsive = ['p8','p38','p39','p40'],
                      hic_res=2000,
                      substitute_nans=False) :
    expr_binned = np.zeros (nsites)
    n_reporters = np.zeros (nsites)
    for r in reporters :
        prom = r['prom']
        if prom != 'p0' and prom in responsive : 
            pos = r['pos']
            i = pos/hic_res
            try :
                expr_binned [i] += r['nexp']
                n_reporters [i] += 1.
            except IndexError :
                continue
    # calculate the average per bin
    with np.errstate (invalid='ignore') :
        expr_binned/=n_reporters
    if substitute_nans :
        expr_binned [np.isnan (expr_binned)] = 0.
    return expr_binned
