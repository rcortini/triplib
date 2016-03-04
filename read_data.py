import numpy as np
import os

# module-wide variables
datadir = os.getenv ("HOME") + "/work/data/"
hic_file_normalized = datadir + "chr2L_normalized.mtx"
hic_file_filled = datadir + "chr2L_filled.mtx"
gene_file = datadir + "Kc_exp_color.bed"
colors_file = datadir + "colors.bed"
reporter_file = datadir + "allprom.txt"

# load hi-c data
def load_hic (chromosome, normalized=False, filled=True) :
    hic_file = datadir + "chr" + chromosome
    if filled :
        return np.loadtxt (hic_file + "_filled.mtx")
    else :
        if normalized :
            return np.loadtxt (hic_file + "_normalized.mtx")
        else :
            return np.loadtxt (hic_file + ".mtx")


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
    types = ['S30','S2','S2','i8','f','S4','i2']
    for i in range (len (types), len (keys)) :
        types.append ('f')
    reporter_dtype = {'names' : keys, 'formats' : types}
    return np.genfromtxt (reporter_file, dtype=np.dtype (reporter_dtype), skip_header=6)

# get hic-binned average reporter expression
def load_expr_binned (reporters, nsites,
                      responsive = ['pI','pII','pIII','pIV'],
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
