import numpy as np
from .read_data import *
from .data_process import log_message

class Chromosome :
    def __init__(self,name,full_init=True,normalized=True,colors=None,reporters=None,genes=None) :
        self.name = name
        self.colors = None
        self.reporters = None
        self.genes = None
        self.H = None
        # load the number of bins in the Hi-C matrix and zerorows
        self.N = chromosome_N (name)
        self.zerorows = chromosome_zerorows (name)
        # if user wishes, we init all
        self.init_all (colors,reporters,genes,full_init=full_init,normalized=normalized)
    def init_colors (self,colors) :
        self.colors = np.array ([c for c in colors if c['chr']==self.name])
    def init_reporters (self,reporters) :
        self.reporters = np.array ([r for r in reporters if r['chr']==self.name])
        self.expr_binned = load_expr_binned (self.reporters,self.N)
    def init_genes (self,genes) :
        self.genes = np.array ([g for g in genes if g['chr']==self.name])
    def init_hic (self,normalized=True) :
        self.H = load_hic (self.name,normalized=normalized)
    def init_all (self,colors,reporters,genes,full_init=False,normalized=True) :
        self.init_colors (colors)
        self.init_reporters (reporters)
        self.init_genes (genes)
        if full_init :
            self.init_hic (normalized=normalized)

def load_all_chromosomes (full_init=True,normalized=True) :
    names = ['2L','2R','3L','3R','X']
    log_message ("load_all_chromosomes", "Loading reporters")
    reporters = load_reporters ()
    genes = load_genes ()
    colors = load_colors ()
    chromosomes = []
    for name in names :
        log_message ("load_all_chromosomes", "Loading chromosome %s"%name)
        chromosomes.append (Chromosome (name,
                                        full_init=full_init,
                                        normalized=normalized,
                                        colors=colors,
                                        reporters=reporters,
                                        genes=genes))
    del reporters
    del genes
    del colors
    return chromosomes
