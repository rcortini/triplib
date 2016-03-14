import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from .data_process import model_r2

def visualize_model_results (population, expr_binned, reporters,
                   n1=1000, n2=2000, hic_res=2000, n_exclude=0, title=None) :
    fig = plt.figure (figsize=(20,10))
    gs = gridspec.GridSpec (3,2)
    ax = plt.subplot(gs[:,0])
    ax.scatter (expr_binned [n_exclude:],np.log2(population))
    ax.axhline (0,-15.,15.,linestyle='--',color='r',linewidth=3)
    ax.axvline (0,-2.5,2.0,linestyle='--',color='r',linewidth=3)
    ax.set_xlabel ("Reporter expression (log)", fontsize=24)
    ax.set_ylabel ("Average visits (log)", fontsize=24)
    r2 = model_r2 (population, expr_binned [n_exclude:])
    xmin, xmax = ax.get_xaxis().get_view_interval()
    ymin, ymax = ax.get_yaxis().get_view_interval()
    deltax = xmax-xmin
    deltay = ymax-ymin
    ratio=0.25
    ax.text (xmax-ratio*deltax,ymax-ratio*deltay, "$r^2 = %.2f$"%r2,
            fontsize=24)
    # total expression
    ax2 = plt.subplot (gs[0,1])
    ax2.plot (np.log2(population))
    # zoom to region
    x = np.arange (n1,n2)
    ax3 = plt.subplot (gs[1,1])
    ax3.bar (x, np.log2(population[n1:n2]))
    ax3.get_yaxis().tick_left()
    ax3.set_ylabel ("Average visits (log)", fontsize=16)
    ax3.get_xaxis().tick_bottom()
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['bottom'].set_visible(False)
    ax3.set_xlim (n1,n2)
    # reporter expression
    ax4 = plt.subplot (gs[2,1],sharex=ax3)
    myreporters = np.array ([r for r in reporters
                             if (r['pos']-n_exclude)>n1*hic_res
                             and (r['pos']-n_exclude)<n2*hic_res])
    ax4.bar ((myreporters['pos']-n_exclude)/hic_res, myreporters['nexp'])
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.get_yaxis().tick_left()
    ax4.get_xaxis().tick_bottom()
    ax4.set_ylabel ("Reporter expression (log)", fontsize=16)
    # set title
    if title :
        fig.suptitle (title,fontsize=32)
    plt.show ()
    return fig
