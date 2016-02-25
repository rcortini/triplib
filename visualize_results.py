import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from .sims import model_r2

def visualize_model_results (population, expr_binned) :
    fig = plt.figure (figsize=(8,10))
    gs = gridspec.GridSpec (2,1,height_ratios=[10,2],hspace=0.2)
    ax = plt.subplot(gs[0,0])
    ax.scatter (expr_binned,np.log2(population))
    ax.axhline (0,-15.,15.,linestyle='--',color='r',linewidth=3)
    ax.axvline (0,-2.5,2.0,linestyle='--',color='r',linewidth=3)
    ax.set_xlabel ("Reporter expression (log)", fontsize=24)
    ax.set_ylabel ("Average visits (log)", fontsize=24)
    r2 = model_r2 (population, expr_binned)
    xmin, xmax = ax.get_xaxis().get_view_interval()
    ymin, ymax = ax.get_yaxis().get_view_interval()
    deltax = xmax-xmin
    deltay = ymax-ymin
    ratio=0.25
    ax.text (xmax-ratio*deltax,ymax-ratio*deltay, "$r^2 = %.2f$"%r2,
            fontsize=24)
    ax2 = plt.subplot (gs[1,0])
    ax2.plot (np.log2(population))
    plt.show ()
