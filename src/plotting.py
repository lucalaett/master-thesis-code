
import matplotlib.pyplot as plt

    
def axes_plot(axes, x, y, label=None, xscale='linear', yscale='linear', xlabel=None, ylabel=None, p_read=1.5):
    """
    Plots datapoints with preset color, marker and fontsizes

    Parameters [all float, all in GeV^x]
    ----------
    axes : axes to plot on
    x, y : datapoint coordinates
    label, xscale, yscale, xlabel, ylabel : standard matplotlib.plot inputs
    p_read : factor to change fontsizes
    """
    axes.plot(x, y, color='black', linestyle='None', marker='d', markeredgecolor='black', markersize=3, label=label)
    axes.set_xscale(xscale)
    axes.set_yscale(yscale)
    axes.tick_params(axis='both', labelsize=8*p_read)
    plt.grid(axis='y', color='gainsboro') 
    plt.xlabel(xlabel, fontsize=10*p_read)
    plt.ylabel(ylabel, fontsize=10*p_read)
    plt.tight_layout()