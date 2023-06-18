from sherpa.astro.ui import *
from matplotlib import pyplot as plt

def plot_range(low, high, title=None, units='energy'):

    filt = get_filter()
    ignore()
    notice(low, high)
    plot(title=title, units=units)
    ignore()
    notice(filt)

def plot(title=None, units='energy'):
    set_ylog()
    if units == 'energy':
        set_xlog()
    plot_fit_ratio()
    if title is not None:
        ax1, ax2 = plt.gcf().axes
        ax1.set_title(title)
    plt.yscale('linear')
    return

    set_ylog()

    if title is not None:
        set_plot_title(title)
    if units == 'energy':
        log_scale(X_AXIS)

    current_plot('plot2')
    if units == 'energy':
        log_scale(X_AXIS)

    #plot_fit_delchi()
    #current_plot('plot2')
    #lin_scale(Y_AXIS)
