from matplotlib import pyplot as mpl
import numpy as np

def relabel(tick):
    if tick < 0:
        return "%d bp" % tick
    elif tick == 0:
        return "-0 bp/0%"
    elif tick < 100:
        return "%d%%" % tick
    elif tick == 100:
        return '100%/+0 bp'
    else:
        return "+%d bp" % (tick - 100)

def plot_averaged_genes(upstream, gene, downstream):
    """ Plots the coverage around genes.

    Here, the upstream and downstream regions are assumed to be in base pairs,
    and the gene itself is assumed to cover 100% of the gene length, averaged
    over all the genes included
    """

    retval = []

    fig = mpl.gcf()
    ax = fig.gca()
    retval.append(ax.plot(np.arange(-len(upstream), 0),
                           upstream,
                           label='Upstream'))
    retval.append(ax.plot(np.linspace(0, 100, len(gene)),
                           gene,
                           label='Gene'))
    retval.append(ax.plot(np.arange(100, len(downstream) + 100),
                           downstream,
                           label='Downstream'))
    ymin, ymax = ax.get_ybound()
    retval.append(ax.vlines([0, 100], ymin, ymax , linestyles='dashed'))
    ticks = list(ax.get_xticks())
    ticks.extend([0,50,100]) # Ensure we have at least the full gene labelled
    ticks = np.unique(ticks)
    ticklabels = [relabel(tick) for tick in ticks]
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)
    mpl.show()
    return retval


