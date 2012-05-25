from matplotlib import pyplot as mpl
import numpy as np

def relabel(tick, bp_equiv):
    if tick < 0:
        return "%d bp" % tick
    elif tick == 0:
        return "-0 bp/0%"
    elif tick < bp_equiv:
        return "%d%%" % (float(tick)/bp_equiv * 100)
    elif tick == bp_equiv:
        return '100%/+0 bp'
    else:
        return "+%d bp" % (tick - bp_equiv)

def plot_averaged_genes(upstream, gene, downstream, bp_equiv = 200):
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
    retval.append(ax.plot(np.linspace(0, bp_equiv, len(gene)),
                           gene,
                           label='Gene'))
    retval.append(ax.plot(np.arange(bp_equiv, len(downstream) + bp_equiv),
                           downstream,
                           label='Downstream'))
    ymin, ymax = ax.get_ybound()
    retval.append(ax.vlines([0, bp_equiv], ymin, ymax , linestyles='dashed'))
    ticks = list(ax.get_xticks())
    ticks.extend([0,bp_equiv/2,bp_equiv]) # Ensure we have at least the full gene labelled
    ticks = np.unique(ticks)
    ticklabels = [relabel(tick, bp_equiv) for tick in ticks]
    ax.set_xticks(ticks)
    ax.set_xticklabels(ticklabels)
    mpl.show()
    return retval


