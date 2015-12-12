import matplotlib.pyplot as plt
import numpy as np

def plot_convergence(block, labels=None, relative=True, savefig=True, figname=None):
    """
    Plots convergence of the chi^2 and of individual parameters.
    :param block:
    :param labels:
    :param relative:
    :return:
    """
    # define the color
    color = 0.1 +0.9*np.random.random(3)

    nrow, ncol = np.shape(block)
    # normalize with the best value
    # if relative is p[assed
    if relative:
        for i in range(0, ncol):
            block[:,i] = block[:, i]/block[-1, i]

    # start a new figure
    fig = plt.figure(dpi=100, figsize=(15, 10))
    ax = fig.add_subplot(111)

    # plot convergence
    for i in range(0, ncol):
        if labels is not None:
            ax.plot(block[:,i], '-', color=color, label=labels[i])
        else:
            ax.plot(block[:,i], '-', color=color)

    # save the plot
    if savefig == True:
        if figname is None:
            figname = 'convergence.png'

    plt.savefig(figname)






