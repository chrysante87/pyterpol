import matplotlib.pyplot as plt
import numpy as np
from pyterpol.synthetic.auxiliary import read_text_file

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

def read_fitlog(f):
    """
    Reads the fitting log and stores it within a dictionary.
    :param f:
    :return:
    """
    # read the file
    lines = read_text_file(f)

    # key counter and ouput dictionary
    fitlog = {}
    hkcounter = 0

    # define header keys
    head_keys = ['name', 'component', 'group']
    for l in lines:
        d = l.split()
        for hk in head_keys:
            if l.find(hk) > -1:
                hkcounter += 1
                # groups are integers of course
                if hk == 'group':
                    d = map(int, d[2:])
                else:
                    d = d[2:]
                break

        # append the header info
        fitlog[hk] = d
        # once we read all data, we end
        if hkcounter == 3:
            break
    # append data
    fitlog['data'] = np.loadtxt(f)

    return fitlog












