import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from pyterpol.synthetic.auxiliary import read_text_file

def plot_convergence(block, labels=None, relative=True, savefig=True, figname=None):
    """
    Plots convergence of the chi^2 and of individual parameters.
    :param block:
    :param labels:
    :param relative:
    :param savefig
    :param figname
    :return:
    """

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
        # define the color
        color = 0.1 +0.9*np.random.random(3)

        if labels is not None:
            ax.plot(block[:,i], '-', color=color, label=labels[i])
        else:
            ax.plot(block[:,i], '-', color=color)
        ax.set_xlabel('Iteration number')
        ax.set_ylabel('Relative value.')
        ax.legend(loc=1, fontsize=10)

    # save the plot
    if savefig == True:
        if figname is None:
            figname = 'convergence.png'

    plt.savefig(figname)

def plot_chi2_map(x, y, nbin=10, labels=None, savefig=True, figname=None):
    """
    Plots a covariance map.
    :param x parameter values
    :param y parameter values
    :param nbin number of bins in a histogram
    :param labels
    :param savefig
    :param figname
    :return:
    """
    # if user did not pass the labels
    if labels == None:
        labels = ['x', 'y']

    # set up the figure
    fig = plt.figure(figsize=(10,10), dpi=100)
    var_axes = [(2, 2, 1), (2, 2, 4)]
    var_data = [x, y]

    # firs the plot of the variance
    for i in range(0, 2):
        ax = fig.add_subplot(var_axes[i], aspect=1.0)

        # plot the histogram
        n, bins, patches = ax.hist(var_data[i], normed=True, label=labels[0])
        x_g = np.linspace(bins.min(), bins.max(), 50)

        # plot the gaussian 'fit'
        mean = var_data[i].mean()
        var = var_data[i].std(ddof=1)
        g = norm(loc=mean, scale=var)
        ax.plot(x_g, g.pdf(x_g), 'r-')

        # labeling
        ax.set_xlabel(labels[0])
        ax.set_ylabel('$n_i/N$')
        ax.title('$\sigma_%s$ = %.3f' % (labels[0], var))

    # plot the 2d chi2 map
    ax = fig.add_subplot(223, aspect=1.0)
    ax.hist2d(x, y, nbin, normed=True)

    # compute the correlation
    cov = ((x-x.mean())*(y-y.mean())).mean()
    cor = cov/(x.std(ddof=1)*y.std(ddof=1))

    # labelling
    ax.set_xlabel(labels[0])
    ax.set_ylabel(labels[1])
    ax.set_title(r'$\rho(%s, %s) = %.3f$' % (labels[0], labels[1], cor))

    # save the figure
    if savefig:
        if figname is None:
            figname = 'covariance' + '_'.join(labels) + 'png'

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
        # print d
        for hk in head_keys:
            if l.find(hk) > -1:
                # groups are integers of course
                if hk == 'group':
                    d[2:] = map(int, d[2:])
                else:
                    d[2:] = d[2:]

                # append the header info
                fitlog[hk] = d[2:]
                hkcounter += 1
                break

        # once we read all data, we end
        if hkcounter == 3:
            break

    # print fitlog
    # append data
    fitlog['data'] = np.loadtxt(f)

    return fitlog


















