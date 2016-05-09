import copy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import numpy as np
from scipy.stats import norm
from pyterpol.synthetic.auxiliary import read_text_file


def get_walker(db, nchain, nwalker, niter):
    """
    Retrieves a walker from the chain.
    :param db:
    :param nchain:
    :param nwalker:
    :param niter:
    :return:
    """
    rows = np.arange(niter)
    rows = nchain + nwalker*rows
    return db[rows]

def plot_walkers_for_one_param(db, ipar, nwalker, niter, ax):
    """
    :param db:
    :param ipar:
    :param nwalker:
    :param niter:
    :param ax:
    :return:
    """
    # set teh iterations
    iters = np.arange(niter)
    # plot each walker
    for i in range(0, nwalker):
        w = get_walker(db, i, nwalker, niter)
        ax.plot(iters, w[:,ipar], '-')

def plot_walkers(block, niter, nwalker, indices=None, labels=None, savefig=True, figname=None):
    """
    :param block:
    :param indices:
    :param niter:
    :param nwalker:
    :param labels:
    :param savefig:
    :param figname:
    :return:
    """

    if figname is not None:
        savefig = True

    # define which parameters are plotted
    if indices is None:
        indices = np.arange(len(block[0]))
    npar = len(indices)

    # definethe plotting grid
    ncol = 3
    nrow = npar / ncol
    if npar % ncol > 0:
        nrow += 1

    # create the grid and the figure
    gs1 = gs.GridSpec(nrow, ncol, hspace=0.2, wspace=0.4)
    fig = plt.figure(figsize=(4*ncol, 3*nrow), dpi=100)

    # plot each figure
    for j, ind in enumerate(indices):

        # set label
        if labels is None:
            label = 'p' + str(ind).zfill(2)
        else:
            label = labels[j]

        # set the position
        icol = j % ncol
        irow = j / ncol
        ax = fig.add_subplot(gs1[irow, icol])

        # plot the walkers
        plot_walkers_for_one_param(block, ind, nwalker, niter, ax)
        ax.set_xlabel('Iteration number', fontsize=8)
        ax.set_ylabel(label, fontsize=8)

    # save the figure
    if savefig:
        if figname is None:
            figname = 'mcmc_convergence.png'

        # plt.tight_layout()
        plt.savefig(figname)


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
        rel_block = copy.deepcopy(block)
        for i in range(0, ncol):
            rel_block[:,i] = block[:, i]/block[-1, i]

    # start a new figure
    fig = plt.figure(dpi=100, figsize=(15, 10))
    ax = fig.add_subplot(111)

    # plot convergence
    for i in range(0, ncol):
        # define the color
        color = 0.1 +0.9*np.random.random(3)

        if labels is not None:
            ax.plot(rel_block[:,i], '-', color=color, label=labels[i])
        else:
            ax.plot(rel_block[:,i], '-', color=color)
        ax.set_xlabel('Iteration number')
        ax.set_ylabel('Relative value.')
        ax.legend(loc=1, fontsize=10)

    # save the plot
    if savefig == True:
        if figname is None:
            figname = 'convergence.png'

    plt.savefig(figname)

    # try to produce another kind of plot
    if ncol % 2 > 0:
        nfigrow = ncol / 2 + 1
    else:
        nfigrow = ncol / 2

    # setup the grid
    gs1 = gs.GridSpec(nfigrow, 2, hspace=0.5)

    # setup the figure
    fig2 = plt.figure(dpi=100, figsize=(10, 3*nfigrow))

    # plot convergence of individual parameters
    for i in range(0, ncol):
        ax = fig2.add_subplot(gs1[i/2, i%2])
        ax.set_xlabel('Iteration number')
        ax.set_ylabel('Value')
        ax.set_ylabel(labels[i], fontsize=8)
        ax.plot(block[:, i], 'k-', label=labels[i])
        # ax.legend(loc=1)

    # save the figure
    fig2.savefig('convergence_2.png')
    plt.close()



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
    fs=8
    # if user did not pass the labels
    if labels == None:
        labels = ['x', 'y']

    # set up the figure
    fig = plt.figure(figsize=(10,10), dpi=100)
    var_axes = [221, 224]
    var_data = [x, y]

    # firs the plot of the variance
    for i in range(0, 2):
        ax = fig.add_subplot(var_axes[i])

        # plot the histogram
        n, bins, patches = ax.hist(var_data[i], nbin, normed=True, label=labels[0])
        x_g = np.linspace(bins.min(), bins.max(), 50)

        # plot the gaussian 'fit'
        mean = var_data[i].mean()
        var = var_data[i].std(ddof=1)
        g = norm(loc=mean, scale=var)
        ax.plot(x_g, g.pdf(x_g), 'r-')

        # labeling
        ax.set_xlabel(labels[i], fontsize=8)
        ax.set_ylabel('$n_i/N$', fontsize=8)
        ax.set_title(r'$\sigma$_%s=%.3f' % (labels[i], var), fontsize=8)

    # plot the 2d chi2 map
    ax = fig.add_subplot(223)
    ax.hist2d(x, y, nbin, normed=True)

    # compute the correlation
    cov = ((x-x.mean())*(y-y.mean())).mean()
    cor = cov/(x.std(ddof=1)*y.std(ddof=1))

    # labelling
    ax.set_xlabel(labels[0], fontsize=8)
    ax.set_ylabel(labels[1], fontsize=8)
    ax.set_title(r'$\rho$(%s, %s) = %.3f' % (labels[0], labels[1], cor), fontsize=8)

    # save the figure
    if savefig:
        if figname is None:
            figname = '_'.join(labels) + '.png'

        plt.savefig(figname)
        plt.close()


def plot_variance(x, nbin=10, label=None, savefig=True, figname=None):
    """
    Plots a covariance map.
    :param x parameter values
    :param nbin number of bins in a histogram
    :param labels
    :param savefig
    :param figname
    :return:
    """
    fs=8
    # if user did not pass the labels
    if label is None:
        label = 'x'

    # set up the figure
    fig = plt.figure(figsize=(6,6), dpi=100)

    # firs the plot of the variance
    ax = fig.add_subplot(111)

    # plot the histogram
    n, bins, patches = ax.hist(x, nbin, normed=True, label=label)
    x_g = np.linspace(bins.min(), bins.max(), 50)

    # plot the gaussian 'fit'
    mean = x.mean()
    var = x.std(ddof=1)
    g = norm(loc=mean, scale=var)
    ax.plot(x_g, g.pdf(x_g), 'r-')

    # labeling
    ax.set_xlabel(label, fontsize=8)
    ax.set_ylabel('$n_i/N$', fontsize=8)
    ax.set_title(r'$\sigma$_%s=%.3f' % (label, var), fontsize=8)
    ax.legend(fontsize=8)

    # save the figure
    if savefig:
        if figname is None:
            figname = label + '.png'
        else:
            figname += label+'.png'

        plt.savefig(figname)
        plt.close()

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

def read_mc_chain(f):
    """
    Reads the mcmc chain created with emcee
    :param f: chain_file
    :return:
    """

    # read the file
    lines = read_text_file(f)

    # key counter and ouput dictionary
    chainlog = {}
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
                chainlog[hk] = d[2:]
                hkcounter += 1
                break

        # once we read all data, we end
        if hkcounter == 3:
            break

    # load the file
    d = np.loadtxt(f)

    # get fit properties
    nwalkers = int(np.max(d[:, 0])) + 1
    niter = len(d[:, 0]) / nwalkers
    npars = len(d[0]) - 2

    # remove the first column with numbering
    chainlog['data'] = d[:, 1:]
    return chainlog, nwalkers, niter, npars























