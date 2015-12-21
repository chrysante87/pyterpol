"""
V746Cas - fitting of a observed spectra.
This example also show, ho we can proceed if
we want to fit parameters step by step.
"""

import pyterpol
import numpy as np
import matplotlib.pyplot as plt

def inspect_spectra(f):
    """
    Plots all spectra.
    :param f:
    :return:
    """
    ifile = open(f, 'r')
    slist = ifile.readlines()
    ifile.close()

    for rec in slist:
        ifile = open(rec.rstrip('\n'), 'r')
        ifile.readline()

        x, y = np.loadtxt(ifile, unpack=True, usecols=[0,1])

        plt.plot(x, y, '-')
    plt.show()

def read_obs_from_list(f):
    """
    Create a list of observations.
    :param f:
    :return:
    """
    ifile = open(f, 'r')
    slist = ifile.readlines()
    ifile.close()

    obs = []
    for i,rec in enumerate(slist[:]):
        o = pyterpol.ObservedSpectrum(filename=rec.rstrip('\n'), group=dict(rv=i))
        o.get_sigma_from_fft(store=True)
        obs.append(o)

    return obs

def setup_interface_more_obs():

    # have a look at one observation
    ol = pyterpol.ObservedList()
    obs = read_obs_from_list('spec.lis')
    ol.add_observations(obs)

    # define a starlist
    sl = pyterpol.StarList()
    sl.add_component(component='primary', teff=15800., logg=4.13, vrot=148., z=1.9, lr=1.0)

    # define regions
    rl = pyterpol.RegionList()
    rl.add_region(wmin=6340, wmax=6410, groups=dict(lr=0))
    rl.add_region(wmin=6540, wmax=6595, groups=dict(lr=0))
    rl.add_region(wmin=6670, wmax=6685, groups=dict(lr=0))

    # create interfaces
    itf = pyterpol.Interface(sl=sl, rl=rl, ol=ol, debug=False)

    # set fit order = 2 to do it fast
    itf.set_grid_properties(order=2, step=0.05)
    itf.setup()

    #review
    print itf

    # save the session
    itf.save('initial.itf')

def optimize_rv(session0, session1):
    """
    Optimizes RV.
    :return:
    """

    # setup the spectra
    itf = pyterpol.Interface.load(session0)

    # set parameters
    itf.set_parameter(parname='rv', fitted=True, vmin=-60., vmax=60.)
    itf.choose_fitter('nlopt_nelder_mead', ftol=1e-6)

    print itf

    # run fit
    itf.run_fit()

    # write result
    itf.write_fitted_parameters(outputname=session1.split('.')[0]+'.dat')

    # plot every comparison
    itf.plot_all_comparisons(figname='rvfit')

    # save the fit
    itf.save(session1)

def optimize_all(session0, session1):
    """
    Optimizes all parameters
    :return:
    """

    # setup the spectra
    itf = pyterpol.Interface.load(session0)
    itf.set_parameter(parname='rv', fitted=True, vmin=-60., vmax=60.)
    itf.set_parameter(parname='teff', fitted=True, vmin=15000., vmax=18000.)
    itf.set_parameter(parname='logg', fitted=True, vmin=3.8, vmax=4.75)
    itf.set_parameter(parname='vrot', fitted=True, vmin=100., vmax=180.)
    itf.set_parameter(parname='z', fitted=True, vmin=1.2, vmax=2.0)
    itf.choose_fitter('nlopt_nelder_mead', ftol=1e-5)

    # run fit
    itf.run_fit()

    # write result
    itf.write_fitted_parameters(outputname=session1.split('.')[0]+'.dat')

    # plot every comparison
    itf.plot_all_comparisons(figname='allfit')

    # save the fit
    itf.save(session1)

# inspect_spectra('spec.lis')
setup_interface_more_obs()
# optimize_rv('initial.itf', 'rvfit.itf')
optimize_all('initial.itf', 'nmallfit.itf')

# lets do some plotting
itf  = pyterpol.Interface.load('nmallfit.itf')
# itf.plot_convergence(figname='chi2.png')
# itf.plot_convergence(parameter='all', figname='convergence_params.png')
# itf.plot_covariances(nbin=50, parameters=['z', 'logg', 'teff', 'vrot'])
# itf.plot_variances(nbin=30, parameters=['rv'])
# itf.write_fitted_parameters(outputname='trial.res')

itf.set_error(error=10)
fitparams = itf.get_fitted_parameters(attribute='value')
itf.fitter.run_mcmc(itf.compute_chi2, 'chain.dat', fitparams, 4*len(fitparams), 200)

# ol = pyterpol.ObservedList()
# ol.load('initial.itf')
# print ol