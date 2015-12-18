"""
V746Cas - fitting of a observed spectra.
This example also show, ho we can proceed if
we want to fit parameters step by step.
"""

import pyterpol
import matplotlib.pyplot as plt

def inspect_spectra(f):
    ifile = open(f, 'r')
    slist = ifile.readlines()
    ifile.close()

    for rec in slist:
        ifile = open(rec.rstrip('\n'), 'r')
        ifile.readline()

        x, y = np.loadtxt(ifile, unpack=True, usecols=[0,1])

        plt.plot(x, y, '-')
        plt.show()

# Setting up interface is something, that should be kept
# separated from the fitting, because consequetive fits
# change teh initial settings.
def setup_interface_single_obs():

    # have a look at one observation
    ol = pyterpol.ObservedList()
    obs = pyterpol.ObservedSpectrum(filename='v7c00001.asc', group=dict(rv=0), error=0.01)

    # two methods how to estimate the error
    print obs.get_sigma_from_continuum(cmin=6665, cmax=6670)
    print obs.get_sigma_from_fft()
    ol.add_observations([obs])

    # define a starlist
    sl = pyterpol.StarList()
    sl.add_component(component='primary', teff=17000., logg=4.0, vrot=180., z=1.0, lr=1.0)

    # define regions
    rl = pyterpol.RegionList()
    rl.add_region(wmin=6340, wmax=6410, groups=dict(lr=0))
    rl.add_region(wmin=6520, wmax=6610, groups=dict(lr=0))
    rl.add_region(wmin=6665, wmax=6690, groups=dict(lr=0))

    # create interfaces
    itf = pyterpol.Interface(sl=sl, rl=rl, ol=ol)

    # set fit order = 2 to do it fast
    itf.set_grid_properties(order=2)
    itf.setup()

    # review the result - one rv group, one lr group
    print itf

    # plot comparison
    itf.plot_all_comparisons(figname='teff17000')

    # try different temperatures - this way we can easilyt review
    # several comparisons
    itf.set_parameter(parname='teff', value=25000.)
    itf.populate_comparisons()
    itf.plot_all_comparisons(figname='teff25000')

    itf.set_parameter(parname='teff', value=13000.)
    itf.populate_comparisons()
    itf.plot_all_comparisons(figname='teff13000')
    itf.save('initial.itf')


# if we want to fit interactively, parameter by parameter
# it is easier to use the save/load mechanism
# itf = pyterpol.Interface.load('tefffit.itf')
# itf = pyterpol.Interface.load('vrotfit.itf')
itf = pyterpol.Interface.load('loggfit.itf')

# choose a fitter
itf.choose_fitter('nlopt_nelder_mead', ftol=1e-6)

# change another parameter
# itf.set_parameter(parname='vrot', vmin=120., vmax=200., fitted=True)
# itf.set_parameter(parname='logg', vmin=3.5, vmax=4.5, fitted=True)
itf.set_parameter(parname='z', vmin=0.5, vmax=2.0, fitted=True)
itf.set_parameter(parname='teff', vmin=15000., fitted=True)
itf.run_fit()

# get the result
# itf.plot_all_comparisons(figname='vrotfit')
itf.plot_all_comparisons(figname='zfit')

# itf.write_fitted_parameters(outputname='iter03.res')
itf.write_fitted_parameters(outputname='iter05.res')

# save a new Interface
# itf.save('vrotfit.itf')
itf.save('zfit.itf')


