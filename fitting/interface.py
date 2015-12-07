# -*- coding: utf-8 -*-
import copy
import warnings
import numpy as np
import matplotlib.pyplot as plt
from pyterpol.synthetic.makespectrum import SyntheticGrid
from pyterpol.synthetic.makespectrum import SyntheticSpectrum
from pyterpol.observed.observations import ObservedSpectrum
from pyterpol.fitting.parameter import Parameter
from pyterpol.fitting.parameter import parameter_definitions
from pyterpol.fitting.fitter import Fitter
from pyterpol.synthetic.auxiliary import generate_least_number
from pyterpol.synthetic.auxiliary import keys_to_lowercase
from pyterpol.synthetic.auxiliary import parlist_to_list
from pyterpol.synthetic.auxiliary import sum_dict_keys
from pyterpol.synthetic.auxiliary import ZERO_TOLERANCE


# TODO Go through similarities in the classes and write a parent class
# TODO to get rid of the redundant code.

# repeat userwarnings
warnings.simplefilter('always', UserWarning)

class Interface(object):
    """
    """
    def __init__(self, sl=None, rl=None, ol=None, fitter=None, debug=False, adaptive_resolution=True):
        """
        :param sl: StarList type
        :param rl: RegionList type
        :param ol: ObservedList type
        :param fitter
        :param debug
        :param adaptive resolution - this (sounds better than it actually is)
                just means that resolution of the grid is set to twice
                the resolution of the spectrum with highest resolution
        :return:
        """

        self.sl = sl
        self.rl = rl
        self.ol = ol
        self.synthetics = {}
        self.grids = {}
        self.fitter = fitter

        # debug mode
        self.debug = debug

        # define empty comparison list
        self.comparisonList = None

        # parameters that cannot be obatined through interpolation
        self._not_given_by_grid = ['lr', 'rv', 'vrot']

        # relation between rv_groups and regions
        self.rel_rvgroup_region = {}

        # properties of synthetic spectra
        self._synthetic_spectrum_kwargs = {}

        # properties of grids
        self._grid_kwargs = {}

        # initialization of various boolean variables
        self.grid_properties_passed = False
        self.fit_is_running = False
        self.adaptive_resolution = adaptive_resolution

    def __str__(self):
        """
        String representation of the class
        :return:
        """
        string = ""
        for attr, name in zip(['sl', 'rl', 'ol', 'fitter'], ['StarList', 'RegionList', 'ObservedList', 'Fitter']):
            string += '%s%s\n' % (name[:len(name)/2].rjust(50, '='), name[len(name)/2:].ljust(50, '='))
            string += str(getattr(self, attr))
        # too much information
        # for key in self.grids.keys():
        #     string += 'Grid for region: %s.\n' % key
        #     string += str(self.grids[key])
        string += ''.ljust(100,'=')

        return string

    def accept_fit(self):
        """
        Propagates the fitting result to the class.
        :return:
        """

        # this should be done more carefully
        # print self.fitter.result
        final_pars = self.fitter.result

        # list fitted parameters
        fitparams = self.get_fitted_parameters()

        # updates the parameters with the result
        for i in range(0, len(final_pars)):
            fitparams[i]['value'] = final_pars[i]

        # update the fitter with new initial parameters
        self.fitter.par0 = copy.deepcopy(final_pars)



    def add_comparison(self, region=None, parameters={}, observed=None, synthetic={}, groups={}):
        """
        :param region the name of the corresponding region
        :param parameters a dictionary of the parameters required for the synthetic
                            spectrum
        :param observed the observed spectrum
        :param synthetic empty dictionaryf for the synthetic spectra
        :param groups

        Add a record to the comparisonList
        :return: None

        """

        if self.debug:
            print 'Settting comparison for region: %s \n groups: %s. \n parameters: %s' % \
                  (str(region), str(groups), str(parameters))


        if self.comparisonList is None:
            raise Exception('The comparisonList has not been defined yet. Use Inteface.get_comparison for that.')
        else:
            self.comparisonList.append(dict(region=region,
                                        parameters = parameters,
                                        observed = observed,
                                        groups = groups,
                                        synthetic = {x: None for x in parameters.keys()},
                                        chi2 = 0.0
                                     ))
    def compute_chi2(self, pars, l=None, verbose=False):
        """
        :param pars:
        :param l
        :param verbose
        :return: chi square
        """
        if l is None:
            l = self.comparisonList

        # propagate the parameters to the
        # parameterlist and update it
        self.propagate_and_update_parameters(l, pars)

        # reads out the chi_2 from individual spectra
        chi2 = self.read_chi2_from_comparisons(l, verbose)

        # store chi_square within the fitter
        if verbose:
            chi2 = chi2[0]
            chi2_detailed = chi2[1]
        else:
            chi2_detailed = []

        # if we are fitting we store the info on the parameters
        if self.fit_is_running:
            # print dict(parameters=pars, chi2=chi2, detailed=chi2_detailed)
            self.fitter.append_iteration(dict(parameters=pars, chi2=chi2, detailed=chi2_detailed))

        # if self.debug:
            # print 'Computed model: %s chi2: %s' % (str(pars), str(chi2))
        print 'Computed model: %s chi2: %s' % (str(pars), str(chi2))

        return chi2

    def clear_all(self):
        """
        Clears the class.
        :return:
        """
        self.comparisonList = None
        self.grids = {}
        self.ol = None
        self.rl = None
        self.sl = None
        self.fitter = None
        self.synthetics = {}
        self._grid_kwargs = {}
        self._synthetic_spectrum_kwargs = {}
        self.rel_rvgroup_region = {}
        self.grid_properties_passed = False

    def choose_fitter(self, *args, **kwargs):
        """
        Just wrapper for the Fitter.choose_fitter method
        see parameter descriptio there.
        :param args:
        :param kwargs:
        :return:
        """
        # fitter is rather simple, so if there is none set, we set an empty
        # one
        # print kwargs
        if self.fitter is None:
            self.fitter = Fitter(debug=self.debug)
        self.fitter.choose_fitter(*args, **kwargs)

    def extract_parameters(self, l, attr='value'):
        """
        Converts a list of parameter class to a
        dictionary.
        :param l
        :param attr
        :return:
        """
        params = {par['name']: par['value'] for par in l}
        return params

    def get_comparisons(self, verbose=False, **kwargs):
        """
        Narrows down the number of comparisons.
        :param verbose
        :param kwargs
        :return:
        """

        clist = []
        indices = []
        keys = kwargs.keys()

        for i in range(0, len(self.comparisonList)):
            # print i, clist
            include = True
            for key in keys:
                # what if the key lies
                if key in self.comparisonList[i]['groups'].keys() \
                        and (kwargs[key] != self.comparisonList[i]['groups'][key]):
                    # print key, kwargs[key], self.comparisonList[i]['groups'][key]
                    include = False
                    break
                if hasattr(self.comparisonList[i]['observed'], key) and \
                             self.comparisonList[i]['observed'].key != kwargs[key]:
                    include = False
                    break
                if key == 'region' and self.comparisonList[i]['region'] != kwargs[key]:
                    include = False
                    break
            # if it survived all tests it is included
            if include:
                clist.append(self.comparisonList[i])
                indices.append(i)

        # if we want to get indices of the found in the original array
        if verbose:
            return clist, indices
        else:
            return clist

    def get_fitted_parameters(self):
        """
        lists all fitted parameters
        :return:
        """

        return self.sl.get_fitted_parameters()

    def get_observed_spectrum(self, filename=None):
        """
        Returns
        :return:
        """
        return self.ol.get_spectra(filename=filename)[0]

    def list_comparisons(self, l=None):
        """
        This function displays all comparisons.
        :return: string
        """
        if l is None:
            l = self.comparisonList

        string = ''
        for i,rec in enumerate(l):
            string += "========================= Comparison %s =========================\n" % str(i).zfill(3)
            reg = rec['region']
            # list region
            string += 'region: %s:(%s,%s)\n' % (reg, str(self.rl.mainList[reg]['wmin']),
                                                str(self.rl.mainList[reg]['wmax']))

            # list observed spectrum
            if rec['observed'] is not None:
                string += "observed: %s\n" % rec['observed'].filename
            else:
                string += "observed: NONE\n"

            # lists all parameters
            for c in rec['parameters'].keys():
                string += 'component: %s ' % c
                # print rec['parameters'][c]
                for par in rec['parameters'][c]:
                    string += "%s: %s " % (par['name'], str(par['value']))
                string += '\n'

            # list all groups
            string += 'groups: %s\n' % str(rec['groups'])
            string += 'chi2: %s\n' % str(rec['chi2'])
        string += "==================================================================\n"

        return string

    def populate_comparisons(self, l=None, demand_errors=False):
        """
        Creates a synthetic spectrum for every record in
        the comparisonList.
        :param l
        :param demand_errors
        :return:
        """
        if l is None:
            l = self.comparisonList
        # go over ech comparison in the list
        for rec in l:

            # print rec
            # get the region
            region = rec['region']
            wmin = self.rl.mainList[region]['wmin']
            wmax = self.rl.mainList[region]['wmax']

            # go over each component
            for c in rec['parameters'].keys():
                pars = self.extract_parameters(rec['parameters'][c])

                # use only those parameters that are not constrained with the grid
                pars = {x:pars[x] for x in pars.keys() if x in self._not_given_by_grid}
                # print pars

                # populate with the intensity vector of each component
                # print rec['observed']
                if rec['observed'] is not None:
                    try:
                        wave, intens, error = rec['observed'].get_spectrum(wmin, wmax)
                    except:

                        # if we are fitting we will demand errors
                        if demand_errors:
                            raise ValueError('It is not allowed to call chi-square without having'
                                             ' uncertainties set. SET THE ERRORS FFS!')
                        else:
                            wave = rec['observed'].get_spectrum(wmin, wmax)[0]
                            error = None

                    korelmode = rec['observed'].korel
                else:
                    wave = np.arange(wmin, wmax, 0.01)
                    korelmode = False

                # print korelmode, pars, wave.min(), wave.max()

                rec['synthetic'][c] = self.synthetics[region][c].get_spectrum(wave=wave,
                                                                              only_intensity=True,
                                                                              korel=korelmode,
                                                                              **pars)

            # it is mandatory to provide errors for
            # computation of the chi2
            if error is not None:
                # sum component spectra
                for i,c in enumerate(rec['synthetic'].keys()):
                    if i == 0:
                        syn = rec['synthetic'][c].copy()
                    else:
                        syn = syn + rec['synthetic'][c]

            # setup the chi2
            rec['chi2'] = np.sum(((intens - syn)/error)**2)

    def plot_all_comparisons(self, l=None, figname=None):
        """
        Creates a plot of all setup comparisons.
        :param l
        :param figname
        :return: None
        """
        if l is None:
            l = self.comparisonList
        if len(l) == 0:
            raise ValueError('The comparison list is empty. Did you run interface.setup() and interface.populate()?')
        for i in range(0, len(l)):
            self.plot_comparison_by_index(i, l=l, savefig=True, figname=figname)

    def plot_comparison_by_index(self, index, l=None, savefig=False, figname=None):
        """
        :param index
        :param l
        :param savefig
        :param figname
        :return:
        """

        # the comparison
        if l is None:
            cpr = self.comparisonList[index]
        else:
            cpr = l[index]

        # boundaries
        reg = cpr['region']
        wmin = self.rl.mainList[reg]['wmin']
        wmax = self.rl.mainList[reg]['wmax']

        # merge the spectra
        if any([cpr['synthetic'][key] == None for key in cpr['synthetic'].keys()]):
            raise ValueError('The synthetic spectra are not computed. Did you run Interface.populate_comparisons()?')
        si = sum_dict_keys(cpr['synthetic'])

        # names
        if cpr['observed'] is not None:
            obsname = cpr['observed'].filename
        else:
            obsname = 'NONE'
        synname = ''
        for c in cpr['parameters']:
            synname += 'Component: %s ' % c
            synname += str(self.extract_parameters(cpr['parameters'][c])) + '\n'

        if cpr['observed'] is not None:
            try:
                w, oi, ei = cpr['observed'].get_spectrum(wmin, wmax)
            except:
                w, oi, = cpr['observed'].get_spectrum(wmin, wmax)
                ei = np.zeros(len(w))
                warnings.warn('Your data observed spectrum: %s has not errors attached!' )
        else:
            w = np.linspace(wmin, wmax, len(si))

        if figname is None:
            figname = "_".join([obsname, 'wmin', str(int(wmin)), 'wmax', str(int(wmax))]) +'.png'
        else:
            figname = "_".join([figname, obsname, 'wmin', str(int(wmin)), 'wmax', str(int(wmax))]) +'.png'

        if self.debug:
            print "Plotting comparison: observed: %s" % obsname
            print "Plotting comparison: synthetics: %s" % synname

        # do the plot
        fig = plt.figure(figsize=(16, 10), dpi=100)
        ax = fig.add_subplot(211)

        if cpr['observed'] is not None:
            ax.errorbar(w, oi, yerr=ei, fmt='-', color='k', label=obsname)
        ax.plot(w, si, 'r-', label=synname)
        ax.set_xlim(wmin, wmax)
        ax.set_ylim(0.95*si.min(), 1.05*si.max())
        ax.set_xlabel('$\lambda(\AA)$')
        ax.set_ylabel('$F_{\lambda}$(rel.)')
        ax.legend(fontsize=8, loc=3)

        if cpr['observed'] is not None:
            ax = fig.add_subplot(212)
            resid = oi-si
            ax.plot(w, resid, 'y', label='residuals')
            ax.set_xlabel('$\lambda(\AA)$')
            ax.set_ylabel('$F_{\lambda}$(rel.)')
            ax.set_xlim(wmin, wmax)
            ax.set_ylim(0.95*resid.min(), 1.05*resid.max())
            ax.legend(fontsize=8, loc=3)

        # save the figure
        if savefig:
            plt.savefig(figname)

    def propagate_and_update_parameters(self, l, pars):
        """
        :param l
        :param pars
        :return:
        """

        # parameters are passed by reference, so
        # this should also change the starlist
        # and corresponding
        fitpars = self.sl.get_fitted_parameters()

        if len(pars) != len(fitpars):
            raise ValueError('Length of the vector passed with the fitting environment does '
                             'mot match length of the parameters marked as fitted.')

        for i,v in enumerate(pars):
            # print fitpars[i]['value'], v
            fitpars[i]['value'] = v


        # we have to recompute the synthetic spectra
        # if one grid parameter was passed
        # first check for which parameters
        # the grid parameters are fitted
        components_to_update = []
        for c in self.sl.fitted_types.keys():
            # print self.sl.fitted_types, components_to_update
            for rec in self.sl.fitted_types[c]:
                # print rec
                if rec not in self._not_given_by_grid:
                    components_to_update.append(c)

        # update the syntrhetic spectra
        if len(components_to_update) > 0:
            self.ready_synthetic_spectra(complist=components_to_update)

        # populate the comparison
        self.populate_comparisons(l=l, demand_errors=True)


    def ready_synthetic_spectra(self, complist=[]):
        """
        Readies the synthetic spectra for each region.
        :param complist
        :return:
        """
        # if there is no list of components
        # for which to set the synthetic
        # parameters
        if len(complist) == 0:
            complist = self.sl._registered_components

        for reg in self.rl._registered_regions:
            # add the region to synthetics
            if reg not in self.synthetics.keys():
                self.synthetics[reg] = dict()

            # wavelength_boundaries
            wmin = self.rl.mainList[reg]['wmin']
            wmax = self.rl.mainList[reg]['wmax']

            # get all parameters for a given region
            reg_groups = self.rl.mainList[reg]['groups'][0]
            reg_groups = {x: reg_groups[x] for x in reg_groups.keys() \
                          if x not in self._not_given_by_grid}
            grid_pars = [x for x in self.sl.get_physical_parameters() \
                         if x not in self._not_given_by_grid]

            # print grid_pars, reg_groups
            # setup default groups - ie zero
            for par in grid_pars:
                if par not in reg_groups.keys():
                    reg_groups[par] = 0

            # get list of Parameters
            parlist = self.sl.get_parameter(**reg_groups)

            for c in complist:
                # convert Parameter list to dictionary
                params = self.extract_parameters(parlist[c])
                # print params

                # padding has to be relatively large, since
                # we do not know what the rvs will be
                self.synthetics[reg][c] = self.grids[reg].get_synthetic_spectrum(params,
                                                                                 np.array([wmin, wmax]),
                                                                                 **self._synthetic_spectrum_kwargs)
    def read_chi2_from_comparisons(self, l=None, verbose=False):
        """
        Reads the chi-squares from the list.
        :param l:
        :return:
        """

        # work with the min comparisonList if no other
        # is provided
        if l is None:
            l = self.comparisonList

        #
        chi2 = 0.0
        if verbose:
            chi2_detailed = []

        # read out the chi squares
        for i in range(0, len(l)):
            chi2 += l[i]['chi2']

            # if verbosity is desired a detailed chi-square
            # info on each region is returned
            if verbose:
                self.chi2_detailed.append(dict(chi2=l[i]['chi2'],
                                          region=self.rl.mainList[l[i]['region']],
                                          rv_group=l[i]['groups']['rv']))
        if verbose:
            return chi2, chi2_detailed
        else:
            return chi2

    def ready_comparisons(self):
        """
        This function creates a dictionary, which is one of the
        cornerstones of the class. It creates a list of all
        combinations of the parameters.
        :return:
        """

        # start a list of comparisons that will
        # be carried out with the given dataset
        self.comparisonList = []

        # go region by region
        for reg in self.rl.mainList.keys():
        # reg = self.rl.mainList.keys()[0]
            # fitted region
            wmin = self.rl.mainList[reg]['wmin']
            wmax = self.rl.mainList[reg]['wmax']

            # region-dfined groups
            reg_groups = copy.deepcopy(self.rl.mainList[reg]['groups'][0])
            phys_pars = [x for x in self.sl.get_physical_parameters() if x not in ['rv']]
            # print reg, phys_pars, reg_groups

            # if the group is not defined, it is zero
            for par in phys_pars:
                if par not in reg_groups.keys():
                    reg_groups[par] = 0

            # print reg_groups

            # extract all parameters
            reg_pars = self.sl.get_parameter(**reg_groups)

            # create a list of unique rv groups
            rv_groups = self.sl.get_defined_groups(parameter='rv')
            rv_groups = [rv_groups[key]['rv'] for key in rv_groups.keys()]
            # print rv_groups
            temp = []
            for row in rv_groups:
                temp.extend(row)
            rv_groups = np.unique(temp)
            # print rv_groups
            # rv_groups = np.unique(np.ravel([rv_groups[key]['rv'] for key in rv_groups.keys()])).tolist()
            # print 'rv_groups: %s' % str(rv_groups)

            for rv_group in rv_groups:
                # append rv_group to groups
                all_groups = copy.deepcopy(reg_groups)
                all_groups['rv'] = rv_group

                # append rv parameter to the remaining parameters
                rv_pars = self.sl.get_parameter(rv=rv_group)
                all_pars = copy.deepcopy(reg_pars)
                for c in rv_pars.keys():
                    all_pars[c].extend(rv_pars[c])

                if self.ol is not None:

                    if rv_group not in self.rel_rvgroup_region[reg]:
                        continue

                    # the wmin wmax is used to check again that
                    # we are in the correct region.
                    obs = self.ol.get_spectra(wmin=wmin, wmax=wmax, rv=rv_group)
                    # print len(obs)
                    if len(obs) == 0:
                        continue
                    # print obs, rv_group
                    # print rv_group, obs
                    # c = obs.component

                    # in case of korel spectrum we compare only one component
                    # if c != 'all':
                    #     all_pars = {c: all_pars[c]}
                else:
                    obs = [None]

                # add the comparison for each observed spectrum
                # because in an unlikely event, when we fit the
                # same RVs for several spectra
                for o in obs:
                    # setup component
                    if o is None:
                        c = 'all'
                    else:
                        c = o.component
                    if c != 'all':
                        temp_all_pars = {c: all_pars[c]}
                    else:
                        temp_all_pars = all_pars

                    self.add_comparison(region=reg,
                                        parameters=temp_all_pars,
                                        groups = all_groups,
                                        observed=o,
                                        )

    def remove_parameter(self, component, parameter, group):
        """
        :param component: component for which the parameter is deleted
        :param parameter:deleted paramer
        :return:
        """

        self.sl.remove_parameter(component, parameter, group)

    def run_fit(self, l=None, verbose=False):
        """
        Starts the fitting
        :param l:
        :param verbose:
        :return:
        """
        # this starts recording of each iteration chi2
        self.fit_is_running = True

        # runs the fitting
        self.fitter(self.compute_chi2,l,verbose)

        # writes the remaining iterations within the file
        self.fitter.flush_iters()

        # turn of the fitting
        self.fit_is_running = False

    def setup(self):
        """
        This function probes the observed and
        region list and propagates group definitions
        from them to the starlist.
        :return:
        """
        # first setup region groups
        if self.rl is not None:
            region_groups = self.rl.get_region_groups()
            self.sl.set_groups(region_groups)
        else:
            self.rl = RegionList(debug=self.debug)
            self.rl.get_regions_from_obs(copy.deepcopy(self.ol.observedSpectraList['spectrum']))

            # TODO setting up the region <-> rv relation better - this is a quick fix
            # TODO and unlikely a robust one
            self.rel_rvgroup_region = {reg:[0] for reg in self.rl._registered_regions}
            # print self.rel_rvgroup_region

            region_groups = self.rl.get_region_groups()
            self.sl.set_groups(region_groups)

        # setup radial velocity groups
        if self.ol is not None:
            self.setup_rv_groups()
        else:
            warnings.warn('There are no data attached, so all regions are set to '
                          'have the same radial velocity. Each component can have'
                          'different velocity of course.')

        # setup grids
        if not self.grid_properties_passed:
            # if the settings of the grids were not passed,
            # we go for default ones
            warnings.warn('Using default grid and synthetic spectra settings.')

            # choose the resolution of the grid to 2*R_best_observed
            if self.adaptive_resolution:
                step = self.ol.get_resolution()

            if self.debug:
                print "The step size of the grid is: %s Angstrom." % str(step)

            self.set_grid_properties(step=step/2.)

        self._setup_grids()

        # print self.grids['region00']

        # create the basic interpolated spectra
        self.ready_synthetic_spectra()

        # prepare list of comparisons
        self.ready_comparisons()

        # setup fitter
        self.fitter = Fitter(debug=self.debug)

    def set_grid_properties(self, **kwargs):
        """
        :param kwargs: parameters for the SyntheticGrid.__init__
        :param kwargs: padding - number of spectra to use for
                padding of synthetic spectra
        :param kwargs: nspectra - maximal number of spectra
                for interpolation
        :return:
        """
        if kwargs is not None:
            # setup how the grids are cretaed
            self._grid_kwargs = dict(mode = kwargs.get('mode', 'default'),
                                    debug = kwargs.get('debug', self.debug))

            # setup keyword for individual spectra
            self._synthetic_spectrum_kwargs = dict(padding = kwargs.get('padding', 20),
                                                   order = kwargs.get('order', 4),
                                                   step = kwargs.get('step', 0.01))

        self.grid_properties_passed = True

    def set_parameter(self, component='all', parname=None, group='all', **kwargs):
        """
        :param component:
        :param parname
        :param group:
        :param kwargs: keywords to be set up for each parameter
        :return:
        """

        # check the results
        if parname is None:
            print "I cannot adjust parameter: %s." % str(parname)
        if len(kwargs.keys()) == 0:
            return

        # setup the components
        if component == 'all':
            component = self.sl._registered_components
        else:
            component = [component]

        # create a list of unique groups if all are needed
        if group is 'all':
            groups = []
            dict_groups = self.sl.get_defined_groups(parameter=parname)
            for c in dict_groups.keys():
                groups.extend(dict_groups[c][parname])
            groups = np.unique(groups)
        else:
            groups = [group]

        # propagate to the star
        for c in component:
            for g in groups:
                self.sl.set_parameter(parname, c, g, **kwargs)


    def _setup_grids(self):
        """
        Initializes grid of synthetic spectra for each region -
        i.e. there is no point in calling the function without
        having the regions set up.

        :params kwargs -see pyterpol.
        :return:
        """
        for reg in self.rl.mainList.keys():
            self.grids[reg] = SyntheticGrid(**self._grid_kwargs)

    def setup_rv_groups(self):
        """
        Setting up the rv_groups is a pain..
        :return:
        """
        # TODO Can this be done better?????
        # empty array for components where cloning
        # was performed - to get rid of the first
        # group
        cloned_comps = []
        registered_groups = []

        # dictionary for newly registered groups
        # this is necessary in case we do not
        # the newly registered groups have to
        # be assigned back to the spectra
        # otherwise we would not know which rv
        # belongs to which spectrum
        new_groups = dict()

        # get wavelength boundaries of defined regions
        wmins, wmaxs, regs = self.rl.get_wavelengths(verbose=True)

        # this dictionary is needed to have
        # unambiguous relationship between
        # rv_group, spectrum and region
        reg2rv = {x: [] for x in regs}

        # for every region we have a look if we have some datas
        for wmin, wmax, reg in zip(wmins, wmaxs, regs):

            # query spectra for each region
            observed_spectra = self.ol.get_spectra(wmin=wmin, wmax=wmax)

            for i, spectrum in enumerate(observed_spectra):

                # read out properties of spectra
                component = spectrum.component
                rv_group = spectrum.group['rv']

                # readout groups that were already defined for all components
                def_groups = self.sl.get_defined_groups(component='all', parameter='rv')['all']['rv']
                # print spectrum, wmin, wmax, component, rv_group, def_groups

                # We define group for our observation
                if rv_group is None:
                    gn = generate_least_number(def_groups)
                    reg2rv[reg].append(gn)

                    # save the newly registered group
                    if spectrum.filename not in new_groups.keys():
                        new_groups[spectrum.filename] = []
                    new_groups[spectrum.filename].append(gn)

                elif rv_group not in def_groups:
                    gn = rv_group
                    reg2rv[reg].append(rv_group)

                # if the group is defined we only need to
                # add it among the user defined one, so it
                # so it is not deleted later
                elif rv_group in def_groups:
                    registered_groups.append(rv_group)
                    reg2rv[reg].append(rv_group)
                    continue

                # attachs new parameter to the StarList
                # print component, gn
                self.sl.clone_parameter(component, 'rv', group=gn)

                if component not in cloned_comps:
                    if component == 'all':
                        cloned_comps.extend(self.sl.get_components())
                    else:
                        cloned_comps.append(component)
                    registered_groups.append(gn)

        # print registered_groups, cloned_comps
        # remove the default groups
        for c in cloned_comps:
            gref = self.sl.componentList[c]['rv'][0]['group']
            if gref not in registered_groups:
                self.remove_parameter(c, 'rv', gref)

        # print new_groups
        # back register the group numbers to the observed spectra
        for filename in new_groups.keys():
            self.ol.set_spectrum(filename=filename, group={'rv': new_groups[filename]})

        # finalize the list of rv_groups for each region
        self.rel_rvgroup_region = {x: np.unique(reg2rv[x]).tolist() for x in reg2rv.keys()}

    def verify(self):
        pass

    def verify_before_fitting(self):
        pass

    def write_synthetic_spectra(self, component=None, region=None, outputname=None):
        """
        Writes the synthetic spectra obtained through the fitting.
        :param component
        :param region
        :param outputname
        :return:
        """

        # set defaults for component
        if component is None:
            components = self.sl._registered_components
        if isinstance(component, str):
            components = [component]

        # set defaults for region
        if region is None:
            regions = self.rl._registered_regions
        if isinstance(region, str):
            regions = [region]

        for r in regions:
            for c in components:
                pass


class List(object):
    """
    Future parent class for all the lists, which are dictionaries... :-)
    """
    def __init__(self, l=None, debug=None):
        """
        :param l: the list stored within the class
        :param debug: debugmode on/off
        :return:
        """
        # list
        if l is not None:
            self.mainList = l
        else:
            self.mainList = {}

            # setup debug mode
            self.debug = debug

    def clear_all(self):
        """
        Clears the list
        :return: None
        """

        self.mainList = {}


class ObservedList(object):
    """
    A helper class which groups all observed spectra and
    prepares necessary parameters for fitting.
    """
    def __init__(self, observedSpectraList=None, debug=False):
        """
        :param observedSpectraList: this should not be used in general, this creates the class
                    assuming that we are passin the self.observedSpectraList,
                    this shoudl not be used probably
        :param debug: debug mode
        :return:
        """

        # dictionary containing all observed spectra, apart from that
        # it also carries information on. A group fro radial velocities
        # has to be always set, because we intend to fit spectra acquired
        # on different times.
        self.observedSpectraList = dict(spectrum=[], group=dict(), properties=dict())
        self.groupValues = dict()

        # list of properties
        self._property_list = ['component', 'filename', 'hasErrors', 'korel', 'loaded', 'wmin', 'wmax']
        # self._queriables = copy.deepcopy(self._property_list).extend(['group'])

        # although wmin, wmax can be queried, it is treated separately from the remaining
        # parameters, because it cannot be tested on equality
        self._queriables = [x for x in self._property_list if x not in ['wmin', 'wmax']]
        self._queriable_floats = ['wmin', 'wmax']

        # initialize with empty lists
        self.observedSpectraList['properties'] = {key: [] for key in self._property_list}

        # debug
        self.debug = debug

        if observedSpectraList is not None:
            self.observedSpectraList = observedSpectraList
            self.read_groups()
            self.read_properties()
            self.groupValues = self.get_defined_groups()

    def __len__(self):
        """
        Returns number of attached observed spectra.
        """
        return len(self.observedSpectraList['spectrum'])

    def __str__(self):
        """
        String method for the class
        :return:
         string.. string representation of teh class
        """
        string = 'List of all attached spectra:\n'
        for i, spectrum in enumerate(self.observedSpectraList['spectrum']):
            string += str(spectrum)
        return string

    def add_one_observation(self, update=True, **kwargs):
        """
        Adds observation to the list.
        :param update - update the observed spectra list
        :param kwargs
            see class ObservedSpectrum (observations module) for details.
        """
        # adds the spectrum and loads it
        if self.debug:
            kwargs['debug'] = True

        obs = ObservedSpectrum(**kwargs)
        self.observedSpectraList['spectrum'].append(obs)

        if self.debug:
            print "Adding spectrum: %s" % (str(obs))

        # builds the observedSpectraList dictionary
        if update:
            self.read_groups()
            self.read_properties()
            self.groupValues = self.get_defined_groups()

    def add_observations(self, spec_list, update=True):
        """
        :param spec_list: list of dictionaries - key words are
                the same as for ObservedSpectrum class constructor
        :param update: whether to update the dictionary
                with the properties of the observed spectra
        """
        # attachs the spectra
        for rec in spec_list:
            self.add_one_observation(update=False, **rec)

        # builds the observedSpectraList dictionary
        if update:
            self.read_groups()
            self.read_properties()
            self.groupValues = self.get_defined_groups()

    def clear_all(self):
        """
        Clears all spectra.
        """
        self.__init__()

    def get_data_groups(self, components):
        """
        Returns a dictionary, containing a record
        on defined group for each component.
        :param components: a list of queried components
        :return:
        """
        groups = dict()
        for component in components:
            osl = self.get_spectra(verbose=True, component=component)

            if self.debug:
                print 'Queried observed spectra: %s for component: %s.' % (str(osl), component)

            if len(osl) > 0:
                groups[component] = ObservedList(observedSpectraList=osl).get_defined_groups()

        return groups

    def get_defined_groups(self, component=None):
        """
        Reads all groups and values that are set
        for the spectra in the list.
        :param component
        :return dictionary of defined group for all/given component
        """
        if component is 'all':
            component = None

        # empty dicitonary for the values
        groups = dict()

        # go through ech spectrum and store defined values
        for spectrum in self.observedSpectraList['spectrum']:

            # select component
            if component is not None and spectrum.component != component:
                continue

            for key in spectrum.group.keys():
                if key not in groups.keys():
                    groups[key] = []
                if isinstance(spectrum.group[key], (list, tuple)):
                    groups[key].extend(spectrum.group[key])
                else:
                    groups[key].append(spectrum.group[key])

        # only unique values are needed
        for key in groups.keys():
            groups[key] = np.unique(groups[key]).tolist()

        return groups

    def get_resolution(self, verbose=False):
        """
        Reads resoolution for each spectrum
        :param verbose
        :return:
        """
        # create a list of resolutions
        resolutions = np.zeros(len(self))
        for i in  range(0, len(self)):
            resolutions[i] = self.observedSpectraList['spectrum'][i].step

        # if verbose is set returns resolution for each spectrum
        if verbose:
            return resolutions

        # or just the maximum value
        else:
            return np.max(resolutions)

    def get_spectra(self, verbose=False, **kwargs):
        """
        :param kwargs.. properties of ObservedSpectrum,
          that we want to return. This function does not
          search the individual spectra, but the dictionary
          observedSpectraList.
        :param verbose return the whole bserved spectra list
          stub

          In general this could be - wmin, wmax, group,
          component etc..
        :return:
          speclist = all spectra that have the queried properties
        """

        # First of all check that all passed arguments are
        # either defined among queriables or is in groups
        for key in kwargs.keys():
            # print key, self._queriables
            if (key not in self._queriables) & (key not in self._queriable_floats):
                if key not in self.groupValues.keys():
                    raise KeyError('Keyword %s is not defined. This either means, that it was not set up for '
                                   'the observed spectra, or is an attribute of Observed spectrum, but is not '
                                   'defined among queriables, or is wrong.' % key)

        # create a copy of the spectralist
        osl = copy.deepcopy(self.observedSpectraList)

        # debug string
        dbg_string = 'Queried: '

        # reduce the list
        for key in kwargs.keys():

            # find all matching for a given key-word
            keytest = key.lower()

            # these can be tested on equality as strings
            if keytest in self._queriables:
                # print osl['properties'][keytest], kwargs[key]
                vind = np.where(np.array(osl['properties'][keytest], dtype=str) == str(kwargs[key]))
            elif keytest == 'component':
                vind = np.where((np.array(osl['properties'][keytest], dtype=str) == str(kwargs[key])) or
                                (np.array(osl['properties'][keytest], dtype=str) == 'all'))[0]
            # that cannot be tested on equality
            elif keytest == 'wmin':
                vind = np.where(np.array(osl['properties'][keytest]) <= kwargs[key])[0]
            elif keytest == 'wmax':
                # print osl['properties'][keytest], kwargs[key]
                vind = np.where(np.array(osl['properties'][keytest]) >= kwargs[key])[0]

            # those that are defined in groups
            elif keytest in osl['group'].keys():
                vind = []
                # print osl['spectrum'], len(osl['spectrum'])
                # print len(osl), keytest, kwargs[key], osl['group'][keytest]
                for i in range(0, len(osl['spectrum'])):
                    if isinstance(osl['group'][keytest][i], (tuple, list)):
                        if kwargs[key] in osl['group'][keytest][i]:
                            vind.append(i)
                    else:
                        if kwargs[key] == osl['group'][keytest][i]:
                            vind.append(i)
                vind = np.array(vind)

            # print keytest, vind

            # TODO Improve this, so the warning is not issued,
            # TODO when the observed spectra are lusted under all.

            if len(vind) == 0:
                warnings.warn('No spectrum matching %s: %s was found in the '
                              'list of observed spectra:\n%sDo not panic, it can '
                              'still be listed among \'all\'.' % (key, str(kwargs[key]), str(self)))
                return []

            if self.debug:
                dbg_string += '%s: %s ' % (key, str(kwargs[key]))
                print "%s.. %s spectra remain." % (dbg_string, str(len(vind)))

            # extract them from the list
            for dic in osl.keys():
                # if the key refers to a dictionary
                if isinstance(osl[dic], dict):
                    for sub_key in osl[dic].keys():
                        osl[dic][sub_key] = (np.array(osl[dic][sub_key])[vind]).tolist()

                # if it refers to a list or array
                else:
                    osl[dic] = (np.array(osl[dic])[vind]).tolist()

        # simple output, just spectra
        if not verbose:
            return osl['spectrum']
        # otherwise the whole remnant of the
        # observed spectra list is returned
        else:
            return osl

    def read_groups(self):
        """
        Updates the dictionary observedSpectraList with group
        records for every single observations and creates
        the dictionary groupValues which contains lists of
        all defined groups for every parameter.

        For parameters != 'rv':
            If at least one spectrum has a group assigned
            it is automatically assumed, that it does not
            belong among the remaining ones. This means
            that all remaining spectra are assigned their
            own group.

        For parameters == 'rv':
            Each spectrum is assigned unique RV group,
            unless this is overriden by the user by setting
            them up. This comes natural, since we are
            likely to fit spectra from different times,
            regions, where slight shifts in rv are
            very likely.
        """
        # TODO Maybe I should assign groups with respect
        # TODO to the components.

        # this stores the last group number for
        # parameter 'rv'
        gn_rv = None

        # First go through each spectrum to see, which
        # groups were defined by user
        groups = self.get_defined_groups()
        # print groups

        # check that rv has been setup - mandatory, because each observed spectrum
        # is assigned its own rv_group
        if 'rv' not in groups.keys():
            groups['rv'] = []

        # assign empty group arrays
        for key in groups.keys():
            self.observedSpectraList['group'][key] = np.zeros(len(self)).astype('int16').tolist()

        # Assigning groups to every spectrum
        for i, spectrum in enumerate(self.observedSpectraList['spectrum']):
            for key in groups.keys():

                # If not user defined the maximal possible
                # group is assigned
                if key is not 'rv':
                    gn = spectrum.get_group(key)
                    def_groups = groups[key]

                    # if spectrum has no group, but some groups have been defined,
                    # the group is assigned to the least number not in defuined groups
                    if gn is None and len(def_groups) > 0:
                        gn = 0
                        while gn in def_groups:
                            gn += 1

                    # if no group is defined for all spectra, start with zero
                    elif gn is None and len(def_groups) == 0:
                        gn = 0

                    # store the groupnumber
                    self.observedSpectraList['group'][key][i] = gn

                else:

                    gn = spectrum.get_group(key)
                    if gn is None:
                        self.observedSpectraList['group'][key][i] = None
                    else:
                        self.observedSpectraList['group'][key][i] = gn

        # propagate the groups back to spectra
        self._set_groups_to_spectra()

    def read_properties(self):
        """
        Goes through the attached spectra and reads
        stores them within the observedSpectraList
        dictionary.
        """
        # initialize with empty lists
        for key in self._property_list:
            self.observedSpectraList['properties'][key] = np.empty(len(self), dtype=object)

        # fill the dictionary
        for i, spectrum in enumerate(self.observedSpectraList['spectrum']):
            for key in self._property_list:
                self.observedSpectraList['properties'][key][i] = getattr(spectrum, key)

    def set_spectrum(self, filename=None, **kwargs):
        """
        Sets spectrum to a given value.
        :param filename
        :param kwargs:
        :return:
        """
        # print kwargs
        for i in range(0, len(self)):
            if self.observedSpectraList['spectrum'][i].filename == filename:
                for key in kwargs.keys():
                    setattr(self.observedSpectraList['spectrum'][i], key, kwargs[key])
        self.read_groups()
        self.groupValues = self.get_defined_groups()

    def _set_groups_to_spectra(self):
        """
        Propagates groups, which are set in observedSpectraList,
        in in dividual spectra.
        """
        for i in range(0, len(self.observedSpectraList['spectrum'])):
            group = {key: self.observedSpectraList['group'][key][i] for key in self.observedSpectraList['group'].keys()}
            self.observedSpectraList['spectrum'][i].set_group(group)


class RegionList(List):
    """
    """
    def __init__(self, **kwargs):
        """
        Class constructor
        :return:None
        """

        # setup the parent class
        super(RegionList, self).__init__(**kwargs)

        # registered keywords
        self._registered_records = ['components', 'groups','wmin', 'wmax']

        # if not given along the class a blank one is created
        if len(self.mainList.keys()) < 1:
            self.mainList = {}
            self._registered_regions = []
            self._user_defined_groups = {}
        else:
            self._registered_regions = self.get_registered_regions()

    def __str__(self):
        """
        String representation of the class.
        :return: string
        """

        string = ''

        # go over regions
        for key0 in self.mainList.keys():
            # region properties
            string += "Region name: %s: (wmin, wmax) = (%s, %s):\n" % (key0, str(self.mainList[key0]['wmin']),
                                                                  str(self.mainList[key0]['wmax']))
            # componentn properties
            for i in range(0, len(self.mainList[key0]['components'])):
                string += "%s: %s " % ('component', str(self.mainList[key0]['components'][i]))
                string += "%s: %s " % ('groups', str(self.mainList[key0]['groups'][i]))
                string += '\n'
        return string

    def add_region(self, component='all', identification=None, wmin=None, wmax=None, groups=None):
        """
        :param component: component for whichg the region apply
        :param identification
        :param wmin: minimal wavelength
        :param wmax: maximal wavelength
        :param groups: group numbers for this region
        :return: None
        """

        # if we are crazy and want to set this up
        # either by wavelength or by identification
        if (wmin is None or wmax is None) and identification is None:
            raise ValueError('Boundaries are not set properly: (wmin,wmax)= (%s, %s)' % (str(wmin), str(wmax)))
        else:
            if (wmin >= wmax) and identification not in self._registered_regions:
                raise ValueError('wmin is greater than wmax: %s > %s '
                                 'or the region: %s is not registered.' % (str(wmin), str(wmax), identification))

        # convert component/group/identification keys to lowercase
        if groups is not None:
            groups = keys_to_lowercase(groups)
        else:
            groups = {}
        component = component.lower()
        ident = identification
        if ident is not None:
            ident = ident.lower()

        # maybe the region has been already defined
        if ident in self.mainList.keys():
            region = ident
        elif ident is None:
            region = self.get_region(wmin, wmax)
        else:
            region = None

        # if there is a region exists and the component is all,
        # there is no point to attach it
        # print region, component
        if (region != None) and (component == 'all'):
            warnings.warn('The region already exists as region: %s -> doing nothing.' % region)
            return

        # if it is not empty
        if region is not None:

            if self.debug:
                print "Adding component: %s to region: %s" % (component, region)

            # check that the component ws not set earlier
            if self.has_component(region, component):
                warnings.warn('The component: %s is already set for region: %s. -> doing nothing.'
                            % (component, region))
                return

            # get lr from the region first record
            # print groups, self.mainList[region]['groups']
            groups['lr'] = self.mainList[region]['groups'][0]['lr']

            self.read_user_defined_groups(groups)

            # store everything apart from the wmin, wmax
            self.mainList[region]['groups'].append(groups)
            self.mainList[region]['components'].append(component)

            # readout user-defined groups
            self.read_user_defined_groups(groups)
        else:

            # print 'tu'
            # setup identification for
            if ident is None:
                ident = 'region' + str(len(self._registered_regions)).zfill(2)

            if self.debug:
                print "Creating new region: %s." % (ident)

            # register the new region
            self.mainList[ident] =  dict(wmin=wmin, wmax=wmax, components=[component], groups=[])
            self._registered_regions.append(ident)

            # if the luminosity group is not defined
            if 'lr' not in groups.keys() or len(groups['lr']) == 0:
                all_groups = self.get_defined_groups()
                if 'lr' in all_groups.keys():
                    def_groups = all_groups['lr']
                else:
                    def_groups = []
                gn = 0
                while gn in def_groups:
                    gn += 1
                groups['lr'] = gn

            # add groups to the list
            self.mainList[ident]['groups'].append(groups)
            # readout user-defined groups
            self.read_user_defined_groups(groups)

        self.setup_undefined_groups()

    def clear_all(self):
        """
        Clears the class.
        :return:
        """

        super(RegionList, self).clear_all()
        self._registered_regions = []
        self._user_defined_groups = {}

    def get_defined_groups(self):
        """
        Returns plain list of all defined groups regardless of their components.
        :return: list of defined groups
        """
        groups = {}
        for reg in self._registered_regions:
            for rec in self.mainList[reg]['groups']:
                for key in rec.keys():
                    if key not in groups.keys():
                        groups[key] = [rec[key]]
                    else:
                        if rec[key] not in groups[key]:
                            groups[key].append(rec[key])

        return groups

    def get_region(self, wmin, wmax):
        """
        Checks that a region with this wavelength range
        does not exist.
        :param wmin
        :param wmax
        :return:
        """

        for region in self.mainList:
            if (abs(self.mainList[region]['wmin'] - wmin) < ZERO_TOLERANCE) & \
               (abs(self.mainList[region]['wmax'] - wmax) < ZERO_TOLERANCE):
                return region
        return None

    def get_region_groups(self):
        """
        A dictionary of groups defined for regions component by component.
        :return: dictionary containing records on groups
                which can be directly passed to type StarList
                through set_groups
        """
        groups = {}

        # go over each region
        for reg in self.mainList.keys():
            for i in range(0, len(self.mainList[reg]['components'])):
                component = self.mainList[reg]['components'][i]
                comp_groups = self.mainList[reg]['groups'][i]

                # setup component
                if component not in groups.keys():
                    groups[component] = {}

                # setup keys
                for key in comp_groups.keys():
                    if key not in groups[component].keys():
                        groups[component][key] = [comp_groups[key]]
                    else:
                        if comp_groups[key] not in groups[component][key]:
                            groups[component][key].append(comp_groups[key])
        return groups



    def get_registered_regions(self):
        """
        Returns an array of registered regions.
        :return:
        """
        return self.mainList.keys()

    def get_wavelengths(self, verbose=False):
        """
        Returns registered wavelengths
        :param verbose
        :return: wmins, wmaxs = arrays of minimal/maximal wavelength for each region
        """
        wmins = []
        wmaxs = []
        regs = []
        for reg in self.mainList.keys():
            wmins.append(self.mainList[reg]['wmin'])
            wmaxs.append(self.mainList[reg]['wmax'])
            regs.append(reg)

        if verbose:
            return wmins, wmaxs, regs
        else:
            return wmins, wmaxs

    def get_regions_from_obs(self, ol, append=False):
        """
        Reads the region from a list of observations. In general this
        function should not be used for fitting, because it
        makes no sense to fit the whole spectrum.

        :param ol: list of ObservedSpectrum
        :param append are we appending to existing list?
        :return: list of unique limits
        """
        if len(ol) == 0:
            raise ValueError('Cannot setup regions from observed spectra, because'
                             ' their list is empty!')

        # clear the regions if needed
        if not append:
            self.clear_all()

        # empty arrays for limits
        limits = {}
        # the rounding is there get over stupid problems with float precision
        for obs in ol:
            component = obs.component
            if component not in limits:
                limits[component] = [[], []]

            limits[component][0].append(np.ceil(obs.wmin))
            limits[component][1].append(np.floor(obs.wmax))

            # get only unique values
            for i in range(0,2):
                limits[component][i] = np.unique(limits[component][i])

        # check that something funny did not happen
        for component in limits.keys():
            if len(limits[component][0]) != len(limits[component][1]):
                raise ValueError('The limits were not read out correctly from observed spectra.')

            # setup the regions
            for i in range(0, len(limits[component][0])):
                self.add_region(component=component,
                                wmin=limits[component][0][i],
                                wmax=limits[component][1][i])

        return limits

    def has_component(self, region, component):
        """
        Checks that certain component was attached for a given
        region.
        :param region:
        :param component:
        :return: bool has/has_not the component
        """

        for regcomp  in self.mainList[region]['components']:
            if (regcomp == component) or (regcomp == 'all'):
                return True
        return False

    def read_user_defined_groups(self, groups):
        """
        When adding new region, all user defined groups
        are read out to properly set the default groups
        :param groups groups to be read
        :return: None
        """
        for key in groups.keys():
            if key not in self._user_defined_groups.keys():
                self._user_defined_groups[key] = [groups[key]]
            else:
                if groups[key] not in self._user_defined_groups[key]:
                    self._user_defined_groups[key].append(groups[key])

    def setup_undefined_groups(self):
        """
        User can be a bit lazy. If we split some parameter
        into more groups, we can only set group for few
        and the remaining dataset gets a default one.

        This nonetheless has to be run after all
        regions were attached. If we do this
        earlier, we will get into serious problems.
        :return:
        """
        # defined groups
        groups = self.get_defined_groups()

        # setup default group numbers for region->component
        # with unset group
        for region in self._registered_regions:
            for i, comp_group in enumerate(self.mainList[region]['groups']):

                # go over each defined group
                for key in groups.keys():
                    # if the key is unset for the component
                    # we have to assign some. This must
                    # not be one of the user-defined.
                    # That is why we maintain dictionary
                    # of user defined groups.
                    if key not in comp_group.keys():
                        gn = 0
                        while gn in self._user_defined_groups[key]:
                            gn += 1
                        self.mainList[region]['groups'][i][key] = gn


class StarList(object):
    """
    """
    def __init__(self, debug=False):
        """
        """

        # set up debug mode
        self.debug = debug

        # define empty list of components
        self.componentList = {}

        # array storing registered components
        self._registered_components = []

        # defined groups
        self.groups = {}

    def __len__(self):
        """
        Returns number of parameters.
        :return: l
        """
        pass


    def __str__(self):
        """
        :return: string = string represantation of the class
        """
        string = ''
        for component in self.componentList.keys():
            string += "Component: %s\n" % (component)
            for parkey in self.componentList[component].keys():
                for par in self.componentList[component][parkey]:
                    string += str(par)

        return string

    def add_component(self, component=None, groups={}, **kwargs):
        """
        Setups a component - if no kwargs are given,
        all parameters from the parameter_definitions
        are taken.

        If one wants to not-include a parameter,
        params = None, has to be passed. If one
        wants to add a parameter, that is not
        defined in parameter definitions, just
        pass parameter + value.

        :param component: Registration string of the component
                if None is given, it is registred as 'componentXX'
        :param group: group set to all parameters of a component
        :param kwargs:
        :return:
        """

        # setup name of the component and create a record within
        # component list
        if component is None:
            component = 'component'+str(len(self._registered_components))

        # register he component
        self._registered_components.append(component)

        # the parameters will be stored in a dictionary
        self.componentList[component] = dict()
        pd = copy.deepcopy(parameter_definitions)

        # setup groups for default parameters
        for key in groups.keys():
            if key in pd.keys():
                pd[key]['group'] = groups[key]

        # process the keyword-arguments
        for key in kwargs.keys():
            keytest = key.lower()
            # if we pass par + value, it is just stored
            if keytest in pd.keys() and kwargs[key] is not None:
                self.componentList[component][keytest] = []
                self.componentList[component][keytest].append(Parameter(**pd[key]))
                self.componentList[component][keytest][-1]['value'] = kwargs[key]
            elif kwargs[key] is None:
                warnings.warn('The parameter %s is set to %s. Therefore it is not '
                              'included into component parameters.' % (key, str(kwargs[key])))
            elif keytest not in pd.keys() and kwargs[key] is not None:

                #set up group
                if keytest in groups.keys():
                    group = groups[keytest]
                self.componentList[component][keytest] = []
                self.componentList[component][keytest].append(Parameter(name=key, value=kwargs[key], group=group))
                self.componentList[component][keytest][-1].set_empty()
                warnings.warn('The parameter %s: %s is not set among the '
                              'parameter definitions. Therefore you should pay '
                              'attention to ist settings.')

        # pass all unset parameters in definitions
        for key in pd.keys():
            if key not in self.componentList[component].keys():
                self.componentList[component][key] = []
                self.componentList[component][key].append(Parameter(**pd[key]))

        # readout the groups
        self.read_groups()
        self.get_fitted_types()

    def add_parameter_to_component(self, component, p=None, **kwargs):
        """
        Adds a parameter to a specific component.
        :param component: component for which we want to add a parameter
        :param p: assigning directly the Parameter type
        :param kwargs: see Parameter class for description
        :return:
        """
        if p is None:
            self.componentList[component][kwargs['name']] = []
            self.componentList[component][kwargs['name']].append(Parameter(**kwargs))
        else:
            # print p['name']
            self.componentList[component][p['name']].append(copy.deepcopy(p))

        # redefine groups
        self.read_groups()
        self.get_fitted_types()

    def add_parameter_to_all(self, **kwargs):
        """
        Adds a parameter to all components
        :param kwargs: see Parameter class
        :return: None
        """
        for component in self._registered_components:
            self.add_parameter_to_component(component, **kwargs)

    def clear(self):
        """
        Clears the component list
        :return: None
        """
        self.componentList = {}
        self._registered_components = []

    # def clone_parameter(self, component, parameter, index=0, all=False, **kwargs):
    #     """
    #     Clones a parameter and stores it for a given component.
    #     This function will be primarily used to clone parameters
    #     to acount for different groups.
    #
    #     :param component: component for which we want to clone the parameter
    #     :param parameter: the cloned parameter
    #     :param index : the specific cloned parameter
    #     :param kwargs: values we want to change for the parameter
    #     :return: clone type_Parameter - the cloned parameter
    #     """
    #     # in case we pass
    #     if component.lower() == 'all':
    #         all = True
    #         component = self._registered_components[0]
    #
    #     # copy the parameter
    #     clone = copy.deepcopy(self.componentList[component][parameter][index])
    #
    #     # adjust its values
    #     for key in kwargs.keys():
    #         keytest = key.lower()
    #         clone[keytest] = kwargs[key]
    #
    #     # append the new component to the componentlist
    #     if all:
    #         self.add_parameter_to_all(p=clone)
    #     else:
    #         self.add_parameter_to_component(component, p=clone)
    #
    #     return clone

    def clone_parameter(self, component, parameter, index=0, all=False, **kwargs):
        """
        Clones a parameter and stores it for a given component.
        This function will be primarily used to clone parameters
        to acount for different groups.

        :param component: component for which we want to clone the parameter
        :param parameter: the cloned parameter
        :param index : the specific cloned parameter
        :param kwargs: values we want to change for the parameter
        :return: clone type_Parameter - the cloned parameter
        """
        # in case we pass
        if component.lower() == 'all':
            components = self._registered_components
        else:
            components = [component]

        clones = []
        # go over each component
        for component in components:

            # copy the parameter
            clone = copy.deepcopy(self.componentList[component][parameter][index])
            clones.append(clone)

            # adjust its values
            for key in kwargs.keys():
                keytest = key.lower()
                clone[keytest] = kwargs[key]

            # append the new component to the componentlist
            self.add_parameter_to_component(component, p=clone)

        return clones

    def delete_hollow_groups(self):
        """
        Goes through parameters and deletes those that
        are set to None.
        :return: None
        """

        for component in self._registered_components:
            for parkey in self.componentList[component].keys():
                i = 0
                while(i < len(self.componentList[component][parkey])):

                    # if the parameter group is not, it is deleted
                    if self.componentList[component][parkey][i]['group'] is None:
                        del self.componentList[component][parkey][i]
                    else:
                        i+=1

    def delete_duplicities(self):
        """
        Delete duplicities in groups.
        :return: None
        """
        for component in self._registered_components:
            # groups can a have to be the same for two components ofc,
            def_groups = []
            for parkey in self.componentList[component].keys():
                i = 0
                while (i < len(self.componentList[component][parkey])):
                    if self.componentList[component][parkey][i]['group'] not in def_groups:
                        def_groups.append(self.componentList[component][parkey][i]['group'])
                        i += 1
                    # if the parameter with the group has been already defined, delete it
                    else:
                        del self.componentList[component][parkey][i]

    def get_common_groups(self):
        """
        Returns a dictionary of groups shared by all components.
        :return: com_groups
        """
        # get the keys of physical parameters
        parkeys = self.get_physical_parameters()

        # get the groups
        com_groups = {}
        for key in parkeys:
            com_groups[key] = []

            # define teh reference component
            comp0 = self._registered_components[0]

            # groups are always common for one parameter
            if len(self._registered_components) < 2:
                is_common = True

            # go over each group of
            for i in range(0, len(self.componentList[comp0][key])):
                refpar = self.componentList[comp0][key][i]
                # print refpar

                # at the beginning
                for component in self._registered_components[1:]:
                    is_common = False
                    for j, par in enumerate(self.componentList[component][key]):
                        # print par
                        if refpar['group'] == par['group']:
                            is_common = True
                            break
                    if not is_common:
                        break
                if is_common:
                    com_groups[key].append(refpar['group'])

        return com_groups

    def get_components(self):
        """
        Returns list of all defined components.
        :return:
        """
        return copy.deepcopy(self._registered_components)

    def get_defined_groups(self, component=None, parameter=None):
        """:
        :param component: starlist component
        :param parameter: physical parameter
        :return: dictionary of groups
        """
        groups = {}

        # setup parameters
        if parameter is None:
            parameters = self.get_physical_parameters()
        else:
            parameters = [parameter]

        # setup components
        if component is None or component == 'all':
            components = self.get_components()
        else:
            components = [component]

        # go over the registered componentss
        for comp in components:
            groups[comp]= {}

            # go over passed parameters
            for param in parameters:
                groups[comp][param] = []
                for regparam in self.componentList[comp][param]:
                    if regparam.name == param:
                        groups[comp][param].append(regparam.group)

        # merge groups if component was 'all'
        if component == 'all':
            for p in parameters:
                groups[component]={}
                temp = []
                for c in components:
                    temp.extend(groups[c][p])
                groups[component][p] =  np.unique(temp).tolist()

        return groups

    def get_fitted_parameters(self):
        """
        Returns a list of fitted parameters wrapped within the Parameter class ofc.
        :return:
        """
        fit_pars = []
        # go over all parameters and components
        for c in self._registered_components:
            for parname in self.get_physical_parameters():
                for par in self.componentList[c][parname]:
                    if par['fitted']:
                        fit_pars.append(par)

        return fit_pars

    def get_fitted_types(self):
        """
        Stores a dictionary of fitted types for
        each component in the class. This should
        be updated whenever a parameter is changed.
        :return:
        """

        fitted_types = {}

        # go over each component
        for c in self.componentList.keys():
            fitted_types[c] = []

            # go over each parameter type
            for parname in self.componentList[c]:
                # print c, parname
                # and finaly over each parameter
                for par in self.componentList[c][parname]:
                    # print par['fitted'], parname
                    if parname not in fitted_types[c]:
                        if par['fitted']:
                            fitted_types[c].append(parname)
                    else:
                        break

        # print fitted_types
        self.fitted_types = fitted_types


    def get_index(self, component, parameter, group):
        """
        Returns index of a component/parameter/group.
        :param component:
        :param parameter:
        :param group:
        :return:
        """

        for i, par in enumerate(self.componentList[component][parameter]):
            if par['group'] == group:
                return i

        warnings.warn('Component: %s Parameter: %s Group: :s'
                      ' not found.' % (component, parameter, group))
        return None

    def get_parameter(self, **kwargs):
        """
        Returns all parameters, which have certain group.
        :param kwargs:
        :return:
        """
        pars = {x: [] for x in self._registered_components}
        for key in kwargs.keys():
            for c in self._registered_components:
                for i, par in enumerate(self.componentList[c][key]):
                    # print i, par
                    if par.group == kwargs[key]:
                        pars[c].append(self.componentList[c][key][i])

        return pars



    def get_physical_parameters(self):
        """
        Reads physical parameters from the starlist.
        :return:
        """
        component = self._registered_components[0]
        return self.componentList[component].keys()

    def read_groups(self):
        """
        Reads all groups from the defined components. This
        is then compared to the list obtained from observations
        and defined regions,
        :return:
        """

        for component in self.componentList.keys():
            self.groups[component] = dict()
            for key in self.componentList[component].keys():
                self.groups[component][key]=[]
                for par in self.componentList[component][key]:
                    self.groups[component][key].append(par['group'])

    def remove_parameter(self, component, parameter, group):
        """
        :param component: component for which the parameter is deleted
        :param parameter:deleted paramer
        :return:
        """
        index = self.get_index(component, parameter, group)
        del self.componentList[component][parameter][index]


    def set_groups(self, groups, overwrite=False):
        """
        Sets up groups - this function is designed to
        use output from ObservedList.get_groups().
        It is assumed that the structure is following:
        dict(component01=dict(par1=[], par2=[]), component2=..)

        This function should be used to primarily
        used to assign rv_groups, where cloning
        is necessary to not to get crazy.

        This function merges groups defined
        in the type and the one passed. In general
        we should not be able to do this.

        :return: None
        """

        for component in groups.keys():
            for parkey in groups[component].keys():

                # bool variable for case, when we want to completely overwrite
                # previous settings
                first_in_list = True

                for group in groups[component][parkey]:
                    # setting group for all components
                    if component.lower() == 'all':
                        for one_comp in self._registered_components:
                            # print one_comp, parkey, self.groups
                            if group not in self.groups[one_comp][parkey]:
                                warnings.warn("Group %s: %s previously undefined."
                                              "Adding to the remaining groups." % (parkey, str(group)))
                                # print one_comp, parkey, group
                                self.clone_parameter(one_comp, parkey, group=group)

                                # deletes all previous groups
                                if overwrite and first_in_list:
                                    while len(self.groups[one_comp][parkey]) > 1:
                                        del self.groups[one_comp][parkey][0]
                                    first_in_list = False

                    # if we are setting group only for one component
                    else:
                        if group not in self.groups[component][parkey]:
                            warnings.warn("Group %s: %s previously undefined."
                                                  "Adding to the remaining groups." % (parkey, str(group)))
                            self.clone_parameter(component, parkey, group=group)

                            # deletes all previous groups
                            if overwrite and first_in_list:
                                while len(self.groups[one_comp][parkey]) > 1:
                                    del self.groups[one_comp][parkey][0]
                                first_in_list = False


    def set_parameter(self, name, component, group, **kwargs):
        """
        Sets values defined in kwargs for a parameter
        of a given component and group.
        :param name:
        :param component:
        :param group:
        :param kwargs
        :return:
        """
        name = name.lower()
        if name not in self.get_physical_parameters():
            raise Exception("Parameter: %s unknown." % name)
        elif component not in self._registered_components:
            # print self._registered_components, component
            raise Exception("Component: %s unknown" % component)
        else:
            for i, par in enumerate(self.componentList[component][name]):
                if par['name']== name and par['group'] == group:
                    for key in kwargs.keys():
                        keytest = key.lower()
                        self.componentList[component][name][i][keytest] = kwargs[key]
        # print self
        # update the list of fitted types
        self.get_fitted_types()


class SyntheticList(List):
    """
    List of resulting synthetic spectra.
    """
    def __init__(self, **kwargs):

        # initialize the parent
        super(SyntheticList, self).__init__(**kwargs)













