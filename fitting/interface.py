# -*- coding: utf-8 -*-
import copy
import warnings
import numpy as np
from pyterpol.observed.observations import ObservedSpectrum
from pyterpol.fitting.parameter import Parameter
from pyterpol.fitting.parameter import parameter_definitions
from pyterpol.synthetic.auxiliary import keys_to_lowercase
from pyterpol.synthetic.auxiliary import generate_least_number

# TODO Go through similarities in the classes and write a parent class
# TODO to get rid of the redundant code.

# tolerance on float comparison
floatToler = 1e-6

# repeat userwarnings
warnings.simplefilter('always', UserWarning)

class Interface(object):
    """
    """
    def __init__(self, sl=None, rl=None, ol=None, debug=False):
        """
        :param sl: StarList type
        :param rl: RegionList type
        :param ol: ObservedList type
        :return:
        """

        self.sl = sl
        self.rl = rl
        self.ol = ol

        # debug mode
        self.debug = debug

    def __str__(self):
        """
        String representation of the class
        :return:
        """
        string = ""
        for attr in ['sl', 'rl', 'ol']:
            string += str(getattr(self, attr))
        return string

    def clear_all(self):
        """
        Clears the class.
        :return:
        """

        self.sl = None
        self.rl = None
        self.ol = None

    def setup_groups(self):
        """
        This function probes the observed and
        region list and propagates group definitions
        from them to the starlist.
        :return:
        """

        # # get registered components
        # components = copy.deepcopy(self.sl._registered_components)
        # components.append('all')
        #
        # # read the groups from observed data
        # # and region definitions
        # if self.ol is not None:
        #     groups_data = self.ol.get_data_groups(components)
        # if self.rl is not None:
        #     groups_regs = self.rl.get_region_groups()
        #
        # if self.debug:
        #     print "Reading groups: %s from data." % str(groups_data)
        #     print "Reading groups: %s from regions." % str(groups_regs)
        #
        # for groups in [groups_data, groups_regs]:
        #     self.sl.set_groups(groups)

        # first setup region groups
        if self.rl is not None:
            region_groups = self.rl.get_region_groups()
            self.sl.set_groups(region_groups)
        else:
            raise ValueError('Cannot setup groups without the RegionList attached to the interface.')

        # setup radial velocity groups
        if self.ol is not None:
            self.setup_rv_groups()
        else:
            warnings.warn('There are no data attached, so all regions are set to '
                          'have the same radial velocity. Each component can have'
                          'different velocity of course.')

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
        wmins, wmaxs = self.rl.get_wavelengths()

        # for every region we have a look if we have some data
        for wmin, wmax in zip(wmins, wmaxs):

            # query spectra for each region
            observed_spectra = self.ol.get_spectra(wmin=wmin, wmax=wmax)

            for i, spectrum in enumerate(observed_spectra):

                # read out properties of spectra
                component = spectrum.component
                rv_group = spectrum.group['rv']

                # readout groups that were already defined
                def_groups = self.sl.get_defined_groups(component=component, parameter='rv')[component]['rv']
                # print spectrum, wmin, wmax, component, rv_group, def_groups

                # We define group for our observation
                if rv_group is None:
                    gn = generate_least_number(def_groups)

                    # save the newly registered group
                    if spectrum.filename not in new_groups.keys():
                        new_groups[spectrum.filename] = []
                    new_groups[spectrum.filename].append(gn)

                elif rv_group not in def_groups:
                    gn = rv_group

                # if the group is defined we only need to
                # add it among the user defined one, so it
                # so it is not deleted later
                elif rv_group in def_groups:
                    registered_groups.append(rv_group)
                    continue

                # the group number is assigned to the spectrum
                # NOT A GOOD IDEA
                # self.ol.observedSpectraList['spectrum'][i].group['rv'] = gn
                # self.ol.observedSpectraList['group']['rv'][i] = gn

                # attachs new parameter to the StarList
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

        print new_groups
        # back register the group numbers to the observed spectra
        for filename in new_groups.keys():
            self.ol.set_spectrum(filename=filename, group={'rv': new_groups[filename]})



    def get_combinations(self):
        """
        This function creates a dictionary, which is one of the
        cornerstones of the class. It creates a list of all
        combinations of the parameters.
        :return:
        """

        # empty dictionary for the combinations
        combos = dict(synthetic=[], observed=[], parameters=[], groups=[])

        # walk over the dictionaries to get the combinations
        common_groups = self.sl.get_common_groups()

    def remove_parameter(self, component, parameter, group):
        """
        :param component: component for which the parameter is deleted
        :param parameter:deleted paramer
        :return:
        """

        self.sl.remove_parameter(component, parameter, group)

    def verify(self):
        pass

    def verify_before_fitting(self):
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
        :param osl: this should not be used in general, this creates the class
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
                vind = np.where(np.array(osl['properties'][keytest], dtype=str) == str(kwargs[key]))[0]

            # that cannot be tested on equality
            elif keytest == 'wmin':
                vind = np.where(np.array(osl['properties'][keytest]) <= kwargs[key])[0]
            elif keytest == 'wmax':
                # print osl['properties'][keytest], kwargs[key]
                vind = np.where(np.array(osl['properties'][keytest]) >= kwargs[key])[0]

            # those that are defined in groups
            elif keytest in osl['group'].keys():
                # vind = []
                vind = np.where(osl['group'][keytest] == kwargs[key])[0]

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
                            gn+=1

                    # if no group is defined for all spectra, start with zero
                    elif gn is None and len(def_groups) == 0:
                        gn = 0

                    # store the groupnumber
                    self.observedSpectraList['group'][key][i] = gn

                else:

                    gn = spectrum.get_group(key)
                    # def_groups = groups[key]

                    # try to set groups component-wise
                    # component = spectrum.component
                    # comp_groups = self.get_defined_groups(component=component)
                    # # print comp_groups
                    # if 'rv' not in comp_groups.keys():
                    #     comp_groups['rv'] = []
                    # def_groups = comp_groups[key]
                    # # print def_groups
                    #
                    # # if spectrum has no group
                    # if gn is None and len(def_groups) > 0:
                    #     if gn_rv is None:
                    #         gn = 0
                    #     else:
                    #         gn = gn_rv+1
                    #     while gn in def_groups:
                    #         gn+=1
                    #
                    # # if no group is defined for all spectra, start with zero
                    # elif gn is None and len(def_groups) == 0:
                    #     if gn_rv is None:
                    #         gn = 0
                    #     else:
                    #         gn = gn_rv + 1
                    # print key, gn, def_groups
                    if gn is None:
                        self.observedSpectraList['group'][key][i] = None
                    else:
                        self.observedSpectraList['group'][key][i] = gn

                    # gn_rv = gn

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
        :param kwargs:
        :return:
        """
        print kwargs
        for i in range(0, len(self)):
            if self.observedSpectraList['spectrum'][i].filename == filename:
                for key in kwargs.keys():
                    setattr(self.observedSpectraList['spectrum'][i], key, kwargs[key])
        self.read_groups()

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

        string =''

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
            if (abs(self.mainList[region]['wmin'] - wmin) < floatToler) & \
               (abs(self.mainList[region]['wmax'] - wmax) < floatToler):
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

    def get_wavelengths(self):
        """
        Returns registered wavelengths
        :return: wmins, wmaxs = arrays of minimal/maximal wavelength for each region
        """
        wmins = []
        wmaxs = []
        for reg in self.mainList.keys():
            wmins.append(self.mainList[reg]['wmin'])
            wmaxs.append(self.mainList[reg]['wmax'])

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

    def add_parameter_to_component(self, component, p=None, **kwargs):
        """
        Adds a parameter to a specific component.
        :param component: component for which we want to add a parameter
        :param p: assigning directly the Parameter type
        :param kwargs: see Parameter class for description
        :return:
        """
        if p is None:
            self.componentList[component][kwargs['name']]= []
            self.componentList[component][kwargs['name']].append(Parameter(**kwargs))
        else:
            # print p['name']
            self.componentList[component][p['name']].append(copy.deepcopy(p))

        # redefine groups
        self.read_groups()

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
            all = True
            component = self._registered_components[0]

        # copy the parameter
        clone = copy.deepcopy(self.componentList[component][parameter][index])

        # adjust its values
        for key in kwargs.keys():
            keytest = key.lower()
            clone[keytest] = kwargs[key]

        # append the new component to the componentlist
        if all:
            self.add_parameter_to_all(p=clone)
        else:
            self.add_parameter_to_component(component,p=clone)

        return clone

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

    def get_index(self, component, parameter, group):
        """
        Returns index of a component/parameter/group.
        :param component:
        :param parameter:
        :param group:
        :return:
        """

        for par in self.componentList[component][parameter]:
            if par['group'] == group:
                return par['group']

        warnings.warn('Component: %s Parameter: %s Group: :s'
                      ' not found.' % (component, parameter, group))
        return None


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
                first_in_list=True

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


class SyntheticList(List):
    """
    List of resulting synthetic spectra.
    """
    def __init__(self, **kwargs):

        # initialize the parent
        super(SyntheticList, self).__init__(**kwargs)













