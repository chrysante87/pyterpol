# -*- coding: utf-8 -*-
import copy
import warnings
import numpy as np
from pyterpol.observed.observations import ObservedSpectrum
from pyterpol.fitting.parameter import Parameter
from pyterpol.fitting.parameter import parameter_definitions

# TODO Go through similarities in the classes and write a parent class
# TODO to get rid of the redundant code.

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
        self._property_list = ['loaded', 'hasErrors', 'component', 'wmin', 'wmax', 'korel']
        # self._queriables = copy.deepcopy(self._property_list).extend(['group'])

        # although wmin, wmax can be queried, it is treated separately from the remaining
        # parameters, because it cannot be tested on equality
        self._queriables = [x for x in self._property_list if x not in ['wmin', 'wmax']];
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

    def get_groups_for_components(self, components):
        """
        Returns a dictionary, containing a record
        on defined group for each component.
        :param components: a list of queried components
        :return:
        """
        groups = dict()
        for component in components:
            osl = self.get_spectra(verbose=True, component=component)
            if len(osl) > 0:
                groups[component] = ObservedList(observedSpectraList=osl).get_defined_groups()

        return groups

    def get_defined_groups(self):
        """
        Reads all groups and values that are set
        for the spectra in the list.
        param:
        OUTPUT:
            groups    dictionary containing all (so far)
                      defined groups (parameters+values).
        """
        # empty dicitonary for the values
        groups = dict(rv=[])

        # go through ech spectrum and store defined values
        for spectrum in self.observedSpectraList['spectrum']:
            # print spectrum
            for key in spectrum.group.keys():
                if key not in groups.keys():
                    groups[key] = []
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
                vind = np.where(osl['group'][keytest] == kwargs[key])[0]

            # print keytest, vind

            if len(vind) == 0:
                warnings.warn('No spectrum matching %s: %s was found in the '
                              'list of observed spectra:\n%s' % (key, str(kwargs[key]), str(self)))
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
        # this stores the last group number for
        # parameter 'rv'
        gn_rv = None

        # First go through each spectrum to see, which
        # groups were defined by user
        groups = self.get_defined_groups()

        # assign empty group arrays
        for key in groups.keys():
            self.observedSpectraList['group'][key] = np.zeros(len(self)).astype('int16')

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
                    def_groups = groups[key]

                    # if spectrum has no group, but some have been defined, set the group = max(number)+1
                    if gn is None and len(def_groups) > 0:
                        if gn_rv is None:
                            gn = 0
                        else:
                            gn = gn_rv+1
                        while gn in def_groups:
                            gn+=1

                    # if no group is defined for all spectra, start with zero
                    elif gn is None and len(def_groups) == 0:
                        if gn_rv is None:
                            gn = 0
                        else:
                            gn = gn_rv + 1
                    # print key, gn, def_groups
                    self.observedSpectraList['group'][key][i] = gn
                    gn_rv = gn

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

    def _set_groups_to_spectra(self):
        """
        Propagates groups, which are set in observedSpectraList,
        in in dividual spectra.
        """
        for i in range(0, len(self.observedSpectraList['spectrum'])):
            group = {key: self.observedSpectraList['group'][key][i] for key in self.observedSpectraList['group'].keys()}
            self.observedSpectraList['spectrum'][i].set_group(group)


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
        self._registered_records = ['wmin', 'wmax', 'group']

        #
        if len(self.mainList.keys()) < 1:
            self.mainList = {'all':{'lr':[]}}

    def __str__(self):
        """
        String representation of the class.
        :return: string
        """

        string =''
        for key0 in self.mainList.keys():
            string += "%s:\n" % (key0)
            for key1 in self.mainList[key0].keys():
                for rec in self.mainList[key0][key1]:
                    string += "%s\n" % str(rec)

        return string


    def add_region(self, wmin=None, wmax=None, group=None):
        """
        :param wmin: minimal wavelength
        :param wmax: maximal wavelength
        :param group: group number for this region
        :return: None
        """

        # if we are crazy and want to set this up
        if wmin is None or wmax is None:
            warnings.warn('One of the region boundaries is set to None. '
                          'dno not forget to set this up later!')
        else:
            if wmin >= wmax:
                raise ValueError('wmin is greater than wmax: %s > %s.' % (str(wmin), str(wmax)))

        # automatic assignment of a group - EACH AUTOMATICALLY
        # ASSIGNED GROUP IS NEW
        if group is None:
            def_groups = self.get_defined_groups()
            if len(def_groups) == 0:
                group = 0
            else:
                group = 0
                while group in def_groups:
                    group += 1

        # append the region
        self.mainList['all']['lr'].append(dict(wmin=wmin, wmax=wmax, group=group))

    def clear_all(self):
        """
        Clears the class.
        :return:
        """

        super(RegionList, self).clear_all()
        self.mainList = {'all':{'lr':[]}}

    def get_defined_groups(self):
        """
        Returns plain list of all defined groups.
        :return: list of defined groups
        """
        groups = []
        for rec in self.mainList['all']['lr']:
            if rec['group'] not in groups:
                groups.append(rec['group'])

        return groups

    def get_region_groups(self):
        """
        A dictionary of groups defined for regions.
        :return: dictionary containing records on groups
                which can be directly passed to type StarList
                through set_groups
        """
        groups = self.get_defined_groups()
        return dict(all=dict(lr=groups))


    def get_regions_from_obs(self, ol, append=False):
        """
        Reads the region from a list of observations. In general this
        function should not be used for fitting, because it
        makes no sense to fit the whole spectrum.

        :param ol: list of ObservedSpectrum
        :param append are we appending to existing list?
        :return: list of unique limits
        """

        # empty arrays for limits
        limits = [[], []]

        # the rounding is there get over stupid problems with float precision
        for obs in ol:
            limits[0].append(np.ceil(obs.wmin))
            limits[1].append(np.floor(obs.wmax))

        # get only unique values
        for i in range(0,2):
            limits[i] = np.unique(limits[i])

        # check that something funny did not happen
        if len(limits[0]) != len(limits[1]):
            raise ValueError('The limits were not read out correctly from observed spectra.s')

        # clear the regions
        if not append:
            self.clear_all()

        # setup the regions
        for i in range(0, len(limits[0])):
            self.add_region(wmin=limits[0][i], wmax=limits[1][i], group=i)

        return limits
















