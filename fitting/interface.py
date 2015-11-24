# -*- coding: utf-8 -*-
import copy
import warnings
import numpy as np
from pyterpol.observed.observations import ObservedSpectrum


class ObservedList(object):
    """
    A helper class which groups all observed spectra and
    prepares necessary parameters for fitting.
    """
    def __init__(self, debug=False):
        """
        Setups the class
        """
        # dictionary containing all observed spectra, apart from that
        # it also carries information on. A group fro radial velocities
        # has to be always set, because we intend to fit spectra acquired
        # on different times.
        self.observedSpectraList = dict(spectrum=[], group=dict(), properties=dict())
        self.groupValues = dict()

        # list of properties
        self._property_list = ['loaded', 'hasErrors', 'component', 'wmin', 'wmax']
        self._queriables = self._property_list.copy().extend(['group'])

        # initialize with empty lists
        self.observedSpectraList['properties'] = {key:[] for key in self._property_list}

        # debug
        self.debug = debug

    def __len__(self):
        """
        :return:
            l = number of the spectra in the list
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
        INPUT:
            see class ObservedSpectrum (observations module) for details.
        """
        # adds the spectrum and loads it
        if self.debug:
            kwargs['debug']=True

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
            self.add_one_observation(update=False,**rec)

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

    def get_defined_groups(self):
        """
        Reads all groups and values that are set
        for the spectra in the list.
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
        :param kwargs:.. properties of ObservedSpectrum,
          that we want to return. This function does not
          search the individual spectra, but the dictionary
          observedSpectraList.

          In general this could be - wmin, wmax, group,
          component etc..
        :return:
          speclist = all spectra that have the queried properties
        """

        # First of all check that all passed arguments are
        # either defined among queriables or is in groups
        for key in kwargs.keys():
            if key not in self._queriables:
                if key not in self.groupValues.keys():
                    raise KeyError('Keyword %s is not defined. This either means, that it was not set up for '
                                   'the observed spectra, or is an attribute of Observed spectrum, but is not '
                                   'defined among queriables, or is wrong.' % (key))

        # create a copy of the spectralist
        osl = copy.deepcopy(self.observedSpectraList)

        # debug string
        dbg_string = 'Queried: '

        # reduce the list
        for key in kwargs.keys():
            # find all matching for a given key-word
            if key in self._queriables:
                vind = np.where(osl['properties'][key] == kwargs[key])[0]
            else:
                vind = np.where(osl['properties'][key] == kwargs[key])[0]

            if len(vind) == 0:
                warnings.warn('No spectrum matching %s:%s was found in the '
                              'list of observed spectra' % (key, str(kwarg[key])))
                return []

            if self.debug:
                dbg_string += '%s: %s ' % (key, str(kwargs[key]))

            # extract them from the list
            for dic in ['groups', 'properties']:
                for sub_key in osl[dic].keys():
                    osl[dic][sub_key] = osl[dic][sub_key][vind]

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
        for i,spectrum in enumerate(self.observedSpectraList['spectrum']):
            for key in groups.keys():

                # If not user defined the maximal possible
                # group is assigned
                if key is not 'rv':
                    gn = spectrum.get_group(key)
                    def_groups = groups[key]

                    # if spectrum has no group, but some have been defined, set the group = max(number)+`
                    if gn is None and len(def_groups) > 0:
                        gn = max(def_groups) + 1

                    # if no group is defined for all spectra, start with zero
                    elif gn is None and len(def_groups) == 0:
                        gn = 0
                    # print key, gn
                    # store the groupnumber
                    self.observedSpectraList['group'][key][i] = gn

                else:
                    gn = spectrum.get_group(key)
                    def_groups = groups[key]

                    # if spectrum has no group, but some have been defined, set the group = max(number)+1
                    if gn is None and len(def_groups) > 0:
                        if gn_rv is None:
                            gn = max(def_groups) + 1
                        else:
                            gn = gn_rv + 1

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
            group = {key:self.observedSpectraList['group'][key][i] for key in self.observedSpectraList['group'].keys()}
            self.observedSpectraList['spectrum'][i].set_group(group)
















