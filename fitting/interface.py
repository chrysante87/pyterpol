# -*- coding: utf-8 -*-
import numpy as np
from pyterpol.observed.observations import ObservedSpectrum


class ObservedList(object):
    """
    A helper class which groups all observed spectra and
    prepares necessary parameters for fitting.
    """
    def __init__(self):
        """
        Setups the class
        """
        # dictionary containing all observed spectra, apart from that
        # it also carries information on. A group fro radial velocities
        # has to be always set, because we intend to fit spectra acquired
        # on different times.
        self.observedSpectraList = dict(spectrum=[], group=dict())
        self.groupValues = dict()

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

    def add_observation(self, **kwargs):
        """
        Adds observation to the list.
        INPUT:
            see class ObservedSpectrum (observations module) for details.
        """
        # adds the spectrum and loads it
        # print kwargs
        obs = ObservedSpectrum(**kwargs)
        # print obs
        self.observedSpectraList['spectrum'].append(obs)

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

    def set_groups(self):
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

    def _set_groups_to_spectra(self):
        """
        Propagates groups, which are set in observedSpectraList,
        in in dividual spectra.
        """
        for i in range(0, len(self.observedSpectraList['spectrum'])):
            group = {key:self.observedSpectraList['group'][key][i] for key in self.observedSpectraList['group'].keys()}
            self.observedSpectraList['spectrum'][i].set_group(group)
















