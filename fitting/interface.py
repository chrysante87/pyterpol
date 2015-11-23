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

    def add_observation(self, **kwargs):
        """
        Adds observation to the list.
        INPUT:
            see class ObservedSpectrum (observations module) for details.
        """
        # adds the spectrum and loads it
        print kwargs
        obs = ObservedSpectrum(**kwargs)
        print obs
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
        groups = dict()

        # go through ech spectrum and store defined values
        for spectrum in self.observedSpectraList['spectrum']:
            print spectrum
            for key in spectrum.group.keys():
                if key not in groups.keys():
                    groups[key] = []
                groups[key].append(spectrum.group[key])

        # only unique values are needed
        for key in groups.keys():
            groups[key] = np.unique(groups[key]).tolist()

        return groups

    def read_groups(self):
        """
        Goes through each record within.
        extracts available groups. If the group is
        not set, it is automatically set to 0 for
        anything apart from RV.

        For RV each additional spectrum is given a
        unique group, unless it is a KOREL spectrum.
        In that case, it assigns the same group
        to all spectra covering the same wavelength
        interval.
        """
        # First go through each spectrum to see, which
        # groups were defined.
        group_keys = self.observedSpectraList['group'].keys()
        for spectrum in self.observedSpectraList['spectrum']:
            for key in spectrum.group.keys():

                # if there is a key, which is not defined, append it
                if key not in group_keys:
                    group_keys.append(key)

        # initialize variables information on the spectra
        for key in group_keys:
            if key not in self.observedSpectraList['group'].keys():
                self.observedSpectraList['group'][key] = []

        # Initial RV group number

        # Assigning the values
        for spectrum in self.observedSpectraList:
            for key in group_keys:
                # for non-RVs it is simple, just store user-defined value, or zero
                if key is not 'rv':
                    temp_gn = spectrum.get_group(key)
                    # if it does not belong to any group, its number is zero
                    if temp_gn is None:
                        groupNumber = 0
                    else:
                        groupNumber = temp_gn
                    self.observedSpectraList['group'][key].append(groupNumber)
                else:
                    # For observed spectra it is simple -
                    if not spectrum.korel:
                        temp_gn = spectrum.get_group(key)











