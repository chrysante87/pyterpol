import os
import sys
import copy
import warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.constants import c
from auxiliary import is_within_interval
from auxiliary import instrumental_broadening
from auxiliary import interpolate_spec
from auxiliary import interpolate_block_faster
from auxiliary import read_text_file
from auxiliary import rotate_spectrum
from auxiliary import shift_spectrum
from auxiliary import ZERO_TOLERANCE
from defaults import default_grid_order
from defaults import gridDirectory
from defaults import grid_files
from defaults import gridListFile
from defaults import ABS_default_grid_order
from defaults import ABS_gridDirectory
from defaults import ABS_grid_files
from defaults import ABS_gridListFile

# CONSTANTS

class SyntheticSpectrum:
    def __init__(self, f=None, wave=None, intens=None, do_not_load=False, **props):
        """
        Reads the synthetic spectrum and its properties.
        input:
          f.. file with the spectrum
          wave.. wavelength vector of the synthetic spectrum
          intens.. intensity vector of the synthetic spectrum
          do_not_load.. switch for cases, when we want to build the
              class but do not want to load the spectrum.
        **props.. properties of the spectrum, in the correct type
        """

        # reads the spectrum
        if f is not None:
            # from file
            # this delay reading of the data
            self.filename = f
            if not do_not_load:
                self.loaded = True
                self.wave, self.intens = np.loadtxt(f, unpack=True, usecols=[0, 1])
                self.measure_spectrum()
            else:
                self.loaded = False
        else:
            # wavelengths and intensities are given
            self.wave = wave
            self.intens = intens
            self.measure_spectrum()
            self.loaded = True
            self.filename = None

        # setups properties of the synthetic spectrum
        self.properties = []
        for key in props.keys():
            setattr(self, key.lower(), props[key])
            self.properties.append(key.lower())

    def __getitem__(self, key):
        """
        Returns an attribute of the synthtic spectrum.
        Works only with properties.
        input:
           key.. searched attribute
        output:
           prop.. value of the attributr, if not present, False
        """

        if not hasattr(self, key):
            return False
        else:
            return getattr(self, key)

    def __setitem__(self, key, value):
        """
        Changes physical attribute. If it does not
        exist, exception is raised.
        input:
            key..   the attribute to be changes
            value.. nbew value of the attribute
        """

        if not hasattr(self, key):
            raise AttributeError('The atribute %s does not exist.' % key)
        else:
            setattr(self, key, value)

    def __str__(self):
        """
        String representation.
        """
        string = ""
        # if taken from file, prints its name
        if self.filename is not None:
            string = string + "filename:%s " % (self.filename)
            string = string + "loaded:%s " % (str(self.loaded))
        # prints properties of the spectrum
        for prop in self.properties:
            string = string + "%s:%s " % (prop, str(self[prop]))

        # get the wavelength boundaries
        if self.loaded:
            string += "(wmin, wmax): (%s, %s)" % (str(self.wmin), str(self.wmax))
            string = string + '\n'

        return string

    def check_boundaries(self, wmin, wmax):
        """
        Checks that the given wavelengths do not
        overlap the sythetic spectra.
        intput:
            wmin = minimal wavelength
            wmax = maximal wavelength
        """

        # lets have a special case, where the boundaries are None
        if wmin is None:
            wmin = self.wmin
        if wmax is None:
            wmax = self.wmax

        if (wmin - (self.wmin - self.step) < ZERO_TOLERANCE) | \
           (wmax - (self.wmax + self.step) > ZERO_TOLERANCE):
            return False
        else:
            return True

    def keys(self):
        """
        Returns a list of properties.
        """
        return self.properties

    def load_spectrum(self, f=None):
        """
        Loads the spectrum and stores it within the type.
        input:
            f.. filename
        """
        if f is not None:
            self.filename = f

        # check if a binary representation exists -- then load it
        binary_file = self.filename + '.npz'
        if os.path.isfile(binary_file):
            npz = np.load(binary_file, mmap_mode='r')
            self.wave = npz['arr_0']
            self.intens = npz['arr_1']

        # otherwise, load ascii (very slow!) and save it as binary
        else:
            self.wave, self.intens = np.loadtxt(self.filename, unpack=True, usecols=[0, 1])
            print("Saving binary file: " + str(binary_file))
            np.savez(binary_file, self.wave, self.intens)

        # measures the spectrum and marks it as loaded
        self.measure_spectrum()
        self.loaded = True

    def measure_spectrum(self):
        """
        Stores maximal, minimal wavelength and step within the type.
        """

        # saves properties of synthetic
        # spectra - min, max, step
        self.wmin = self.wave.min()
        self.wmax = self.wave.max()
        self.step = self.wave[1] - self.wave[0]

    def pad_continuum(self, wave, intens, bumpsize):
        """
        Pads synthetic spectrum with continua at
        each end.
        input:
          wave, intens.. the input spectrum
        output:
          bump_wave, bump_intens.. the  1-padded spectrum
        """
        # gets properties of teh spectrum
        w0 = wave[0]
        wn = wave[-1]
        step = wave[1] - wave[0]

        # left bump
        l_bump_wave = np.arange(w0 - bumpsize, w0, step)

        # left bump
        r_bump_wave = np.arange(wn + step, wn + bumpsize, step)

        # continuum - just ones
        # l_cont = np.ones(len(l_bump_wave))
        # r_cont = np.ones(len(r_bump_wave))

        # continuum - just ones
        l_cont = 1.0 - np.linspace(0, bumpsize, len(l_bump_wave)) * (1.0 - intens[0]) / bumpsize
        r_cont = intens[-1] + np.linspace(0, bumpsize, len(r_bump_wave)) * (1.0 - intens[-1]) / bumpsize

        # cretes empty arrays
        total_length = len(l_bump_wave) + len(wave) + len(r_bump_wave)
        bump_wave = np.zeros(total_length)
        bump_intens = np.zeros(total_length)

        # copy the bumpers and the spectra
        imin = 0
        imax = 0
        for w, c in zip([l_bump_wave, wave, r_bump_wave], [l_cont, intens, r_cont]):
            imax += len(w)
            bump_wave[imin:imax] = w
            bump_intens[imin:imax] = c
            imin = imax

        return bump_wave, bump_intens

    def get_spectrum(self, wave=None, rv=None, vrot=None, lr=1.0, korel=False,
                     only_intensity=False, wmin=None, wmax=None, keep=False,
                     fwhm=None):
        """
        Return the sythetic spectrum stored within the class. If
        a set of wavelengths is provided, an interpolated spectrum
        is returned.
        input:
          optional:
          wave.. array of deesired wavelengths
          rv.. radila velocity in km/s
          vrot.. projected rotational velocity in km/s
          only_intensity.. returns intensity only
          :param korel
          :param keep
        output:
          wave, intens.. synthetic spectrum
        """
        # checks that we do not pass negative values
        if vrot is not None and vrot < 0.0:
            warnings.warn('vrot cannot be negative! Setting to zero!')
            vrot = 0.0

        if wave is None:
            # for some reason we want to work with the
            # whole spectrum
            # print wmin, wmax
            if wmin is not None and wmax is not None:
                wave, intens = self.select_interval(wmin, wmax)
            else:
                wave = self.wave
                intens = self.intens
            syn_wave = wave.copy()

            # adds the instrumental broadening
            if fwhm is not None and fwhm > ZERO_TOLERANCE:
                intens = instrumental_broadening(syn_wave, intens, width=fwhm)

            if vrot is not None and vrot > ZERO_TOLERANCE:
                # rotates the spectrum
                # print vrot
                intens = rotate_spectrum(syn_wave, intens, vrot)

            if rv is not None and abs(rv) > 0.0:
                # if we want to shift it, we need to pad it,
                # so it does not have to extrapolate
                w0min = wave.min()
                w0max = wave.max()
                mins = np.array([w0min, w0max])
                WAVE_BUMP = np.ceil(np.max(np.absolute(mins * (1 + 1000 * rv / c.value) - mins)))
                syn_wave, intens = self.pad_continuum(syn_wave, intens, WAVE_BUMP)

                # shift it in RV
                syn_wave = shift_spectrum(syn_wave, rv)

            # the spectrum is shrinked
            if lr is not None and abs(lr - 1.0) > ZERO_TOLERANCE:
                intens = intens*lr

            if np.any([x != None for x in [rv, vrot]]):
                # interpolates back
                intens = interpolate_spec(syn_wave, intens, wave)
        else:
            # we are interpolating, so
            # we check boundaries and we
            # also add some more points
            # at each end of the spectrum
            # because we might want to
            # operate with it a bit

            # usually, if it does not fit,
            # we can take a longer spectrum,
            # so there is no point in padding the
            # spectrum # the extension is
            w0min = wave.min()
            w0max = wave.max()

            # the velocity shift rounded up
            mins = np.array([w0min, w0max])

            # Securing additional points on spectrum sides
            # has sense only if we plan to shift it in RV
            if rv is not None and abs(rv) > ZERO_TOLERANCE:
                WAVE_BUMP = np.ceil(np.max(np.absolute(mins * (1 + 1000 * rv / c.value) - mins)))
            else:
                WAVE_BUMP = 0.0

            wmin = w0min - WAVE_BUMP
            wmax = w0max + WAVE_BUMP
            # print wmin, wmax, self.wmin, self.wmax
            if not self.check_boundaries(wmin, wmax):
                warnings.warn('Synthetic spectra do not cover the whole wavelength region' \
                              ' extrapolation has to be employed and THAT IS DANGEROUS! Note that' \
                              ' each spectrum is extended by %f Angstrom at each side.' % (WAVE_BUMP))

            # the part of the spectrum is selected
            # there is no point in working with the
            # whole dataset
            # print wmin, wmax
            syn_wave, intens = self.select_interval(wmin, wmax)

            # adds the instrumental broadening
            if fwhm is not None and fwhm > ZERO_TOLERANCE:
                # intens = instrumental_broadening(syn_wave, intens, width=fwhm)
                intens = instrumental_broadening(syn_wave, intens, width=fwhm)

            # rotates the spectrum
            if vrot is not None and vrot > ZERO_TOLERANCE:
                # intens = rotate_spectrum(syn_wave, intens, vrot)
                # print syn_wave
                intens, syn_wave = rotate_spectrum(syn_wave, intens, vrot, interpolate_back=False)

            # adjusts the spectrum for the radial velocity
            if rv is not None and abs(rv) > ZERO_TOLERANCE:
                syn_wave = shift_spectrum(syn_wave, rv)

            # the spectrum is shrinked
            if lr is not None and abs(lr - 1.0) > ZERO_TOLERANCE:
                intens = intens*lr

            # interpolates to the user specified wavelengths
            intens = interpolate_spec(syn_wave, intens, wave)

        # if we want to extract the spectra in KOREL format
        if korel:
            intens = 1.0 - (lr - intens)

        # if we want to update the class with what
        # we computed
        if keep:
            # update the size of the spectrum
            self.intens = intens
            self.wave = wave
            self.measure_spectrum()

            #update its parameters
            for attr, val in zip(['rv', 'vrot', 'lr', 'korel'], [rv, vrot, lr, korel]):
                if val is not None:
                    setattr(self, attr, val)
                    self.properties.append(attr)
            return

        if only_intensity:
            return intens
        else:
            return wave, intens

    def get_size(self):
        """
        Gets the size of the spectrum i.e. wmin, wmax and step.
        output:
            props.. dictionary with records 'wmin', 'wmax', 'step'
        """

        if self.loaded:
            # guarantees fresh result
            self.measure_spectrum()
            return self.wmin, self.wmax, self.step
        else:
            raise Exception('Spectrum has not been loaded yet.')

    def get_properties(self):
        """
        Returns dictionary with the physical properties of the
        synthetic spectrum.
        Output:
        props.. physical properties of the sythetic spectrum
        """
        # return dictionary with the physical properties
        props = {key: self[key] for key in self.properties}
        return props

    def plot(self, ax=None, savefig=False, figname=None, **kwargs):
        """
        :param figname
        :param savefig
        :param ax: AxesSubplot
        :param kwargs:
        :return:
        """
        w = self.wave
        i = self.intens
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        props = str({prop: self[prop] for prop in self.properties})
        ax.plot(w, i, label=props, **kwargs)
        ax.set_xlim(self.wmin, self.wmax)
        ax.set_ylim(0.95*i.min(), 1.05*i.max())
        ax.set_xlabel('$\lambda(\AA)$')
        ax.set_ylabel('$F_{\lambda}$(rel.)')
        ax.legend(fontsize=10)

        # save the figure
        if savefig:
            if figname is None:
                figname = []
                for key in self.properties:
                    # print key, self.properties
                    figname.extend([key, str(self[key])])
                figname.extend(['wmin', str(self.wmin)])
                figname.extend(['wmax', str(self.wmax)])
                figname = '_'.join(figname) + '.png'
            # save the plot
            plt.savefig(figname)

    def select_interval(self, wmin, wmax):
        """
        Selects a spectral interval from the
        synthetic spectrum.
        :param wmin minimal wavelength
        :param wmax maximal wavelength
        :return wave wavelength vector
        :return intens iuntensity vector
        """

        ind = np.where((self.wave >= wmin) & (self.wave <= wmax))[0]
        wave = self.wave[ind]
        intens = self.intens[ind]

        return wave, intens

    def set_linear_wavelength(self, wmin, wmax, step):
        """
        In case we want to attach linear wavelengths.
        :param wmin
        :param wmax
        :param step
        """

        self.wave = np.arange(wmin, wmax + step / 2., step)

    def truncate_spectrum(self, wmin=None, wmax=None):
        """
        Truncates the spectrum.
        input:
            wmin, wmax.. boundaries in wavelength
        """
        if self.loaded is False:
            raise Exception('The spectrum was not loaded.')
        else:
            # if there is a boundary missing
            # if both, we are waisting compupter time
            if wmin is None:
                wmin = self.wave.min()
            elif wmax is None:
                wmax = self.wave.max()

            # checks that the interval (wmin, wmax) lies within the
            # spectrum. If not exception is raised
            if (self.wave.min() > wmin) | (self.wave.max() < wmax):
                raise ValueError('The spectrum %s does not cover the whole spectral region <%s,%s>.' % \
                                 (str(self).rstrip('\n'), str(wmin), str(wmax)))

            # does the trunctation
            ind = np.where((self.wave >= wmin) & (self.wave <= wmax))[0]
            # ind = np.where(((self.wave - wmin) >= -ZERO_TOLERANCE) & ((self.wave - wmax) <= ZERO_TOLERANCE))[0]
            self.wave = self.wave[ind]
            self.intens = self.intens[ind]

    def write_spectrum(self, filename='synspec.dat', fmt='%12.6f %12.8e', **kwargs):
        """
        Writes the current synthetic spectrum.
        :return:
        """
        header = str(self.get_properties())
        np.savetxt(filename, np.column_stack([self.wave, self.intens]), fmt=fmt, header=header)


class SyntheticGrid:
    def __init__(self, mode='default', flux_type='relative', debug=False):
        """
        Setup the grid.
        input:
            mode..
        """

        # import from defaults
        self.default_grid_order = default_grid_order
        self.gridDirectory = gridDirectory
        self.grid_files = grid_files
        self.gridListFile = gridListFile

        # Table containing list of the SyntheticSpectrum types
        self.SyntheticSpectraList = []

        # Table containing all eligible values
        self.parameterList = []
        self.columns = []

        # grid preference order
        self.gridOrder = None

        # reads default grids
        if mode.lower() != 'custom':
            self.setup_defaults(mode, flux_type)

        # updates debug mode
        self.debug = debug

        # Initializes array in which the wavelength
        # vector is stored
        self.wave = None

    def __str__(self):
        """
        String representation.
        """
        string = "List of spectra:\n"
        for rec in self.SyntheticSpectraList:
            string = string + str(rec)
        return string

    def check_properties(self, **kwargs):
        """
        Checks that at least some spectra have
        """

    def clear_all(self):
        """
        Empties SyntheticSpectraList
        """
        self.SyntheticSpectraList = []

    def deselect_exact(self, l, **props):
        """
        Deletes all cases where we interpolate at exact value.
        :param l: list of spectra selected with self.select_parameters
        :param props: values at which we interpolate
        :return:
        """

        l = np.array(l)
        # print np.shape(l)[-1]
        keys = props.keys()

        # go column by column
        for i in range(0, np.shape(l)[-1]):
            v = props[keys[i]]
            # print v, np.unique(l[:,i])

            # deselects extact matches
            if np.any(abs(np.unique(l[:,i]) - v) < ZERO_TOLERANCE):
                ind = np.where(abs(l[:,i] - v) < ZERO_TOLERANCE)
                l = l[ind]

        return l

    def get_synthetic_spectrum(self, params, wave, order=2, step=0.01, padding=20.0):
        """
        Method which computes the interpolated spectrum
        and wraps it within SyntheticSpectrum class. This
        function should be accessed by the user.
        input:
            params.. dictionary containing values at which
                 we want to interpolate
            order..  number of spectra at which we are going
                 to interpolate, i.e. the order of the
                 fit is k = order-1 for order < 4 and
                 k = 3 for order > 4.
            wave..   wavelength vector for which the synthetic
                 spectrum should be created.
        """
        if isinstance(wave, (list, tuple)):
            wave = np.array(wave)

        # sets up the equidistant wavelength vector
        wmin = wave.min() - padding
        wmax = wave.max() + padding

        # print wmin, wmax, step, order, padding

        # overwrite the wave vector for
        self.set_wavelength_vector(wmin, wmax, step)

        # first of all we need to get list of parameters,
        # which the program will interpolate in
        parlist, vals, keys = self.select_and_verify_parameters(order=order, **params)

        # second creates a list of the spectra used for interpolation
        spectra = self.get_spectra_for_interpolation(parlist, keys, step=step,
                                                     wmin=wmin, wmax=wmax)

        # interpolates the spectra
        intens = self.interpolate_spectra(parlist, spectra, vals)

        # wrapps the interpolated synthetic spectrum within SyntheticSpectrum class
        spectrum = SyntheticSpectrum(wave=self.wave.copy(), intens=intens, **params)

        return spectrum

    def get_all(self, **kwargs):
        """
        Returns all spectra having a certain property, or
        the property is equal to some value. One can list
        all spectra that have a certain value.
        input:
            kwargs.. a dictionary of property = value
        """
        # just in case we got empty dictionary
        if len(kwargs.keys()) == 0:
            return self.SyntheticSpectraList.copy()

        # goes through each stored synthetic spectrum
        spectra = []
        for rec in self.SyntheticSpectraList:
            # goes through passed kwargs
            for i, key in enumerate(kwargs.keys()):
                # print rec[key]
                # print kwargs[key]
                # keys agrees
                if key.lower() in rec.keys():
                    # values are the same
                    if (abs(kwargs[key] - rec[key.lower()]) < ZERO_TOLERANCE):
                        # and we are checking the last record
                        if i == (len(kwargs.keys()) - 1):
                            spectra.append(rec)
                    else:
                        break

        if len(spectra) == 0:
            warnings.warn("No eligible spectrum was found! You probably got out of the grid.")
        return spectra

    def get_available_values(self, prop, **constraints):
        """
        Lists all available properties in a sorted array.
        input:
          prop.. the searched property
          grid.. list of SyntheticSpectrum types -
             in case we want to search a narrowed
             down list
        """
        # lets say we want to contraint the grid
        if len(constraints.keys()) == 0:
            grid = self.SyntheticSpectraList
        else:
            grid = self.get_all(**constraints)

        # returns all eligible values
        values = []
        prop = prop.lower()
        for rec in grid:
            if prop in rec.keys() and rec[prop] not in values:
                values.append(rec[prop])

        return np.sort(values)

    def get_available_values_fast(self, prop, **constraints):
        """
        Returns possible values of a parameter.
        """

        parLis = np.array(self.parameterList)
        # print parLis

        for key in constraints.keys():
            # value
            v = constraints[key]

            # constraining column
            col = self.columns.index(key)
            ind = np.where(abs(parLis[:, col] - v) < ZERO_TOLERANCE)[0]

            # narrow down the list
            parLis = parLis[ind]

        col = self.columns.index(prop.lower())

        # return sorted list
        return sorted(set(parLis[:, col]))

    def get_spectra_for_interpolation(self, parList, header, step=0.01, wmin=None, wmax=None):
        """
        Creates a list of spectra - physical spectra - this
        will require some more tweaking. I would like to account
        for:
        1) spectrum is/is not loaded.
        2) spectrum should be somehow trucated
           there is no need to load all - maybe it
           will even be better to load only small.
           This will be controlled through this method
           which loads the data. Also some spectra
           can be unloaded after several iterations.
        input:
        output:
        """
        # empty list for the spectra
        syntheticSpectra = []

        # switch for spectra truncation
        if (wmin is None) and (wmax is None):
            truncateSpectrum = False
        else:
            truncateSpectrum = True

        # go through each row of the parameter list
        for i, row in enumerate(parList):
            # retrieve spectrum
            props = {prop: row[j] for j, prop in enumerate(header)}
            spectrum = self.get_all(**props)

            # if there are two or more spectra for the
            # same temperature
            if len(spectrum) > 1:
                spectrum = self.resolve_degeneracy(spectrum)
            else:
                spectrum = spectrum[0]

            # load the spectrum if not
            # already loaded
            if not spectrum.loaded:
                if self.debug:
                    print "Loading spectrum: %s" % (str(spectrum).rstrip('\n'))
                else:
                    print "Loading spectrum: %s" % (str(spectrum).rstrip('\n'))

                spectrum.load_spectrum()

                # check that the synthetic spectrum has sufficient size
                spectrum.check_boundaries(wmin, wmax)

                # truncates the loaded spectrum
                if truncateSpectrum:
                    if self.debug:
                        print "Truncating spectrum to: (%f,%f)" % (wmin, wmax)
                    spectrum.truncate_spectrum(wmin, wmax)
            else:
                if self.debug:
                    print "Spectrum loaded: %s" % (str(spectrum).rstrip('\n'))

            # We have to be sure that the spectra aren't off
            # each other by less than one step
            swmin, swmax, sstep = spectrum.get_size()
            if np.any(np.abs([swmin - wmin, swmax - wmax, sstep - step]) > ZERO_TOLERANCE):
                if self.debug:
                    print "Spectrum %s does not have the wavelength scale (wmin, wmax,step)=(%s, %s, %s)" % \
                          (str(spectrum).rstrip('\n'), str(wmin), str(wmax), str(step))

                # if they do not agree - we have to interpolate
                # it is cruacial that all spectra have the same
                # wavelength scale
                if self.wave == None:
                    wave = np.arange(wmin, wmax + step / 2., step)
                else:
                    wave = self.wave

                # interpolate the spectrum to the wavelength scale
                intens = spectrum.get_spectrum(wave=wave, only_intensity=True)

            else:
                if self.debug:
                    print "Wavelenght scale of spectrum: %s is (wmin, wmax,step)=(%s, %s, %s)." % \
                          (str(spectrum).rstrip('\n'), str(wmin), str(wmax), str(step))
                # read out the intensities
                intens = spectrum.get_spectrum(only_intensity=True)

            # print len(intens)
            # append spectrum to the list
            syntheticSpectra.append(intens)

        return syntheticSpectra

    def interpolate_spectra(self, parList, synspectra, parameters):
        """
        Interpolates in all parameters.
        input:
            parlist.. list generated with select parameters method
            synspectra.. list generated with the get_spectra_for_interpolation
                 method
            parameters.. list of parameter values in which we interpolate
                 the order must be the same as in case of parlist
                 this is guaranteed by the ouput of select_and_verify_parameters method
        output:
            intens.. the resulting array of intensities
        """

        # convert to arrays, easier to handle
        plist = np.array(parList)
        syns = np.array(synspectra)
        ncol = len(plist[0])
        pars = parameters
        while ncol > 0:

            # extract new value
            xnew = pars[ncol - 1]

            # print xnew

            new_plist = []
            new_syns = []

            # take the first row
            j = 0
            while j < len(plist):
                row = plist[j]

                # narrow it down - all values
                # that have the first ncol-1
                # values the same are chosen
                t_plist = plist.copy()
                t_syns = syns.copy()
                for i in range(ncol - 1):
                    ind = np.where(abs(t_plist[:, i] - row[i]) < ZERO_TOLERANCE)[0]
                    t_plist = t_plist[ind]
                    t_syns = t_syns[ind]

                # if there is really nothing to interpolate in
                # the one value is copied and we proceed to next
                # step
                if len(t_plist) == 1:

                    if self.debug:
                        print "Skipping interpolation in %s - there is only one spectrum for values %s." % \
                              (str(xnew), str(t_plist[:, :ncol - 1]))

                    intens = t_syns[0]
                    new_plist.append(row[:ncol - 1])
                    new_syns.append(intens)
                    j += len(ind)
                    continue

                # sort according to the last columns
                ind = np.argsort(t_plist[:, ncol - 1])

                # extract the abscissa
                x = t_plist[ind, ncol - 1]
                t_syns = t_syns[ind]

                if self.debug:
                    print "Interpolating in vector: %s at value %s." % (str(x), xnew)

                # everything is sorted, we can interpolate
                # unless our value is exact ofc.
                intens = interpolate_block_faster(x, t_syns, xnew)

                # add it to new plists and syns
                new_plist.append(row[:ncol - 1])
                new_syns.append(intens)
                j += len(ind)

            syns = np.array(new_syns)
            plist = np.array(new_plist)
            ncol = len(plist[0])

        return syns[0]

    @staticmethod
    def list_modes():
        """
        This method lists available modes for the SyntheticGrid.
        :return:
        """
        # go over differents modes
        string = 'List of registered modes and their properties follows:\n'
        for i in range(0, len(self.grid_files['identification'])):
            string += ''.ljust(100,'=') + '\n'
            string += 'mode: %s:\n' % self.grid_files['identification'][i]
            string += 'directories: %s \n' % str(self.grid_files['directories'][i])
            string += 'columns: %s\n' % str(self.grid_files['columns'][i])
            string += 'families: %s\n' % str(self.grid_files['families'][i])
        string += ''.ljust(100,'=') + '\n'

        return string

    def narrow_down_grid(self, **kwargs):
        """
        To speed up computations, one can
        narrow down the grid, to certain
        family, or parameter range.
        input:
          One can either fix a parameter:
          par = value
          or fix an interval:
          par = (vmin, vmax)
        output:
          list of synthetic spectra
        """

        # separate fixed from free
        fixed = {}
        free = {}
        for key in kwargs.keys():
            if not isinstance(kwargs[key], (tuple, list)):
                fixed[key] = kwargs[key]
            else:
                free[key] = kwargs[key]

        # first narrow down the fixed ones
        grid = self.get_all(**fixed)

        # if there are no other restrictions -
        # this option is covered with get_all
        # method ofc.
        if len(free.keys()) == 0:
            return grid
        else:
            narrowed_grid = []
            for rec in grid:
                for i, key in enumerate(free.keys()):
                    if key.lower() in rec.keys():
                        if (rec[key.lower()] >= free[key][0]) & \
                                (rec[key.lower()] <= free[key][1]):
                            # all keys must agree
                            if i == len(free.keys()) - 1:
                                narrowed_grid.append(rec)

        # if we narrowed it down to zero
        if len(narrowed_grid) == 0:
            warnings.warn("The narrowed-down grid is empty!!")
        return narrowed_grid

    def read_list_from_file(self, f, columns, directory=None, family=None):
        """
        Reads list of grid spectra from a file.
        input:
          f.. file containing the records
          columns.. column description
          family.. family is provided
          directory.. directory where the files are stored
                  this option should be used in case
                  when the path to the spectra = filename
                  is relative only within the file
        """
        # read from text_file
        if columns is None:
            raise KeyError('Description of the input file was not specified.')

        # There are two mandatory records - 1) path
        # to the sythetic spectrum and 2) family tp
        # which the spectrum belongs. By family I mean
        # a published grid. In this family, there should
        # not exist 2 spectra with the same properties
        # Missing filename will raise an error,
        # missing family will raise warnig, because
        # all read spectra will be assigned the same family
        hasFamily = False
        if family is None:
            for rec in columns:
                if rec.upper() == 'FAMILY':
                    hasFamily = True
                    addFamily = False
        else:
            hasFamily = True
            addFamily = True
        if not hasFamily:
            warnings.warn("The family (aka the grid) of spectra was not specified. Assigning family...")
            families = self.get_available_values('FAMILY')
            family = 'family' + str(len(families))
            addFamily = True

        # WE CHECK THAT THERE IS FILENAME RECORD
        # FOR EACH SPECTRUM - WITHOU THAT WE WONT
        # COMPUTE ANYTHING
        hasFilename = False
        for rec in columns:
            if rec.upper() == 'FILENAME':
                hasFilename = True
        if not hasFilename:
            raise KeyError('Record filename = path to the spectrum is missing in the column description!')

        lines = read_text_file(f)
        # go through file, line by line
        for j, line in enumerate(lines):

            # store one line = one spectrum info
            rec = {}
            data = line.split()

            # make sure, we have description of all
            # properties
            if len(data) > len(columns):
                raise KeyError('Description of some columns is missing.')
            for i, col in enumerate(columns):
                # the name should be the only string
                if col.upper() in ['FAMILY']:
                    rec[col.upper()] = data[i]
                elif col.upper() in ['FILENAME']:
                    rec[col.upper()] = os.path.join(directory, data[i])
                else:
                    rec[col.upper()] = float(data[i])

            # Adds family if needed
            if addFamily:
                rec['FAMILY'] = family

            filename = rec.pop('FILENAME')
            synspec = SyntheticSpectrum(f=filename, do_not_load=True, **rec)

            # Adds the record so the synthetic spectra list
            self.SyntheticSpectraList.append(synspec)

            # Adds the record to the parameterList - so without family
            rec.pop('FAMILY')
            phys_cols = [x.lower() for x in columns if x not in ['FAMILY', 'FILENAME']]
            self.parameterList.append([rec[col.upper()] for col in phys_cols])
            # also stores identification of columns
            if j == 0:
                self.columns = phys_cols

    def resolve_degeneracy(self, speclist):
        """
        If there are more spectra having the same
        parameters, one has to choose one prefered.
        input:
            speclist.. list of SyntheticSpectrum types corresponding to same properties
        """
        # if we did not set up the order -> error
        if self.gridOrder is None:
            raise KeyError('There are same spectra for the same parameters.'
                           ' I think it is because we have more grids, that overlap.'
                           ' You can overcome this by setting gridOrder variable.')

        indices = []
        for i in range(0, len(speclist)):
            # print speclist[i]['family']
            indices.append(self.gridOrder.index(speclist[i]['family']))

        # justr in case there was something peculiar
        if np.any(indices < -1):
            warnings.warn('At least one grid was not found in the gridOrder variable.'
                          ' Verify that the names set in gridOrder agree with family names of spectra.')

        # return spectrum with the smallest index
        return speclist[np.argmin(indices)]

    def select_and_verify_parameters(self, order=2, **props):
        """
        A wrapper to the select_parameters method.
        This method can deal with overlaps of grids.
        But since it does not know the grid apriori,
        it is unable to recognize wrong result.

        This wrapper checks the result.

        input:
          order.. maximal number of interpolated spectra
          props.. dictionary of interpolated parameters
        output:
          parlist.. each row represents one spectrum
                which is needed to interpolate in
                give props
          vals.. values in which we interpolate
          keys.. names of the interpolated
        """

        # all parameters are defined lowercase
        # so we have to convert it
        for key in props.keys():
            v = props.pop(key)
            props[key.lower()] = v

        if self.debug:
            print "In select_and_verify_parameters: order=%i properties:" % (order)
            print str(props)

        # keys and values
        keys = props.keys()
        vals = [props[key] for key in props.keys()]

        # gets the parameter list
        # print order, props
        parlist = self.select_parameters(order=order,**props)
        # print parlist

        # deselect reduntdant spectra
        parlist = self.deselect_exact(parlist, **props)

        if len(parlist) == 0:
            raise Exception('Do %s lie within the grid? I do not think so...' % (str(props)))

        if self.debug:
            print 'Following parameters were chosen with select_parameters method:'
            for row in parlist:
                print row

        # checks the result
        temp = np.array(parlist)
        # print temp, vals
        for i, val in enumerate(vals):
            # print val,  temp[:, i], is_within_interval(val, temp[:, i])
            if not is_within_interval(val, temp[:, i]):
                raise ValueError('Parameters %s lie outside the grid.' % (str(props)))

        return parlist, vals, keys

    def select_parameters(self, values=[], order=2, constraints={}, **props):
        """
        Creates a final list - this is still
        first guess. I think that searching up
        eligible values and spectra can be done
        better.
        input:
          grid - synthetic spectraList, which is searched
          order - how many spectra are used this is
              adjusted dynamically if there are not
              enough values
          constraints - resolve conflicts between grids
          props - properties in which we fit
        output:
          values of spectra for interpolation
        """


        # extract the parameter and its values
        key = props.keys()[0].lower()
        v = props.pop(key)

        # print key, constraints, props
        # list eligible values for a given parameter
        elig_vals = np.array(self.get_available_values_fast(key, **constraints))
        # print key, elig_vals

        # sorts the grid, from nearest to farthest
        ind = np.argsort(abs(elig_vals - v))
        vals = elig_vals[ind]

        # equality check

        # print vals, v, key
        # what if the grid step is inhomogeneous? - actually it is
        # in z - what shall we do, what shall we do???
        if vals[:order].min() > v or vals[:order].max() < v:
            # TODO think of something better than this!!!!!!
            try:
                lower = np.max(vals[np.where(vals - v < ZERO_TOLERANCE)[0]])
                upper = np.min(vals[np.where(vals - v > ZERO_TOLERANCE)[0]])
                vals = np.array([lower, upper])
            except:
                pass
            # print lower, upper, vals


        # checks that there is not equality
        # if np.any(abs(vals - v) < ZERO_TOLERANCE):
        #     ind = np.argmin(abs(vals - v))
        #     vals = [vals[ind]]
        #
        #     if self.debug:
        #         print "%s=%s is precise. Skipping choice of parameters." % (key, str(v))

        # if the eligible values do not surround the parameter
        if not is_within_interval(v, vals):
            return values

        # if there are no other spectra to interpolate in
        if len(props.keys()) == 0:
            for i in range(0, len(vals)):
                row = []
                # append those that are already fixed
                for key in constraints.keys():
                    row.append(constraints[key])

                # append the last parameter
                row.append(vals[i])
                # append the row
                values.append(row)

                # once 'order' spectra are appended, we can
                # end
                if i == order - 1:
                    break
            return values
        else:
            j = 0
            for i in range(0, len(vals)):

                # add a constraint
                constraints[key] = vals[i]

                # recursively calls the function
                values_new = self.select_parameters(values=copy.deepcopy(values), order=order, constraints=constraints,
                                                    **props)

                # some searches are in vain - so we
                # wait until meaningful calls accumulate
                if len(values_new) > len(values):
                    j += 1

                # copis the result, so we can go on
                values = values_new

                # remove constraint
                constraints.pop(key)
                if j == order:
                    break
        return values

    def setup_defaults(self, mode, flux_type):
        """
        Given a key loads a grid stored within
        the directory.
        input:
            mode.. one of the defaults mode OSTAR, BSTAR, POLLUX, AMBRE
                   defaulst = all
            flux_type.. either relative or absolute
        """

        # we do not want to bother with the case
        mode = mode.upper()
        flux_type = flux_type.upper()

        # select the correct type of flux
        # note we cannot overwrite globals, but only class variables
        if flux_type == 'ABSOLUTE':
            self.grid_files = ABS_grid_files
            self.gridDirectory = ABS_gridDirectory
            self.gridListFile = ABS_gridListFile
            self.default_grid_order = ABS_default_grid_order

        # select properties
        ind = self.grid_files['identification'].index(mode)

        if ind < 0:
            raise ValueError('Default settings named %s not found.' % (mode))

        dirs = self.grid_files['directories'][ind]
        cols = self.grid_files['columns'][ind]
        fams = self.grid_files['families'][ind]

        # reads the grid files
        for i, d in enumerate(dirs):
            spectralist = os.path.join(self.gridDirectory, d, self.gridListFile)
            directory = os.path.join(self.gridDirectory, d)
            self.read_list_from_file(spectralist, cols, family=fams[i], directory=directory)

        # also sets the default grid order
        self.set_grid_order(self.default_grid_order)

    def set_mode(self, mode='default'):
        """
        Set different mode.
        :param mode:
        :return:
        """

        debug = self.debug
        self.__init__(mode=mode, debug=debug)

    def set_grid_order(self, arr):
        """
        Sets grid preference.
        input:
            arr = list of spectra grid1 > grid2 > grid3...
        """

        self.gridOrder = arr

    def set_wavelength_vector(self, wmin, wmax, step):
        """
        Store the wavelength vector within the class.
        input:
            wmin.. minimal wavelength
            wmax.. maximal wavelength
            step.. step size in the wavelength
        """
        nstep = int((wmax - wmin)/step)+1
        self.wave = np.linspace(wmin, wmax, nstep)
