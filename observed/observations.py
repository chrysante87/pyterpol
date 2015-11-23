import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep
from scipy.interpolate import splev


class ObservedSpectrum:
    """
    A wrapper class for the observed spectra.
    """
    def __init__(self, wave=None, intens=None, error=None, filename=None, component='ALL', korel=False):
        """
        Setups the class.
        INPUT:
        wave.. 	wavelength vector
        intens..   	intensity vector
        error..  	error vector
        filename..	source spectrum
        component.. component to which the spectrum belongs
        korel.. 	korel mode
        """

        # pass all arguments
        self.wave = wave
        self.intens = intens

        # lets have a look at the errors
        if error is None:
            warnings.warn("I found no array with errorbars of observed intensities. "
                          "Do not forget to assign them later!")
            self.error = None
            self.hasErrors = False

        # Making assumption here, that we are
        # not passing errors, if we are not
        # passing wavelengths and intensities
        elif isinstance(error, (float, int)):
            self.error = np.ones(len(wave)) * error
            self.hasErrors = True
        else:
            self.error = error
            self.hasErrors = True

            # sets that the spectrum is loaded
            if (wave is not None) and (intens is not None):
                self.loaded = True
                self.read_size()

                # check lengths of intens and wave
                self.check_length()
            else:
                self.loaded = False

            # if we provided the filename
            self.filename = filename
            if (not self.loaded) and (self.filename is not None):
                self.read_spectrum_from_file(filename)
            elif (not self.loaded) and (self.filename is None):
                warnings.warn('No spectrum was loaded. This class is kinda useless without a spectrum. '
                            'I hope you know what you are doing.')

            # setup korel and check that it is proper
            self.component = component
            self.korel = korel
            self.check_korel()

    def __str__(self):
        """
        String representation of the class.
        """
        pass

    def check_korel(self):
        """
        If korel is set, component must be set too.
        """
        if (self.korel) and (self.component.upper() == 'ALL'):
            raise ValueError('In the korel regime, each spectrum must be assigned component! '
                             'Currently it is set to %s.' % str(self.component))

    def check_length(self):
        """
        Checks that wavelengths and intensities have the same length.
        """
        if len(self.wave) != len(self.intens):
            raise ValueError('Wavelength vector and intensity vector do not have the same length!')

    def check_loaded(self):
        """
        Checks that spectrum is loaded.
        """
        if not self.loaded:
            raise ValueError('The spectrum is not loaded.')

    def free_spectrum(self):
        """
        Deletes the stored spectrum.
        """
        self.wave = None
        self.intens = None
        self.error = None
        self.loaded = False
        self.hasErrors = False

    def get_boundaries(self):
        """
        Returns the minimal and the maximal wavelength
        of the spectrum.
        """
        self.read_size()
        return self.wmin, self.wmax

    def get_sigma_from_continuum(self, cmin, cmax, store=True):
        """
        Estimates the error of the flux from the scatter in
        continuum.
        INPUT:
          cmin..	minimal continuum wavelength
          cmax..	maximal continuum wavelength
          store..	save the error as an error bar
        OUTPUT:
          stddev	estimated error bar
        """
        # is the spectrum loaded ?
        self.check_loaded()

        # get the part around continue
        intens = self.get_spectrum(wmin=cmin, wmax=cmax)[1]

        # get the scatter
        stddev = intens.std(ddof=1)

        # save it as an error
        if store:
            self.error = stddev*np.ones(len(self.wave))

        return stddev

    def get_sigma_from_fft(self, nlast=20, store=True):
        """
        Estimates the noise using the FFT.
        """

        # check that everything is loaded
        self.check_loaded()
        self.read_size()

        # get the linear scale
        lin_wave = np.linspace(self.wmin, self.wmax, self.npixel)

        # interpolate to linear scale
        tck = splrep(self.wave, self.intens)
        lin_intens = splev(lin_wave, tck)

        # perform the FFT and shift it
        fft_intens = np.fft.fftshift(np.fft.fft(lin_intens))

        # ontly tbe absolute values are interesting
        fft_intens = np.absolute(fft_intens)
        fft_intens[:-nlast] = 0.0 + 0.0j

        fft_intens = np.fft.ifft(fft_intens)
        #print np.absolute(fft_intens.mean())
        plt.plot(fft_intens, 'k-')
        plt.plot(fft_intens[-nlast:], 'ro')
        plt.show()

        #the nlast points are in general scatter - this should be done better...

        stddev = fft_intens[-nlast:].std(ddof=-1) * fft_intens[-nlast:].mean()

        # store the value as an erro if needed
        if store:
            self.error = stddev * np.ones(len(self.wave))

        return stddev

    def get_spectrum(self, wmin=None, wmax=None):
        """
        Returns the spectrum.
        OUPUT:
          self.wave..	wavelengths
          self.intens..	intensities
          self.error..	errors
        """
        if not self.loaded:
            raise Exception('The spectrum %s has not been loaded yet!' % str(self))
        else:
            # the whole spectrum
            if wmin is None and wmax is None:
                return self.wave.copy(), self.intens.copy(), self.error.copy()
            else:
                # corrects boundaries if needed
                if wmin is None:
                    wmin = self.wave.min()
                if wmax is None:
                    wmax = self.wave.max()

                # selects the spectrum part
                ind = np.where((self.wave >= wmin) & (self.wave <= wmax))
                return self.wave[ind].copy(), self.intens[ind].copy(), self.error[ind].copy()

    def get_wavelength(self):
        """
        Returns the wavelength vector.
        OUPUT:
          self.wave..	wavelengths
        """
        if not self.loaded:
            raise Exception('The spectrum %s has not been loaded yet!' % str(self))
        else:
            return self.wave.copy()

    def read_size(self):
        """
        Gets the minimal wavelength, maximal wavelenbgth
        and the mean step. Linearity in wavelength is not
        required.
        """
        if not self.loaded:
            raise Exception('The spectrum %s has not been loaded yet!' % str(self))

        self.wmin = self.wave.min()
        self.wmax = self.wave.max()
        self.npixel = len(self.wave)
        self.step = np.mean(self.wave[1:] - self.wave[:-1])

    def read_spectrum_from_file(self, filename, global_error=None):
        """
        Reads the spectrum from a file. Following format
        is assumed: %f %f %f (wavelength, intensity, error).
        If user does not provide errors, we still attempt
        to load teh spectrum.
        Comments are denoted with '#'.
        INPUT:
             filename..    the file, from which the spectra are loaded.
             global_error..global error applies to every single observation
        """
        try:
            # first we try to load 3 columns, i.e with errors
            self.wave, self.intensity, self.error = np.loadtxt(filename, unpack=True, usecols=[0, 1, 2])
            self.hasErrors = True
        except:

            # we failed, so we attempt to load two columns
            self.wave, self.intensity = np.loadtxt(filename, unpack=True, usecols=[0, 1])

            # error was not set up
            if global_error is None:
                warnings.warn("I found no errorbars of the observed intensities in file: %s! "
                              "I assume they will be provided later. I remember!!" % (filename))
                self.hasErrors = False
            # error was set up
            else:
                self.error = global_error * np.ones(len(self.wave))
                self.hasErrors = True

                # the spectrum is marked as loaded
                self.loaded = True

        # the spectrum is checked
        self.check_length()
        self.read_size()

    def set_error(self, vec_error=None, global_error=None):
        """
        A tool to set the error.
        INPUT:
            vec_error..  vector of errors, one for each intensity point
            glob_error.. one global errror for each point
        """
        if vec_error is not None:
            self.error = vec_error
            self.hasErrors = True
        if global_error is not None:
            self.error = global_error * np.ones(len(self.wave))
            self.hasErrors = True

    def set_spectrum_from_arrays(self, wave, intens, error):
        """
        Stores the spectrum from arrays. It is assumed
        that user also provides error vector.
        INPUT:
            wave..	wavelength array
            intens..	intensity array
            error..	error array
        """
        self.wave = wave
        self.intens = intens
        self.error = error
        self.loaded = True
        self.hasErrors = True

        # checking and reading
        self.check_length()
        self.read_size()
