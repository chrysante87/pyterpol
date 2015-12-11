import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splrep
from scipy.interpolate import splev
# from pyterpol.synthetic.auxiliary import ZERO_TOLERANCE

# repeat userwarnings
warnings.simplefilter('always', UserWarning)


class ObservedSpectrum:
    """
    A wrapper class for the observed spectra.
    """
    def __init__(self, wave=None, intens=None, error=None, filename=None,
                 component='all', korel=False, group=None, debug=False):
        """
        Setups the class.
        INPUT:
            wave.. 	wavelength vector
            intens..   	intensity vector
            error..  	error vector
            filename..	source spectrum
            component.. component to which the spectrum belongs
            korel.. 	korel mode
            group         group in which the spectrum belongs for a given parameter.
                          all parameters within one group have the same value of a
                          a given parameter. Type = dictionary(param=groupnumber)
        """
        # empty arrays, taht will be filled
        # with read_size
        self.wmin = None
        self.wmax = None
        self.step = None
        self.npixel = None

        # pass all arguments
        self.wave = wave
        self.intens = intens

        # lets have a look at the errors
        if error is None:
            warnings.warn("I found no array with errorbars of observed intensities. "
                          "Do not forget to assign them later!")
            self.error = None
            self.global_error = None
            self.hasErrors = False

        # sets that the spectrum is loaded
        if (wave is not None) and (intens is not None):
            self.loaded = True
            self.read_size()

            # check lengths of intens and wave
            self.check_length()

            # set the error
            if isinstance(error, (float, int)) and error is not None:
                self.error = np.ones(len(wave)) * error
                self.hasErrors = True
                self.global_error = error
            elif error is not None:
                self.error = error
                self.hasErrors = True
                self.global_error = None
        else:
            self.loaded = False

        # if we provided the filename
        self.filename = filename
        if (not self.loaded) and (self.filename is not None):
            self.read_spectrum_from_file(filename, global_error=error)
        elif (not self.loaded) and (self.filename is None):
            warnings.warn('No spectrum was loaded. This class is kinda useless without a spectrum. '
                          'I hope you know what you are doing.')

        # setup korel and check that it is proper
        self.component = component
        self.korel = korel
        self.check_korel()

        # setup the group
        self.group = dict()
        if group is not None:
            self.set_group(group)

        # setup debug mode
        self.debug = debug

    def __str__(self):
        """
        String representation of the class.
        """
        string = ''
        for var in ['filename', 'component', 'korel', 'loaded', 'hasErrors', 'global_error', 'group']:
            string += "%s: %s " % (var, str(getattr(self, var)))
        if self.loaded:
            string += "%s: %s " % ('(min, max)', str(self.get_boundaries()))
        string += '\n'
        return string

    def check_korel(self):
        """
        If korel is set, component must be set too.
        """
        if self.korel and str(self.component).lower() == 'all':
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

    def get_group(self, param):
        """
        Checks that group is assigned for a given
        parameter. If not - zero is automatically
        assigned.
        :param: param string representing the queried parameter
        :return the group for the given parameter
        """
        if param.lower() in self.group:
            return self.group[param]
        else:
            return None

    def get_sigma_from_continuum(self, cmin, cmax, store=True):
        """
        Estimates the error of the flux from the scatter in
        continuum.
        :param cmin the minimal continuum value
        :param cmax the maximal continuum value
        :param store save the found error as an error
        :return stddev the standard deviation
        """
        # is the spectrum loaded ?
        self.check_loaded()

        # get the part around continue
        intens = self.get_spectrum(wmin=cmin, wmax=cmax)[1]

        # get the scatter
        stddev = intens.std(ddof=1)

        # save it as an error
        if store:
            self.error = stddev * np.ones(len(self.wave))

        return stddev

    def get_sigma_from_fft(self, nlast=20, store=True):
        """
        Estimates the noise using the FFT.
        :param nlast length opf the FFT spectrum tail used to estimate the scatter
        :param store should we save the standard deviation
        """
        # TODO The function is probably not correct in its current state
        # check that everything is loaded
        self.check_loaded()
        self.read_size()

        # get the linear scale
        lin_wave = np.linspace(self.wmin, self.wmax, self.npixel)

        # interpolate to linear scale
        tck = splrep(self.wave, self.intens)
        lin_intens = splev(lin_wave, tck)

        # perform the FFT and shift it
        fft_intens = np.fft.fft(lin_intens)

        # get absolute values
        abs_fft_intens = np.absolute(fft_intens)

        # get the high frequency tail
        abs_fft_intens = abs_fft_intens[len(abs_fft_intens) / 2 - nlast + 1:len(abs_fft_intens) / 2 + nlast]

        # estimate the error
        stddev = abs_fft_intens.std() * abs_fft_intens.mean()

        # store the value as an erro if needed
        if store:
            self.error = stddev * np.ones(len(self.wave))

        return stddev

    def get_spectrum(self, wmin=None, wmax=None):
        """
        Returns the spectrum with wavelengths wmin -> wmax
        :param wmin minimal wavelength
        :param wmax maximal wavelength
        :return wave, intens. error (optional) - the observed spectrum,
        wavelength, intensity and error (if it is given)
        """
        if not self.loaded:
            raise Exception('The spectrum %s has not been loaded yet!' % str(self))
        else:
            # the whole spectrum
            if wmin is None and wmax is None:
                if self.error is not None:
                    return self.wave.copy(), self.intens.copy(), self.error.copy()
                else:
                    return self.wave.copy(), self.intens.copy()
            else:
                # corrects boundaries if needed
                if wmin is None:
                    wmin = self.wave.min()
                if wmax is None:
                    wmax = self.wave.max()

                # selects the spectrum part
                ind = np.where((self.wave >= wmin) & (self.wave <= wmax))[0]

                if self.error is not None:
                    return self.wave[ind].copy(), self.intens[ind].copy(), self.error[ind].copy()
                else:
                    return self.wave[ind].copy(), self.intens[ind].copy()

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

        props = str({'filename': self.filename})
        ax.plot(w, i, label=props, **kwargs)
        ax.set_xlim(self.wmin, self.wmax)
        ax.set_ylim(0.95*i.min(), 1.05*i.max())
        ax.set_xlabel('$\lambda(\AA)$')
        ax.set_ylabel('$F_{\lambda}$(rel.)')
        ax.legend(fontsize=10)

        # save the figure
        if savefig:
            if figname is None:
                figname = self.filename + '.png'

            # save the plot
            plt.savefig(figname)

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
        :param filename spectrum source file
        :param global_error the error applicable to the spectrum
        :return None
        """
        try:
            # first we try to load 3 columns, i.e with errors
            self.wave, self.intens, self.error = np.loadtxt(filename, unpack=True, usecols=[0, 1, 2])
            self.hasErrors = True
        except:
            # we failed, so we attempt to load two columns
            self.wave, self.intens = np.loadtxt(filename, unpack=True, usecols=[0, 1])

            # error was not set up
            if global_error is None:
                warnings.warn("I found no errorbars of the observed intensities in file: %s! "
                              "I assume they will be provided later. I remember!!" % filename)
                self.hasErrors = False
                self.global_error = None

            # error was set up
            else:
                self.error = global_error * np.ones(len(self.wave))
                self.hasErrors = True
                self.global_error = global_error

        # the spectrum is marked as loaded
        self.loaded = True

        # the spectrum is checked
        self.check_length()
        self.read_size()

    def set_error(self, vec_error=None, global_error=None):
        """
        Sets error to the spectrum..either local or global.
        :param vec_error vector error len(vec_error) = len(spectrum)
        :param global_error int float error applied to the whole spectrum
        """
        if vec_error is not None:
            self.error = vec_error
            if len(vec_error) != len(self.npixel):
                raise ValueError('The lenght of the error vector and the length of the spectrum do not match (%s, %s)'
                                 % (len(vec_error), str(self.npixel)))
            self.hasErrors = True
            self.global_error = None
        if global_error is not None:
            self.error = global_error * np.ones(len(self.wave))
            self.hasErrors = True
            self.global_error = global_error

    def set_group(self, group):
        """
        Sets a group to the spectrum
        :param group a dictionary of pairs parameter + group
        """
        for key in group.keys():
            self.group[key.lower()] = group[key]

    def set_spectrum_from_arrays(self, wave, intens, error):
        """
        Stores the spectrum from arrays. It is assumed
        that user also provides error vector.
        :param wave wavelength vector
        :param intens intensity vector
        :param error eror vector
        """
        self.wave = wave
        self.intens = intens
        self.error = error
        self.loaded = True
        self.hasErrors = True

        # checking and reading
        self.check_length()
        self.read_size()
