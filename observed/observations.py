import warnings
import numpy as np

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
	    warnings.warn("I found no array with errorbars of observed intensities. " \
			  "Do not forget to assign them later!")
	    self.error = None
	    self.hasErrors = False
	
	# Making assumption here, that we are
	# not passing errors, if we are not
	# passing wavelengths and intensities
	elif isinstance(error, (float, int)):
	    self.error = np.ones(len(wave))*error
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
	    raise ValueError('In the korel regime, each spectrum must be assigned component! ' \
			     'Currently it is %s' % str(self.component))
    def check_length(self):
	"""
	Checks that wavelengths and intensities have the same length.
	"""
	if len(self.wave) != len(self.intens):
	    raise ValueError('Wavelength vector and intensity vector do not have the same length!')
	
	  
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
      
    def get_spectrum(self):
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
	    return self.wave.copy(), self.intens.copy(), self.error.copy()
	  
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
	self.step = np.mean(self.wave[1:]-self.wave[:-1])	  
	
	 
    def read_spectrum_from_file(self, filename, global_error=None):
	"""
	Reads the spectrum from a file. Following format
	is assumed: %f %f %f (wavelength, intensity, error).
	If user does not provide errors, we still attempt
	to 
	Comments are denoted with '#'.
	INPUT:
	     filename..    the file, from which the spectra are loaded.
	     global_error..global error applies to every single observation 
	"""
	try:
	    # first we try to load 3 columns, i.e with errors
	    self.wave, self.intensity, self.error = np.loadtxt(filename, unpack=True, usecols=[0,1,2])
	    self.hasErrors = True
	except:
	  
	    # we failed, so we attempt to load two columns
	    self.wave, self.intensity = np.loadtxt(filename, unpack=True, usecols=[0,1])
	    
	    # error was not set up
	    if global_error is None:
	        warnings.warn("I found no errorbars of the observed intensities in file: %s! " \
			      "I assume they will be provided later. I remember!!" % (filename))
		self.hasErrors = False
	    # error was set up
	    else:
		self.error = global_error*np.ones(len(self.wave))
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
	    self.error = global_error*np.ones(len(self.wave))
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
	