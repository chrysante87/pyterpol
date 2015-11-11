import warnings

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
	    component
	    korel
	"""
	
	# pass all arguments
	self.wave = wave
	self.intens = intens
	
	# lets have a look at the errors
	if error is None:
	    warnings.warn("I found no errorbars of the observed intensities in file: %s! " \
			  "I assume they will be provided later. I remember!!" % (filename))
	    self.hasErrors = False
	
	# Making assumption here, that we are
	# not passing errors, if we are not
	# passing wavelengths and intensities
	elif isinstance(error, float, int):
	    self.error = np.ones(len(wave))*error
	    self.hasErrors = True
	else:
	    self.error = error
	    self.hasErrors = True
	    
	
	# sets that the spectrum is loaded
	if (wave is not None) and (intens is not None):
	    self.loaded = True
	else:
	    self.loaded = False
	
	# if we provided the filename
	self.filename = filename
	if (not self.loaded) and (self.filename is not None):
	    self.read_spectrum_from_file(filename)
	
	# setup korel and check that it is proper
	self.component = component
	self.korel = korel
	self.check_korel()
	
    def check_korel(self):
	"""
	If korel is set, component must be set too.
	"""
	if (korel) and (self.component.upper() == 'ALL'):
	    raise ValueError('In the korel regime, each spectrum must be assigned component! ' \
			     'Currently it is %s' % str(self.component))
    def free_spectrum(self):
	"""
	Deletes the stored spectrum.
	"""
	self.wave = None
	self.intens = None
	self.error = None
	self.loaded = False
	 
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
	    self.wave, self.intensity = np.loadtxt(filename, unpack=True, usecols=[0,1,2])
	    
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
	
    def set_error(vec_error=None, glob_error=None):
	"""
	A tool to set the error.
	INPUT:
	    vec_error..  vector of errors, one for each intensity point
	    glob_error.. one global errror for each point
	""" 
	if vec_error is not None:
	    self.error = vec_error
	    self.hasErrors = True
	if glob_error is not None:
	    self.error = glob_error*(len(self.wave))
	    self.hasErrors = True
	