import sys
sys.path.append('/home/jana/work/pyterpol')
import pyterpol
import numpy as np

obsfile = 'ua170034.asc'
w, i = np.loadtxt(obsfile, unpack=True, usecols=[0,1])
e = np.random.random(len(w))*0.01

# No spectrum is passes - only some warnings are issues -checked
#os = pyterpol.ObservedSpectrum()

# spectrum is passed - only warning regarding errors qre issued - checked
#os = pyterpol.ObservedSpectrum(wave=w, intens=i)

# spectrum is passed with filename -- warnings regarding error bars issued
#os = pyterpol.ObservedSpectrum(filename=obsfile)

# errors are also pass - locally and globally - checked
#os = pyterpol.ObservedSpectrum(wave=w, intens=i, error=e)
#print os.error
#os = pyterpol.ObservedSpectrum(wave=w, intens=i, error=0.01)
#print os.error

# errors are passed after creation of grid - checked
#os = pyterpol.ObservedSpectrum(wave=w, intens=i)
#os.set_error(vec_error=e)
#print os.error
#os.set_error(global_error=0.01)
#print os.error

# everything is passed after creation of type - checked
#os = pyterpol.ObservedSpectrum()
#os.set_spectrum_from_arrays(w,i,e)
#print os.error

# check that the spectrum is measured properly - checked
#print w[0], w[-1], os.get_boundaries()

# free the spectrum - checked
#os.free_spectrum()
#print os.wave, os.intens, os.error, os.loaded

# tests some precautions
# what if the intensities and wavelengths do not have the same length -checked
#os.set_spectrum_from_arrays(w[:-3],i,e)

# what if we set up korel, b ut do not set up component - checked
#os = pyterpol.ObservedSpectrum(wave=w, intens=i, error=e, korel=True)

# try to do estimate of the error
os = pyterpol.ObservedSpectrum(wave=w, intens=i, error=e)
err_cont = os.get_sigma_from_continuum(4453, 4459)
print err_cont, os.error
err_fft = os.get_sigma_from_fft(nlast=50)
print err_fft, os.error



