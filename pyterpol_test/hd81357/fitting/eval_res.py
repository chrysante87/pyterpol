import numpy as np
import matplotlib.pyplot as plt
import pyterpol

itf = pyterpol.Interface()
itf.load('hd81357.sav')
itf.populate_comparisons()
# print itf.get_degrees_of_freedom()
# print itf.compute_chi2_treshold()

# try to load the data
data = pyterpol.read_fitlog('fit.log')

# try to plot the convergence
itf.plot_convergence(parameter='all')

# try to plot covariances
itf.plot_covariances(nbin=10, parameters=['teff', 'vrot', 'lr'])




