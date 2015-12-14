"""
This shows how to evaluate the outcome of the fitting
We were fitting the disentangled spectra of the secondary.
"""

import pyterpol

# 1) Load the last session - create an empty Interface
itf = pyterpol.Interface()
# fill it with teh last session
itf.load('hd81357.sav')

# 2) Have a look at the comparisons
itf.plot_all_comparisons(figname='finalfit')

# 3) Export the disentangle spectra
itf.write_synthetic_spectra(outputname='final_spectra')

# 4) Have a look how everything converged
itf.plot_convergence(figname='covergence_hd81357.png')

# 5) Have look at uncertainty of the fit
itf.plot_covariances(nbin=20)





