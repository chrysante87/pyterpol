import pyterpol
pyterpol.Interface.plot_convergence_mcmc('chain.sav.dat', figname='mcmc_convergence.png')
pyterpol.Interface.plot_covariances_mcmc('chain.sav.dat', parameters=['vrot', 'teff', 'logg'], figname='mcmc_correlations.png')
pyterpol.Interface.plot_variances_mcmc('chain.sav.dat', parameters=['rv'], figname='rv_var')
pyterpol.Interface.write_mc_result('chain.sav.dat', outputname='mcmc.res')