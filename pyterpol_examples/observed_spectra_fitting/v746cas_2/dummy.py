import pyterpol
# pyterpol.Interface.plot_convergence_mcmc('chain.sav.dat', figname='mcmc_convergence.png')
# pyterpol.Interface.plot_covariances_mcmc('chain.sav.dat', parameters=['vrot', 'teff', 'logg'], figname='mcmc_correlations.png')
# pyterpol.Interface.plot_variances_mcmc('chain.sav.dat', parameters=['rv'], groups=[0, 1, 19, 20, 21], figname='rv_var')
pyterpol.Interface.write_mc_result('chain.sav.dat', outputname='mcmc.res')