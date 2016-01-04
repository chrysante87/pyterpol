import pyterpol

# check convergence of individual parameters
pyterpol.Interface.plot_convergence_mcmc('chain.dat', figname='mcmc_convergence.png')

# plot covariance of radiative parameters
pyterpol.Interface.plot_covariances_mcmc('chain.dat', parameters=['vrot', 'teff', 'logg'], figname='mcmc_correlations.png')

# plot variance of rvs
pyterpol.Interface.plot_variances_mcmc('chain.dat', parameters=['rv'], figname='rv_var')

# write result
pyterpol.Interface.write_mc_result('chain.dat', outputname='mcmc.res')