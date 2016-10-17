import pyterpol

def main():
    """
    :return:
    """

    ## 1) Generate region
    rl = pyterpol.RegionList()
    rl.add_region(wmin=6337., wmax=6410.0, groups=dict(lr=0))
    rl.add_region(wmin=6530., wmax=6600.0, groups=dict(lr=0))
    rl.add_region(wmin=6660., wmax=6690.0, groups=dict(lr=0))

    ## 2) Load observed data
    # first wrap them into observed spectrum class
    o1 = pyterpol.ObservedSpectrum(filename='blb00001.clean.asc', instrumental_width=0.46, group=dict(rv=0))
    o2 = pyterpol.ObservedSpectrum(filename='blb00002.clean.asc', instrumental_width=0.46, group=dict(rv=1))
    o1.get_sigma_from_continuum(6630., 6640.)
    o2.get_sigma_from_continuum(6630., 6640.)

    # create list of observed spectra
    ol = pyterpol.ObservedList()
    ol.add_observations([o1, o2])

    # alternatively
    ol = pyterpol.ObservedList()
    obs = [dict(filename='blb00001.clean.asc', instrumental_width=0.46, group=dict(rv=0), error=0.01),
           dict(filename='blb00002.clean.asc', instrumental_width=0.46, group=dict(rv=1), error=0.01)]
    ol.add_observations(obs)
        
    ## 3) Generate components
    sl = pyterpol.StarList()
    sl.add_component('primary', teff=16000.,  logg=4.285, lr=1.0, vrot=90., z=1.0)   

    ## 4) construct the interface
    itf = pyterpol.Interface(sl=sl, rl=rl, ol=ol)
    itf.set_grid_properties(order=4, step=0.05)
    itf.setup()

    ## 5) Set parameters for fitting
    itf.set_parameter(parname='teff', vmin=14000., vmax=17000., fitted=True)
    itf.set_parameter(parname='logg', vmin=3.5, vmax=4.5, fitted=True)
    itf.set_parameter(parname='rv', vmin=-30., vmax=0., fitted=True)
    itf.set_parameter(parname='vrot', vmin=70., vmax=110., fitted=True)
    itf.set_parameter(parname='lr', vmin=0.99, vmax=1.01, fitted=True)

    ## 6) Choose fitting environment
    itf.choose_fitter('nlopt_nelder_mead', ftol=1e-4)

    ## 7) Run the fitting
    itf.run_fit()

    # write down rvs
    itf.write_rvs('hd.rvs')

    ## 8) plot and save results
    itf.plot_all_comparisons(figname='final')
    itf.save('test.final.sav')


if __name__ == '__main__':
    main()
