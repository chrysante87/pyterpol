import numpy as np
import pyterpol

def load_observations(f):
    """
    :param f: file
    :return:
    """

    # load the observations
    flist = np.loadtxt(f, usecols=[0], unpack=True, dtype=str)
    hjd = np.loadtxt(f, usecols=[1], unpack=True).tolist()
    hjd[0] = None

    # create list of observations
    obs = []
    for i, sf in enumerate(flist):

        # wrap the spectrum into observed spectrum class
        # o = pyterpol.ObservedSpectrum(filename=sf, group=dict(rv=0))
        # o = pyterpol.ObservedSpectrum(filename=sf, group=dict(rv=0), hjd=hjd[i])
        o = pyterpol.ObservedSpectrum(filename=sf, group=dict(rv=i), hjd=hjd[i])

        # estimate uncertainty from continuum
        o.get_sigma_from_continuum(6620., 6630.)
        obs.append(o)

    # create ObservedList
    ol = pyterpol.ObservedList()
    ol.add_observations(obs)

    return ol, flist, hjd

def main():
    """
    :return:
    """
    # parameters
    niter = 2

    # 1) Generate region
    rl = pyterpol.RegionList()
    rl.add_region(wmin=6337., wmax=6410.0, groups=dict(lr=0))
    rl.add_region(wmin=6530., wmax=6600.0, groups=dict(lr=0))
    rl.add_region(wmin=6660., wmax=6690.0, groups=dict(lr=0))

    # 2) Load observed data
    ol = load_observations('prekor.lst')[0]
        
    ## 3) Generate components
    sl = pyterpol.StarList()
    sl.add_component('primary',   teff=16000.,  logg=4.285, lr=1.0, vrot=90., z=1.0)   

    ## 4) construct the interface
    itf = pyterpol.Interface(sl=sl, rl=rl, ol=ol)
    itf.set_grid_properties(order=4, step=0.05)
    itf.setup()

    print itf

    ## 5) write rvs
    itf.save('test.sav')
    itf.write_rvs('test.rv.dat')
    
    # 6) try to load it 
    itf.load('test.sav')
    
    # 7) and save it again
    itf.save('test2.sav')

if __name__ == '__main__':
    main()
