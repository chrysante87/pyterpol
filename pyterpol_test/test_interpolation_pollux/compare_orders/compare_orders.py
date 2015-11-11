import sys
import numpy as np
import matplotlib.pyplot as plt


speclist = np.loadtxt(sys.argv[1], dtype=str)
plt.figure(figsize=(10,15), dpi=100)
for i,spec in enumerate(speclist):
    new = np.loadtxt(spec)
    
    if i == 0:
	fst = new.copy()
	old = np.loadtxt('LOGG_3.8_Z_1.0_TEFF_12200.old')
	ind = np.where((old[:,0] >= new[:,0].min()) & (old[:,0] <= new[:,0].max()))[0]
	old = old[ind]
    
    color = np.random.random(3)
    
    plt.subplot(311)
    plt.plot(new[:,0], new[:,1]+i*0.1, '-', label=spec[:-4], color=color)
    plt.xlim(new[:,0].min(),new[:,0].max())
    plt.ylim(0.0,2.2)
    plt.legend()
    
    plt.subplot(312)
    plt.plot(new[:,0], fst[:,1]-new[:,1], '-', label=spec[:-4], color=color)
    plt.xlim(new[:,0].min(),new[:,0].max())
    plt.legend()
    
    plt.subplot(313)
    plt.plot(new[:,0], old[:,1]-new[:,1], '-', label=spec[:-4], color=color)
    plt.xlim(new[:,0].min(),new[:,0].max())
    plt.legend()

plt.savefig('comparison.png')
plt.show()