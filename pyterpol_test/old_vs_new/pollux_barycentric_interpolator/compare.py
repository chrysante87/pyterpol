import sys
import numpy as np
import matplotlib.pyplot as plt

speclist = np.loadtxt(sys.argv[1], dtype=str)
plt.figure(figsize=(10,12), dpi=100)
for i,spec in enumerate(speclist):
    new = np.loadtxt(spec)
    old = np.loadtxt(spec[:-3]+'old')
    
    ind = np.where((old[:,0] >= new[:,0].min()) & (old[:,0] <= new[:,0].max()))[0]
    old = old[ind]
    
    #print old
    
    plt.subplot(211)
    plt.plot(new[:,0], new[:,1]+i*0.5, 'k-', label=spec[:-4])
    plt.plot(old[:,0], old[:,1]+i*0.5, 'r-', label=spec[:-4])
    plt.xlim(new[:,0].min(),new[:,0].max())
    plt.ylim(0.0,2.2)
    plt.legend()
    
    plt.subplot(212)
    plt.plot(new[:,0], old[:,1]-new[:,1], 'k-', label=spec[:-4])
    plt.xlim(new[:,0].min(),new[:,0].max())
    plt.legend()

plt.savefig('comparison.png')
plt.show()