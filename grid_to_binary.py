import os 
import argparse
import numpy as np

def main():
    ps = argparse.ArgumentParser()
    ps.add_argument('--remove', action='store_true', default=False, help='Removes ascii files.')
    ps.add_argument('--overwrite', action='store_true', default=False, help='Overwrites binary files -- mandatory for every machine swap. ')
    args = ps.parse_args()
    print args
    
    # get grids directory
    for gdname in ['grids', 'grids_ABS']:
        cwd = os.getcwd()
        gdir = os.path.join(cwd, gdname)
        dl = os.listdir(gdir)
        
        # go through each grid directory
        for direc in dl:
            # path to directory
            path = os.path.join(gdir, direc)
            
            # directories only
            if not os.path.isdir(path):
                continue
         
            # list of spectra
            gl = os.path.join(path, 'gridlist')
            
            # load the list
            synlist = np.loadtxt(gl, dtype=str, unpack=True, usecols=[0])
            
            # transform each spectrum to binary
            for synspec in synlist:
              
                # define name of the binary file
                bin_synspec = synspec + '.npz'
                if os.path.isfile(os.path.join(path, bin_synspec)) and not args.overwrite:
                    print "File: %s exists." % bin_synspec
                    if os.path.isfile(os.path.join(path, synspec)) and args.remove:
                        os.remove(os.path.join(path, synspec))       
                    continue
        
                # load the ascii spectrum and save it as binary file
                w, i = np.loadtxt(os.path.join(path, synspec), unpack=True, usecols=[0, 1])
                np.savez(os.path.join(path, bin_synspec), w, i)
                if os.path.isfile(os.path.join(path, synspec)) and args.remove:
                    os.remove(os.path.join(path, synspec))
            
if __name__ == '__main__':
    main()
    
    