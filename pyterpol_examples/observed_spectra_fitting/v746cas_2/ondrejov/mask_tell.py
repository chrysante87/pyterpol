import sys
import numpy as np

tellbase = [
    [6522., 6525.5],
    [6530., 6538.5],
    [6541.8, 6550.37],
    [6551.75, 6554.9],
    [6557., 6560],
    [6563.6, 6564.8],
    [6568.38, 6576.3],
    [6580.2, 6588.2],
    [6594.2, 6596.],
    [6598.8, 6603.4]
]

def remove_telluric(f):
    """
    Removes intervals defined in tellbase
    :param f:
    :return:
    """

    # load the data
    w,i = np.loadtxt(f, unpack=True, usecols=[0,1])

    #remove wavelength intervals line by line
    for lim in tellbase:
        ind = np.where((w <= lim[0]) | (w >= lim[1]))[0]

        w = w[ind]
        i = i[ind]

    np.savetxt(f, np.column_stack([w,i]), fmt='%12.6f')


def main():
    f = sys.argv[1]
    remove_telluric(f)

if __name__ == '__main__':
    main()




