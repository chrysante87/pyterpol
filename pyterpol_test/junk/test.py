import numpy as np
from scipy.interpolate import splrep, splev, interp1d
import matplotlib.pyplot as plt

x = np.linspace(0, 2*np.pi, 10)
y = np.sin(x)
xnew = np.linspace(0, 2*np.pi, 100)

prec = np.sin(xnew)

tck = splrep(x, y, k=3)
int_splev = splev(xnew, tck)

f = interp1d(x, y, kind='cubic')
int_interp1d = f(xnew)

plt.subplot(211)
plt.plot(x, y, 'ro', label='original')
plt.plot(xnew, prec, 'r-', label='precise')
plt.plot(xnew, int_splev, 'b-', label='splev')
plt.plot(xnew, int_interp1d, 'm-', label='interp1d')
plt.legend()

plt.subplot(212)
plt.plot(xnew, prec-int_splev, 'b-', label='splev')
plt.plot(xnew, prec-int_interp1d, 'm-', label='interp1d')
plt.legend()
plt.show()