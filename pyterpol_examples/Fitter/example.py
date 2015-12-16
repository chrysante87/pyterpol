"""
This function demonstrates usage of the class Fitter.
"""

import pyterpol
import numpy as np
import matplotlib.pyplot as plt
import time
################################## BARE MINIMUM #########################################

# define a function that will be minimized
def func(x):
    """
    Polynomial of order 4
    :param x:
    :param p:
    :return:
    """
    x = x[0]
    return 0.5*x**4 - 2*x**3 - 5*x**2 + 12*x - 2

# create an empty fitter
fitter = pyterpol.Fitter()

# fitter is designed to work with sets of Parameter types.
# so we create one.
par = pyterpol.Parameter(name='x', value=-5., vmin=-100., vmax=100., fitted=True)

# What kind of fitter will we choose... lets have a look at
# the available ones. Note that optional arguments thatr
# control each fit are also listed. For detail have
# a look at the homepage of each environment
print fitter.list_fitters()

# nlopt_nelder_mead is always a good choice - note that
# this function expects a list of parameters not a single
# parameter, so we have to put it into brackets. Ftol
# sets that the fitting will end once the relative change
# of the cost function is less than 1e-6.
fitter.choose_fitter('nlopt_nelder_mead', fitparams=[par], ftol=1e-6)

# check the fitter
print fitter

# we can run the fitting by calling the fitter
t0 = time.time()
fitter(func)
dt1 = time.time() - t0

# have a look at the minimum and the value at minimum
print "func(%s) = %s" % (fitter.result, func(fitter.result))

# lets plot the function to be sure that we are in the global minimum
x = np.linspace(-10, 10, 200)
plt.plot(x, [func([v]) for v in x], 'k-', label='func(x)')
plt.plot(fitter.result, func(fitter.result), 'ro')
plt.ylim(-100, 100)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.savefig('result_nm.png')
plt.close()

# We see that the minimizer failed to find the global minimum
# it is not very unusual if we have similar minima.
# lets choose a more capable fitting environment
fitter.choose_fitter('sp_diff_evol', fitparams=[par])

# run the fitting
t0 = time.time()
fitter(func)
dt2 = time.time() - t0

# have a look at the minimum and the value at minimum
print "func(%s) = %s" % (fitter.result, func(fitter.result))

# lets plot the function to be sure that we are in the global minimum
x = np.linspace(-10, 10, 200)
plt.plot(x, [func([v]) for v in x], 'k-', label='func(x)')
plt.plot(fitter.result, func(fitter.result), 'ro')
plt.ylim(-100, 100)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.savefig('result_de.png')

# The simplex was faster, but converged only
# locally where the differential evolution
# converged correctly at the cost of ten times
# longer computation time
print "T(simplex) = %s" % str(dt1)
print "T(differential_evolution) = %s" % str(dt2)


