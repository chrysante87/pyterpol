from scipy.optimize import fmin
from pyterpol.synthetic.auxiliary import parlist_to_list

fitters = dict(
    np_nelder_mead=dict(par0type='value',
                        optional_kwargs=['xtol', 'ftol', 'maxiter', 'maxfun'],
                        object=fmin)
)

class Fitter(object):
    """
    """
    def __init__(self, fitparams=None, verbose=False, debug=False, output='fit.log'):
        """
        :param fitparams a list of Parameter types
        :param verbose whether to save detailed chi_square information
        :param debug
        :return:
        """

        # pass the parameters
        self.fitparams = fitparams
        self.verbose = verbose
        self.debug = debug

        # empty parameters
        self.fitter = None
        self.fit_kwargs = {}
        self.par0 = []

        # empty list of all trial fits
        self.iters = []

    def __call__(self, func, **args, **kwargs):
        """
        :param func:
        :param args:
        :param kwargs:
        :return:
        """
        # reset the counter
        self.iter_number = 0

        # debug
        if self.debug:
            print "Started fitted with fitting environment: %s\n" \
                  " vector of parameters: %s and optional" \
                  " enviromental parameters: %s." % (self.fittername, str(self.par0), str(self.fit_kwargs))

        if len(self.par0) == 0:
            raise ValueError('No initial vector of parameters (wrapped in Parameter class) was passeed.')

        # run fitting
        self.result = self.fitter(func, self.par0, **self.fit_kwargs)

    def choose_fitter(self, name, fitparams=None, **kwargs):
        """
        Prepares the variables for the fitting
        :param name:
        :param kwargs:
        :return:
        """

        # check the input
        if name.lower() not in fitters:
            raise ValueError('Fitter: %s is unknown.' % name)
        else:
            self.fitter = fitters[name].object
            self.fittername = name
        for key in kwargs.keys():
            if key not in fitters[name]['optional_kwargs']:
                raise KeyError('The parameter: %s is not listed among '
                               'optional_kwargs for fitter: %s. The eligible'
                               'optional_kwargs are: %s'  % (key, name, str(fitters[name]['optional_kwargs'])))
            else:
                self.fit_kwargs[key] = kwargs[key]

        if self.debug:
            print 'Choosing environment: %s\n' \
                  ' environmental parameters: %s.' % (name, str(self.fit_kwargs))

        # if we want to change the fitted parameters
        if fitparams is None:
            fitparams = self.fitparams

        if fitters[name]['par0type'] == 'value':
            self.par0 = parlist_to_list(fitparams, property='value')
        if fitters[name]['par0type'] == 'limit':
            vmins = parlist_to_list(fitparams, property='vmin')
            vmaxs = parlist_to_list(fitparams, property='vmax')
            self.par0 = [[vmin, vmax] for vmin, vmax in zip(vmins, vmaxs)]

        if self.debug:
            print 'Setting initial parameters: %s' % str(self.par0)

    def append_iteration(self, iter):
        """
        Appends each iteration.
        :return:
        """
        self.iter_number += 1
        self.iters.append(iter)














