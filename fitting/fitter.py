import os
import nlopt
import warnings
from scipy.optimize import fmin
from pyterpol.synthetic.auxiliary import parlist_to_list

fitters = dict(
    sp_nelder_mead=dict(par0type='value',
                        optional_kwargs=['xtol', 'ftol', 'maxiter', 'maxfun'],
                        object=fmin,
                        uses_bounds=False
                        ),
    nlopt_nelder_mead=dict(par0type='value',
                           optional_kwargs=['xtol', 'ftol', 'maxiter', 'maxfun'],
                           object = nlopt.opt,
                           environment = nlopt.LN_NELDERMEAD,
                           uses_bounds=True
                           ),
)

class Fitter(object):
    """
    """
    def __init__(self, fitparams=[], verbose=False, debug=False, fitlog='fit.log'):
        """
        :param fitparams a list of Parameter types
        :param verbose whether to save detailed chi_square information
        :param debug
        :return:
        """

        # pass the parameters
        self.fitparams = fitparams
        self.verbose = verbose
        self.fitlog = fitlog
        self.debug = debug

        # empty parameters
        self.fitter = None
        self.fittername = None
        self.fit_kwargs = {}
        self.par0 = []
        self.uses_bounds=False
        self.family=None
        self.vmins = None
        self.vmaxs = None
        self.nlopt_environment = None

        # empty list of all trial fits
        self.iters = []

        # clear the fitting log
        if os.path.isfile(fitlog):
            warnings.warn('A fitlog from previous fitting was found and overwritten..muhahahaha!')
            open(fitlog, 'w')


    def __call__(self, func, *args, **kwargs):
        """
        :param func:
        :param args:
        :param kwargs:
        :return:
        """
        # reset the counter and clear the fitting
        self.iter_number = 0
        self.iters = []

        # debug
        if self.debug:
            print "Started fitted with fitting environment: %s\n" \
                  " vector of parameters: %s and optional" \
                  " enviromental parameters: %s." % (self.fittername, str(self.par0), str(self.fit_kwargs))

        if len(self.par0) == 0:
            raise ValueError('No initial vector of parameters (wrapped in Parameter class) was passeed.')

        # run fitting
        self.result = self.fitter(func, self.par0, args=args, **self.fit_kwargs)

    def __str__(self):
        """
        String representation of the class.
        :return:
        """
        string = 'Initial parameters:'
        string += 'Fitter: %s optional_arguments: %s\n' % (self.fittername, str(self.fit_kwargs))
        for i, par in enumerate(self.fitparams):
            string += "(%s, group): (%s, %s); " % (par['name'], str(self.par0[i]), str(par['group']))
        string += '\n'

        return string


    def choose_fitter(self, name, fitparams=None, **kwargs):
        """
        Prepares the variables for the fitting
        :param name:
        :param kwargs:
        :return:
        """
        # print fitparams

        # check the input
        if name.lower() not in fitters:
            raise ValueError('Fitter: %s is unknown.' % name)
        else:
            self.fitter = fitters[name]['object']
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
        else:
            self.fitparams = fitparams

        # set up initial value
        if fitters[name]['par0type'] == 'value':
            self.par0 = parlist_to_list(fitparams, property='value')
        if fitters[name]['par0type'] == 'limit':
            vmins = parlist_to_list(fitparams, property='vmin')
            vmaxs = parlist_to_list(fitparams, property='vmax')
            self.par0 = [[vmin, vmax] for vmin, vmax in zip(vmins, vmaxs)]

        if self.debug:
            print 'Setting initial parameters: %s' % str(self.par0)

        # checks that there are any fitting boundaries
        if fitters[name]['uses_bounds']:
            self.uses_bounds = True
            self.vmins = vmins
            self.vmaxs = vmaxs

        # set up family
        self.family = name.split('_')[0]

        if self.family == 'nlopt':
            self.nlopt_environment = fitters[name]['environment']
            self.setup_nlopt()


    def append_iteration(self, iter):
        """
        Appends each iteration.
        :return:
        """
        self.iter_number += 1
        self.iters.append(iter)

        # if the number of iterations exceeds a certain number
        # they are written to a file
        if len(self.iters) > 1000:
            self.flush_iters()
            self.iters = []


    def flush_iters(self, f=None):
        """
        Flushes all records within self.iters to a file
        :param f:
        :return:
        """
        if f is None:
            f = self.fitlog

        # create a block of lines
        lines = []
        for row in self.iters:
            line = ''
            for key in row.keys():
                line += "%s: %s " % (key, str(row[key]))
            line += '\n'
            lines.append(line)

        # write the to a file
        ofile = open(f, 'a')
        ofile.writelines(lines)
        ofile.close()

    def setup_nlopt(self):
        """
        Sets up the the NLOPT fitter.
        :return:
        """

        if self.debug:
            print "Setting up NLOPT minimizer."

        # length of the fitted parameters
        n = len(self.fitparams)

        # configures the fitter
        self.fitter = self.fitter(self.nlopt_environment, n)

        # setup parameters for fitting terminatio
        for key in self.fit_kwargs.keys():
            if key == 'xtol':
                self.fitter.set_xtol_rel(self.fit_kwargs[key])
            if key == 'ftol':
                self.fitter.set_ftol_rel(self.fit_kwargs[key])
            if key == 'maxfun':
                self.fitter.set_maxeval(self.fit_kwargs[key])

        # setup boundaries
        self.fitter.set_lower_bounds(self.vmins)
        self.fitter.set_upper_bounds(self.vmaxs)

        # setup initial step
        stepsize = ((np.array(self.vmaxs) - np.array(self.vmins))/2.).tolist()
        self.fitter.set_initial_step(stepsize)



















