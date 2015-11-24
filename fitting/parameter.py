# definition of parameters - here I add parameters which apply for implemented
# grids. To fit additional parameters, one will have to define along with
# addition of new grids, or here..
parameter_definitions=dict(
    teff=dict(name='teff', value=10000., vmin=6000., vmax=50000., unit='K', fitted=False, group=None, typedef=(float)),
    logg=dict(name='logg', value=3.5, vmin=0.0, vmax=5.0, unit='log(g.cm^-2)', fitted=False, group=None, typedef=(float)),
    vrot=dict(name='vrot', value=0.0, vmin=0.0, vmax=500., unit='km.s^-1', fitted=False, group=None, typedef=(float)),
    rv=dict(name='rv', value=0.0, vmin=-1000., vmax=1000., unit='km.s^-1', fitted=False, group=None, typedef=(float)),
    lr=dict(name='lr', value=1.0, vmin=0.0, vmax=1.0, unit='relative', fitted=False, group=None, typedef=(float)),
    z=dict(name='z', value=1.0, vmin=0.0, vmax=2.0, unit='Z_solar', fitted=False, group=None, typedef=(float)),
)

class Parameter(object):
    """
    """
    def __init__(self, name=None, value=None, vmin=None, vmax=None, unit=None, fitted=None, group=None, typedef=None, debug=False):
        """
        :param name: name of the parameter
        :param value: of the parameter
        :param vmin: minimal value
        :param vmax: maximal value
        :param fitted: is optimized
        :param component: component to which it belongs
        :param group: group for which the parameter applies
        :return: None
        """
        # pass all arguments
        self.name = name
        self.value = value
        self.vmin = vmin
        self.vmax = vmax
        self.unit = unit
        self.fitted = fitted
        self.group = group
        self._typedef = typedef
        if typedef:
            for var in [value, vmin, vmax]:
                self.check_type(var)

        # define floats
        self._float_attributes = ['value', 'vmin', 'vmax']

        # define debug_mode
        self.debug = debug

    def __getitem__(self, item):
        """
        :param item: desired attribute of the class
        :return: value of the parameter
        """
        if hasattr(self, item):
            return getattr(self, item)
        else:
            raise AttributeError('The parameter %s has no attribute %s.' % (str(self), item))

    def __setitem__(self, item, value):
        """
        :param item: changed attribute
        :param value: new value
        :return:None
        """

        if hasattr(self, item):
            # if type was passed, than check them all the time they are changed
            if item in self._float_attributes and self._typedef is not None:
                self.check_type(value)
            setattr(self, item, value)
        else:
            raise AttributeError('The parameter %s has no attribute %s.' % (str(self), item))

    def __str__(self):
        """
        :return: string represantation of the class
        """
        string = ''
        for var in ['name', 'value', 'vmin', 'vmax', 'fitted', 'group', '_typedef']:
            string += "%s: %s " % (var, str(getattr(self, var)))
        string += '\n'

        return string

    def check_type(self, value):
        """
        :param value: value to be checked
        :return: bool whether the tested value has correct type
        """
        if self._typedef is None:
            raise ValueError('The type of parameter %s has not been set.' % str(self))
        else:
            if not isinstance(value, self._typedef):
                raise TypeError('The passed value has not correct type.'
                                ' Correct type is %s.' % str(self._typedef))

    def get_range(self):
        """
        Returns boundaries of the fitting.
        :return: list, boundaries
        """
        return [self['vmin'], self['vmax']]

    def get_error(self, relative=None):
        """
        :param relative: relative error desired by the user
        :return: error - uncertainty of the parameter.,
        """
        # if relative error is not give, the uncertainty is
        # taken from boundaries
        if relative is None:
            error = (self['vmin'] + self['vmax']) / 2.

        # otherwise it is a fraction of the value
        else:
            error = relative * self['value']

        return error








