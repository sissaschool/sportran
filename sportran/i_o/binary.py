# -*- coding: utf-8 -*-

from sportran.utils.attributedict import AttributeDict
from sportran.current import Current
import pickle

__all__ = ['SportranBinaryFile', 'SportranInput', 'SportranOutput', 'SportranSettings']


def multi_pickle_dump(self, filename, **objs):
    # dump a dictionary of all the passed objects
    pickle.dump(objs, open(filename, 'wb'))


def multi_pickle_load(self, filename):
    # return a dictionary of all the stored objects
    return pickle.load(open(filename, 'rb'))


class SportranInput(AttributeDict):

    pass


class SportranOutput(AttributeDict):

    pass


class SportranSettings(AttributeDict):

    pass


class SportranBinaryFile():
    """
    A SporTran Binary File object, that stores information about a SporTran calculation, its inputs, data, and outputs.
    The data stored can be used to restore a calculation that was previously run.

    Attributes:
        input_parameters    a SportranInput object, the input parameters of the calculation
        settings            a SportranSettings object, the calculation settings
        current             a Current object, the original current
        current_resampled   a Current object, the current obtained by resampling current
        output_results      a SportranOutput object, the calculation outputs

    Methods:
        :classmethod: load(filename)
        :method: dump(filename)
    """

    _storable_attrs = ('input_parameters', 'settings', 'current', 'current_resampled', 'output_results')

    def __init__(self, **kwargs):
        for key in self._storable_attrs:
            self.__setattr__(key, kwargs.pop(key, None))
        if kwargs:
            raise ValueError('Keys {} are not valid.'.format(kwargs.keys()))

    @classmethod
    def load(cls, filename):
        """Load SporTran data from a binary file."""
        return cls(**multi_pickle_load(filename))

    def dump(self, filename):
        """Dump SporTran data to a binary file."""
        multi_pickle_dump(filename, **{key: self.key for key in self._storable_attrs})
        print(f'Binary file "{filename}" successfully created, with the following data:')
        print([key for key in self._storable_attrs if key is not None])

    @property
    def input_parameters(self):
        return self._input_parameters

    @input_parameters.setter
    def input_parameters(self, value):
        if value is None:
            self._input_parameters = None
        elif isinstance(value, SportranInput):
            self._input_parameters = value
        else:
            raise TypeError()

    @property
    def settings(self):
        return self._settings

    @settings.setter
    def settings(self, value):
        if value is None:
            self._settings = None
        elif isinstance(value, SportranSettings):
            self._settings = value
        else:
            raise TypeError()

    @property
    def current(self):
        return self._current

    @current.setter
    def current(self, value):
        if value is None:
            self._current = None
        elif isinstance(value, Current):
            self._current = value
        else:
            raise TypeError()

    @property
    def current_resampled(self):
        return self._current_resampled

    @current_resampled.setter
    def current_resampled(self, value):
        if value is None:
            self._current_resampled = None
        elif isinstance(value, Current):
            self._current_resampled = value
        else:
            raise TypeError()

    @property
    def output_results(self):
        return self._output_results

    @output_results.setter
    def output_results(self, value):
        if value is None:
            self._output_results = None
        elif isinstance(value, SportranOutput):
            self._output_results = value
        else:
            raise TypeError()
