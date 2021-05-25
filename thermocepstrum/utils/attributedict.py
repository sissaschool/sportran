# -*- coding: utf-8 -*-


class AttributeDict(dict):   # pylint: disable=too-many-instance-attributes
    """
    This class internally stores values in a dictionary, but exposes
    the keys also as attributes, i.e. asking for attrdict.key
    will return the value of attrdict['key'] and so on.

    Raises an AttributeError if the key does not exist, when called as an attribute,
    while the usual KeyError if the key does not exist and the dictionary syntax is
    used.
    """

    def __init__(self, dictionary=None):
        """Recursively turn the `dict` and all its nested dictionaries into `AttributeDict` instance."""
        super().__init__()
        if dictionary is None:
            dictionary = {}

        for key, value in dictionary.items():
            if isinstance(value, Mapping):
                self[key] = AttributeDict(value)
            else:
                self[key] = value

    def __repr__(self):
        """Representation of the object."""
        return f'{self.__class__.__name__}({dict.__repr__(self)})'

    def __getattr__(self, attr):
        """Read a key as an attribute.

        :raises AttributeError: if the attribute does not correspond to an existing key.
        """
        try:
            return self[attr]
        except KeyError:
            errmsg = f"'{self.__class__.__name__}' object has no attribute '{attr}'"
            raise AttributeError(errmsg)

    def __setattr__(self, attr, value):
        """Set a key as an attribute."""
        try:
            self[attr] = value
        except KeyError:
            raise AttributeError(
                f"AttributeError: '{attr}' is not a valid attribute of the object '{self.__class__.__name__}'")

    def __delattr__(self, attr):
        """Delete a key as an attribute.

        :raises AttributeError: if the attribute does not correspond to an existing key.
        """
        try:
            del self[attr]
        except KeyError:
            errmsg = f"'{self.__class__.__name__}' object has no attribute '{attr}'"
            raise AttributeError(errmsg)

    def __deepcopy__(self, memo=None):
        """Deep copy."""
        from copy import deepcopy

        if memo is None:
            memo = {}
        retval = deepcopy(dict(self))
        return self.__class__(retval)

    def __getstate__(self):
        """Needed for pickling this class."""
        return self.__dict__.copy()

    def __setstate__(self, dictionary):
        """Needed for pickling this class."""
        self.__dict__.update(dictionary)

    def __dir__(self):
        return self.keys()
