# -*- coding: utf-8 -*-

from functools import wraps


def add_method(cls):
    """A decorator to dynamically add a method to a class."""

    def decorator(func):

        @wraps(func)
        def wrapper(self, *args, **kwargs):
            return func(self, *args, **kwargs)

        setattr(cls, func.__name__, wrapper)
        # Note we are not binding func, but wrapper which accepts self but does exactly the same as func
        return func   # returning func means func can still be used normally

    return decorator


## Example:
## let A be a class
##
## Non-decorator way (note the function must accept self)
## def foo(self):
##     print('hello world!')
## setattr(A, 'foo', foo)
#
## def bar(self, s):
##     print(f'Message: {s}')
## setattr(A, 'bar', bar)
#
## Decorator can be written to take normal functions and make them methods
#@add_method(A)
#def foo():
#    print('hello world!')
#
#@add_method(A)
#def bar(s):
#    print(f'Message: {s}')
