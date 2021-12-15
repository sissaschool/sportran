# -*- coding: utf-8 -*-
"""
==========================================
    Sportran graphic user interface
==========================================
--------------------------------------------
    Utils file
--------------------------------------------

This file contains utility functions.
"""


class PrintMethod:
    """
    PrintMethod is used to define and change
    the output method.
    """

    _print_func = None
    _METHOD = 'bash'

    @classmethod
    def write_log(cls, s, s2='', *args, **kwargs):
        if cls._METHOD == 'bash':
            print(s, *args, **kwargs)
        elif cls._METHOD == 'other':
            if cls._print_func:
                cls._print_func(s, s2, *args, **kwargs)
            else:
                print(s, s2, *args, **kwargs)
        else:
            print(s, s2, *args, **kwargs)

    @classmethod
    def set_func(cls, func):
        """
        This function sets the output function.
        :param func: The output function to be used.
                      It must have at least one input string.

                      Will be called in this way
                      _print_func(s, s2, *args, **kwargs)
        """
        cls._print_func = func

    @classmethod
    def set_method(cls, method):
        """
        This function sets the output method.
        Two output methods are allowed: bash and other.
        The default value is bash that uses the normal print() function.

        :param method: the method used to output the data.
            bash - prints the output on the console using the print built-in function.
            other -  allows you to define a custom output function using set_func()
        """

        if method is 'bash' or method is 'other':
            cls._METHOD = method
        else:
            raise ValueError('Invalid input method')
