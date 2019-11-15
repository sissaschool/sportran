from __future__ import print_function


class PrintMethod:

    _print_func = None
    _METHOD = 'bash'

    def __init__(self):
        pass

    @classmethod
    def write_log(cls, s, *args, **kwargs):
        if cls._METHOD == 'bash':
            print(s, *args, **kwargs)
        elif cls._METHOD == 'other':
            if cls._print_func:
                cls._print_func(s, *args, **kwargs)
            else:
                print(s, *args, **kwargs)
        else:
            print(s, *args, **kwargs)

    @classmethod
    def set_func(cls, func):
        cls._print_func = func

    @classmethod
    def set_method(cls, method):
        cls._METHOD = method
