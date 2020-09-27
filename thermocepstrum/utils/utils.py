# -*- coding: utf-8 -*-


class PrintMethod:

    _print_func = None
    _METHOD = 'bash'

    def __init__(self):
        pass

    @classmethod
    def write_log(cls, *args, **kwargs):
        if cls._METHOD == 'bash':
            print(*args, **kwargs)
        elif cls._METHOD == 'other':
            if cls._print_func:
                cls._print_func(*args, **kwargs)
            else:
                print(*args, **kwargs)
        else:
            print(*args, **kwargs)

    @classmethod
    def set_func(cls, func):
        cls._print_func = func

    @classmethod
    def set_method(cls, method):
        cls._METHOD = method
