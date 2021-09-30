# -*- coding: utf-8 -*-


class PrintMethod:

    _print_func = None
    _METHOD = ['bash']

    def __init__(self):
        pass

    @classmethod
    def open_file(cls, fname):
        cls.lfile = open(fname, 'w')

    @classmethod
    def close_file(cls):
        cls.lfile.close()

    @classmethod
    def write_log(cls, *args, **kwargs):
        if 'bash' in cls._METHOD:
            print(*args, **kwargs)
        if 'file' in cls._METHOD:
            s = ''
            for a in args:
                s += str(a)
            cls.lfile.write(s+'\n')
        if 'other' in cls._METHOD:
            if cls._print_func:
                cls._print_func(*args, **kwargs)
            else:
                print(*args, **kwargs)

    @classmethod
    def set_func(cls, func):
        cls._print_func = func

    @classmethod
    def append_method(cls, method):
        cls._METHOD.append(method)

    @classmethod
    def set_method(cls, method):
        cls._METHOD = [method]
