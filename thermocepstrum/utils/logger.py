# -*- coding: utf-8 -*-


class PrintMethod:
    """
    This class is intended to manage all the messages that should be put in a file,
    in the stdout with print, or in any other place. Is has to be used as a global thing.
    It is just a routing point. In the whole library please use this.
    By default it is equivalent to a call to print()
    """

    _print_func = None   #: print function that can be called by :func:`write_log`
    _METHOD = ['bash']   #: specify the list of methods to be called by :func:`write_log`

    def __init__(self):
        pass

    @classmethod
    def open_file(cls, fname):
        """opens the log file that will be used globally"""
        cls.lfile = open(fname, 'w')

    @classmethod
    def close_file(cls):
        cls.lfile.close()

    @classmethod
    def write_log(cls, *args, **kwargs):
        """
        Calls all the methods added by :func:`append_method`
        """
        if 'bash' in cls._METHOD:
            print(*args, **kwargs)
        if 'file' in cls._METHOD:
            s = ''
            for a in args:
                s += str(a)
            cls.lfile.write(s + '\n')
        if 'other' in cls._METHOD:
            if cls._print_func:
                cls._print_func(*args, **kwargs)
            else:
                print(*args, **kwargs)

    @classmethod
    def set_func(cls, func):
        """Set the function to call when :func:`write_log` is called.
        The function is called only if the 'other' method is setted by :func:`append_method`
        or :func:`set_method`"""
        cls._print_func = func

    @classmethod
    def append_method(cls, method):
        """append the method to the list.
           :param method: the method to be added to the method list, can be any of 'bash', 'file' or 'other'
           :type method: str
           if 'file' remember to call :func:`open_file` somewhere
        """
        cls._METHOD.append(method)

    @classmethod
    def set_method(cls, method):
        """Removes all method and set only the provided one"""
        cls._METHOD = [method]
