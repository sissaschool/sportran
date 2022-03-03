# -*- coding: utf-8 -*-

from sportran_gui.utils.custom_widgets import *
from sportran_gui.core import control_unit as cu

INDENT = 0


def print_name(func):

    def inner(*args, **kwargs):
        global INDENT
        print('{}BEGIN {}'.format(' ' * INDENT, func.__name__))
        INDENT = INDENT + 1
        func(*args, **kwargs)
        INDENT = INDENT - 1
        print('{}END   {}'.format(' ' * INDENT, func.__name__))

    return inner


class PStarSelector(Frame):

    def __init__(self, parent, main):
        Frame.__init__(self, parent)

        self.main = main

        self.next_frame = None
        self.prev_frame = None

        self.parent = parent
        self.main_frame_scroll = ScrollFrame(self, self)
        self.main_frame = self.main_frame_scroll.viewPort

        #self.main_frame.grid(column=0, row=0, sticky='nsew')

        sections = Frame(self.main_frame)
        sections.grid(row=0, column=0, sticky='nsew', padx=20, columnspan=2)

        self.graph = GraphWidget(sections, sections, size=(7, 4), toolbar=True)

        self.container_frame = Frame(self.main_frame)
        self.container_frame.grid(row=1, column=0, sticky='nsew')

        variable_frame = Frame(self.container_frame, bd=1, relief=SOLID)
        variable_frame.pack(side=TOP, anchor='w', padx=20, fill='x', expand=1, pady=20)

        self.fstar_label = Label(variable_frame, text='', font=('Arial 24'))
        self.fstar_label.grid(row=0, column=0, sticky='w')
        self.kmin_label = Label(variable_frame, text='', font=('Arial 24'))
        self.kmin_label.grid(row=1, column=0, sticky='w')

        value_frame = Frame(self.container_frame)
        value_frame.pack(side=TOP, anchor='w', padx=20, fill=BOTH, expand=1)

        Label(value_frame, text=LANGUAGES[settings.LANGUAGE]['stp5'],
              font='Arial 12 bold').grid(row=0, column=0, sticky='w')
        ttk.Separator(value_frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we', pady=10, columnspan=5)

        Label(value_frame, text=LANGUAGES[settings.LANGUAGE]['pstar_cmt'],
              font='Arial 12').grid(row=2, column=0, sticky='w')
        Label(value_frame, text='P*: ', font='Arial 12').grid(row=3, column=0, sticky='w')
        self.value_entry = Spinbox(value_frame, bd=1, relief=SOLID, increment=1)
        self.value_entry.grid(row=3, column=1, sticky='w')

        self.increment = IntVar()
        Label(value_frame, text=LANGUAGES[settings.LANGUAGE]['inc'],
              font='Arial 12').grid(row=4, column=0, sticky='w', pady=10)

        rdbt_frame = Frame(value_frame)
        rdbt_frame.grid(row=4, column=1, sticky='w')

        Radiobutton(rdbt_frame, text='1', font='Arial 11 bold', variable=self.increment, value=1,
                    command=self._change_increment).pack(side=LEFT)
        Radiobutton(rdbt_frame, text='10', font='Arial 11 bold', variable=self.increment, value=10,
                    command=self._change_increment).pack(side=LEFT)
        Radiobutton(rdbt_frame, text='100', font='Arial 11 bold', variable=self.increment, value=100,
                    command=self._change_increment).pack(side=LEFT)

        Button(value_frame, text=LANGUAGES[settings.LANGUAGE]['recalculate'], font='Arial 12 bold', bd=1, relief=SOLID,
               command=self._recalc, width=20).grid(row=2, column=2, sticky='wens', rowspan=2, padx=50)

        value_frame.columnconfigure(0, weight=1, minsize=110)
        value_frame.columnconfigure(1, weight=1, minsize=150)
        value_frame.columnconfigure(2, weight=1, minsize=1)

        self.main_frame.rowconfigure(1, weight=1)
        button_frame = Frame(self.main_frame)
        button_frame.grid(row=2, column=0, sticky='w', padx=10, pady=20)

        back_button = Button(button_frame, text=LANGUAGES[settings.LANGUAGE]['back'], bd=1, relief=SOLID,
                             command=lambda: self.back(), width=10)
        back_button.grid(row=0, column=0, sticky='we', padx=5)

        new_a = Button(button_frame, text=LANGUAGES[settings.LANGUAGE]['new_a'], bd=1, relief=SOLID,
                       command=lambda: cu.new(self.main.root), width=10)
        new_a.grid(row=0, column=1, sticky='we', padx=5)

        self.main_frame.columnconfigure(0, weight=1, minsize=500)

        self.setted = False

        if cu.info:
            cu.update_info(cu.info)

    def set_prev_frame(self, frame):
        self.prev_frame = frame

    def back(self):
        cu.data.pstar = int(self.value_entry.get())
        if self.prev_frame:
            self.main.show_frame(self.prev_frame)
        else:
            raise ValueError('Prev frame isn\'t defined')

    def _get_pstar(self, aic_type='aic', cutoffK=None):
        cu.data.xf.cepstral_analysis(aic_type=aic_type, manual_cutoffK=((cutoffK - 1) if cutoffK is not None else None))

    def _pstar(self):
        self.value_entry.config(from_=2, to=cu.data.xf.NFREQS)

        if cu.data.xf.cepf:
            self.value_entry.delete(0, END)
            self.value_entry.insert(0, (cu.data.xf.cepf.cutoffK + 1))
            self.fstar_label.config(text='f*: {:.3f}    P*: {}'.format(cu.data.fstar, cu.data.xf.cepf.cutoffK + 1))
            self.kmin_label.config(text=u'\u03f0: {:18f} +/- {:8f}  {}\n'.format(cu.data.xf.kappa, cu.data.xf.kappa_std,
                                                                                 cu.data.xf._KAPPA_SI_UNITS))

    def _change_increment(self):
        self.value_entry.config(increment=int(self.increment.get()))

    def _recalc(self):
        if self.value_entry.get():
            kmin_c = int(self.value_entry.get())
        else:
            kmin_c = 0
        self._get_pstar(aic_type='aic', cutoffK=kmin_c)
        self.graph.add_graph(cu.gm.plot_cepstral_spectrum, 'cepstral', mode='linear', current=cu.data.xf,
                             kappa_units=True)
        xf = cu.data.xf
        self.graph.update_cut()
        cu.data.xf = xf
        self._pstar()

    def _setup_pstar(self):
        cu.data.xf.cepstral_analysis(aic_type='aic', manual_cutoffK=None)
        self._pstar()
        cu.data.pstar = int(self.value_entry.get())

    def _draw_graph(self):
        self.graph.graph.clear()
        self.graph.show(cu.gm.plot_periodogram, current=cu.data.j, mode='linear', kappa_units=True)
        cu.data.xf = cu.data.j.resample(fstar_THz=cu.data.fstar, PSD_FILTER_W=cu.data.psd_filter_width, plot=False)
        self.graph.add_graph(cu.gm.plot_resample, 'resample', xf=cu.data.xf, mode='linear', x=cu.data.j,
                             PSD_FILTER_W=cu.data.psd_filter_width)

    def recalculate(self):
        self._setup_pstar()
        self._draw_graph()
        self._recalc()

        if cu.info:
            cu.update_info(cu.info)

    def update_data(self):
        pass

    def update(self):
        super().update()

        self.setted = cu.data.recalc_pstar

        if self.setted:
            self.recalculate()

        if cu.info:
            cu.update_info(cu.info)

        self._recalc()
        self.graph.update_cut()

        self.main_frame_scroll.update_view()
