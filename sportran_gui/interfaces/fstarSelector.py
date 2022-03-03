# -*- coding: utf-8 -*-

from sportran_gui.utils.custom_widgets import *


class FStarSelector(Frame):

    def __init__(self, parent, main):
        Frame.__init__(self, parent)

        self.main = main

        self.next_frame = None
        self.prev_frame = None

        self.parent = parent
        self.main_frame_scroll = ScrollFrame(self, self)
        self.main_frame = self.main_frame_scroll.viewPort

        self.filter_width_value = DoubleVar()
        # self.main_frame.grid(column=0, row=0, sticky='nsew')

        self.sections = Frame(self.main_frame, pady=20)
        self.sections.grid(row=0, column=0, sticky='nsew')

        self.graph = GraphWidget(self.sections, self.sections, size=(7, 4), toolbar=True)

        self.slider_locked = False

        slider_frame = Frame(self.sections)
        slider_frame.pack(side=TOP, anchor='w', padx=20, fill=BOTH)

        Label(slider_frame, text=LANGUAGES[settings.LANGUAGE]['fslide'],
              font='Arial 12').grid(row=0, column=1, sticky='w', padx=20)

        self.slider = ttk.Scale(slider_frame, from_=0, to_=0.1)
        self.slider.grid(row=1, column=1, sticky='we', columnspan=1, padx=20, pady=5)
        slider_frame.columnconfigure(0, weight=2)
        slider_frame.columnconfigure(1, weight=11)
        slider_frame.columnconfigure(2, weight=2)

        slider_options_frame = Frame(slider_frame)
        slider_options_frame.grid(row=2, column=1, sticky='w', padx=20, pady=2)

        lock_slider = Button(slider_options_frame, command=lambda: self._lock_unlock_slider(), bd=1, relief=SOLID)
        lock_slider.grid(row=0, column=0, padx=2, sticky='w')

        self.change_view_button = Button(slider_options_frame, text='Zoom-in', command=lambda: self._change_view())
        self.change_view_button.grid(row=0, column=1, sticky='w', padx=2)

        value_frame = Frame(self.sections)
        value_frame.pack(side=TOP, pady=10, padx=20, fill=BOTH, expand=1)

        Label(value_frame, text=LANGUAGES[settings.LANGUAGE]['stp4'],
              font='Arial 12 bold').grid(row=0, column=0, sticky='w')

        ttk.Separator(value_frame, orient=HORIZONTAL).grid(row=0, column=0, sticky='we', columnspan=4, pady=5, padx=20)

        Label(value_frame, text=LANGUAGES[settings.LANGUAGE]['slct_v'] + ': ')\
            .grid(row=1, column=0, sticky='w', pady=4, padx=20)

        self.value_entry = Entry(value_frame, bd=1, relief=SOLID)
        self.value_entry.grid(row=1, column=1, sticky='we', padx=20)

        self.graph.attach_entry(self.value_entry)
        self.graph.attach_slider(self.slider)

        value_frame.columnconfigure(1, weight=1, minsize=150)
        value_frame.columnconfigure(2, weight=1, minsize=10)
        value_frame.columnconfigure(3, weight=1, minsize=300)
        value_frame.columnconfigure(4, weight=1, minsize=150)

        Label(value_frame, text=LANGUAGES[settings.LANGUAGE]['fl_w']+': ')\
            .grid(row=2, column=0, sticky='w', padx=20)
        self.filter_width = Spinbox(value_frame, from_=0.1, to=10, increment=0.1, bd=1, relief=SOLID,
                                    textvariable=self.filter_width_value)
        self.filter_width.grid(row=2, column=1, sticky='we', pady=10, padx=20)

        self.fstar_screen = Label(value_frame, text='F*: ', font='Arial 14 bold', width=20, bd=1, relief=SOLID)
        self.fstar_screen.grid(row=1, column=3, sticky='we', padx=50)

        Button(value_frame, text=LANGUAGES[settings.LANGUAGE]['resample'], font='Arial 12 bold', bd=1, relief=SOLID,
               command=self.resample, width=20).grid(row=2, column=3, sticky='wens', rowspan=1, padx=50)

        value_frame.rowconfigure(3, weight=1)
        button_frame = Frame(self.main_frame)
        button_frame.grid(row=3, column=0, padx=20, pady=10, sticky='w')

        back_button = Button(button_frame, text=LANGUAGES[settings.LANGUAGE]['back'], bd=1, relief=SOLID,
                             command=lambda: self.back(), width=10)
        back_button.grid(row=0, column=0, sticky='we', padx=5)

        next_button = Button(button_frame, text=LANGUAGES[settings.LANGUAGE]['next'], bd=1, relief=SOLID,
                             command=lambda: self.next(), width=10)
        next_button.grid(row=0, column=1, sticky='we', padx=5)

        self.main_frame.rowconfigure(0, weight=1)
        self.main_frame.columnconfigure(0, weight=1)   # , minsize=720)

        self.setted = False

        if cu.info:
            cu.update_info(cu.info)

    def _lock_unlock_slider(self, force=False):
        if self.slider_locked or not force:
            self.slider_locked = False
            self.slider.state(['!disabled'])
            self.value_entry.config(state=NORMAL)
        else:
            self.slider_locked = True
            self.slider.state(['disabled'])
            self.value_entry.config(state=DISABLED)

    def _change_view(self):
        if self.graph.cut_line > 0:
            if self.graph.show_selected_area:
                self.graph.show_selected_area = False
                self.change_view_button.config(text='Zoom-in')
            else:
                self.graph.show_selected_area = True
                self.change_view_button.config(text=LANGUAGES[settings.LANGUAGE]['reset_view'])

            self.graph.change_view()
            self.graph.update_cut()

    def resample(self):
        cu.data.fstar = float(self.value_entry.get())
        cu.data.psd_filter_width = float(self.filter_width_value.get())
        if cu.data.fstar > 0:
            if cu.data.changes:
                cu.data.xf = cu.data.j.resample(fstar_THz=cu.data.fstar, PSD_FILTER_W=cu.data.psd_filter_width,
                                                plot=False)
            self.graph.add_graph(cu.gm.plot_resample, 'resample', xf=cu.data.xf, mode='linear', x=cu.data.j,
                                 PSD_FILTER_W=cu.data.psd_filter_width)
            self.graph.update_cut()

            self.graph.cut_line = cu.data.xf.Nyquist_f_THz
            self.fstar_screen.config(text='F*: {}'.format(round(cu.data.xf.Nyquist_f_THz, 3)))
            cu.data.changes = False
        else:
            msg.showwarning('Value error', 'F* must be greater than zero')
            return False

        if self.graph.show_selected_area:
            self.graph.show_selected_area = True
            self.graph.change_view()

        self.update()
        return True

    def set_next_frame(self, frame):
        self.next_frame = frame

    def set_prev_frame(self, frame):
        self.prev_frame = frame

    def back(self):
        cu.data.psd_filter_width = self.filter_width_value.get()
        cu.data.fstar = float(self.value_entry.get())
        cu.Data.loaded = True
        if self.prev_frame:
            self.main.show_frame(self.prev_frame)
        else:
            raise ValueError('Prev frame isn\'t defined')

    def next(self):
        if (self.resample()):
            #cu.data.fstar = cu.data.xf.Nyquist_f_THz

            if self.next_frame:
                self.main.show_frame(self.next_frame)
            else:
                raise ValueError('Next frame isn\'t defined')

    def recalculate(self, slider_config=None):
        cu.data.psd_filter_width = float(self.filter_width_value.get())
        self.graph.other_graph.clear()
        self.graph.graph.clear()
        self.graph.show(cu.gm.plot_periodogram, mode='linear', current=cu.data.j, PSD_FILTER_W=cu.data.psd_filter_width,
                        slider_config=slider_config, kappa_units=True)
        if float(self.value_entry.get()):
            print('resampled')
            self.resample()
        if cu.info:
            cu.update_info(cu.info)

    def update_data(self):
        self.filter_width_value.set(cu.data.psd_filter_width)

    def update(self):
        super().update()
        cu.data.fstar = float(self.value_entry.get())

        if cu.data.first_fstar and '_FSTAR' in cu.data.jdata.keys():
            cu.data.fstar = cu.data.jdata['_FSTAR']
            self.value_entry.delete(0, END)
            self.value_entry.insert(0, cu.data.fstar)
            cu.log.write_log(LANGUAGES[settings.LANGUAGE]['fs_loaded'].format(cu.data.jdata['_FSTAR']))
            self.recalculate(slider_config=cu.data.fstar)
        elif cu.data.changes:
            self.recalculate()

        if cu.info:
            cu.update_info(cu.info)

        self.main_frame_scroll.update_view()
