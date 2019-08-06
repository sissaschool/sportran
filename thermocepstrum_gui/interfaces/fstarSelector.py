from tkinter import messagebox as msg
from thermocepstrum_gui.utils.custom_widgets import *
from thermocepstrum_gui.core.control_unit import log


class FStarSelector(Frame):

    def __init__(self, parent, main):
        Frame.__init__(self, parent)

        self.main = main

        self.next_frame = None
        self.prev_frame = None

        self.parent = parent
        main_frame = self

        main_frame.grid(column=0, row=0, sticky='nswe')

        sections = Frame(main_frame)
        sections.grid(row=0, column=0, sticky='nswe')

        self.graph = GraphWidget(sections, sections, size=(7, 4), toolbar=True)
        self.graph.pack(side=TOP, anchor='w', padx=10, fill=BOTH, expand=1)

        self.slider_locked = False

        slider_frame = Frame(sections)
        slider_frame.pack(side=TOP, anchor='w', padx=80, fill=BOTH, expand=1)

        self.slider = ttk.Scale(slider_frame, from_=0, to_=0.1)
        self.slider.grid(row=0, column=0, sticky='we')
        slider_frame.columnconfigure(0, weight=9)

        lock_slider = Button(slider_frame, command=lambda: self._lock_unlock_slider(),
                             bd=1, relief=SOLID)
        lock_slider.grid(row=0, column=1, padx=2)
        slider_frame.columnconfigure(1, weight=1)
        self.graph.attach_slider(self.slider)

        self.change_view_button = Button(slider_frame, text='Zoom-in', command=lambda: self._change_view())
        self.change_view_button.grid(row=0, column=2, sticky='w')
        slider_frame.columnconfigure(2, weight=1)

        value_frame = Frame(sections)
        value_frame.pack(side=LEFT, pady=10, padx=20, fill=BOTH, expand=1)

        Label(value_frame, text='Selected value:').grid(row=0, column=0, sticky='w', pady=4)
        self.value_entry = Entry(value_frame, bd=1, relief=SOLID)
        self.value_entry.grid(row=0, column=1, sticky='w')
        self.graph.attach_entry(self.value_entry)

        Label(value_frame, text='Filter width:').grid(row=1, column=0, sticky='w')
        self.filter_width = Spinbox(value_frame, from_=0.1, to=10, increment=0.1, bd=1, relief=SOLID)
        self.filter_width.grid(row=1, column=1, sticky='w', pady=10)

        Button(value_frame, text='Resample', bd=1, relief=SOLID,
               command=self.resample).grid(row=2, column=0, sticky='w')

        button_frame = Frame(value_frame)
        button_frame.grid(row=3, column=0)

        back_button = Button(button_frame, text='Back', bd=1, relief=SOLID, command=lambda: self.back())
        back_button.grid(row=0, column=0, sticky='w', padx=5)

        next_button = Button(button_frame, text='Next', bd=1, relief=SOLID, command=lambda: self.next())
        next_button.grid(row=0, column=1, sticky='w', padx=5)

        self.info_section = Frame(main_frame)
        self.info_section.grid(row=0, column=1, sticky='nswe', pady=20)

        main_frame.columnconfigure(0, weight=3)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(0, weight=1)
        self.logs = None
        self.info = None

        self._init_output_frame()

        # StatusFrame(parent, controller)

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
                self.change_view_button.config(text='Reset view')

            self.graph.change_view()
            self.graph.update_cut()

    def _init_output_frame(self):
        if TopBar.show_logs.get():
            if not self.logs:
                self.logs = TextWidget(self.parent, self.info_section, 'Logs', 15, 45)
        else:
            if self.logs:
                self._del_out_frames()

        if TopBar.show_info.get():
            if not self.info:
                self.info = TextWidget(self.parent, self.info_section, 'Info', 10, 45)
        else:
            if self.info:
                self._del_out_frames()

    def _del_out_frames(self):
        for el in self.info_section.winfo_children():
            el.destroy()
        log.set_func(None)
        self.logs = None
        self.info = None
        self.update()

    def resample(self):
        cu.Data.fstar = float(self.value_entry.get())
        filter_width = float(self.filter_width.get())
        cu.Data.psd_filter_width = filter_width

        if cu.Data.fstar > 0:
            self.graph.add_graph(cu.gm.resample_current, 'resample', x=cu.Data.j, fstar_THz=cu.Data.fstar,
                                 PSD_FILTER_W=cu.Data.psd_filter_width)
            self.graph.update_cut()

        else:
            msg.showwarning('Value error', 'F* must be greater than zero')

        self.graph.cut_line = cu.Data.xf.Nyquist_f_THz

        if self.graph.show_selected_area:
            self.graph.show_selected_area = True
            self.graph.change_view()
        self.update()

    def set_next_frame(self, frame):
        self.next_frame = frame

    def set_prev_frame(self, frame):
        self.prev_frame = frame

    def back(self):
        response = msg.askyesnocancel('Back to file manager?',
                                      "Save changes?\nIf reopen the same file "
                                      "\nthe values that you chosed will not be deleted!")

        log.set_func(None)
        if response:
            cu.Data.fstar = float(self.value_entry.get())
            cu.Data.loaded = True
            if self.next_frame:
                self.main.show_frame(self.next_frame)
            else:
                raise ValueError('Prev frame isn\'t defined')

        elif not response:
            cu.Data.fstar = 0.0
            cu.Data.loaded = False
            cu.Data.temperature = 0.0
            cu.Data.volume = 0.0
            cu.Data.DT_FS = 0.0
            cu.Data.psd_filter_width = 0.1

            self.graph.other_graph.clear()
            self.graph.graph.clear()
            self.graph.cut_line = 0
            if self.next_frame:
                self.main.show_frame(self.next_frame)
            else:
                raise ValueError('Prev frame isn\'t defined')
        else:
            pass

    def next(self):
        self.resample()
        cu.Data.fstar = cu.Data.xf.Nyquist_f_THz
        if self.next_frame:
            self.main.show_frame(self.next_frame)
        else:
            raise ValueError('Next frame isn\'t defined')

    def update(self):
        super().update()

        self.graph.show(cu.gm.GUI_plot_periodogram, x=cu.Data.j, PSD_FILTER_W=cu.Data.psd_filter_width)
        self._init_output_frame()
        if self.info:
            cu.update_info(self.info)
        if self.logs:
            log.set_func(self.logs.write)
