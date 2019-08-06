from thermocepstrum_gui.utils.custom_widgets import *
from thermocepstrum_gui.core import control_unit as cu
from thermocepstrum_gui.core.control_unit import log


class PStarSelector(Frame):

    def __init__(self, parent, main):
        Frame.__init__(self, parent)

        self.main = main

        self.next_frame = None
        self.prev_frame = None

        self.parent = parent

        sections = Frame(self)
        sections.pack(side=LEFT, anchor='n', fill=BOTH, expand=1)
        self.graph = GraphWidget(parent, sections, size=(7, 4), toolbar=True)

        variable_frame = Frame(sections, bd=1, relief=SOLID)
        variable_frame.pack(side=TOP, pady=3, fill=BOTH, expand=1)

        Label(variable_frame, text='f*: ', font='Arial 12 bold').grid(row=0, column=0)
        self.fstar_label = Label(variable_frame, text='')
        self.fstar_label.grid(row=0, column=1)
        Label(variable_frame, text='P*: ', font='Arial 12 bold').grid(row=0, column=2)
        self.pstar_label = Label(variable_frame, text='')
        self.pstar_label.grid(row=0, column=3)
        Label(variable_frame, text='\u03f0:', font='Arial 12 bold').grid(row=0, column=4)
        self.kmin_label = Label(variable_frame, text='')
        self.kmin_label.grid(row=0, column=5)
        value_frame = Frame(sections)
        value_frame.pack(side=TOP)

        Label(value_frame, text='P*').pack(side=TOP, pady=10)
        self.value_entry = Spinbox(value_frame, bd=1, relief=SOLID, increment=1)
        self.value_entry.pack()

        self.increment = IntVar()
        Label(value_frame, text='Increment by').pack(side=TOP, pady=4)
        Radiobutton(value_frame, text='1', variable=self.increment, value=1,
                    command=self._change_increment).pack(side=LEFT, anchor='n', padx=2)
        Radiobutton(value_frame, text='10', variable=self.increment, value=10,
                    command=self._change_increment).pack(side=LEFT, anchor='n', padx=2)
        Radiobutton(value_frame, text='100', variable=self.increment, value=100,
                    command=self._change_increment).pack(side=LEFT, anchor='n', padx=2)

        button_frame = Frame(sections)
        button_frame.pack(pady=20)

        back_button = Button(button_frame, text='Back', bd=1, relief=SOLID, command=lambda: self.back())
        back_button.grid(row=0, column=0, sticky='w', padx=5)

        next_button = Button(button_frame, text='Recalculate', bd=1, relief=SOLID, command=self._reload)
        next_button.grid(row=0, column=1, sticky='w', padx=5)

        self.info_section = Frame(self)
        self.info_section.pack(side=RIGHT, anchor='n', pady=30, padx=20, fill='x', expand=True)

        self.logs = None
        self.info = None

        self._init_output_frame()

        self.setted = False

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
        self.logs = None
        self.info = None
        self.update()

    def set_prev_frame(self, frame):
        self.prev_frame = frame

    def back(self):
        if self.next_frame:
            self.main.show_frame(self.next_frame)
        else:
            raise ValueError('Prev frame isn\'t defined')

    def _get_pstar(self, aic_type='aic', Kmin_corrfactor=1.0):
        cu.Data.xf.cepstral_analysis(aic_type=aic_type, K_PSD=Kmin_corrfactor-1)

    def _pstar(self):
        self.value_entry.config(from_=2, to=cu.Data.xf.Nfreqs)
        self.value_entry.delete(0, END)
        self.value_entry.insert(0, (cu.Data.xf.dct.aic_Kmin+1))

        self.fstar_label.config(text='{:4f}'.format(cu.Data.fstar))
        self.pstar_label.config(text=f'{cu.Data.xf.dct.aic_Kmin+1}')
        self.kmin_label.config(text='{:18f} +/- {:10f} W/mK'.format(cu.Data.xf.kappa_Kmin, cu.Data.xf.kappa_Kmin_std))

    def _change_increment(self):
        self.value_entry.config(increment=int(self.increment.get()))

    def _reload(self):
        self._get_pstar(aic_type='aic', Kmin_corrfactor=int(self.value_entry.get()))
        self._get_pstar(aic_type='aic', Kmin_corrfactor=int(self.value_entry.get()))
        self.graph.add_graph(cu.gm.plot_cepstral_spectrum, 'cepstral', x=cu.Data.xf)
        self.graph.update_cut()

    def _setup_pstar(self):
        cu.Data.xf.cepstral_analysis(aic_type='aic', K_PSD=None)
        self._pstar()

    def update(self):
        super().update()

        if cu.Data.fstar == cu.Data.old_fstar:
            self.setted = True
        else:
            self.setted = False
            cu.Data.old_fstar = cu.Data.fstar

        if not self.setted:
            self.setted = True
            self._setup_pstar()

        self.graph.show(cu.gm.GUI_plot_periodogram, x=cu.Data.j)
        self.graph.add_graph(cu.gm.resample_current, 'resample', x=cu.Data.j, fstar_THz=cu.Data.fstar,
                             PSD_FILTER_W=cu.Data.psd_filter_width)

        self._init_output_frame()
        if self.info:
            cu.update_info(self.info)
        if self.logs:
            log.set_func(self.logs.write)
        self._reload()
        self.graph.update_cut()
