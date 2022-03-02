# -*- coding: utf-8 -*-

from tkinter.ttk import Separator
from tkinter import messagebox as msg
from sportran_gui.utils.custom_widgets import *
from sportran_gui.assets import LANGUAGES


class OtherVariables(Frame):

    def __init__(self, parent, main):
        Frame.__init__(self, parent)

        self.main = main
        self.next_frame = None
        self.prev_frame = None

        self.main_frame = ScrollFrame(self, self)
        self.main_frame.grid(row=0, column=0, sticky='nswe')

        self.temperature_value = DoubleVar(value=0.1)
        self.volume_value = DoubleVar(value=0.1)
        self.DT_FS_value = DoubleVar(value=0.1)
        self.filter_width_value = DoubleVar(value=0.1)

        variable_frame = Frame(self.main_frame.viewPort)

        variable_frame.grid(column=0, row=0, sticky='nswe', padx=20, pady=5)

        Label(variable_frame, text=LANGUAGES[settings.LANGUAGE]['stp3'],
              font='Arial 14 bold').grid(row=0, column=0, pady=2, sticky='w')

        variable_frame.columnconfigure(0, weight=0, minsize=150)
        variable_frame.columnconfigure(1, weight=1, minsize=200)

        Label(variable_frame, text=LANGUAGES[settings.LANGUAGE]['e_v'],
              font='Arial 11').grid(row=1, column=0, pady=20, sticky='w')
        Separator(variable_frame, orient=HORIZONTAL).grid(row=2, sticky='we', columnspan=3)

        Label(variable_frame, text=LANGUAGES[settings.LANGUAGE]['tmp'] + ' (K):').grid(row=3, column=0, sticky='w')
        self.temperature_entry = Spinbox(variable_frame, from_=0.1, to=10000, increment=0.1, bd=1, relief=SOLID,
                                         textvariable=self.temperature_value)
        self.temperature_entry.grid(row=3, column=1, padx=2, sticky='w', pady=10)
        self.temp_advertise = Label(variable_frame, text='', font='Arial 10')
        self.temp_advertise.grid(row=3, column=2, sticky='w')

        Label(variable_frame,
              text=LANGUAGES[settings.LANGUAGE]['volume'] + ' (Angstrom^3):').grid(row=4, column=0, sticky='w')
        self.volume_entry = Spinbox(variable_frame, from_=0.1, to=10000, increment=0.1, bd=1, relief=SOLID,
                                    textvariable=self.volume_value)
        self.volume_entry.grid(row=4, column=1, padx=2, sticky='w')
        self.volume_advertise = Label(variable_frame, text='', font='Arial 10')
        self.volume_advertise.grid(row=4, column=2, sticky='w')

        Label(variable_frame,
              text=LANGUAGES[settings.LANGUAGE]['timestep'] + ' (fs):').grid(row=5, column=0, sticky='w')
        self.DT_FS_entry = Spinbox(variable_frame, from_=0.1, to=10000, increment=0.1, bd=1, relief=SOLID,
                                   textvariable=self.DT_FS_value)
        self.DT_FS_entry.grid(row=5, column=1, padx=2, sticky='w', pady=10)
        self.DT_FS_advertise = Label(variable_frame, text='', font='Arial 10')
        self.DT_FS_advertise.grid(row=5, column=2, sticky='w')

        Label(variable_frame, text=LANGUAGES[settings.LANGUAGE]['fl_v'],
              font='Arial 11').grid(row=6, column=0, pady=20, sticky='w')
        Separator(variable_frame, orient=HORIZONTAL).grid(row=7, sticky='we', columnspan=3)

        Label(variable_frame, text=LANGUAGES[settings.LANGUAGE]['fl_w'] + ' (THz):').grid(row=8, column=0, sticky='w')
        self.filter_width_entry = Spinbox(variable_frame, from_=0.1, to=10.0, increment=0.1, bd=1, relief=SOLID,
                                          textvariable=self.filter_width_value)
        self.filter_width_entry.grid(row=8, column=1, padx=2, sticky='w', pady=10)

        variable_frame.rowconfigure(9, weight=1)
        button_frame = Frame(variable_frame)
        button_frame.grid(row=9, column=0, sticky='ws', pady=20)

        Button(button_frame, text=LANGUAGES[settings.LANGUAGE]['back'], bd=1, relief=SOLID, font='Arial 12',
               command=lambda: self.back(), width=10).grid(row=0, column=0, sticky='we')

        Button(button_frame, text=LANGUAGES[settings.LANGUAGE]['next'], bd=1, relief=SOLID, font='Arial 12',
               command=lambda: self.next(), width=10).grid(row=0, column=1, padx=5, sticky='we')

        self.update_data()

    def set_next_frame(self, frame):
        self.next_frame = frame

    def set_prev_frame(self, frame):
        self.prev_frame = frame

    def get_entry_data(self):
        if self.filter_width_value.get():
            cu.data.psd_filter_width = float(self.filter_width_value.get())
        if self.temperature_value.get():
            cu.data.temperature = float(self.temperature_value.get())
        if self.volume_value.get():
            cu.data.volume = float(self.volume_value.get())
        if self.DT_FS_value.get():
            cu.data.DT_FS = float(self.DT_FS_value.get())

    def next(self):
        self.get_entry_data()

        er = False
        msgs = []
        if cu.data.temperature:
            if cu.data.temperature < 0 and 'Temperature' not in cu.data.description:
                msgs.append(LANGUAGES[settings.LANGUAGE]['temp_low'])
                er = True
        else:
            msgs.append(LANGUAGES[settings.LANGUAGE]['temp_void'])

        if cu.data.volume:
            if cu.data.volume < 0:
                msgs.append(LANGUAGES[settings.LANGUAGE]['vol_low'])
                er = True
        else:
            msgs.append(LANGUAGES[settings.LANGUAGE]['vol_void'])
            er = True

        if cu.data.DT_FS:
            if cu.data.DT_FS < 0:
                msgs.append(LANGUAGES[settings.LANGUAGE]['DT_low'])
                er = True
        else:
            msgs.append(LANGUAGES[settings.LANGUAGE]['DT_void'])
            er = True

        if cu.data.psd_filter_width:
            if cu.data.psd_filter_width < 0:
                msgs.append(LANGUAGES[settings.LANGUAGE]['fw_low'])
                er = True
        else:
            msgs.append(LANGUAGES[settings.LANGUAGE]['fw_void'])
            er = True

        if not er:
            cu.load_data(cu.data.CURRENT_FILE, cu.data.inputformat, _selected_keys=cu.data.keys,
                         _descriptions=cu.data.description, temperature=cu.data.temperature, units=cu.data.units,
                         volume=cu.data.volume, psd_filter_w=cu.data.psd_filter_width, DT_FS=cu.data.DT_FS)

            if self.next_frame:
                self.main.show_frame(self.next_frame)
            else:
                raise ValueError('Next frame isn\'t defined')
        else:
            ermsg = '\n'.join(mser for mser in msgs)
            msg.showerror(LANGUAGES[settings.LANGUAGE]['value_error'], ermsg)

    def back(self):
        self.get_entry_data()

        if self.prev_frame:
            self.main.show_frame(self.prev_frame)
        else:
            raise ValueError('Prev frame isn\'t defined')

    def update_data(self):
        self.temperature_value.set(cu.data.temperature)
        self.DT_FS_value.set(cu.data.DT_FS)
        self.volume_value.set(cu.data.volume)
        self.filter_width_value.set(cu.data.psd_filter_width)

    def update(self):
        super().update()

        self.temperature_entry.config(state=NORMAL)
        self.DT_FS_entry.config(state=NORMAL)
        self.volume_entry.config(state=NORMAL)

        if cu.Data.options[3] in cu.data.description:
            self.temp_advertise.config(text=LANGUAGES[settings.LANGUAGE]['automatic_T'], fg='red')
            if cu.data.inputformat == 'dict':   #only for dict I have the data here
                temp = cu.data.jdata[cu.data.keys[cu.data.description.index(cu.Data.options[3])]]
                if type(temp) == float:
                    cu.data.temperature = temp
                else:
                    cu.data.temperature = temp.mean()
                    cu.data.temperature_std = temp.std()
            self.temperature_entry.config(state=DISABLED)
        else:
            self.temperature_entry.config(state=NORMAL)
            self.temp_advertise.config(text='')
        if cu.Data.options[4] in cu.data.description:
            self.volume_advertise.config(text=LANGUAGES[settings.LANGUAGE]['automatic_V'], fg='red')
            if cu.data.inputformat == 'dict':   #only for dict I have the data here
                vol = cu.data.jdata[cu.data.keys[cu.data.description.index(cu.Data.options[4])]]
                if type(vol) == float:
                    cu.data.volume = vol
                else:
                    cu.data.volume = vol.mean()
            else:
                raise RuntimeError('NOT IMPLEMENTED')
            self.volume_entry.config(state=DISABLED)
        else:
            self.volume_advertise.config(text='')
            self.volume_entry.config(state=NORMAL)
        if cu.Data.options[5] in cu.data.description:
            self.DT_FS_advertise.config(text=LANGUAGES[settings.LANGUAGE]['automatic_DT'], fg='red')
            if cu.data.inputformat == 'dict':   #only for dict I have the data here
                cu.data.DT_FS = cu.data.jdata[cu.data.keys[cu.data.description.index(cu.Data.options[5])]]
            else:
                raise RuntimeError('NOT IMPLEMENTED')
            self.DT_FS_entry.config(state=DISABLED)
        else:
            self.DT_FS_advertise.config(text='')
            self.DT_FS_entry.config(state=NORMAL)

        self.update_data()

        self.main_frame.viewPort.columnconfigure(0, weight=1)
        self.main_frame.viewPort.rowconfigure(0, weight=1)
        self.main_frame.update_view()
