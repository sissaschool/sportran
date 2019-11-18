from tkinter.ttk import Separator
from tkinter import messagebox as msg
from thermocepstrum_gui.utils.custom_widgets import *
from thermocepstrum_gui.assets import LANGUAGES


class OtherVariables(Frame):

    def __init__(self, parent, main):
        Frame.__init__(self, parent)

        self.main = main
        self.next_frame = None
        self.prev_frame = None

        self.main_frame = ScrollFrame(self, self)
        self.main_frame.grid(row=0, column=0, sticky='nswe')

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
        self.temperature_entry = Spinbox(variable_frame, from_=0, to=100000, increment=0.1, bd=1, relief=SOLID)
        self.temperature_entry.grid(row=3, column=1, padx=2, sticky='w', pady=10)
        self.temp_advertise = Label(variable_frame, text='', font='Arial 10')
        self.temp_advertise.grid(row=3, column=2, sticky='w')

        Label(variable_frame,
              text=LANGUAGES[settings.LANGUAGE]['volume'] + ' (angstrom^3):').grid(row=4, column=0, sticky='w')
        self.volume_entry = Spinbox(variable_frame, increment=1, bd=1, relief=SOLID)
        self.volume_entry.grid(row=4, column=1, padx=2, sticky='w')

        Label(variable_frame,
              text=LANGUAGES[settings.LANGUAGE]['timestep'] + ' (fs):').grid(row=5, column=0, sticky='w')
        self.DT_FS_entry = Spinbox(variable_frame, increment=0.1, bd=1, relief=SOLID)
        self.DT_FS_entry.grid(row=5, column=1, padx=2, sticky='w', pady=10)

        Label(variable_frame, text=LANGUAGES[settings.LANGUAGE]['fl_v'],
              font='Arial 11').grid(row=6, column=0, pady=20, sticky='w')
        Separator(variable_frame, orient=HORIZONTAL).grid(row=7, sticky='we', columnspan=3)

        Label(variable_frame, text=LANGUAGES[settings.LANGUAGE]['fl_w'] + ' (THz):').grid(row=8, column=0, sticky='w')
        self.filter_width_entry = Spinbox(variable_frame, from_=0.1, to=10.0, increment=0.1, bd=1, relief=SOLID)
        self.filter_width_entry.grid(row=8, column=1, padx=2, sticky='w', pady=10)

        variable_frame.rowconfigure(9, weight=1)
        button_frame = Frame(variable_frame)
        button_frame.grid(row=9, column=0, sticky='ws', pady=20)

        Button(button_frame, text=LANGUAGES[settings.LANGUAGE]['back'], bd=1, relief=SOLID, font='Arial 12',
               command=lambda: self.back(), width=10).grid(row=0, column=0, sticky='we')

        Button(button_frame, text=LANGUAGES[settings.LANGUAGE]['next'], bd=1, relief=SOLID, font='Arial 12',
               command=lambda: self.next(), width=10).grid(row=0, column=1, padx=5, sticky='we')

    def set_next_frame(self, frame):
        self.next_frame = frame

    def set_prev_frame(self, frame):
        self.prev_frame = frame

    def next(self):
        psd_filter_width = None
        temperature = None
        volume = None
        DT_FS = None

        if self.filter_width_entry.get():
            psd_filter_width = float(self.filter_width_entry.get())
        if self.temperature_entry.get():
            temperature = float(self.temperature_entry.get())
        if self.volume_entry.get():
            volume = float(self.volume_entry.get())
        if self.DT_FS_entry.get():
            DT_FS = float(self.DT_FS_entry.get())

        er = False
        msgs = []
        if temperature:
            if temperature < 0 and 'Temperature' not in cu.data.description:
                msgs.append(LANGUAGES[settings.LANGUAGE]['temp_low'])
                er = True
        else:
            msgs.append(LANGUAGES[settings.LANGUAGE]['temp_void'])

        if volume:
            if volume < 0:
                msgs.append(LANGUAGES[settings.LANGUAGE]['vol_low'])
                er = True
        else:
            msgs.append(LANGUAGES[settings.LANGUAGE]['vol_void'])
            er = True

        if DT_FS:
            if DT_FS < 0:
                msgs.append(LANGUAGES[settings.LANGUAGE]['DT_low'])
                er = True
        else:
            msgs.append(LANGUAGES[settings.LANGUAGE]['DT_void'])
            er = True

        if psd_filter_width:
            if psd_filter_width < 0:
                msgs.append(LANGUAGES[settings.LANGUAGE]['fw_low'])
                er = True
        else:
            msgs.append(LANGUAGES[settings.LANGUAGE]['fw_void'])
            er = True

        if not er:
            cu.load_data(cu.data.CURRENT_FILE, cu.data.inputformat, _selected_keys=cu.data.keys,
                         _descriptions=cu.data.description, temperature=temperature, units=cu.data.units, volume=volume,
                         psd_filter_w=cu.data.psd_filter_width, DT_FS=DT_FS)

            if self.next_frame:
                self.main.show_frame(self.next_frame)
            else:
                raise ValueError('Next frame isn\'t defined')
        else:
            ermsg = '\n'.join(mser for mser in msgs)
            msg.showerror(LANGUAGES[settings.LANGUAGE]['value_error'], ermsg)

    def back(self):
        if self.filter_width_entry.get():
            cu.data.psd_filter_width = float(self.filter_width_entry.get())
        if self.temperature_entry.get():
            cu.data.temperature = float(self.temperature_entry.get())
        if self.volume_entry.get():
            cu.data.volume = float(self.volume_entry.get())
        if self.DT_FS_entry.get():
            cu.data.DT_FS = float(self.DT_FS_entry.get())

        if self.prev_frame:
            self.main.show_frame(self.prev_frame)
        else:
            raise ValueError('Prev frame isn\'t defined')

    def update(self):
        super().update()

        self.temperature_entry.config(state=NORMAL)
        self.temperature_entry.config(value=cu.data.temperature)
        self.DT_FS_entry.config(state=NORMAL)
        self.volume_entry.config(state=NORMAL)

        if cu.Data.options[3] in cu.data.description:
            self.temperature_entry.config(state=DISABLED)
            self.temp_advertise.config(text=LANGUAGES[settings.LANGUAGE]['automatic_T'], fg='red')
        else:
            self.temperature_entry.config(state=NORMAL)
            self.temp_advertise.config(text='')

        if cu.data.inputformat == 'dict':
            if True:   #not cu.Data.loaded:
                if cu.Data.options[3] in cu.data.description:
                    temp = cu.data.jdata[cu.data.keys[cu.data.description.index(cu.Data.options[3])]]
                    if type(temp) == float:
                        self.temperature_entry.config(value=temp)
                        cu.data.temperature = temp
                    else:
                        cu.data.temperature = temp.mean()
                        cu.data.temperature_std = temp.std()
                        self.temperature_entry.config(value=cu.data.temperature)
                if cu.Data.options[4] in cu.data.description:
                    vol = cu.data.jdata[cu.data.keys[cu.data.description.index(cu.Data.options[4])]]
                    if type(vol) == float:
                        cu.data.volume = vol
                    else:
                        cu.data.volume = vol.mean()
                    self.volume_entry.config(state=DISABLED)
                if cu.Data.options[4] in cu.data.description:
                    cu.data.DT_FS = cu.data.jdata[cu.data.keys[cu.data.description.index(cu.Data.options[5])]]
                    self.DT_FS_entry.config(state=DISABLED)

        self.volume_entry.config(value=cu.data.volume)
        self.DT_FS_entry.config(value=cu.data.DT_FS)
        self.filter_width_entry.config(value=cu.data.psd_filter_width)
        self.main_frame.viewPort.columnconfigure(0, weight=1)
        self.main_frame.viewPort.rowconfigure(0, weight=1)
        self.main_frame.update_view()
