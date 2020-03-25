from tkinter import messagebox as msg
from thermocepstrum_gui.utils.custom_widgets import *
from thermocepstrum import HeatCurrent

import traceback


class HeaderSelector(Frame):

    def __init__(self, parent, main):
        Frame.__init__(self, parent)

        self.main = main
        self.next_frame = None
        self.prev_frame = None

        self.main_frame = ScrollFrame(self, self)

        header_frame = self.main_frame.viewPort
        #header_frame.grid(row=0, column=0, sticky='nswe', padx=20, pady=5)

        Label(header_frame, text=LANGUAGES[settings.LANGUAGE]['stp2'],
              font='Arial 14 bold').grid(row=0, column=0, sticky='w')

        definitions_frame = Frame(header_frame)
        definitions_frame.grid(row=1, column=1, sticky='wn', padx=20)
        Label(definitions_frame, text=LANGUAGES[settings.LANGUAGE]['df1'], wraplengt=500, justify=LEFT)\
            .grid(row=0, column=0, sticky='wn')
        Label(definitions_frame, text=LANGUAGES[settings.LANGUAGE]['df2'], wraplengt=500, justify=LEFT)\
            .grid(row=1, column=0, sticky='wn')
        # todo: put definition
        Label(definitions_frame, text=LANGUAGES[settings.LANGUAGE]['df3'], wraplengt=500, justify=LEFT)\
            .grid(row=2, column=0, sticky='wn')
        Label(definitions_frame, text=LANGUAGES[settings.LANGUAGE]['df4'], wraplengt=500, justify=LEFT)\
            .grid(row=3, column=0, sticky='wn')

        definitions_frame.columnconfigure(0, weight=1)
        for i in range(0, 3):
            definitions_frame.rowconfigure(i, weight=1)

        header_list_frame = Frame(header_frame)
        header_list_frame.grid(row=1, column=0, sticky='nswe', pady=10)

        self.scrollable_header_list = ScrollFrame(header_list_frame, header_list_frame, bd=1)
        self.check_list = CheckList(self.scrollable_header_list.viewPort, self.scrollable_header_list.viewPort,
                                    start_row=1)

        Label(header_frame, text=LANGUAGES[settings.LANGUAGE]['select_units'], font='Arial 14 bold') \
            .grid(row=2, column=0, sticky='w', pady=10)
        self.units_selector_frame = Frame(header_frame)
        self.units_selector = ttk.Combobox(self.units_selector_frame, values=HeatCurrent.get_units_list(),
                                           state='readonly')
        self.units_selector.current(0)
        Label(self.units_selector_frame,
              text=LANGUAGES[settings.LANGUAGE]['units']).grid(row=0, column=0, sticky='we', pady=2)
        self.units_selector.grid(row=0, column=1, sticky='we', pady=2)
        self.units_selector_frame.grid(row=3, column=0, sticky='nswe', pady=2)

        button_frame = Frame(header_frame)
        button_frame.grid(row=4, column=0, sticky='w')

        header_frame.rowconfigure(0, weight=1, minsize=50)
        header_frame.rowconfigure(1, weight=10, minsize=150)
        header_frame.rowconfigure(2, weight=1, minsize=50)
        header_frame.rowconfigure(3, weight=10, minsize=50)
        header_frame.rowconfigure(4, weight=1, minsize=50)
        header_frame.columnconfigure(0, weight=1, minsize=300)
        header_frame.columnconfigure(1, weight=1, minsize=500)

        Button(button_frame, text=LANGUAGES[settings.LANGUAGE]['back'], bd=1, relief=SOLID, font='Arial 12',
               command=lambda: self.back(), width=10).grid(row=0, column=0, sticky='we')

        Button(button_frame, text=LANGUAGES[settings.LANGUAGE]['next'], bd=1, relief=SOLID, font='Arial 12',
               command=lambda: self.next(), width=10).grid(row=0, column=1, padx=5, sticky='we')

    def set_next_frame(self, frame):
        self.next_frame = frame

    def set_prev_frame(self, frame):
        self.prev_frame = frame

    def next(self):
        keys, description = self.check_list.get_list()
        if cu.Data.options[1] in description:   #energy current
            if description.count(cu.Data.options[1]) == 1:
                if description.count(cu.Data.options[3]) <= 1:   #temperature
                    if description.count(cu.Data.options[4]) <= 1:
                        if description.count(cu.Data.options[5]) <= 1:
                            cu.data.keys = keys
                            cu.data.description = description
                            cu.Data.loaded = True
                            cu.data.units = self.units_selector.get()
                            if self.next_frame:
                                self.main.show_frame(self.next_frame)
                            else:
                                raise ValueError('Next frame isn\'t defined')
                        else:
                            msg.showerror(LANGUAGES[settings.LANGUAGE]['value_error'],
                                          LANGUAGES[settings.LANGUAGE]['only_one_DT'])
                    else:
                        msg.showerror(LANGUAGES[settings.LANGUAGE]['value_error'],
                                      LANGUAGES[settings.LANGUAGE]['only_one_Volume'])
                else:
                    msg.showerror(LANGUAGES[settings.LANGUAGE]['value_error'],
                                  LANGUAGES[settings.LANGUAGE]['only_one_T'])
            else:
                msg.showerror(LANGUAGES[settings.LANGUAGE]['value_error'], LANGUAGES[settings.LANGUAGE]['only_one_E'])
        else:
            msg.showerror(LANGUAGES[settings.LANGUAGE]['no_key'], LANGUAGES[settings.LANGUAGE]['no_key_t'])

    def back(self):
        keys, description = self.check_list.get_list()
        cu.data.keys = keys
        cu.data.description = description

        if self.prev_frame:
            self.main.show_frame(self.prev_frame)
        else:
            raise ValueError('Prev frame isn\'t defined')

    def update(self):
        super().update()
        self.main_frame.update_view()

        try:
            print('cu.Data.loaded={}'.format(cu.Data.loaded))
            if not cu.Data.loaded:
                keys = cu.load_keys(cu.data.CURRENT_FILE)
            else:
                keys = {key: '' for key in cu.data.keys}

            self.check_list.set_list(keys)

            #try to set units (if given)
            try:
                self.units_selector.current(HeatCurrent.get_units_list().index(cu.data.jdata['_UNITS']))
                cu.log.write_log(LANGUAGES[settings.LANGUAGE]['units_loaded'].format(cu.data.jdata['_UNITS']))
            except:
                pass

            if cu.Data.loaded:
                for i, check in enumerate(self.check_list.controller.winfo_children()):
                    check.winfo_children()[1].current(cu.Data.options.index(cu.data.description[i]))
        except Exception as e:
            cu.Data.loaded = False
            cu.log.write_log(str(e))
            traceback.print_exc()
            msg.showerror(LANGUAGES[settings.LANGUAGE]['read_error'], LANGUAGES[settings.LANGUAGE]['read_error_t'])
            if self.prev_frame:
                self.main.show_frame(self.prev_frame)
            else:
                raise ValueError('Prev frame isn\'t defined')
