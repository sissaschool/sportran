'''
--------------------------------------------
    Thermocepstrum graphic user interface

    Version: 0.0.1
    Release state: Beta
    Last update: 26/06/2019

    Developer: Sebastiano Bisacchi
--------------------------------------------

This file contains the GUI of the Thermocepstrum project developed at SISSA
'''
# todo: Put an accurate description of the project?

import os

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.patches as patches

import numpy as np

from tkinter import *
from tkinter import ttk
from tkinter.ttk import Separator, Progressbar
from tkinter import messagebox as msg
import tkinter.filedialog as dialog
from tkinter.font import Font

from core import settings
import core.control_unit as cu


# Main app
class ThermocepstrumGUI(Tk):
    open_windows = []
    frames = []
    frame = None

    def __init__(self, version, dev_state, last_release, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)

        # Define some info
        self.version = version
        self.dev_state = dev_state
        self.last_release = last_release
        # Add the main window to the open windows
        ThermocepstrumGUI.open_windows.insert(0, self)

        # Configure the main window
        self.title('Thermocepstrum')
        self.geometry('{}x{}+{}+{}'.format(settings.X_SIZE,
                                           settings.Y_SIZE,
                                           settings.X_SPACING,
                                           settings.Y_SPACING))

        self.configure(bg=settings.BG_COLOR)
        self.resizable(settings.X_RESIZE, settings.Y_RESIZE)
        self.protocol('WM_DELETE_WINDOW', func=lambda: cu.secure_exit(ThermocepstrumGUI))

        # Creating the main frame
        container = Frame(self)
        container.pack(side=TOP, fill=BOTH, expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # Setting up multiple window system

        ThermocepstrumGUI.frames = {}

        for F in (FileManager, Cutter, PStar, Output):
            ThermocepstrumGUI.frame = F(container, self)
            ThermocepstrumGUI.frame.configure(bg=settings.BG_COLOR)

            ThermocepstrumGUI.frames[F] = ThermocepstrumGUI.frame
            ThermocepstrumGUI.frame.grid(row=0, column=0, sticky='nsew')

        self.show_frame(FileManager)

    @staticmethod
    def show_frame(frame):
        ThermocepstrumGUI.frame = ThermocepstrumGUI.frames[frame]
        ThermocepstrumGUI.frame.tkraise()
        ThermocepstrumGUI.frame.update()


class TopBar(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        # Setup the top menu
        top_menu = Menu(self)
        controller.configure(menu=top_menu)

        # Create the file section of the top menu
        file_menu = Menu(top_menu, tearoff=False)
        top_menu.add_cascade(label='File', menu=file_menu)

        file_menu.add_command(label='Import data')
        file_menu.add_command(label='Export data')
        file_menu.add_separator()
        file_menu.add_command(label='Preferences')
        file_menu.add_separator()
        file_menu.add_command(label='Exit', command=lambda: cu.secure_exit(ThermocepstrumGUI))

        # Create the tool section of the top menu
        # tool_menu = Menu(top_menu, tearoff=False)
        # top_menu.add_cascade(label='Tools', menu=tool_menu)

        # Create the info section of the top menu
        file_menu = Menu(top_menu, tearoff=False)
        top_menu.add_cascade(label='Info', menu=file_menu)

        file_menu.add_command(label='Version')
        file_menu.add_separator()
        file_menu.add_command(label='Developers')
        file_menu.add_command(label='Contacts')
        file_menu.add_command(label='About')
        file_menu.add_separator()
        file_menu.add_command(label='Help')


class StatusFrame(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        status_frame = Frame(controller)
        status_frame.pack(fill='x', side=BOTTOM)

        self.status = Label(status_frame, text=('Status: ' + settings.STATUS_NOW))
        self.status.pack(side=LEFT, padx=4, pady=2)


class GraphWidget(Frame):

    def __init__(self, parent, controller, size=(4, 4), type=111, toolbar=False):
        Frame.__init__(self, parent)

        self.title = ''
        self.size = size
        self.type = type

        self.cut_line = 0
        self.max_x = 1
        self.max_y = 1

        self.f = Figure(figsize=self.size, dpi=100)
        self.graph = self.f.add_subplot(self.type)

        self.canvas = FigureCanvasTkAgg(self.f, controller)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, anchor='w', padx=10)

        self.func = None

        self.other_graph = []

        self.slider = None
        self.entry = None
        self.line2D = None
        self.plot_call = None
        self.plot_call_kwargs = None

        if toolbar:
            toolbar = NavigationToolbar2Tk(self.canvas, controller)
            toolbar.update()
            self.canvas._tkcanvas.pack(side=TOP)

    def set_title(self, title):
        self.f.suptitle(title)

    def get_max_x(self):
        if self.graph:
            return self.graph.get_xlim()[1]
        else:
            return 1.0

    def get_max_y(self):
        if self.graph:
            return self.graph.get_ylim()[1]
        else:
            return 1.0

    def show(self, func, **kwargs):
        self.func = func
        # todo: add kwargs
        cu.set_graph(self.graph, func, **kwargs)
        self.max_x = self.get_max_x()
        self.max_y = self.get_max_y()
        if self.slider:
            self.slider.config(to_=self.max_x)
        self.update_cut()

    def add_graph(self, func, name, **kwargs):
        exist = False
        pos = 0
        for p, n in enumerate(self.other_graph):
            if n[0] == name:
                exist = True
                pos = p
                break
        if not exist:
            self.other_graph.append([name, func, kwargs])
        else:
            del self.other_graph[pos]
            self.other_graph.append([name, func, kwargs])

    def update_cut(self):
        if self.graph:
            if self.entry:
                self.entry.delete(0, END)
                self.entry.insert(0, self.cut_line)

            self.graph.clear()
            cu.set_graph(self.graph, self.func, x=cu.Data.j, PSD_FILTER_W=cu.Data.psd_filter_width)
            for graph in self.other_graph:
                cu.set_graph(self.graph, graph[1], **graph[2])

            rect = patches.Rectangle((0, 0), self.cut_line, self.max_y, linewidth=0, facecolor=(0.1, 0.2, 0.5, 0.3))
            self.graph.plot([self.cut_line, self.cut_line], [0, self.max_y])
            self.graph.add_patch(rect)
        self.canvas.draw()

    def get_graph(self):
        return self.graph

    def set_plot_call(self, f, **f_args):
        self.plot_call = f
        self.plot_call_kwargs = f_args

    def attach_slider(self, slider):
        self.slider = slider
        self.slider.config(command=self._on_slider_change, to_=self.get_max_x())

    def attach_entry(self, entry):
        self.entry = entry
        self.entry.bind('<Key-Return>', self._on_entry_change)
        self.entry.delete(0, END)
        self.entry.insert(0, self.cut_line)

    def _on_entry_change(self, ev):
        self.cut_line = float(self.entry.get())
        self.update_cut()

    def _on_slider_change(self, ev):
        self.cut_line = self.slider.get()
        self.update_cut()


class TextWidget(Frame):

    def __init__(self, parent, controller, title, height):
        Frame.__init__(self, parent, controller)

        text_frame = LabelFrame(controller, text=title, bg=settings.BG_COLOR, bd=1, relief=SOLID)
        text_frame.pack(side=TOP)

        self.text_box = Text(text_frame, height=height, bd=0)
        self.text_box.pack(padx=4, pady=4)

    def write(self, text):
        self.text_box.insert(INSERT, str(text)+'\n')


class CheckList(Frame):

    def __init__(self, parent, controller, check_list=dict()):
        Frame.__init__(self, parent, controller)

        self.controller = controller
        if list:
            for row, el in enumerate(list(check_list.keys())):
                chk = ttk.Checkbutton(self.controller, text=el)
                chk.grid(row=row, column=0)
                chk.state(['!selected'])

    def set_list(self, check_list):
        self.clear_list()
        for row, el in enumerate(list(check_list.keys())):
                chk = ttk.Checkbutton(self.controller, text=el)
                chk.grid(row=row, column=0)
                #chk.deselect()

    def clear_list(self):
        for el in self.controller.winfo_children():
            el.destroy()

    def get_list(self):
        check = []

        for el in self.controller.winfo_children():
            if el.instate(['selected']):
                check.append(el['text'])

        return check


class FileManager(Frame):
    # todo: add a function to update the file manager
    SortDir = True

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        TopBar(parent, controller)

        file_manager = Frame(self, width=400)
        file_manager.pack(fill=BOTH, padx=100, pady=30)

        prev_frame = Frame(self, width=400, height=10)
        prev_frame.pack(fill=BOTH, padx=100, pady=20)

        self.preview = Text(prev_frame, bd=1, relief=SOLID, height=10)
        self.preview.pack(fill=BOTH, side=TOP)

        settings_frame = Frame(self, bg=settings.BG_COLOR, width=400)
        settings_frame.pack(fill=BOTH, padx=100)

        selection_frame = Frame(settings_frame, bg=settings.BG_COLOR)
        selection_frame.pack(fill=BOTH, padx=100)

        Label(selection_frame, text='Selected: ', bg=settings.BG_COLOR).grid(row=0, column=0)

        self.selected = Entry(selection_frame, width=80, relief=SOLID, bd=1)
        self.selected.grid(row=0, column=1, ipadx=1, ipady=1, sticky='w')

        self.find_button = Button(selection_frame, text='...', relief=SOLID, bd=1,
                                  command=lambda: self._select_file_with_manager())
        self.find_button.grid(row=0, column=2, padx=4, sticky='w')

        Label(selection_frame, text='Input format: ', bg=settings.BG_COLOR).grid(row=0, column=3, padx=5)
        self.input_selector = ttk.Combobox(selection_frame, values=["table", "dict", "lammps"], state='readonly')
        self.input_selector.current(0)
        self.input_selector.grid(row=0, column=4, sticky='w')

        Label(selection_frame, text='Filter width: ', bg=settings.BG_COLOR).grid(row=1, column=0)
        self.filter_width_entry = Spinbox(selection_frame, from_=0.1, to=10.0, increment=0.1, bd=1, relief=SOLID)
        self.filter_width_entry.grid(row=1, column=1, padx=2, sticky='w', pady=10)

        Label(selection_frame, text='Keys ', bg=settings.BG_COLOR).grid(row=2, column=0)
        check_frame = Frame(selection_frame, bg=settings.BG_COLOR, width=100)
        check_frame.grid(row=3, column=0, pady=10)

        self.check_list = CheckList(self, check_frame)

        enviroment_settings = Frame(selection_frame, bg=settings.BG_COLOR, width=200)
        enviroment_settings.grid(row=3, column=1)

        Label(enviroment_settings, text='Temperature: ', bg=settings.BG_COLOR).grid(row=0, column=0)
        self.temperature_entry = Spinbox(enviroment_settings, from_=0, to=100000, increment=0.1, bd=1, relief=SOLID)
        self.temperature_entry.grid(row=0, column=1, padx=2, sticky='w', pady=10)

        Label(enviroment_settings, text='Volume: ', bg=settings.BG_COLOR).grid(row=1, column=0)
        self.volume_entry = Entry(enviroment_settings, bd=1, relief=SOLID)
        self.volume_entry.grid(row=1, column=1, padx=2, sticky='w')

        Label(enviroment_settings, text='DT_FS: ', bg=settings.BG_COLOR).grid(row=2, column=0)
        self.DT_FS_entry = Entry(enviroment_settings, bd=1, relief=SOLID)
        self.DT_FS_entry.grid(row=2, column=1, padx=2, sticky='w', pady=10)

        self.auto_flag = Checkbutton(enviroment_settings, text='Auto')
        self.auto_flag.grid(row=3, column=0)

        self.start_button = Button(selection_frame, text='Start analysis', relief=SOLID, bd=1,
                                   command=self._start_analysis)
        self.start_button.grid(row=5, column=0, sticky='w')

        self._start_file_manager(file_manager)
        self._parse_files()

        # StatusFrame(controller, self)

    def _start_file_manager(self, parent):
        inner_frame = parent

        # create the tree and scrollbars
        self.headers = ('File name', 'File type', 'Size')
        self.file_list = ttk.Treeview(inner_frame, columns=self.headers,
                                      show='headings')

        ysb = ttk.Scrollbar(inner_frame, orient=VERTICAL, command=self.file_list.yview)
        self.file_list['yscroll'] = ysb.set

        # add scrollbars to frame
        self.file_list.grid(row=0, column=0, sticky=NSEW)
        ysb.grid(row=0, column=1, sticky='ns')

        self.file_list.bind('<<TreeviewSelect>>', self._select_file)
        self.file_list.bind('<Double-1>', self._start_analysis)
        # set frame resize priorities
        inner_frame.rowconfigure(0, weight=1)
        inner_frame.columnconfigure(0, weight=1)

    def _parse_files(self):

        files = os.listdir(settings.DATA_PATH)

        self.loaded_files = []

        # Get the info of each file
        for file in files:
            file_name, file_type = file.split('.')
            file_size = cu.get_file_size(os.path.join(settings.DATA_PATH, file))


            self.loaded_files.append((file_name, file_type, file_size))

        # Set column header
        for header in self.headers:
            self.file_list.heading(header, text=header.title(),
                                   command=lambda c=header: self._column_sort(c, FileManager.SortDir))

        # Populate the file manager
        for file in self.loaded_files:
            self.file_list.insert('', 'end', values=file)

            # and adjust column widths if necessary
            for index, value in enumerate(file):
                iwidth = Font().measure(value)
                if self.file_list.column(self.headers[index], 'width') < iwidth:
                    self.file_list.column(self.headers[index], width=iwidth)

    def _column_sort(self, column, descending=False):

        files = [(self.file_list.set(child, column), child) for child in self.file_list.get_children('')]

        files.sort(reverse=descending)
        for index, file in enumerate(files):
            self.file_list.move(file[1], '', index)

        # reverse sort
        FileManager.SortDir = not descending

    def _select_file(self, event):
        '''
        This function is called when a file in the listbox is selected.
        This function set the value of the entry to the path of the file.
        '''
        self.selected.delete(0, END)
        name = '.'.join(el for el in self.file_list.item(self.file_list.selection())['values'][:2])
        path = os.path.join(settings.DATA_PATH, name)
        self.selected.insert(0, path)

        with open(path, 'r') as file:
            lines = file.readlines()[0:settings.PREVIEW_LINES]
            self.preview.delete('1.0', END)
            self.preview.insert('1.0', lines)

        keys = cu.load_keys(path)
        self.check_list.set_list(keys)

    def _select_file_with_manager(self):
        '''
        This function allow the user to search in a more accurately way the file
        by using the OS manager.
        '''
        path = dialog.askopenfile(initialdir="/",
                                  title="Select file",
                                  filetypes=(("all files", "*.*"), ))

        # self.selected.delete(0, END)
        self.selected.insert(INSERT, path.name)

    def _start_analysis(self, ev=None):
        if self.selected.get():
            if os.path.exists(self.selected.get()):
                if self.selected.get().split('.')[-1] in settings.FILE_EXTENSIONS:
                    cu.CURRENT_FILE = self.selected.get()
                    # load_process = LoadingWindow(cu.get_file_size(cu.CURRENT_FILE))
                    if not cu.Data.loaded:
                        psd_filter_w = float(self.filter_width_entry.get())
                        cu.Data.psd_filter_width = psd_filter_w
                        keys = self.check_list.get_list()

                        if keys:

                            # if self.auto_flag.instate(['selected']):
                            #     temperature = None
                            #     volume = None
                            # else:
                            #     pass

                            temperature = float(self.temperature_entry.get())
                            volume = float(self.volume_entry.get())
                            DT_FS = float(self.DT_FS_entry.get())

                            if temperature > 0 and volume > 0 and DT_FS > 0:
                                cu.load_data(self.selected.get(),
                                    self.input_selector.get(),
                                    keys,
                                    temperature=temperature,
                                    units='metal',
                                    volume=volume,
                                    psd_filter_w=psd_filter_w,
                                    DT_FS=DT_FS,
                                    logs=ThermocepstrumGUI.frames[Cutter].logs)

                                ThermocepstrumGUI.show_frame(Cutter)
                                ThermocepstrumGUI.frame.update()
                            else:
                                msg.showerror('Value error', 'Temperature, volume and DT_FS can\'t be less than 0')
                        else:
                            msg.showerror('No keys selected', 'You must select almost one header key!')
                else:
                    msg.showerror('Invalid format!', 'The file that you have selected has an invalid format!')
            else:
                msg.showerror('File doesn\'t exists!', 'The file that you have selected doesn\'t exists!')
        else:
            msg.showerror('No file selected!', 'You must select a data file!')


class Cutter(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        TopBar(parent, controller)

        sections = Frame(self, bg=settings.BG_COLOR)
        sections.pack(side=LEFT, anchor='n')
        self.graph = GraphWidget(parent, sections, size=(7, 4), toolbar=True)

        self.slider_locked = False

        slider_frame = Frame(sections, bg=settings.BG_COLOR)
        slider_frame.pack(side=TOP, anchor='w', padx=115)

        self.slider = ttk.Scale(slider_frame, from_=0, to_=0.1, length=520)
        self.slider.grid(row=0, column=0)

        lock_slider = Button(slider_frame, command=lambda: self._lock_unlock_slider(),
                             bg=settings.BG_COLOR, bd=1, relief=SOLID)
        lock_slider.grid(row=0, column=1, padx=2)
        self.graph.attach_slider(self.slider)

        value_frame = Frame(sections, bg=settings.BG_COLOR)
        value_frame.pack(side=TOP)

        Label(value_frame, text='Selected value:', bg=settings.BG_COLOR).pack(side=TOP, pady=10)
        self.value_entry = Entry(value_frame, bd=1, relief=SOLID)
        self.value_entry.pack(side=LEFT)
        self.graph.attach_entry(self.value_entry)

        self.filter_width = Spinbox(value_frame, from_=0.1, to=10, increment=0.1, bd=1, relief=SOLID)
        self.filter_width.pack(side=LEFT)

        resample_button = Button(value_frame, text='Resample', bd=1, relief=SOLID,
                                 command=self.resample).pack(side=RIGHT, padx=10)
        button_frame = Frame(sections, bg=settings.BG_COLOR)
        button_frame.pack(pady=20)

        back_button = Button(button_frame, text='Back', bd=1, relief=SOLID, command=lambda: self.back())
        back_button.grid(row=0, column=0, sticky='w', padx=5)

        next_button = Button(button_frame, text='Next', bd=1, relief=SOLID, command=lambda: self.next())
        next_button.grid(row=0, column=1, sticky='w', padx=5)

        info_section = Frame(self, bg=settings.BG_COLOR)
        info_section.pack(side=RIGHT, anchor='n', pady=30, padx=20, fill='x', expand=True)


        self.logs = TextWidget(parent, info_section, 'Logs', 15)
        self.info = TextWidget(parent, info_section, 'Info', 10)
        # StatusFrame(parent, controller)

    def _lock_unlock_slider(self):
        if self.slider_locked:
            self.slider_locked = False
            self.slider.state(['!disabled'])
            self.value_entry.config(state=NORMAL)
        else:
            self.slider_locked = True
            self.slider.state(['disabled'])
            self.value_entry.config(state=DISABLED)

    def resample(self):
        cu.Data.fstar = float(self.value_entry.get())
        filter_width = float(self.filter_width.get())
        cu.Data.psd_filter_width=filter_width

        if cu.Data.fstar > 0:
            self.graph.add_graph(cu.gm.resample_current, 'resample', x=cu.Data.j, fstar_THz=cu.Data.fstar,
                                 PSD_FILTER_W=filter_width)
            self.graph.update_cut()
        else:
            msg.showwarning('Value error', 'F* must be greater than zero')

    def back(self):
        response = msg.askyesnocancel('Back to file manager?', "Save changes?\nIf reopen the same file \nthe values that you chosed will not be deleted!")

        if response:
            cu.Data.fstar = float(self.value_entry.get())
            cu.Data.loaded = True
            ThermocepstrumGUI.show_frame(FileManager)
        elif response == False:
            cu.Data.fstar = 0.0
            cu.Data.loaded = False
            ThermocepstrumGUI.show_frame(FileManager)
        else:
            pass

    def next(self):
        cu.Data.fstar = float(self.value_entry.get())
        ThermocepstrumGUI.show_frame(PStar)

    def update(self):
        self.graph.show(cu.gm.GUI_plot_periodogram, x=cu.Data.j, PSD_FILTER_W=cu.Data.psd_filter_width)
        # self.graph.cut_line = cu.Data.fstar


class PStar(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        TopBar(parent, controller)

        sections = Frame(self, bg=settings.BG_COLOR)
        sections.pack(side=LEFT, anchor='n')
        self.graph = GraphWidget(parent, sections, size=(7, 4), toolbar=True)

        value_frame = Frame(sections, bg=settings.BG_COLOR)
        value_frame.pack(side=TOP)

        Label(value_frame, text='Correction factor', bg=settings.BG_COLOR).pack(side=TOP, pady=10)
        self.value_entry = Spinbox(value_frame, bd=1, relief=SOLID, increment=1)
        self.value_entry.pack()

        self.increment = IntVar()
        Label(value_frame, text='Increment by', bg=settings.BG_COLOR).pack(side=TOP, pady=4)
        Radiobutton(value_frame, text='1', variable=self.increment, value=1,
                    bg=settings.BG_COLOR, command=self._change_increment).pack(side=LEFT, anchor='n', padx=2)
        Radiobutton(value_frame, text='10', variable=self.increment, value=10,
                    bg=settings.BG_COLOR, command=self._change_increment).pack(side=LEFT, anchor='n', padx=2)
        Radiobutton(value_frame, text='100', variable=self.increment, value=100,
                    bg=settings.BG_COLOR, command=self._change_increment).pack(side=LEFT, anchor='n', padx=2)

        button_frame = Frame(sections, bg=settings.BG_COLOR)
        button_frame.pack(pady=20)

        back_button = Button(button_frame, text='Back', bd=1, relief=SOLID, command=lambda: self.back())
        back_button.grid(row=0, column=0, sticky='w', padx=5)

        next_button = Button(button_frame, text='Next', bd=1, relief=SOLID, command=self._reload)
        next_button.grid(row=0, column=1, sticky='w', padx=5)

        info_section = Frame(self, bg=settings.BG_COLOR)
        info_section.pack(side=RIGHT, anchor='n', pady=30, padx=20, fill='x', expand=True)

        self.logs = TextWidget(parent, info_section, 'Logs', 15)
        self.info = TextWidget(parent, info_section, 'Info', 10)

    def back(self):
        ThermocepstrumGUI.show_frame(Cutter)

    def _get_pstar(self, aic_type='aic', Kmin_corrfactor=1.0):
        cu.Data.xf.cepstral_analysis(aic_type=aic_type, Kmin_corrfactor=Kmin_corrfactor)

    def _corr_factor(self):
        self.value_entry.config(from_=1.0, to=cu.Data.xf.Nfreqs)
        self.value_entry.delete(0, END)
        self.value_entry.insert(0, 1)

    def _change_increment(self):
        self.value_entry.config(increment=int(self.increment.get()))

    def _reload(self):
        self._get_pstar(aic_type='aic', Kmin_corrfactor=int(self.value_entry.get()))
        self.graph.add_graph(cu.gm.plot_cepstral_spectrum, 'cepstral', x=cu.Data.xf)
        self.graph.update_cut()

    def update(self):
        self._corr_factor()
        self.graph.show(cu.gm.GUI_plot_periodogram, x=cu.Data.j)
        self.graph.add_graph(cu.gm.resample_current, 'resample', x=cu.Data.j, fstar_THz=cu.Data.fstar,
                             PSD_FILTER_W=cu.Data.psd_filter_width)
        self._reload()
        self.graph.update_cut()


class Output(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        TopBar(parent, controller)


class LoadingWindow(Tk):

    def __init__(self, size, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)

        ThermocepstrumGUI.open_windows.insert(0, self)
        self.protocol('WM_DELETE_WINDOW', func=lambda: self.kill())
        frame = Frame(self, bg=settings.BG_COLOR)
        frame.pack(fill=BOTH, expand=1)

        Label(frame, text='Loading data').pack(side=TOP)
        Label(frame, text=f'{size} to load').pack(side=TOP)

        self.mainloop()

    def kill(self):
        index = ThermocepstrumGUI.open_windows.index(self)
        del ThermocepstrumGUI.open_windows[index]
        self.destroy()


def run():
    # Load data
    cu.load_path()
    # Run GUI
    app = ThermocepstrumGUI(version='0.0.1', dev_state='beta', last_release='dd/mm/yyyy')
    app.mainloop()


if __name__ == '__main__':
    run()
