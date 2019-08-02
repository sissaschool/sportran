# -*- coding: future_fstrings -*-

'''
--------------------------------------------
    Thermocepstrum graphic user interface

    Version: 0.0.1
    Release state: Beta
    Last update: 08/2019

    Main Developer:     Sebastiano Bisacchi
    Other developer:    Riccardo   Bertossa
--------------------------------------------

This file contains the GUI of the Thermocepstrum project developed at SISSA
'''
# todo: Put an accurate description of the project?

import os

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import NavigationToolbar2
from matplotlib.figure import Figure
import matplotlib.patches as patches

from tkinter import *
from tkinter import ttk
from tkinter.ttk import Separator
from tkinter import messagebox as msg
import tkinter.filedialog as dialog
from tkinter.font import Font

from thermocepstrum_gui.core import settings
import thermocepstrum_gui.core.control_unit as cu
try:
    from thermocepstrum.utils.utils import PrintMethod
except ImportError:
    from thermocepstrum_gui.utils.utils import PrintMethod

# Verify that thermocepstrum is installed
try:
    import thermocepstrum
except ImportError:
    raise ImportError('Couldn\'t find thermocepstrum')

# Init print method
log = PrintMethod()


# Main app
class ThermocepstrumGUI(Tk):
    open_windows = []
    frames = []
    frame = None
    root = None
    container = None

    def __init__(self, version, dev_state, last_release, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)

        ThermocepstrumGUI.root = self
        self.option_add("*Font", "{} {}".format(settings.FONT, settings.FONT_SIZE))
        self.option_add("*Background", "{}".format(settings.BG_COLOR))
        self.option_add("*selectBackground", "light blue")
        self.option_add("*selectForeground", "black")

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

        self.resizable(settings.X_RESIZE, settings.Y_RESIZE)
        self.protocol('WM_DELETE_WINDOW', func=lambda: cu.secure_exit(ThermocepstrumGUI))

        # Creating the main frame
        container = Frame(self)
        ThermocepstrumGUI.container=container
        container.grid(row=0, column=0, sticky='nsew')
        #container=self
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        #container.pack(side=TOP, fill=BOTH, expand=True)
        container.grid_rowconfigure(0, weight=10)
        container.grid_columnconfigure(0, weight=10)

        # Setting up multiple window system

        ThermocepstrumGUI.frames = {}

        self.topbar = TopBar(self,self)

        ThermocepstrumGUI.root.grid_propagate(True)
        for F in ( HeaderSelector, OtherVariables, FStarSelector, PStarSelector, FileManager):
            ThermocepstrumGUI.frame = F(container, self)
            #ThermocepstrumGUI.frame.grid(row=0, column=0, sticky='nsew')
            #ThermocepstrumGUI.frame.grid_forget()

            ThermocepstrumGUI.frames[F] = ThermocepstrumGUI.frame

        self.show_frame(FileManager)

    @staticmethod
    def show_frame(frame):
        #ThermocepstrumGUI.frame.grid_forget()
        #ThermocepstrumGUI.container.grid_forget()
        #ThermocepstrumGUI.container.grid(row=0, column=0, sticky='nsew')
        ThermocepstrumGUI.container.grid_rowconfigure(0, weight=0)
        ThermocepstrumGUI.container.grid_columnconfigure(0, weight=0)

        #ThermocepstrumGUI.container.grid(row=0, column=0, sticky='nsew')
        ThermocepstrumGUI.frame = ThermocepstrumGUI.frames[frame]
        #ThermocepstrumGUI.frame.__init__(ThermocepstrumGUI.container,ThermocepstrumGUI.container)
        ThermocepstrumGUI.frame.grid(row=0, column=0, sticky='nsew')
        ThermocepstrumGUI.container.grid_rowconfigure(0, weight=1)
        ThermocepstrumGUI.container.grid_columnconfigure(0, weight=1)
        ThermocepstrumGUI.frame.tkraise()
        ThermocepstrumGUI.frame.update()

class TopBar(Frame):

    show_logs = None
    show_info = None

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

#        self.parent = parent

        TopBar.show_logs = BooleanVar()
        TopBar.show_info = BooleanVar()

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

        # Create the view section of the top menu
        view_menu = Menu(top_menu, tearoff=False)
        top_menu.add_cascade(label='View', menu=view_menu)
        view_menu.add_checkbutton(label='Show logs', variable=TopBar.show_logs, onvalue=1, offvalue=0, command=self._update_window)
        view_menu.add_checkbutton(label='Show info', variable=TopBar.show_info, onvalue=1, offvalue=0, command=self._update_window)

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

    def _update_window(self):
        ThermocepstrumGUI.frame.update()


class StatusFrame(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        status_frame = Frame(controller)
        status_frame.pack(fill='x', side=BOTTOM)

        self.status = Label(status_frame, text=('Status: ' + settings.STATUS_NOW))
        self.status.pack(side=LEFT, padx=4, pady=2)


class GraphToolbar(NavigationToolbar2):

    def __init__(self, parent, controller):
        self.toolitems = (
            ('Home', 'Reset the view', 'home', 'home'),
            ('Back', 'Go back to the previous view', 'back', 'back'),
            ('Forward', 'Go back to the last view', 'forward', 'forward'),
            (None, None, None, None),
            ('Pan', 'Use to navigate in the graph', 'move', 'pan'),
            ('Zoom', 'Use to zoom a section of the graph', 'zoom_to_rect', 'zoom'),
            (None, None, None, None),
            (None, None, None, None),
            ('Save', 'sollemnes in futurum', 'filesave', 'save_figure'),
        )

        NavigationToolbar2.__init__(self, parent, controller)


class GraphWidget(Frame):

    def __init__(self, parent, controller, size=(4, 4), type=111, toolbar=False):
        Frame.__init__(self, parent)

        self.title = ''
        self.size = size
        self.type = type

        self.cut_line = 0
        self.max_x = 1
        self.max_y = 1
        self.new_view_x = 1

        self.f = Figure(figsize=self.size, dpi=100)
        self.graph = self.f.add_subplot(self.type)

        self.canvas = FigureCanvasTkAgg(self.f, controller)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, anchor='w', padx=10, fill=BOTH, expand=1)

        self.func = None

        self.other_graph = []

        self.slider = None
        self.entry = None
        self.show_selected_area = False
        self.plot_call = None
        self.plot_call_kwargs = None

        if toolbar:
            toolbar = NavigationToolbar2Tk(self.canvas, controller)
            toolbar.pack(side=TOP, pady=10, padx=50, fill=BOTH, expand=1)
            toolbar.update()
            self.canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

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
        cu.set_graph(self.graph, func, **kwargs)
        self.max_x = self.get_max_x()
        self.max_y = self.get_max_y()
        if self.slider:
            if self.show_selected_area:
                self.change_view()
            else:
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
            self.graph.set_ylim([0, self.max_y])
            if self.show_selected_area:
                self.graph.set_xlim([0, self.new_view_x])
        self.canvas.draw()

    def get_graph(self):
        return self.graph

    def set_plot_call(self, f, **f_args):
        self.plot_call = f
        self.plot_call_kwargs = f_args

    def attach_slider(self, slider):
        self.slider = slider
        self.slider.config(command=self._on_slider_change, to_=self.get_max_x())

    def change_view(self):
        if self.show_selected_area:
            self.new_view_x = self.cut_line
            self.slider.config(to_=self.new_view_x)
        else:
            self.slider.config(to_=self.max_x)

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

    def __init__(self, parent, controller, title, height, width):
        Frame.__init__(self, parent, controller)

        text_frame = LabelFrame(controller, text=title, bd=1, relief=SOLID)
        text_frame.pack(side=TOP)

        self.text_box = Text(text_frame, height=height, width=width, bd=0)
        self.text_box.pack(padx=4, pady=4)
        self.text_box.config(state=DISABLED)

    def clear(self):
        self.text_box.config(state=NORMAL)
        self.text_box.delete('0.1', END)
        self.text_box.config(state=DISABLED)

    def write(self, text):
        self.text_box.config(state=NORMAL)
        self.text_box.insert(INSERT, str(text)+'\n')
        self.text_box.config(state=DISABLED)


class CheckList(Frame):

    def __init__(self, parent, controller, check_list=dict()):
        Frame.__init__(self, parent, controller)

        self.controller = controller
        self.combo_func = None
        if list:
            self.set_list(check_list)

    def set_list(self, check_list):
        self.clear_list()
        for row, el in enumerate(list(check_list.keys())):
                frame = Frame(self.controller)
                frame.grid(row=row, column=0, sticky='we', pady=2)
                Label(frame, text=el, font="{} 12 bold".format(settings.FONT)).grid(row=0, column=0)
                cmb = ttk.Combobox(frame, values=["None", "Energy current", "Other current", "Temperature"], state='readonly', width=12)
                cmb.bind('<<ComboboxSelected>>', self.combo_func)
                cmb.current(0)
                cmb.grid(row=0, column=1, sticky='e')

    def clear_list(self):
        for el in self.controller.winfo_children():
            el.destroy()

    def get_list(self):
        check = []
        combo = []

        for el in self.controller.winfo_children():
            header = el.winfo_children()[0]
            cmb = el.winfo_children()[1]

            if cmb.get() is not None:
                check.append(header['text'])
                combo.append(cmb.get())

        return check, combo

    def attach_function_on_combo(self, func):
        self.combo_func = func


class ScrollFrame(Frame):

    def __init__(self, parent, controller, width=0, height=0, bd=0):
        Frame.__init__(self, parent)

        bgcol= settings.BG_COLOR
        bgcol2=settings.BG_COLOR

        if width or height:
            self.canvas = Canvas(controller, bd=bd, relief=SOLID, highlightthickness=0, bg=bgcol, width=width, height=height)
        else:
            self.canvas = Canvas(controller, bd=bd, relief=SOLID, highlightthickness=0, bg=bgcol)
        self.viewPort = Frame(self.canvas, background=bgcol2)
        self.vsb = Scrollbar(controller, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)

        #self.vsb.pack(side="right", fill="y")
        #self.canvas.pack(fill=BOTH, expand=True)
        self.vsb.grid(row=0,column=1,sticky='nse')
        self.canvas.grid(row=0, column=0, sticky='nsew')
        controller.rowconfigure(0,weight=10)
        controller.columnconfigure(0,weight=10)
        controller.columnconfigure(1,weight=0)
        #self.viewPort.grid(row=0, column=0, sticky='nsew')
        self.canvas.create_window(0, 0, window=self.viewPort, tags="self.viewPort")
        #self.viewPort.pack(fill=BOTH, side='left', expand=True)
        self.viewPort.grid(row=0, column=0, sticky='nsew')
        self.canvas.grid_rowconfigure(0, weight=1)
        self.canvas.grid_columnconfigure(0, weight=1)
        #self.viewPort.grid_rowconfigure(0, weight=1)
        #self.viewPort.grid_columnconfigure(0, weight=1)


        self.viewPort.bind("<Configure>", self.on_frame_configure)

        self.viewPort.bind("<MouseWheel>", self._on_mousewheel)
        self.viewPort.bind("<Button-4>", self._on_mousewheel)
        self.viewPort.bind("<Button-5>", self._on_mousewheel)

    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(-1 * int((event.delta / 120)), "units")

    def on_frame_configure(self, event):
        '''
        Reset the scroll region to encompass the inner frame
        '''

        self.canvas.configure(scrollregion=self.canvas.bbox('all'))


class FileManager(Frame):
    SortDir = True

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)
        #TopBar(parent, controller)

        self.main_frame = ScrollFrame(self, self)

        file_manager = Frame(self.main_frame.viewPort)
        #file_manager.pack(fill=BOTH, expand=True, padx=100, pady=30)
        file_manager.grid(column=0,row=0,sticky='nswe')


        prev_frame = Frame(self.main_frame.viewPort, height=10)
        #prev_frame.pack(fill=BOTH, expand=True, padx=100, pady=20)
        prev_frame.grid(column=0,row=1,sticky='nswe')

        Label(prev_frame, text='File preview', font='Arial 11 bold').pack(side=TOP, anchor='w',fill=BOTH,expand=True)
        self.preview = Text(prev_frame, bd=1, relief=SOLID, height=10)
        self.preview.pack(fill=BOTH, expand=True, side=BOTTOM)
        self.preview.config(state=NORMAL)

        #settings_frame = Frame(self.main_frame.viewPort)
        #settings_frame.pack(fill=BOTH, expand=True, padx=100)
        #settings_frame.grid(column=0,row=2,sticky='nswe')

        selection_frame = Frame(self.main_frame.viewPort)
        #selection_frame.pack(fill=BOTH, expand=True)
        selection_frame.grid(column=0,row=2,sticky='nswe')

        Label(selection_frame, text='Selected: ').grid(row=0, column=0, sticky='e')

        self.selected = Entry(selection_frame, width=60, relief=SOLID, bd=1)
        self.selected.grid(row=0, column=1, ipadx=1, ipady=1, sticky='we')

        self.find_button = Button(selection_frame, text='...', relief=SOLID, bd=1,
                                  command=lambda: self._select_file_with_manager())
        self.find_button.grid(row=0, column=2, padx=4, sticky='we')

        Label(selection_frame, text='Input format: ').grid(row=0, column=3, padx=5, sticky='we')
        self.input_selector = ttk.Combobox(selection_frame, values=["table", "dict", "lammps"], state='readonly',width=10)
        self.input_selector.current(0)
        self.input_selector.grid(row=0, column=4, sticky='w')


        self.next_button = Button(selection_frame, text='Next', relief=SOLID, bd=1,
                                   command=lambda: self.next())
        self.next_button.grid(row=1, column=0, sticky='w', pady=10)
        selection_frame.rowconfigure(0,weight=1)
        selection_frame.rowconfigure(1,weight=1)

        self._start_file_manager(file_manager)
        self._parse_files()

        for i,w,m in zip(range(0,4),[2,6,1,2,2],[70,70,25,90,50]):
            selection_frame.columnconfigure(i,weight=w,minsize=m)
        for i, w, m in zip(range(0, 2), [7, 7, 4], [70, 70, 40]):
            self.main_frame.viewPort.rowconfigure(i, weight=w, minsize=m)
        self.main_frame.viewPort.columnconfigure(0, weight=1, minsize=300)



    def _start_file_manager(self, parent):
        inner_frame = parent

        Label(inner_frame, text='Select a file', font='Arial 11 bold').grid(row=0, column=0, sticky='nswe')
        # create the tree and scrollbars
        self.headers = ('File name', 'File type', 'Size')
        self.file_list = ttk.Treeview(inner_frame, columns=self.headers,
                                      show='headings')

        ysb = ttk.Scrollbar(inner_frame, orient=VERTICAL, command=self.file_list.yview)
        self.file_list['yscroll'] = ysb.set

        # add scrollbars to frame
        self.file_list.grid(row=1, column=0, sticky='nswe')
        ysb.grid(row=1, column=1, sticky='nswe')

        self.file_list.bind('<<TreeviewSelect>>', self._select_file)
        self.file_list.bind('<Double-1>', self.next)

        # set frame resize priorities
        inner_frame.rowconfigure(0, weight=1,minsize=15)
        inner_frame.rowconfigure(1, weight=4,minsize=40)
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

        self.load_file_settings(path)

    def _select_file_with_manager(self):
        '''
        This function allow the user to search in a more accurately way the file
        by using the OS manager.
        '''
        path = dialog.askopenfile(initialdir=os.getcwd(),
                                  title="Select file",
                                  filetypes=(("all files", "*.*"), ))

        # self.selected.delete(0, END)

        for item in self.file_list.get_children():
            self.file_list.delete(item)

        settings.DATA_PATH = os.path.dirname(path.name)
        self._parse_files()

        self.selected.insert(INSERT, path.name)
        self.load_file_settings(path.name)
    
    def load_file_settings(self, path):
        with open(path, 'r') as file:
            lines = file.readlines()[0:settings.PREVIEW_LINES]
            # todo: remove {}

            self.preview.config(state=NORMAL)
            self.preview.delete('1.0', END)
            self.preview.insert('1.0', lines)
            self.preview.config(state=DISABLED)

    def next(self, ev=None):
        if self.selected.get():
            if os.path.exists(self.selected.get()):
                if self.selected.get().split('.')[-1] in settings.FILE_EXTENSIONS:
                    cu.Data.CURRENT_FILE = self.selected.get()
                    cu.Data.inputformat = self.input_selector.get()

                    ThermocepstrumGUI.show_frame(HeaderSelector)
                else:
                    msg.showerror('Invalid format!', 'The file that you have selected has an invalid format!')
            else:
                msg.showerror('File doesn\'t exists!', 'The file that you have selected doesn\'t exists!')
        else:
            msg.showerror('No file selected!', 'You must select a data file!')

    def update(self):
        self.main_frame.canvas.yview_moveto(0)


        super().update()


class HeaderSelector(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        #TopBar(parent, controller)
        self.main_frame = ScrollFrame(self, self)

        header_frame = self.main_frame.viewPort #Frame(self.main_frame.viewPort,bg='#FF0000')
        #header_frame.pack(fill=BOTH, expand=True, padx=20, pady=20)
        header_frame.grid(row=0, column=0, sticky='nswe')
        #self.main_frame.viewPort.grid_rowconfigure(0, weight=1)
        #self.main_frame.viewPort.grid_columnconfigure(0, weight=1)

        Label(header_frame, text='Define the use of the headers').grid(row=0, column=0, sticky='w')
        definitions_frame = Frame(header_frame)
        definitions_frame.grid(row=1, column=1, sticky='nswe', padx=20)
        Label(definitions_frame, text='None: the header will not be used').grid(row=0, column=0, sticky='wns')
        Label(definitions_frame, text='Temperature: the header that will be used to calculate the temperature').grid(row=1, column=0, sticky='wns')
        # todo: put definition
        Label(definitions_frame, text='Energy current: put definition ').grid(row=2, column=0, sticky='wns')
        Label(definitions_frame, text='Other current: put definition').grid(row=3, column=0, sticky='wns')
        definitions_frame.columnconfigure(0,weight=1)
        for i in range(0,3):
            definitions_frame.rowconfigure(i,weight=1)

        header_list_frame = Frame(header_frame)
        header_list_frame.grid(row=1, column=0, sticky='nswe', pady=10)

        scrollable_header_list = ScrollFrame(self.main_frame, header_list_frame, bd=1)
        self.check_list = CheckList(scrollable_header_list, scrollable_header_list.viewPort)

        button_frame = Frame(header_frame)
        button_frame.grid(row=2, column=0, sticky='w')


        header_frame.rowconfigure(0,weight=1)
        header_frame.rowconfigure(1,weight=10)
        header_frame.rowconfigure(2,weight=1)
        header_frame.columnconfigure(0,weight=1)
        header_frame.columnconfigure(1,weight=1)

        Button(button_frame, text='Back', bd=1, relief=SOLID, font='Arial 12',
               command=lambda: self.back()).grid(row=0, column=0)

        Button(button_frame, text='Next', bd=1, relief=SOLID, font='Arial 12',
               command=lambda: self.next()).grid(row=0, column=1, padx=5)

    def next(self):
        keys, description = self.check_list.get_list()
        if 'Energy current' in description:
            if description.count('Energy current') == 1:
                if description.count('Temperature') <= 1:
                    cu.Data.keys = keys
                    cu.Data.description = description

                    ThermocepstrumGUI.show_frame(OtherVariables)
                else:
                    msg.showerror('Value error', 'You can\'t assign more than one time the value "Temperature"')
            else:
                msg.showerror('Value error', 'You must assign only one "Energy current" value')
        else:
            msg.showerror('No keys selected', 'You must select almost one header key!')

    def back(self):
        cu.Data.keys = None
        cu.Data.description = None

        ThermocepstrumGUI.show_frame(FileManager)

    def update(self):
        keys = cu.load_keys(cu.Data.CURRENT_FILE)
        self.check_list.set_list(keys)
        super().update()


class OtherVariables(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        #TopBar(parent, controller)

        self.main_frame = ScrollFrame(controller, self)

        variable_frame = Frame(self.main_frame.viewPort)
        #variable_frame.pack(fill=BOTH, expand=True, padx=100)
        variable_frame.grid(column=0,row=0,sticky='nswe')

        Label(variable_frame, text='Set variables', font='Arial 11 bold').grid(row=0, column=0, pady=2, sticky='w')

        Label(variable_frame, text='Environment variables', font='Arial 11').grid(row=1, column=0, pady=20)
        Separator(variable_frame, orient=HORIZONTAL).grid(row=2, sticky='we', columnspan=3)

        Label(variable_frame, text='Temperature: ').grid(row=3, column=0, sticky='w')
        self.temperature_entry = Spinbox(variable_frame, from_=0, to=100000, increment=0.1, bd=1, relief=SOLID)
        self.temperature_entry.grid(row=3, column=1, padx=2, sticky='w', pady=10)
        self.temp_advertise = Label(variable_frame, text='', font='Arial 10')
        self.temp_advertise.grid(row=3, column=2, sticky='w')

        Label(variable_frame, text='Volume: ').grid(row=4, column=0, sticky='w')
        self.volume_entry = Entry(variable_frame, bd=1, relief=SOLID)
        self.volume_entry.grid(row=4, column=1, padx=2, sticky='w')

        Label(variable_frame, text='DT_FS: ').grid(row=5, column=0, sticky='w')
        self.DT_FS_entry = Entry(variable_frame, bd=1, relief=SOLID)
        self.DT_FS_entry.grid(row=5, column=1, padx=2, sticky='w', pady=10)

        Label(variable_frame, text='Filter variables', font='Arial 11').grid(row=6, column=0, pady=20, sticky='w')
        Separator(variable_frame, orient=HORIZONTAL).grid(row=7, sticky='we', columnspan=3)

        Label(variable_frame, text='Filter width: ').grid(row=8, column=0, sticky='w')
        self.filter_width_entry = Spinbox(variable_frame, from_=0.1, to=10.0, increment=0.1, bd=1, relief=SOLID)
        self.filter_width_entry.grid(row=8, column=1, padx=2, sticky='w', pady=10)

        button_frame = Frame(variable_frame)
        button_frame.grid(row=9, column=0, sticky='ws')

        Button(button_frame, text='Back', bd=1, relief=SOLID, font='Arial 12',
               command=lambda: self.back()).grid(row=0, column=0)

        Button(button_frame, text='Next', bd=1, relief=SOLID, font='Arial 12',
               command=lambda: self.next()).grid(row=0, column=1, padx=5)

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
            if temperature < 0 and not 'Temperature' in cu.Data.description:
                msgs.append('Temperature can\'t be less than 0')
                er = True
        else:
            msgs.append('Temperature can\'t be void')

        if volume:
            if volume < 0:
                msgs.append('Volume can\'t be less than 0')
                er = True
        else:
            msgs.append('Volume can\'t be void')
            er = True

        if DT_FS:
            if DT_FS < 0:
                msgs.append('DT_FS can\'t be less than 0')
                er = True
        else:
            msgs.append('DT_FS can\'t be void')
            er = True

        if psd_filter_width:
            if psd_filter_width < 0:
                msgs.append('Filter width can\'t be less than 0')
                er = True
        else:
            msgs.append('Filter width can\'t be void')
            er = True

        if not er:
            cu.Data.psd_filter_width = psd_filter_width
            cu.load_data(cu.Data.CURRENT_FILE,
                         cu.Data.inputformat,
                         cu.Data.keys,
                         descriptions=cu.Data.description,
                         temperature=temperature,
                         units='metal',
                         volume=volume,
                         psd_filter_w=cu.Data.psd_filter_width,
                         DT_FS=DT_FS)

            ThermocepstrumGUI.show_frame(FStarSelector)
        else:
            ermsg = '\n'.join(mser for mser in msgs)
            msg.showerror('Value error', ermsg)

    def back(self):
        cu.Data.psd_filter_width = 0.1

        ThermocepstrumGUI.show_frame(HeaderSelector)

    def update(self):

        self.temperature_entry.config(state=NORMAL)
        self.temperature_entry.config(value=cu.Data.temperature)

        if 'Temperature' in cu.Data.description:
            self.temperature_entry.config(state=DISABLED)
            self.temp_advertise.config(text='The temperature will be automatically calculated')
        else:
            self.temperature_entry.config(state=NORMAL)
            self.temp_advertise.config(text='')

        self.volume_entry.delete(0, END)
        self.volume_entry.insert(0, cu.Data.volume)
        self.DT_FS_entry.delete(0, END)
        self.DT_FS_entry.insert(0, cu.Data.DT_FS)
        self.filter_width_entry.config(value=cu.Data.psd_filter_width)
        self.main_frame.viewPort.columnconfigure(0,weight=1)
        self.main_frame.viewPort.rowconfigure(0,weight=1)
        super().update()


class FStarSelector(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        #TopBar(parent, controller)

        self.parent = parent
        main_frame = self #Frame(self)
        #main_frame.pack(expand=True, fill=BOTH)
        main_frame.grid(column=0,row=0,sticky='nswe')

        sections = Frame(main_frame)
        sections.grid(row=0, column=0, sticky='nswe')

        self.graph = GraphWidget(sections, sections, size=(7, 4), toolbar=True)
        self.graph.pack(side=TOP, anchor='w', padx=10, fill=BOTH, expand=1)

        self.slider_locked = False

        slider_frame = Frame(sections)
        slider_frame.pack(side=TOP, anchor='w', padx=80, fill=BOTH, expand=1)

        self.slider = ttk.Scale(slider_frame, from_=0, to_=0.1)#, length=520)
        self.slider.grid(row=0, column=0,sticky='we')
        slider_frame.columnconfigure(0,weight=9)

        lock_slider = Button(slider_frame, command=lambda: self._lock_unlock_slider(),
                             bd=1, relief=SOLID)
        lock_slider.grid(row=0, column=1, padx=2)
        slider_frame.columnconfigure(1,weight=1)
        self.graph.attach_slider(self.slider)

        self.change_view_button = Button(slider_frame, text='Zoom-in', command=lambda: self._change_view())
        self.change_view_button.grid(row=0, column=2, sticky='w')
        slider_frame.columnconfigure(2,weight=1)

        value_frame = Frame(sections)
        value_frame.pack(side=LEFT, pady=10, padx=20,fill=BOTH, expand=1)

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

        main_frame.columnconfigure(0,weight=3)
        main_frame.columnconfigure(1,weight=1)
        main_frame.rowconfigure(0,weight=1)
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

    def back(self):
        response = msg.askyesnocancel('Back to file manager?', "Save changes?\nIf reopen the same file \nthe values that you chosed will not be deleted!")

        log.set_func(None)
        if response:
            cu.Data.fstar = float(self.value_entry.get())
            cu.Data.loaded = True
            ThermocepstrumGUI.show_frame(FileManager)
        elif response == False:
            cu.Data.fstar = 0.0
            cu.Data.loaded = False
            cu.Data.temperature = 0.0
            cu.Data.volume = 0.0
            cu.Data.DT_FS = 0.0
            cu.Data.psd_filter_width = 0.1

            self.graph.other_graph.clear()
            self.graph.graph.clear()
            self.graph.cut_line = 0
            ThermocepstrumGUI.show_frame(FileManager)
        else:
            pass

    def next(self):
        self.resample()
        cu.Data.fstar = cu.Data.xf.Nyquist_f_THz
        ThermocepstrumGUI.show_frame(PStarSelector)

    def update(self):
        self.graph.show(cu.gm.GUI_plot_periodogram, x=cu.Data.j, PSD_FILTER_W=cu.Data.psd_filter_width)
        self._init_output_frame()
        if self.info:
            cu.update_info(self.info)
        if self.logs:
            log.set_func(self.logs.write)
        # self.graph.cut_line = cu.Data.fstar
        super().update()


class PStarSelector(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        self.parent = parent

        #TopBar(parent, controller)

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

    def back(self):
        ThermocepstrumGUI.show_frame(FStarSelector)

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
        super().update()


class LoadingWindow(Tk):

    def __init__(self, size, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)

        ThermocepstrumGUI.open_windows.insert(0, self)
        self.protocol('WM_DELETE_WINDOW', func=lambda: self.kill())
        frame = Frame(self)
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
    log.set_method('other')
    run()
