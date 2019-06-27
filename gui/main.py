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
import glob
import matplotlib.pyplot as plt
from tkinter import *
from tkinter import ttk
from tkinter.ttk import Separator, Progressbar
from tkinter import messagebox as msg
import tkinter.filedialog as dialog
from tkinter.font import Font
from core import settings
import core.control_unit as cu
import core.gui_functions as guif


# Main app
class ThermocepstrumGUI(Tk):
    open_windows = []
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

        self.frames = {}

        for F in (FileManager, Cutter, PStar, Output):
            ThermocepstrumGUI.frame = F(container, self)
            ThermocepstrumGUI.frame.configure(bg=settings.BG_COLOR)

            self.frames[F] = ThermocepstrumGUI.frame
            ThermocepstrumGUI.frame.grid(row=0, column=0, sticky='nsew')

        self.show_frame(FileManager)

    def show_frame(self, frame):
        ThermocepstrumGUI.frame = self.frames[frame]
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
        status_frame.pack(fill='x')

        self.status = Label(status_frame, text=('Status: ' + settings.STATUS_NOW))
        self.status.pack(side=LEFT, padx=4, pady=2)


class FileManager(Frame):

    SortDir = True

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        TopBar(parent, controller)

        file_manager = Frame(self, width=400)
        file_manager.pack(fill=BOTH, padx=100, pady=50)

        selection_frame = Frame(self, bg=settings.BG_COLOR)
        selection_frame.pack(fill=BOTH, padx=100)

        Label(selection_frame, text='Selected: ', bg=settings.BG_COLOR).grid(row=0, column=0)

        self.selected_label = Entry(selection_frame, width=80, relief=SOLID, bd=1)
        self.selected_label.grid(row=0, column=1, ipadx=1, ipady=1)

        self.find_button = Button(selection_frame, text='...', relief=SOLID, bd=1,
                                  command=lambda: self._select_file_with_manager())
        self.find_button.grid(row=0, column=2, padx=4)

        self.start_button = Button(selection_frame, text='Start analysis', relief=SOLID, bd=1)
        self.start_button.grid(row=1, column=0, pady=20)

        self._start_file_manager(file_manager)
        self._parse_files()

        StatusFrame(parent, controller)

    def _start_file_manager(self, parent):
        inner_frame = ttk.Frame(parent)
        inner_frame.pack(side=TOP, fill=BOTH, expand=Y)

        # create the tree and scrollbars
        self.headers = ('File name', 'File type', 'Size')
        self.file_list = ttk.Treeview(columns=self.headers,
                                      show='headings')

        ysb = ttk.Scrollbar(orient=VERTICAL, command=self.file_list.yview)
        self.file_list['yscroll'] = ysb.set

        # add scrollbars to frame
        self.file_list.grid(row=0, column=0, sticky=NSEW, in_=inner_frame)
        ysb.grid(row=0, column=1, sticky='ns', in_=inner_frame)

        self.file_list.bind('<<TreeviewSelect>>', self._select_file)
        # set frame resize priorities
        inner_frame.rowconfigure(0, weight=1)
        inner_frame.columnconfigure(0, weight=1)

    def _parse_files(self):

        files = os.listdir(settings.DATA_PATH)

        self.loaded_files = []

        # Get the info of each file
        for file in files:
            file_name, file_type = file.split('.')
            file_size = os.path.getsize(os.path.join(settings.DATA_PATH, file))
            if file_size >= 1000000:
                file_size //= 1000000
                self.loaded_files.append((file_name, file_type, f"{file_size} MB"))
            else:
                file_size //= 1000
                self.loaded_files.append((file_name, file_type, f"{file_size} KB"))

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
        self.selected_label.delete(0, END)
        name = '.'.join(el for el in self.file_list.item(self.file_list.selection())['values'][:2])
        path = os.path.join(settings.DATA_PATH, name)
        self.selected_label.insert(0, path)

    def _select_file_with_manager(self):
        '''
        This function allow the user to search in a more accurately way the file
        by using the OS manager.
        '''
        path = dialog.askopenfile(initialdir="/",
                                  title="Select file",
                                  filetypes=(("all files", "*.*"), ))

        self.selected_label.delete(0, END)
        self.selected_label.insert(0, path.name)


class Cutter(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        TopBar(parent, controller)


class PStar(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        TopBar(parent, controller)


class Output(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        TopBar(parent, controller)


def run():
    # Load data
    cu.load_path()
    # Run GUI
    app = ThermocepstrumGUI(version='0.0.1', dev_state='beta', last_release='dd/mm/yyyy')
    app.mainloop()


if __name__ == '__main__':
    run()
