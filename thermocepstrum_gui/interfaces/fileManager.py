"""
--------------------------------------------
    Thermocepstrum graphic user interface

    Interface file - file manager interface
--------------------------------------------

This file contains the file manager interface and functionality.
"""

import os

from tkinter import messagebox as msg
import tkinter.filedialog as fdialog
from tkinter.font import Font

from thermocepstrum_gui.utils.custom_widgets import *
from thermocepstrum_gui.assets import LANGUAGES


class FileManager(Frame):
    """
    This file manager is used to select the data files and see
    their preview.
    """
    SortDir = True

    def __init__(self, parent, main):
        Frame.__init__(self, parent)

        # Setup some variables to deal with
        # the other interfaces and the window
        self.main = main
        self.next_frame = None

        # Create the main frame that will contains all the
        # widgets of the file manager.
        self.main_frame = ScrollFrame(self, self)

        # Setup the file manager
        file_manager = Frame(self.main_frame.viewPort)
        file_manager.grid(column=0, row=0, sticky='nswe', padx=20, pady=5)

        # Setup the preview screen
        prev_frame = Frame(self.main_frame.viewPort, height=10)

        prev_frame.grid(column=0, row=1, sticky='nswe', padx=20)

        # Setup some widgets
        Label(prev_frame, text=LANGUAGES[settings.LANGUAGE]['f_prv'],
              font='Arial 11 bold').pack(side=TOP, anchor='w', fill=BOTH, expand=True)
        self.preview = Text(prev_frame, bd=1, relief=SOLID, height=10)
        self.preview.pack(fill=BOTH, expand=True, side=BOTTOM)
        self.preview.config(state=NORMAL)

        selection_frame = Frame(self.main_frame.viewPort)
        selection_frame.grid(column=0, row=2, sticky='nswe', padx=20, pady=5)

        Label(selection_frame, text=LANGUAGES[settings.LANGUAGE]['slct']).grid(row=0, column=0, sticky='w')

        self.selected = Entry(selection_frame, width=60, relief=SOLID, bd=1)
        self.selected.grid(row=0, column=1, ipadx=1, ipady=1, sticky='we')

        self.find_button = Button(selection_frame, text='...', relief=SOLID, bd=1,
                                  command=lambda: self._select_file_with_manager())
        self.find_button.grid(row=0, column=2, padx=4, sticky='we')

        Label(selection_frame, text=LANGUAGES[settings.LANGUAGE]['in_frm']).grid(row=0, column=3, padx=5, sticky='we')
        self.input_selector = ttk.Combobox(selection_frame, values=['table', 'dict', 'lammps'], state='readonly',
                                           width=10)
        self.input_selector.current(0)
        self.input_selector.grid(row=0, column=4, sticky='w')
        self.next_button = Button(selection_frame, text=LANGUAGES[settings.LANGUAGE]['next'], relief=SOLID, bd=1,
                                  command=lambda: self.next(), width=10)
        self.next_button.grid(row=1, column=0, pady=20, sticky='w')
        selection_frame.rowconfigure(0, weight=1)
        selection_frame.rowconfigure(1, weight=1)

        # Start the file manager and init
        self._start_file_manager(file_manager)
        self._parse_files()

        # Generate the layout
        for i, w, m in zip(range(0, 4), [2, 6, 1, 2, 2], [70, 70, 25, 90, 50]):
            selection_frame.columnconfigure(i, weight=w, minsize=m)
        for i, w, m in zip(range(0, 2), [7, 7, 4], [70, 70, 40]):
            self.main_frame.viewPort.rowconfigure(i, weight=w, minsize=m)
        self.main_frame.viewPort.columnconfigure(0, weight=1, minsize=300)

    def _start_file_manager(self, parent):
        """
        This function attach to the parent frame
        the file manager widget and setup his layout.

        :param parent: The frame to attach the file manager.
        """
        inner_frame = parent

        Label(inner_frame, text=LANGUAGES[settings.LANGUAGE]['stp1'],
              font='Arial 11 bold').grid(row=0, column=0, sticky='nswe')

        # create the tree and scrollbars
        self.headers = (LANGUAGES[settings.LANGUAGE]['f_nm'], LANGUAGES[settings.LANGUAGE]['f_tp'],
                        LANGUAGES[settings.LANGUAGE]['size'])
        self.file_list = ttk.Treeview(inner_frame, columns=self.headers, show='headings')

        ysb = ttk.Scrollbar(inner_frame, orient=VERTICAL, command=self.file_list.yview)
        self.file_list['yscroll'] = ysb.set

        # add scrollbars to frame
        self.file_list.grid(row=1, column=0, sticky='nswe')
        ysb.grid(row=1, column=1, sticky='nswe')

        self.file_list.bind('<<TreeviewSelect>>', self._select_file)
        self.file_list.bind('<Double-1>', self.next)

        # set frame resize priorities
        inner_frame.rowconfigure(0, weight=1, minsize=15)
        inner_frame.rowconfigure(1, weight=4, minsize=40)
        inner_frame.columnconfigure(0, weight=1)

    def _parse_files(self):
        """
        This function is used to parse the files in
        the selected directory and display them in the file manager.
        """
        global settings
        try:
            files = os.listdir(settings.DATA_PATH)
        except:
            settings.DATA_PATH = './'
            files = os.listdir(settings.DATA_PATH)

        self.loaded_files = []

        # Get the info of each file
        for file in files:
            file_component = file.split('.')
            if len(file_component) == 2:
                file_name, file_type = file.split('.')
                file_size = cu.get_file_size(os.path.join(settings.DATA_PATH, file))
                self.loaded_files.append((file_name, file_type, file_size))
            else:
                continue

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
        """
        This function sort the files contained in the file
        manger in descending/ascending order.
        :param column: the column to sort.
        :param descending: the way to sort the column.
        """

        # Take the files from the file manager
        files = [(self.file_list.set(child, column), child) for child in self.file_list.get_children('')]

        # Sort the files
        files.sort(reverse=descending)
        for index, file in enumerate(files):
            self.file_list.move(file[1], '', index)

        # reverse sort
        FileManager.SortDir = not descending

    def _select_file(self, event):
        """
        This function is called when a file in the listbox is selected.
        This function set the value of the entry to the path of the file.
        """
        self.selected.delete(0, END)
        name = '.'.join(el for el in self.file_list.item(self.file_list.selection())['values'][:2])
        path = os.path.join(settings.DATA_PATH, name)
        self.selected.insert(0, path)

        self.load_file_settings(path)

    def _select_file_with_manager(self):
        """
        This function allow the user to search in a more accurately way the file
        by using the OS manager.
        """
        path = fdialog.askopenfile(initialdir=os.getcwd(), title='Select file', filetypes=(('all files', '*.*'),))
        if path == None:
            return
        if path.name:
            for item in self.file_list.get_children():
                self.file_list.delete(item)

            settings.DATA_PATH = os.path.dirname(path.name)
            self._parse_files()
            self.selected.delete(0, END)
            self.selected.insert(INSERT, path.name)
            self.load_file_settings(path.name)

    def load_file_settings(self, path):
        """
        This function load the preview of the file and
        display it in the preview screen.

        :param path: The location of the file.
        """

        # Read the file, only if format is not dict
        if self.input_selector.get() != 'dict':
            try:
                with open(path, 'r') as file:
                    lines = file.readlines()[0:int(settings.PREVIEW_LINES)]

                    # Clean the file
                    prev = []
                    for line in lines:
                        prev.append(line.replace('{', '').replace('}', ''))

                    # Enable/disable a spacing between the lines of the preview
                    prev_spacing = False
                    if prev_spacing:
                        schr = '\n'
                    else:
                        schr = ''

                    self.preview.config(state=NORMAL)
                    self.preview.delete('1.0', END)
                    self.preview.insert('1.0', schr.join(prev))
                    self.preview.config(state=DISABLED)
            except Exception as e:
                cu.log.write_log(str(e))
                msg.showinfo(LANGUAGES[settings.LANGUAGE]['display_error'],
                             LANGUAGES[settings.LANGUAGE]['display_error_t'])
                set_to_dict = msg.askyesno(LANGUAGES[settings.LANGUAGE]['display_error'],
                                           LANGUAGES[settings.LANGUAGE]['swich_to_dict_t'])
                if set_to_dict == True:
                    self.input_selector.current(1)

    def set_next_frame(self, frame):
        """
        Set the next frame to be displayed.
        :param frame: the next frame to be displayed.
        """
        self.next_frame = frame

    def next(self, ev=None):
        """
        This function verify the data and
        show the next interface.
        """
        if self.selected.get():
            if os.path.exists(self.selected.get()):
                if self.selected.get().split('.')[-1] in settings.FILE_EXTENSIONS:
                    cu.data.CURRENT_FILE = self.selected.get()

                    # Show the next interface
                    if self.next_frame:
                        if not cu.Data.loaded:
                            #cu.data = cu.Data()
                            cu.data.CURRENT_FILE = self.selected.get()
                            cu.data.inputformat = self.input_selector.get()

                        self.main.show_frame(self.next_frame)
                    else:
                        raise ValueError('Next frame isn\'t defined')
                else:
                    msg.showerror(LANGUAGES[settings.LANGUAGE]['invalid_format'],
                                  LANGUAGES[settings.LANGUAGE]['invalid_format_t'])
            else:
                msg.showerror(LANGUAGES[settings.LANGUAGE]['file_not_exist'],
                              LANGUAGES[settings.LANGUAGE]['file_not_exist_t'])
        else:
            msg.showerror(LANGUAGES[settings.LANGUAGE]['no_file'], LANGUAGES[settings.LANGUAGE]['no_file_t'])

    def update(self):
        """
        This function update the GUI.
        """
        super().update()
        self.main_frame.canvas.yview_moveto(0)
        self.main_frame.update_view()
