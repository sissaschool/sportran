# -*- coding: utf-8 -*-
"""
==========================================
    Sportran graphic user interface
==========================================
--------------------------------------------
    Custom widgets file
--------------------------------------------

This file contains the tools and widgets that the GUI use.
"""

from tkinter import *
from tkinter import ttk
from .tk_html_widgets import HTMLLabel, HTMLScrolledText
import tkinter.filedialog as fdialog
from tkinter import messagebox as msg
import os

import matplotlib
try:
    matplotlib.use('TkAgg')
except:
    print('Error: cannot load backend TkAgg. Are you inside a graphical session? Try to change terminal!')
    raise
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import NavigationToolbar2
from matplotlib.figure import Figure
import matplotlib.patches as patches

from sportran_gui.core import control_unit as cu
from sportran_gui.core import settings
from sportran_gui.assets import ICON, METADATA, LANGUAGES, README_MD, README_GUI_MD, dev_state
import webbrowser
import markdown2


class TopBar(Frame):
    """
    This widget is the top-bar.
    The top-bar contains some cascade menus with some functions.

    TopBar is a subclass of Frame that is the root frame.

    :param parent: the main frame
    :param controller: the frame where the widget is displayed
    :param main: the main window class
    """

    show_logs = None
    show_info = None

    def __init__(self, parent, controller, main):
        Frame.__init__(self, parent)

        self.main = main

        # Setup the top menu
        top_menu = Menu(self)
        controller.configure(menu=top_menu)

        # Create the file section of the top menu
        file_menu = Menu(top_menu, tearoff=False)
        top_menu.add_cascade(label='File', menu=file_menu)

        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['new_a'],
                              command=lambda: cu.new(main))   # Starts a new analysis
        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['export_data'],
                              command=lambda: self._exportData())   # Let you export the data
        file_menu.add_separator()
        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['settings'],
                              command=lambda: run_new_window(main.root, Settings, main))   # Opens the settings window
        file_menu.add_separator()
        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['exit'],
                              command=lambda: cu.secure_exit(main))   # Close the software

        # Create the view section of the top menu
        view_menu = Menu(top_menu, tearoff=False)
        top_menu.add_cascade(label=LANGUAGES[settings.LANGUAGE]['view'], menu=view_menu)
        view_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['s_log'],
                              command=lambda: run_new_window(main.root, Logs, main))   # Opens the logs window
        view_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['s_inf'],
                              command=lambda: run_new_window(main.root, Info, main))   # Opens the outputs window

        # Create the info section of the top menu
        file_menu = Menu(top_menu, tearoff=False)
        top_menu.add_cascade(label=LANGUAGES[settings.LANGUAGE]['info'], menu=file_menu)   # Opens the info menu

        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['version'],
                              command=lambda: run_new_window(main.root, Version, main))   # Opens the version menu
        file_menu.add_separator()
        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['developers'],
                              command=lambda: run_new_window(main.root, Developers, main))   # Opens the developers menu
        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['contacts'],
                              command=lambda: run_new_window(main.root, Contacts, main))   # Opens the contacts menu
        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['about'],
                              command=lambda: run_new_window(main.root, About, main))   # Opens the about menu
        file_menu.add_separator()
        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['help'],
                              command=lambda: run_new_window(main.root, Help, main))   # Opens the help menu

    def _exportData(self):
        """
        This function is used to choose the path to export the data.
        """

        file = fdialog.asksaveasfilename(initialdir=os.getcwd())
        if file:
            if cu.export_data(file):
                msg.showinfo(LANGUAGES[settings.LANGUAGE]['export_success'],
                             LANGUAGES[settings.LANGUAGE]['export_success_t'].format(file))
            else:
                msg.showerror(LANGUAGES[settings.LANGUAGE]['export_fail'],
                              LANGUAGES[settings.LANGUAGE]['export_fail_t'].format(file))

    def _update_window(self):
        self.main.frame.update()


# TODO: add a status bar
# class StatusFrame(Frame):
#
#     def __init__(self, parent, controller):
#         Frame.__init__(self, parent)
#
#         status_frame = Frame(controller)
#         status_frame.pack(fill='x', side=BOTTOM)
#
#         self.status = Label(status_frame, text=('Status: ' + settings.STATUS_NOW))
#         self.status.pack(side=LEFT, padx=4, pady=2)


class GraphWidget(Frame):
    """
    This widget generate graphs and can be connected to a slider
    to select a value.

    GraphWidget is a subclass of Frame that is the root frame.

    :param parent: the main frame
    :param controller: the frame where the widget is displayed
    :param size: the size of the graph
    :param type_: the type of the graph
    :param toolbar: set to true to display the graph toolbar

    """

    def __init__(self, parent, controller, size=(4, 4), type_=111, toolbar=False):
        Frame.__init__(self, parent)

        # Defines some settings for the graph
        self.title = ''
        self.size = size
        self.type_ = type_
        self.mode = 'linear'

        # Defines the position of the cut line,
        # This variable is used only if the widget is attached to a slider
        self.cut_line = 0

        # Defines graph limits
        self.max_x = 1
        self.max_y = 1
        self.new_view_x = 1

        # Generate the graph frame and canvas
        self.f = Figure(figsize=self.size, dpi=100)
        self.graph = self.f.add_subplot(self.type_)

        self.canvas = FigureCanvasTkAgg(self.f, controller)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, anchor='w', padx=20, fill=BOTH)
        self.func = None

        # This variable contains other possible graphs to be drawn.
        self.other_graph = []

        # Defines the variables to attach a slider or an entry.
        self.slider = None
        self.entry = None
        self.show_selected_area = False

        self.plot_call = None
        self.plot_call_kwargs = None

        # Display the graph toolbar if selected
        if toolbar:
            toolbar = NavigationToolbar2Tk(self.canvas, controller)
            toolbar.pack(side=TOP, pady=10, padx=50, fill=BOTH)
            toolbar.update()
            self.canvas._tkcanvas.pack(side=TOP, fill=BOTH)

    def set_title(self, title):
        """
        This function allows you to set the title of the graph.
        :param title: the title to be displayed
        """

        self.f.suptitle(title)

    def get_max_x(self):
        """
        This function returns the limit of the x axis.
        :return: the limit of the x axis.
        """

        if self.graph:
            return self.graph.get_xlim()[1]
        else:
            return 1.0

    def get_max_y(self):
        """
        This function returns the limit of the y axis.
        :return: the limit of the y axis.
        """

        if self.graph:
            return self.graph.get_ylim()[1]
        else:
            return 1.0

    def show(self, func, slider_config=None, **kwargs):
        """
        This function sets the main graph.

        :param func: the function to draw the graph
        :param slider_config: the configs of the attacched slider
        """

        self.func = func
        cu.set_graph(self.graph, func, **kwargs)
        self.max_x = self.get_max_x()
        self.max_y = self.get_max_y()
        if self.slider:
            if self.show_selected_area:
                self.change_view()
            else:
                self.slider.config(to_=self.max_x)
            if slider_config is not None:
                self.slider.set(slider_config)
        self.update_cut()

    def add_graph(self, func, name, **kwargs):
        """
        This function adds other graphs to draw in the canvas.

        :param func: the function of the new graph
        :param name: a tag for the new graph
        """

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
        """
        This function draws the graphs and the cut line.

        """

        if self.graph:
            if self.entry:
                self.entry.delete(0, END)
                self.entry.insert(0, self.cut_line)

            self.graph.clear()
            cu.set_graph(self.graph, self.func, mode=self.mode, current=cu.data.j,
                         PSD_FILTER_W=cu.data.psd_filter_width, kappa_units=True)
            for graph in self.other_graph:
                cu.set_graph(self.graph, graph[1], **graph[2])

            rect = patches.Rectangle((0, 0), self.cut_line, self.max_y, linewidth=0, facecolor=(0.1, 0.2, 0.5, 0.3))
            self.graph.plot([self.cut_line, self.cut_line], [0, self.max_y])
            self.graph.add_patch(rect)

            self.graph.set_ylim([0, self.max_y])
            if self.show_selected_area:
                self.graph.set_xlim([0, self.new_view_x])
        self.graph.autoscale(True, axis='y')
        self.canvas.draw()

    def get_graph(self):
        """
        This function returns the canvas of the graph widget.
        """

        return self.graph

    def set_plot_call(self, f, **f_args):
        self.plot_call = f
        self.plot_call_kwargs = f_args

    def attach_slider(self, slider):
        """
        This function attach a slider to the widget
        :param slider: a slider widget
        """

        self.slider = slider
        self.slider.config(command=self._on_slider_change, to_=self.get_max_x())

    def change_view(self):
        """
        This function make a zoom-in on the selected graph region.
        """

        if self.show_selected_area:
            self.new_view_x = self.cut_line
            self.slider.config(to_=self.new_view_x)
        else:
            self.slider.config(to_=self.max_x)

    def attach_entry(self, entry):
        """
        This function attach an entry to the widget
        :param entry: an entry widget
        """

        self.entry = entry
        self.entry.bind('<Key-Return>', self._on_entry_change)
        self.entry.delete(0, END)
        self.entry.insert(0, self.cut_line)

    def _on_entry_change(self, ev):
        """
        This function updates the canvas when the entry change.
        """

        self.cut_line = float(self.entry.get())
        self.update_cut()

    def _on_slider_change(self, ev):
        """
        This function updates the canvas when the slider change
        """

        self.cut_line = self.slider.get()
        self.update_cut()


class TextWidget(Frame):
    """
    This widget is used to display text.

    TextWidget is a subclass of Frame

    :param parent: the main frame
    :param controller: the frame where the widget is displayed
    :param title: the title of the text container
    :param height: the height of the text container
    :param width: the width of the text container
    """

    def __init__(self, parent, controller, title, height, width):
        Frame.__init__(self, parent, controller)

        text_frame = LabelFrame(controller, text=title, bd=1, relief=SOLID)
        text_frame.pack(side=TOP, fill='x', padx=20)

        self.text_box = Text(text_frame, height=height, width=width, bd=0)
        self.text_box.pack(side=TOP, fill=BOTH, expand=1)
        self.text_box.config(state=DISABLED)
        self.text_box.see(END)

    def clear(self):
        """
        This function clear the the container.
        """

        self.text_box.config(state=NORMAL)
        self.text_box.delete('0.1', END)
        self.text_box.config(state=DISABLED)

    def write(self, text, text2='', *args, **kwargs):
        """
        This function display formatted text in the widget
        :param text: a string
        :param text2: a string
        """

        self.text_box.config(state=NORMAL)
        self.text_box.insert(INSERT, str(text) + str(text2) + '\n')
        self.text_box.config(state=DISABLED)
        self.text_box.see(END)


class CheckList(Frame):
    """
    This widget generates a list that contains a selector
    of combo box.

    :param parent: the main frame
    :param controller: the frame where the widget is displayed
    :param check_list: a dictionary that contains the name of the
                            checkbox and the content of the combobox.
    :param start_row: the row in the frame where the elements start to
                       be generated.
    """

    def __init__(self, parent, controller, check_list=dict(), start_row=0):
        Frame.__init__(self, parent, controller)

        self.controller = controller
        self.start_row = start_row
        self.combo_func = None
        if list:
            self.set_list(check_list)

    def set_list(self, check_list):
        """
        This function generate the check list.

        :param check_list: a dictionary that contains the name of the
                            checkbox and the content of the combobox.
        """

        self.clear_list()
        hidden_cont = 0
        for row, el in enumerate(list(check_list.keys())):
            if el[0] == '_':
                hidden_cont = hidden_cont + 1
                continue

            # Generate the element
            frame = Frame(self.controller)
            frame.grid(row=self.start_row + row - hidden_cont, column=0, sticky='we', pady=2)
            Label(frame, text=el, font='{} 12 bold'.format(settings.FONT)).grid(row=0, column=0)
            cmb = ttk.Combobox(frame, values=cu.Data.options, state='readonly', width=12)
            cmb.bind('<<ComboboxSelected>>', self.combo_func)   # cu.data.inputformat == 'dict':

            try:
                cmb.current(
                    cu.Data.options.index(
                        cu.data.jdata['_HEADERS']['description'][cu.data.jdata['_HEADERS']['keys'].index(el)]))
                cu.log.write_log(LANGUAGES[settings.LANGUAGE]['description_loaded'].format(
                    el, cu.data.jdata['_HEADERS']['description'][cu.data.jdata['_HEADERS']['keys'].index(el)]))

            except:   #try to guess
                if el.upper() == 'TEMP' or el.upper() == 'TEMPERATURE':
                    cmb.current(3)
                elif el.upper() == 'VOL_A' or el.upper() == 'VOLUME_A':
                    cmb.current(4)
                elif el.upper() == 'DT_FS':
                    cmb.current(5)
                else:
                    cmb.current(0)
            # else:
            #    cmb.current(0)
            cmb.grid(row=0, column=1, sticky='e')

    def clear_list(self):
        """
        This function remove all the elements from the list.
        """

        for el in self.controller.winfo_children():
            el.destroy()

    def get_list(self):
        """
        This function returns a formatted list of the check list.

        :return: returns a tuple of two lists in the first there are the
                  checked box and in second the content of combo.
        """

        check = []
        combo = []

        for el in self.controller.winfo_children():
            header = el.winfo_children()[0]
            cmb = el.winfo_children()[1]

            if cmb.get() is not None:   # and cmb.get() != 'None':
                check.append(header['text'])
                combo.append(cmb.get())

        # print(check)
        # print(combo)

        return check, combo

    def attach_function_on_combo(self, func):
        """
        Sets a function to the combo box.

        :param func: the function to attach
        """
        self.combo_func = func


class ScrollFrame(Frame):
    """
        This widget generate a scrollable frame.

        ScrollFrame is a subclass of Frame

        :param parent: the main frame
        :param controller: the frame where the widget is displayed
        :param height: the height of the frame, if 0 will be automatically calculated
        :param width: the width of the frame, if 0 will be automatically calculated
        :param bd: the size of the border
        """

    def __init__(self, parent, controller, width=0, height=0, bd=0):
        Frame.__init__(self, parent)

        bgcol = settings.BG_COLOR
        bgcol2 = settings.BG_COLOR

        # todo: should I use self in place of controller???
        # controller=self
        if width or height:
            self.canvas = Canvas(controller, bd=bd, relief=SOLID, highlightthickness=0, bg=bgcol, width=width,
                                 height=height)
        else:
            self.canvas = Canvas(controller, bd=bd, relief=SOLID, highlightthickness=0,
                                 bg=bgcol)   # '#00FF00')# bg=bgcol)

        self.viewPort = Frame(self.canvas, background=bgcol2)   # "#FF0000")#background=bgcol2)
        self.vsb = Scrollbar(controller, orient='vertical', command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)

        self.vsb.grid(row=0, column=1, sticky='nse')
        self.canvas.grid(row=0, column=0, sticky='nwse')
        controller.rowconfigure(0, weight=10)

        controller.columnconfigure(0, weight=10)
        controller.columnconfigure(1, weight=0)

        self.viewPort_id = self.canvas.create_window(0, 0, window=self.viewPort, tags='self.viewPort', anchor='n')
        # self.viewPort.grid(row=0, column=0, sticky='nsew')
        # THE VIEWPORT MUST NOT BE "GRIDDED" OR "PACKED"
        self.canvas.grid_rowconfigure(0, weight=1)
        self.canvas.grid_columnconfigure(0, weight=1)

        self.viewPort.bind('<Configure>', self.on_frame_configure)
        self.canvas.bind('<Configure>', self.on_canvas_configure)

        self.viewPort.bind('<MouseWheel>', self._on_mousewheel)
        self.viewPort.bind('<Button-4>', self._on_mousewheel)
        self.viewPort.bind('<Button-5>', self._on_mousewheel)

    def _on_mousewheel(self, event):
        """
        This function scrolls the frame on mouse scroll
        """

        if event.delta:
            self.canvas.yview_scroll(int(-1 * (event.delta / 120)), 'units')
        else:
            if event.num == 5:
                move = 1
            else:
                move = -1
            self.canvas.yview_scroll(move, 'units')

    def update_view(self):
        """
        This function configure and render the right portion of the frame
        """

        self.canvas.configure(scrollregion=self.canvas.bbox('all'))
        self.canvas.itemconfig(self.viewPort_id, width=self.canvas.winfo_width())

        print('viewPort height: {}\ncanvas height: {}'.format(self.viewPort.winfo_height(), self.canvas.winfo_height()))
        if self.canvas.winfo_height() > self.viewPort.winfo_height():
            self.canvas.itemconfig(self.viewPort_id, height=self.canvas.winfo_height())

    def on_canvas_configure(self, event):
        """
        Reset the scroll region to encompass the inner frame
        """

        self.canvas.itemconfig(self.viewPort_id, width=event.width)
        self.canvas.itemconfig(self.viewPort_id, height=event.height)

    def on_frame_configure(self, event):
        """
        Reset the scroll region to encompass the inner frame
        """

        # print ("on_frame_configure: {} {}".format(self.canvas.winfo_width(),self.viewPort.winfo_width()))
        self.canvas.configure(scrollregion=self.canvas.bbox('all'))
        self.canvas.itemconfig(self.viewPort_id, width=event.width)
        # print ("after on_frame_configure: {} {}".format(self.canvas.winfo_width(),self.viewPort.winfo_width()))


def run_new_window(root, window, main=None, *args, **kwargs):
    """
    This function starts a new window.

    :param root: the main window
    :param window: the kind of window to start
    :param main: the main window class
    :return: True if the window is generated otherwise False
    """

    # Checks if a window of the same type is already open
    for instance_window in main.open_windows:
        if isinstance(instance_window, window):
            return False

    new_window = Toplevel(root)

    # Starts the new window
    window(new_window, main=main, *args, **kwargs)
    return True


class Email:
    """
    This widget generates a clickable email element.

    :param master: the frame where the widget is displayed
    :param email: the email to show
    """

    def __init__(self, master, email):
        self.email = email

        self.email_link = Label(master, text=self.email, fg='blue', cursor='hand2')

        self.email_link.bind('<Button-1>', lambda e: self.callback())

    def grid(self, *args, **kwargs):
        self.email_link.grid(*args, **kwargs)

    def callback(self):
        """
        This function opens the Mail app.
        """

        webbrowser.open('mailto:{}'.format(self.email))


class Link:
    """
    This widget generates a clickable link element.

    :param master: the frame where the widget is displayed
    :param url: the url to the site
    :param text: the text to be displayed over the url, if None is used the url.
    """

    def __init__(self, master, url, text=None):
        self.url = url

        if text:
            self.link = Label(master, text=text, fg='blue', cursor='hand2')
        else:
            self.link = Label(master, text=self.url, fg='blue', cursor='hand2')

        self.link.bind('<Button-1>', lambda e: self.callback())

    def grid(self, *args, **kwargs):
        self.link.grid(*args, **kwargs)

    def callback(self):
        """
        This function opens the link in the browser.
        """

        webbrowser.open(self.url)


# Windows


class Version:
    """
    This window display some info about the version
    of the software.

    :param master: the new process
    :param main: the main window class
    """

    def __init__(self, master, main):

        # Set up the window layout
        self.master = master
        self.main = main

        self.master.resizable(False, False)

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(side=LEFT, padx=20, pady=15)

        # Display info

        Label(self.frame, text='Sportran GUI', font='Arial 12 bold').grid(row=0, column=0, sticky='we', pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        Label(self.frame, text='Sportran {}: {}'.format(LANGUAGES[settings.LANGUAGE]['version'].lower(),
                                                        METADATA['version'])).grid(row=2, column=0, sticky='w', pady=5)

        Label(
            self.frame,
            text='GUI {}: {} ({})'.format(LANGUAGES[settings.LANGUAGE]['version'].lower(), METADATA['gui_version'],
                                          dev_state)).grid(row=3, column=0, sticky='w')

        # Label(self.frame, text=LANGUAGES[settings.LANGUAGE]['last_release'].format(METADATA['release_date'])).grid(
        #     row=4, column=0, sticky='w', pady=5)

        # Display logo
        icon = PhotoImage(data=ICON)

        image = Label(self.master)
        image.image = icon
        image.pack(side=RIGHT, padx=20, pady=15)
        image.config(image=icon)

        self.quitButton = Button(self.frame, text=LANGUAGES[settings.LANGUAGE]['exit'], command=self.close_windows,
                                 width=10, bd=1, relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class Developers:
    """
    This window display the names of the developers.

    :param master: the new process
    :param main: the main window class
    """

    def __init__(self, master, main):

        # Set up the window layout
        self.master = master
        self.main = main

        self.master.geometry('350x190')
        self.master.resizable(False, False)

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        # Display developers name
        Label(self.frame, text=LANGUAGES[settings.LANGUAGE]['developers'],
              font='Arial 12 bold').grid(row=0, column=0, sticky='we', pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        Label(self.frame, text='Loris Ercole').grid(row=2, column=0, sticky='w', pady=5)

        Label(self.frame, text='Riccardo Bertossa').grid(row=3, column=0, sticky='w')

        Label(self.frame, text='Sebastiano Bisacchi').grid(row=4, column=0, sticky='w', pady=5)

        self.frame.columnconfigure(0, weight=1)
        self.quitButton = Button(self.frame, text=LANGUAGES[settings.LANGUAGE]['exit'], command=self.close_windows,
                                 width=10, bd=1, relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class Contacts:
    """
    This window display some contacts.

    :param master: the new process
    :param main: the main window class
    """

    def __init__(self, master, main):

        # Set up the window layout
        self.master = master
        self.main = main

        self.master.geometry('350x190')
        self.master.resizable(False, False)

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        # Display contacts
        Label(self.frame, text=LANGUAGES[settings.LANGUAGE]['contacts'],
              font='Arial 12 bold').grid(row=0, column=0, sticky='we', pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        Email(self.frame, email=METADATA['author_email']).grid(row=2, column=0, pady=5, sticky='w')
        # Email(self.frame, email='riccardomail@mail.com').grid(row=3, column=0, sticky='w')
        # Email(self.frame, email='sebastianobisacchi@outlook.it').grid(row=4, column=0, pady=5, sticky='w')

        self.frame.columnconfigure(0, weight=1)
        self.quitButton = Button(self.frame, text=LANGUAGES[settings.LANGUAGE]['exit'], command=self.close_windows,
                                 width=10, bd=1, relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class About:
    """
    This window display an info page about the software.

    :param master: the new process
    :param main: the main window class
    """

    def __init__(self, master, main):

        # Set up the window layout
        self.master = master
        self.main = main

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        # Display the info page
        Label(self.frame, text=LANGUAGES[settings.LANGUAGE]['about'],
              font='Arial 12 bold').grid(row=0, column=0, sticky='we', pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        html = markdown2.markdown(README_MD)

        html_view = HTMLScrolledText(self.frame, html=html)
        html_view.grid(row=2, column=0, sticky='wens')

        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(2, weight=1)
        self.quitButton = Button(self.frame, text=LANGUAGES[settings.LANGUAGE]['exit'], command=self.close_windows,
                                 width=10, bd=1, relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class Help:
    """
    This window display a help page.

    :param master: the new process
    :param main: the main window class
    """

    def __init__(self, master, main):

        # Set up the window layout
        self.master = master
        self.main = main

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        # Display the help page
        Label(self.frame, text=LANGUAGES[settings.LANGUAGE]['help'],
              font='Arial 12 bold').grid(row=0, column=0, sticky='we', pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        html = markdown2.markdown(README_GUI_MD)

        html_view = HTMLScrolledText(self.frame, html=html)
        html_view.grid(row=2, column=0, sticky='wens')

        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(2, weight=1)
        self.quitButton = Button(self.frame, text=LANGUAGES[settings.LANGUAGE]['exit'], command=self.close_windows,
                                 width=10, bd=1, relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class Info:
    """
    This window display some info about the actual analysis.

    :param master: the new process
    :param main: the main window class
    """

    def __init__(self, master, main):

        # Set up the window layout
        self.master = master
        self.main = main

        self.master.resizable(False, False)

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        # Display the info container
        Label(self.frame, text='Info', font='Arial 12 bold').grid(row=0, column=0, sticky='we', pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')
        info_frame = Frame(self.frame)
        info_frame.grid(row=3, column=0, sticky='nsew')
        self.info = TextWidget(info_frame, info_frame, 'Info', 10, 45)

        cu.info = self.info

        self.frame.columnconfigure(0, weight=1)
        self.quitButton = Button(self.frame, text=LANGUAGES[settings.LANGUAGE]['exit'], command=self.close_windows,
                                 width=10, bd=1, relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):
        cu.info = None
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class Logs:
    """
    This window display the logs of the software.

    :param master: the new process
    :param main: the main window class
    """

    def __init__(self, master, main):

        # Set up the window layout
        self.master = master
        self.main = main

        self.master.resizable(False, False)

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        # Display the logs container
        Label(self.frame, text='Logs', font='Arial 12 bold').grid(row=0, column=0, sticky='we', pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        log_frame = Frame(self.frame)
        log_frame.grid(row=3, column=0, sticky='nsew')
        self.logs = TextWidget(log_frame, log_frame, 'Logs', 20, 40)

        cu.log.set_func(self.logs.write)
        cu.log.set_method('other')
        self.frame.columnconfigure(0, weight=1)
        self.quitButton = Button(self.frame, text=LANGUAGES[settings.LANGUAGE]['exit'], command=self.close_windows,
                                 width=10, bd=1, relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):

        cu.log.set_func(None)
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class Settings:

    def __init__(self, master, main):

        # Set up the window layout
        self.master = master
        self.main = main

        self.master.resizable(False, False)

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        # Display the settings
        Label(self.frame, text='Settings', font='Arial 12 bold').grid(row=0, column=0, sticky='we', pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        tabs = ttk.Notebook(self.frame)
        tabs.grid(row=2, column=0, sticky='nsew')

        # Define the general settings tab
        general_settings = Frame(tabs)

        # Populate with general settings
        self.font_size_var = IntVar(value=settings.FONT_SIZE)
        self.preview_lines_var = IntVar(value=settings.PREVIEW_LINES)
        if settings.LANGUAGE == 'en-EN':
            self.language_var = StringVar(value='English')
        elif settings.LANGUAGE == 'it-IT':
            self.language_var = StringVar(value='Italian')
        else:
            self.language_var = StringVar(value='English')

        cbx_values = ['English', 'Italiano']

        Label(general_settings, text=LANGUAGES[settings.LANGUAGE]['fs']).grid(row=0, column=0, sticky='w')
        font_size = Spinbox(general_settings, from_=11, to=15, textvariable=self.font_size_var)
        font_size.grid(row=0, column=1, sticky='w')   # Font size settings

        Label(general_settings, text=LANGUAGES[settings.LANGUAGE]['pl']).grid(row=1, column=0, sticky='w')
        preview_line = Spinbox(general_settings, from_=1, to=100, textvariable=self.preview_lines_var)
        preview_line.grid(row=1, column=1, sticky='w')   # Preview lines settings

        Label(general_settings, text=LANGUAGES[settings.LANGUAGE]['lang']).grid(row=2, column=0, sticky='w')
        language = ttk.Combobox(general_settings, values=cbx_values, textvariable=self.language_var, state='readonly')
        language.grid(row=2, column=1, sticky='w')   # Language settings

        # Define the paths settings tab
        paths_settings = Frame(tabs)

        # Populate with paths settings
        self.data_path_var = StringVar(value=settings.DATA_PATH)
        self.logs_path_var = StringVar(value=settings.LOG_PATH)
        self.output_path_var = StringVar(value=settings.OUTPUT_PATH)

        Label(paths_settings, text='Data path').grid(row=0, column=0, sticky='w')
        self.data_dir_entry = Entry(paths_settings, textvariable=self.data_path_var)
        self.data_dir_entry.grid(row=0, column=1, sticky='w', padx=10)
        Button(paths_settings, text='...', relief=SOLID, bd=1,
               command=lambda: self.chose_path('dat')).grid(row=0, column=2, sticky='w')   # Data path

        Label(paths_settings, text='Logs path').grid(row=1, column=0, sticky='w')
        self.logs_dir_entry = Entry(paths_settings, textvariable=self.logs_path_var)
        self.logs_dir_entry.grid(row=1, column=1, sticky='w', padx=10)
        Button(paths_settings, text='...', relief=SOLID, bd=1,
               command=lambda: self.chose_path('log')).grid(row=1, column=2, sticky='w')   # Logs path

        Label(paths_settings, text='Outputs path').grid(row=2, column=0, sticky='w')
        self.out_dir_entry = Entry(paths_settings, textvariable=self.output_path_var)
        self.out_dir_entry.grid(row=2, column=1, sticky='w', padx=10)
        Button(paths_settings, text='...', relief=SOLID, bd=1,
               command=lambda: self.chose_path('out')).grid(row=2, column=2, sticky='w')   # Outputs path

        tabs.add(general_settings, text='General')
        tabs.add(paths_settings, text='Paths')
        self.frame.columnconfigure(0, weight=1)
        self.saveButton = Button(self.frame, text=LANGUAGES[settings.LANGUAGE]['save'], command=self.save, width=10,
                                 bd=1, relief=SOLID)
        self.saveButton.grid(row=5, column=1, sticky='w', pady=5)

        self.quitButton = Button(self.frame, text=LANGUAGES[settings.LANGUAGE]['exit'], command=self.close_windows,
                                 width=10, bd=1, relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def chose_path(self, var):
        """
        This function open a file manager to select a new path
        :param var: the variable to store the new path
        """

        path = fdialog.askdirectory()

        if var == 'dat':
            self.data_path_var.set(path)
        elif var == 'log':
            self.logs_path_var.set(path)
        elif var == 'out':
            self.output_path_var.set(path)

    def save(self):
        """
        This function save to the ini file the settings and load
        them in the settings.py file.
        """

        global settings

        settings.FONT_SIZE = int(self.font_size_var.get())
        settings.PREVIEW_LINES = int(self.preview_lines_var.get())
        settings.DATA_PATH = self.data_path_var.get()
        settings.LOG_PATH = self.logs_path_var.get()
        settings.OUTPUT_PATH = self.output_path_var.get()
        if self.language_var.get() == 'English':
            settings.LANGUAGE = 'en-EN'
        if self.language_var.get() == 'Italiano':
            settings.LANGUAGE = 'it-IT'

        with open('thcp.ini', 'w') as settings_file:
            settings_file.write('DP:{}\n'.format(settings.DATA_PATH))
            settings_file.write('LP:{}\n'.format(settings.LOG_PATH))
            settings_file.write('OP:{}\n'.format(settings.OUTPUT_PATH))
            settings_file.write('FS:{}\n'.format(settings.FONT_SIZE))
            settings_file.write('PL:{}\n'.format(settings.PREVIEW_LINES))
            settings_file.write('LANG:{}\n'.format(settings.LANGUAGE))

    def close_windows(self):

        cu.log.set_func(None)
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()
