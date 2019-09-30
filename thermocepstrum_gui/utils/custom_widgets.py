from tkinter import *
from tkinter import ttk
from tk_html_widgets import HTMLLabel, HTMLScrolledText
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

from thermocepstrum_gui.core import control_unit as cu
from thermocepstrum_gui.core import settings
from thermocepstrum_gui.assets import ICON, METADATA, LANGUAGES, README_MD
import webbrowser
import markdown2


class TopBar(Frame):

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

        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['new_a'], command=lambda: cu.new(main))
        file_menu.add_command(label='Export data', command=lambda: self._exportData())
        file_menu.add_separator()
        file_menu.add_command(label='Preferences')
        file_menu.add_separator()
        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['exit'], command=lambda: cu.secure_exit(main))

        # Create the view section of the top menu
        view_menu = Menu(top_menu, tearoff=False)
        top_menu.add_cascade(label=LANGUAGES[settings.LANGUAGE]['view'], menu=view_menu)
        view_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['s_log'],
                              command=lambda: run_new_window(main.root, Logs, main))
        view_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['s_inf'],
                              command=lambda: run_new_window(main.root, Info, main))

        # Create the info section of the top menu
        file_menu = Menu(top_menu, tearoff=False)
        top_menu.add_cascade(label=LANGUAGES[settings.LANGUAGE]['info'], menu=file_menu)

        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['version'],
                              command=lambda: run_new_window(main.root, Version, main))
        file_menu.add_separator()
        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['developers'],
                              command=lambda: run_new_window(main.root, Developers, main))
        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['contacts'],
                              command=lambda: run_new_window(main.root, Contacts, main))
        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['about'],
                              command=lambda: run_new_window(main.root, About, main))
        file_menu.add_separator()
        file_menu.add_command(label=LANGUAGES[settings.LANGUAGE]['help'],
                              command=lambda: run_new_window(main.root, Help, main))

    def _exportData(self):
        file = fdialog.asksaveasfilename(initialdir=os.getcwd())
        if file:
            if cu.export_data(file):
                msg.showinfo('Export success', 'Numpy data exported to "{}"'.format(file))
            else:
                msg.showerror('Export fail', 'Cannot save numpy data to file "{}"'.format(file))

    def _update_window(self):
        self.main.frame.update()


class StatusFrame(Frame):

    def __init__(self, parent, controller):
        Frame.__init__(self, parent)

        status_frame = Frame(controller)
        status_frame.pack(fill='x', side=BOTTOM)

        self.status = Label(status_frame, text=('Status: ' + settings.STATUS_NOW))
        self.status.pack(side=LEFT, padx=4, pady=2)


class GraphWidget(Frame):

    def __init__(self, parent, controller, size=(4, 4), type_=111, toolbar=False):
        Frame.__init__(self, parent)

        self.title = ''
        self.size = size
        self.type_ = type_

        self.cut_line = 0
        self.max_x = 1
        self.max_y = 1
        self.new_view_x = 1

        self.f = Figure(figsize=self.size, dpi=100)
        self.graph = self.f.add_subplot(self.type_)

        self.canvas = FigureCanvasTkAgg(self.f, controller)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=TOP, anchor='w', padx=10, fill=BOTH)
        self.f.subplots_adjust(left=0.03, right=0.89, top=0.95, bottom=0.2)
        self.func = None

        self.other_graph = []

        self.slider = None
        self.entry = None
        self.show_selected_area = False
        self.plot_call = None
        self.plot_call_kwargs = None

        if toolbar:
            toolbar = NavigationToolbar2Tk(self.canvas, controller)
            toolbar.pack(side=TOP, pady=10, padx=50, fill=BOTH)
            toolbar.update()
            self.canvas._tkcanvas.pack(side=TOP, fill=BOTH)

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
            cu.set_graph(self.graph, self.func, x=cu.data.j, PSD_FILTER_W=cu.data.psd_filter_width)
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
        text_frame.pack(side=TOP, fill='x', padx=20)

        self.text_box = Text(text_frame, height=height, width=width, bd=0)
        self.text_box.pack(side=TOP, fill=BOTH, expand=1)
        self.text_box.config(state=DISABLED)
        self.text_box.see(END)

    def clear(self):
        self.text_box.config(state=NORMAL)
        self.text_box.delete('0.1', END)
        self.text_box.config(state=DISABLED)

    def write(self, text, text2='', *args, **kwargs):
        self.text_box.config(state=NORMAL)
        self.text_box.insert(INSERT, str(text) + str(text2) + '\n')
        self.text_box.config(state=DISABLED)
        self.text_box.see(END)


class CheckList(Frame):

    def __init__(self, parent, controller, check_list=dict(), start_row=0):
        Frame.__init__(self, parent, controller)

        self.controller = controller
        self.start_row = start_row
        self.combo_func = None
        if list:
            self.set_list(check_list)

    def set_list(self, check_list):
        self.clear_list()
        for row, el in enumerate(list(check_list.keys())):
            frame = Frame(self.controller)
            frame.grid(row=self.start_row + row, column=0, sticky='we', pady=2)
            Label(frame, text=el, font='{} 12 bold'.format(settings.FONT)).grid(row=0, column=0)
            cmb = ttk.Combobox(frame, values=cu.Data.options, state='readonly', width=12)
            cmb.bind('<<ComboboxSelected>>', self.combo_func)

            if el.upper() == 'TEMP' or el.upper() == 'TEMPERATURE':
                cmb.current(3)
            elif el.upper() == 'VOL_A' or el.upper() == 'VOLUME_A':
                cmb.current(4)
            elif el.upper() == 'DT_FS':
                cmb.current(5)
            else:
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

            if cmb.get() is not None and cmb.get() != 'None':
                check.append(header['text'])
                combo.append(cmb.get())

        return check, combo

    def attach_function_on_combo(self, func):
        self.combo_func = func


class ScrollFrame(Frame):

    def __init__(self, parent, controller, width=0, height=0, bd=0):
        Frame.__init__(self, parent)

        bgcol = settings.BG_COLOR
        bgcol2 = settings.BG_COLOR

        #todo: should I use self in place of controller???
        #controller=self
        if width or height:
            self.canvas = Canvas(controller,
                                 bd=bd,
                                 relief=SOLID,
                                 highlightthickness=0,
                                 bg=bgcol,
                                 width=width,
                                 height=height)
        else:
            self.canvas = Canvas(controller, bd=bd, relief=SOLID, highlightthickness=0, bg=bgcol)
        self.viewPort = Frame(self.canvas, background=bgcol2)
        self.vsb = Scrollbar(controller, orient='vertical', command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)

        self.vsb.grid(row=0, column=1, sticky='nse')
        self.canvas.grid(row=0, column=0, sticky='nsew')
        controller.rowconfigure(0, weight=10)

        controller.columnconfigure(0, weight=10)
        controller.columnconfigure(1, weight=0)

        self.canvas.create_window(0, 0, window=self.viewPort, tags='self.viewPort')
        self.viewPort.grid(row=0, column=0, sticky='nsew')
        self.canvas.grid_rowconfigure(0, weight=1)
        self.canvas.grid_columnconfigure(0, weight=1)

        self.viewPort.bind('<Configure>', self.on_frame_configure)

        self.viewPort.bind('<MouseWheel>', self._on_mousewheel)
        self.viewPort.bind('<Button-4>', self._on_mousewheel)
        self.viewPort.bind('<Button-5>', self._on_mousewheel)

    def _on_mousewheel(self, event):
        self.canvas.yview_scroll(-1 * int((event.delta / 120)), 'units')

    def on_frame_configure(self, event):
        """
        Reset the scroll region to encompass the inner frame
        """

        self.canvas.configure(scrollregion=self.canvas.bbox('all'))


def run_new_window(root, window, main=None, *args, **kwargs):
    new_window = Toplevel(root)
    window(new_window, main=main, *args, **kwargs)


class Email:

    def __init__(self, master, email):
        self.email = email

        self.email_link = Label(master, text=self.email, fg='blue', cursor='hand2')

        self.email_link.bind('<Button-1>', lambda e: self.callback())

    def grid(self, *args, **kwargs):
        self.email_link.grid(*args, **kwargs)

    def callback(self):
        webbrowser.open('mailto:{}'.format(self.email))


class Link:

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
        webbrowser.open(self.url)


class Version:

    def __init__(self, master, main):
        self.master = master
        self.main = main

        self.master.resizable(False, False)

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(side=LEFT, padx=20, pady=15)

        Label(self.frame, text='Thermocepstrum GUI', font='Arial 12 bold').grid(row=0, column=0, sticky='we', pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        Label(self.frame,
              text='Thermocepstrum {}: {}'.format(LANGUAGES[settings.LANGUAGE]['version'].lower(),
                                                  METADATA['version'])).grid(row=2, column=0, sticky='w', pady=5)

        Label(self.frame,
              text='GUI {}: {} ({})'.format(LANGUAGES[settings.LANGUAGE]['version'].lower(), METADATA['gui_version'],
                                            METADATA['dev_state'])).grid(row=3, column=0, sticky='w')

        Label(self.frame, text='Last release: {}'.format(METADATA['release_date'])).grid(row=4,
                                                                                         column=0,
                                                                                         sticky='w',
                                                                                         pady=5)

        icon = PhotoImage(data=ICON)

        image = Label(self.master)
        image.image = icon
        image.pack(side=RIGHT, padx=20, pady=15)
        image.config(image=icon)

        self.quitButton = Button(self.frame,
                                 text=LANGUAGES[settings.LANGUAGE]['exit'],
                                 command=self.close_windows,
                                 width=10,
                                 bd=1,
                                 relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class Developers:

    def __init__(self, master, main):
        self.master = master
        self.main = main

        self.master.geometry('350x190')
        self.master.resizable(False, False)

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        Label(self.frame, text=LANGUAGES[settings.LANGUAGE]['developers'], font='Arial 12 bold').grid(row=0,
                                                                                                      column=0,
                                                                                                      sticky='we',
                                                                                                      pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        Label(self.frame, text='Loris Ercole').grid(row=2, column=0, sticky='w', pady=5)

        Label(self.frame, text='Riccardo Bertossa').grid(row=3, column=0, sticky='w')

        Label(self.frame, text='Sebastiano Bisacchi').grid(row=4, column=0, sticky='w', pady=5)

        self.frame.columnconfigure(0, weight=1)
        self.quitButton = Button(self.frame,
                                 text=LANGUAGES[settings.LANGUAGE]['exit'],
                                 command=self.close_windows,
                                 width=10,
                                 bd=1,
                                 relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class Contacts:

    def __init__(self, master, main):
        self.master = master
        self.main = main

        self.master.geometry('350x190')
        self.master.resizable(False, False)

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        Label(self.frame, text=LANGUAGES[settings.LANGUAGE]['contacts'], font='Arial 12 bold').grid(row=0,
                                                                                                    column=0,
                                                                                                    sticky='we',
                                                                                                    pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        Email(self.frame, email=METADATA['author_email']).grid(row=2, column=0, pady=5, sticky='w')
        # Email(self.frame, email='riccardomail@mail.com').grid(row=3, column=0, sticky='w')
        # Email(self.frame, email='sebastianobisacchi@outlook.it').grid(row=4, column=0, pady=5, sticky='w')

        self.frame.columnconfigure(0, weight=1)
        self.quitButton = Button(self.frame,
                                 text=LANGUAGES[settings.LANGUAGE]['exit'],
                                 command=self.close_windows,
                                 width=10,
                                 bd=1,
                                 relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class About:

    def __init__(self, master, main):
        self.master = master
        self.main = main

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        Label(self.frame, text=LANGUAGES[settings.LANGUAGE]['about'], font='Arial 12 bold').grid(row=0,
                                                                                                 column=0,
                                                                                                 sticky='we',
                                                                                                 pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        html = markdown2.markdown(README_MD)

        html_view = HTMLScrolledText(self.frame, html=html)
        html_view.grid(row=2, column=0, sticky='wens')

        self.frame.columnconfigure(0, weight=1)
        self.frame.rowconfigure(2, weight=1)
        self.quitButton = Button(self.frame,
                                 text=LANGUAGES[settings.LANGUAGE]['exit'],
                                 command=self.close_windows,
                                 width=10,
                                 bd=1,
                                 relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class Help:

    def __init__(self, master, main):
        self.master = master
        self.main = main

        self.master.geometry('350x190')
        self.master.resizable(False, False)

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        Label(self.frame, text=LANGUAGES[settings.LANGUAGE]['help'], font='Arial 12 bold').grid(row=0,
                                                                                                column=0,
                                                                                                sticky='we',
                                                                                                pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        # todo: put description
        Label(self.frame, text='').grid(row=2, column=0, pady=5, sticky='w')

        self.frame.columnconfigure(0, weight=1)
        self.quitButton = Button(self.frame,
                                 text=LANGUAGES[settings.LANGUAGE]['exit'],
                                 command=self.close_windows,
                                 width=10,
                                 bd=1,
                                 relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class Info:

    def __init__(self, master, main):
        self.master = master
        self.main = main

        self.master.resizable(False, False)

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        Label(self.frame, text='Info', font='Arial 12 bold').grid(row=0, column=0, sticky='we', pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')
        info_frame = Frame(self.frame)
        info_frame.grid(row=3, column=0, sticky='nsew')
        self.info = TextWidget(info_frame, info_frame, 'Info', 10, 45)

        cu.info = self.info

        self.frame.columnconfigure(0, weight=1)
        self.quitButton = Button(self.frame,
                                 text=LANGUAGES[settings.LANGUAGE]['exit'],
                                 command=self.close_windows,
                                 width=10,
                                 bd=1,
                                 relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):
        cu.info = None
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()


class Logs:

    def __init__(self, master, main):
        self.master = master
        self.main = main

        self.master.resizable(False, False)

        self.main.open_windows.insert(0, self)
        self.master.protocol('WM_DELETE_WINDOW', func=lambda: self.close_windows())

        self.frame = Frame(self.master)
        self.frame.pack(fill=BOTH, expand=1, padx=20, pady=15)

        Label(self.frame, text='Logs', font='Arial 12 bold').grid(row=0, column=0, sticky='we', pady=5)
        ttk.Separator(self.frame, orient=HORIZONTAL).grid(row=1, column=0, sticky='we')

        log_frame = Frame(self.frame)
        log_frame.grid(row=3, column=0, sticky='nsew')
        self.logs = TextWidget(log_frame, log_frame, 'Logs', 20, 40)

        cu.log.set_func(self.logs.write)
        cu.log.set_method('other')
        self.frame.columnconfigure(0, weight=1)
        self.quitButton = Button(self.frame,
                                 text=LANGUAGES[settings.LANGUAGE]['exit'],
                                 command=self.close_windows,
                                 width=10,
                                 bd=1,
                                 relief=SOLID)
        self.quitButton.grid(row=5, column=0, sticky='w', pady=5)

    def close_windows(self):

        cu.log.set_func(None)
        del self.main.open_windows[self.main.open_windows.index(self)]
        self.master.destroy()
