from tkinter import *
from tkinter import ttk
from thermocepstrum_gui.core import control_unit as cu
from thermocepstrum_gui.core import settings

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.backend_bases import NavigationToolbar2
from matplotlib.figure import Figure
import matplotlib.patches as patches


class TopBar(Frame):

    show_logs = None
    show_info = None

    def __init__(self, parent, controller, main):
        Frame.__init__(self, parent)

        self.main = main
        
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
        file_menu.add_command(label='Exit', command=lambda: cu.secure_exit(main))

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
        self.main.frame.update()


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

        bgcol = settings.BG_COLOR
        bgcol2 = settings.BG_COLOR

        if width or height:
            self.canvas = Canvas(controller, bd=bd, relief=SOLID, highlightthickness=0, bg=bgcol, width=width, height=height)
        else:
            self.canvas = Canvas(controller, bd=bd, relief=SOLID, highlightthickness=0, bg=bgcol)
        self.viewPort = Frame(self.canvas, background=bgcol2)
        self.vsb = Scrollbar(controller, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=self.vsb.set)

        self.vsb.grid(row=0, column=1, sticky='nse')
        self.canvas.grid(row=0, column=0, sticky='nsew')
        controller.rowconfigure(0, weight=10)

        controller.columnconfigure(0, weight=10)
        controller.columnconfigure(1, weight=0)

        self.canvas.create_window(0, 0, window=self.viewPort, tags="self.viewPort")
        self.viewPort.grid(row=0, column=0, sticky='nsew')
        self.canvas.grid_rowconfigure(0, weight=1)
        self.canvas.grid_columnconfigure(0, weight=1)

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