# -*- coding: utf-8 -*-
"""
===========================================
    Sportran graphic user interface
===========================================

    Developers: Sebastiano Bisacchi, Riccardo Bertossa

This file contains the GUI of the Sportran project developed at SISSA.
"""

# todo: Put an accurate description of the project?

from sportran_gui.interfaces import *
from sportran_gui.utils.custom_widgets import *
from sportran_gui.assets import ICON, METADATA, LANGUAGES, dev_state
# Verify that sportran is installed
try:
    import sportran
except ImportError:
    raise ImportError('Couldn\'t find sportran')


# Main app
class SportranGUI(Tk):
    """
    This class is used to initialize all
    the interfaces and to setup the multi frame functionality.

    SportranGUI is a subclass of Tk that is the root window.
    """

    # Class variables to store some main parameters
    open_windows = []
    frames = []
    frame = None
    root = None
    container = None
    home = FileManager

    def __init__(self, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)

        SportranGUI.root = self

        # Setup the default colors and font to use in the interface
        self.option_add('*Font', '{} {}'.format(settings.FONT, settings.FONT_SIZE))
        self.option_add('*Background', '{}'.format(settings.BG_COLOR))
        self.option_add('*selectBackground', 'light blue')
        self.option_add('*selectForeground', 'black')

        self.show_software_info()

        # Add the main window to the open windows
        SportranGUI.open_windows.insert(0, self)

        # Configure the window

        window_icon = PhotoImage(master=self, data=ICON)
        self.iconphoto(True, window_icon)
        #self.tk.call('wm', 'iconphoto', self._w, window_icon)
        #self.iconbitmap(bitmap=window_icon)
        self.title('Sportran')
        self.geometry('{}x{}+{}+{}'.format(settings.X_SIZE, settings.Y_SIZE, settings.X_SPACING, settings.Y_SPACING))

        self.resizable(settings.X_RESIZE, settings.Y_RESIZE)

        # Define the exit function
        self.protocol('WM_DELETE_WINDOW', func=lambda: cu.secure_exit(SportranGUI))

        # Creating the main frame
        container = Frame(self)
        SportranGUI.container = container
        container.grid(row=0, column=0, sticky='nsew')

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        container.grid_rowconfigure(0, weight=10)
        container.grid_columnconfigure(0, weight=10)

        SportranGUI.root.grid_propagate(True)

        # ## Setting up multiple window system ## #
        self.topbar = TopBar(self, self, SportranGUI)
        SportranGUI.frames = {}
        frame_order = [FileManager, HeaderSelector, OtherVariables, FStarSelector, PStarSelector]

        # Load and setup the interfaces
        for n, F in enumerate(frame_order):
            SportranGUI.frame = F(container, SportranGUI)
            try:
                SportranGUI.frame.set_next_frame(frame_order[n + 1])
            except:
                pass
            try:
                SportranGUI.frame.set_prev_frame(frame_order[n - 1])
            except:
                pass

            SportranGUI.frames[F] = SportranGUI.frame

        # Init the main interface
        self.show_frame(SportranGUI.home)

    @staticmethod
    def show_software_info():
        """
        This function displays some software info at the startup.

        The data displayed are took from METADATA
        """
        print('------------------- Sportran GUI -------------------')
        print('')
        print('\t\t\tGUI version: {}'.format(METADATA['gui_version']))
        print('\t\t\tSportran version: {}'.format(METADATA['version']))
        print('\t\t\tDev state: {}'.format(dev_state))
        #print('\t\t\tLast release: {}'.format(METADATA['release_date']))
        print('\t\t\tDevelopers: {}'.format(METADATA['author']))
        print('\t\t\tURL: {}'.format(METADATA['url']))
        print('')
        print('This software is an open-source project licensed under {}'.format(METADATA['license']))
        print(METADATA['credits'])
        print('')
        print(METADATA['description'])   # todo: Add other project infos
        print('----------------------------------------------------------')

    @staticmethod
    def show_frame(frame):
        """
        This function is used to display a frame.

        :param frame: the frame to be displayed. Must be a Frame object.
        """
        SportranGUI.container.grid_rowconfigure(0, weight=0)
        SportranGUI.container.grid_columnconfigure(0, weight=0)

        SportranGUI.frame = SportranGUI.frames[frame]
        SportranGUI.frame.grid(row=0, column=0, sticky='nsew')

        SportranGUI.container.grid_rowconfigure(0, weight=1)
        SportranGUI.container.grid_columnconfigure(0, weight=1)

        SportranGUI.frame.tkraise()
        SportranGUI.frame.update_data()
        SportranGUI.frame.update()


def run():
    """
    This function is called only one time at the
    startup and it load the .ini file and start the
    software.
    """

    # Load data
    cu.load_settings()

    # Start the software
    app = SportranGUI()
    app.mainloop()


if __name__ == '__main__':

    # Set the output method
    cu.log.set_method('other')
    run()
