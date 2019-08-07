# -*- coding: future_fstrings -*-

"""
--------------------------------------------
    Thermocepstrum graphic user interface

    Version: 0.0.2
    Release state: Beta
    Last update: 09/2019

    Main Developer:     Sebastiano Bisacchi
    Other developer:    Riccardo   Bertossa
--------------------------------------------

This file contains the GUI of the Thermocepstrum project developed at SISSA.
"""

# todo: Put an accurate description of the project?

from thermocepstrum_gui.interfaces import *
from thermocepstrum_gui.utils.custom_widgets import *

# Verify that thermocepstrum is installed
try:
    import thermocepstrum
except ImportError:
    raise ImportError('Couldn\'t find thermocepstrum')


# Main app
class ThermocepstrumGUI(Tk):
    """
    This class is used to initialize all
    the interfaces and to setup the multi frame functionality.

    ThermocepstrumGUI is a subclass of Tk that is the root window.
    """

    open_windows = []
    frames = []
    frame = None
    root = None
    container = None

    def __init__(self, version, dev_state, last_release, *args, **kwargs):
        Tk.__init__(self, *args, **kwargs)

        ThermocepstrumGUI.root = self

        # Setup the default colors and font to use in the interface
        self.option_add("*Font", "{} {}".format(settings.FONT, settings.FONT_SIZE))
        self.option_add("*Background", "{}".format(settings.BG_COLOR))
        self.option_add("*selectBackground", "light blue")
        self.option_add("*selectForeground", "black")

        # Define some info
        self.version = version
        self.dev_state = dev_state
        self.last_release = last_release

        self.show_software_info()

        # Add the main window to the open windows
        ThermocepstrumGUI.open_windows.insert(0, self)

        # Configure the window
        self.title('Thermocepstrum')
        self.geometry('{}x{}+{}+{}'.format(settings.X_SIZE,
                                           settings.Y_SIZE,
                                           settings.X_SPACING,
                                           settings.Y_SPACING))

        self.resizable(settings.X_RESIZE, settings.Y_RESIZE)

        # Define the exit function
        self.protocol('WM_DELETE_WINDOW', func=lambda: cu.secure_exit(ThermocepstrumGUI))

        # Creating the main frame
        container = Frame(self)
        ThermocepstrumGUI.container = container
        container.grid(row=0, column=0, sticky='nsew')

        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)

        container.grid_rowconfigure(0, weight=10)
        container.grid_columnconfigure(0, weight=10)

        ThermocepstrumGUI.root.grid_propagate(True)

        # ## Setting up multiple window system ## #
        self.topbar = TopBar(self, self, ThermocepstrumGUI)
        ThermocepstrumGUI.frames = {}
        frame_order = [FileManager, HeaderSelector, OtherVariables, FStarSelector, PStarSelector]

        # Load and setup the interfaces
        for n, F in enumerate(frame_order):
            ThermocepstrumGUI.frame = F(container, ThermocepstrumGUI)
            try:
                ThermocepstrumGUI.frame.set_next_frame(frame_order[n + 1])
            except:
                pass
            try:
                ThermocepstrumGUI.frame.set_prev_frame(frame_order[n - 1])
            except:
                pass

            ThermocepstrumGUI.frames[F] = ThermocepstrumGUI.frame

        # Init the main interface
        self.show_frame(FileManager)

    def show_software_info(self):
        print('------------------- Thermocepstrum GUI -------------------')
        print('')
        print('\t\t\tVersion: {}'.format(self.version))
        print('\t\t\tDev state: {}'.format(self.dev_state))
        print('\t\t\tLast release: {}'.format(self.last_release))
        print('')
        print('This software is an open-source project')
        print('developed at SISSA, Via Bonomea, 265 - 34136 Trieste ITALY')
        print('ADD OTHER INFOS')    # todo: Add other project infos
        print('----------------------------------------------------------')

    @staticmethod
    def show_frame(frame):
        """
        This function is used to display a frame.

        :param frame: the frame to be displayed. Must be a Frame object.
        """
        ThermocepstrumGUI.container.grid_rowconfigure(0, weight=0)
        ThermocepstrumGUI.container.grid_columnconfigure(0, weight=0)

        ThermocepstrumGUI.frame = ThermocepstrumGUI.frames[frame]
        ThermocepstrumGUI.frame.grid(row=0, column=0, sticky='nsew')

        ThermocepstrumGUI.container.grid_rowconfigure(0, weight=1)
        ThermocepstrumGUI.container.grid_columnconfigure(0, weight=1)

        ThermocepstrumGUI.frame.tkraise()
        ThermocepstrumGUI.frame.update()


def run():
    """
    This function is called only one time at the
    startup and it load the .ini file and start the
    software.
    """

    # Load data
    cu.load_path()

    # Start the software
    app = ThermocepstrumGUI(version='0.0.2', dev_state='beta', last_release='dd/mm/yyyy')
    app.mainloop()


if __name__ == '__main__':

    # Set the output method
    cu.log.set_method('other')
    run()
