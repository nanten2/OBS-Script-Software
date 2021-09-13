import platform
import tkinter as tk


def gvar_init():
    """Assign/initialized necessary global variables"""
    # Tkinter root
    global root
    root = tk.Tk()

    # List of tkinter.Toplevel objects
    global tl_windows
    tl_windows = []

    # Platform name
    global curOS
    curOS = platform.system()
