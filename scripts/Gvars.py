import tkinter as tk


def tlwind_init():
    """Create tl_windows, which is to be a list of top-level windows"""
    global tl_windows
    tl_windows = []


def tkroot_init():
    """Create Tkinter root"""
    global root
    root = tk.Tk()