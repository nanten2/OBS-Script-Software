import tkinter as tk

def tlwind_init():
    """tl_windows keeps track of the number of top-level windows"""
    global tl_windows
    tl_windows = []

def tkroot_init():
    global root
    root = tk.Tk()