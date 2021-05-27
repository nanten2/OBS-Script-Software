import tkinter as tk
from tkinter import ttk

from FileEdit import Tabs, Files
from Parameters import Parameters

class InitDialog:
    """Handles the initial menu upon start up."""

    def __init__(self, root):
        """Create a root window as the start up menu.

        Global variable:
            tl_windows : Keeps track of the number of top-level windows
        """
        self.root = root
        InitFrame = tk.Frame(self.root)
        InitFrame.grid(row=0, column=0, sticky="nsew")

        global tl_windows
        tl_windows = []

        self.newButton = ttk.Button(InitFrame, text="New", takefocus=0, command=self.init_new)
        self.newButton.grid(row=0, column=0, sticky="nsew")

    def init_new(self):
        """Create an instance of MainApplication."""
        tl_windows.append(tk.Toplevel(self.root))
        tl_windows[-1].protocol("WM_DELETE_WINDOW", lambda tl_id=tl_windows[-1]: self.close_tl(tl_id))
        MainApplication(tl_windows[-1], "Untitled.obs", first=True)
        self.root.withdraw()

    def close_tl(self, tl_id):
        """Remove top-level window objet id from tl_windows. Destroy root if tl_windows is empty.

        Args:
            tl_id : id of the top-level window object
        """
        tl_windows.remove(tl_id)
        tl_id.destroy()
        if len(tl_windows) == 0:
            self.root.destroy()

class MainApplication:
    def __init__(self, master, title, **kwargs):
        """Create the skeleton frame widgets for population and create instances of the necessary classes."""
        self.master = master
        master.title(title)
        master.grid_rowconfigure(0, weight=1)
        master.grid_columnconfigure([2], weight=1)

        mainFrame = tk.Frame(master)
        mainFrame.grid(row=0, column=0, sticky="nsew")
        ttk.Separator(master, orient="vertical").grid(row=0, column=1, sticky="nsew")
        mainFrame.grid_rowconfigure(0, weight=3)
        mainFrame.grid_columnconfigure(0, weight=3)

        self.mainFrame2 = tk.Frame(master, bg="#eaeaea")
        self.mainFrame2.grid(row=0, column=2, sticky="nsew")
        self.mainFrame2.grid_rowconfigure(0, weight=1)
        self.mainFrame2.grid_columnconfigure(0, weight=1)

        ttk.Style(mainFrame).configure('leftalign.TNotebook', tabposition='nw')
        self.notebook = ttk.Notebook(mainFrame, style="leftalign.TNotebook")
        self.notebook.grid(row=0, column=0, sticky="nsew")

        self.param = Parameters()

        # This list implementation is an artifact of a previous version.
        # self.tabsList is not further appended to after this.
        self.tabsList = []
        self.tabsList.append(0)
        self.tabsList.append(Tabs(master, self.notebook, self.mainFrame2, self.param.paramlist, "Untitled"))

        try:
            self.Files = Files(master, self.notebook, self.mainFrame2, self.tabsList, self.param.paramlist, from_new=kwargs["first"])
        except KeyError:
            self.Files = Files(master, self.notebook, self.mainFrame2, self.tabsList, self.param.paramlist)

        master.update()
        self.tabsList[-1].grph.canvas.config(width=self.tabsList[-1].Frames[4].winfo_width(), height=self.tabsList[-1].Frames[4].winfo_height())

