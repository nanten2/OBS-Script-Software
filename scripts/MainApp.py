import tkinter as tk
from tkinter import ttk

from . import Gvars
from .FileEdit import Files
from .Parameters import Parameters


class InitDialog:
    """Handles the initial menu upon start up."""

    def __init__(self, root):
        """Create a root window as the start-up menu.

        Parameters
        ----------
        root : tkinter.Tk
            tkinter root
        """
        InitFrame = tk.Frame(root)
        InitFrame.grid(row=0, column=0, sticky="nsew")

        #Gvars.tlwind_init()
        #vars.curOS_init()

        newButton = ttk.Button(InitFrame, text="New", takefocus=0, command=self.init_new)
        newButton.grid(row=0, column=0, sticky="nsew")

        print(f"obs_pS, running on Tcl/Tk {tk.TkVersion}")

    def init_new(self):
        """Create an instance of MainApplication."""
        Gvars.tl_windows.append(tk.Toplevel(Gvars.root))
        Gvars.tl_windows[-1].protocol("WM_DELETE_WINDOW", lambda tl_id=Gvars.tl_windows[-1]: self.close_tl(tl_id))
        MainApplication(Gvars.tl_windows[-1], "Untitled.obs")
        Gvars.root.withdraw()

    def close_tl(self, tl_id):
        """Remove top-level window object id from tl_windows. Destroy root if tl_windows is empty.

        Parameters
        ----------
        tl_id : tkinter.Toplevel
            id of the top-level window instance
        """
        Gvars.tl_windows.remove(tl_id)
        tl_id.destroy()
        if len(Gvars.tl_windows) == 0:
            Gvars.root.destroy()


class MainApplication:
    """Handles the initialization of the GUI."""
    def __init__(self, master, title):
        """Create the skeleton frame widgets for population and create instances of the necessary classes.

        Parameters
        ----------
        master : tkinter.Tk
            Tkinter root instance
        title : str
            .obs filename on initialization
        first : bool, optional
        """
        master.title(title)
        master.grid_rowconfigure(0, weight=1)
        master.grid_columnconfigure([2], weight=1)

        # Frame for the left side (OBS tab, FITS tab, ...)
        mainFrame = tk.Frame(master)
        mainFrame.grid(row=0, column=0, sticky="nsew")
        ttk.Separator(master, orient="vertical").grid(row=0, column=1, sticky="nsew")
        mainFrame.grid_rowconfigure(0, weight=3)
        mainFrame.grid_columnconfigure(0, weight=3)

        # Frame for the right side (Canvas, quickoptions, ...)
        mainFrame_g = tk.Frame(master)
        mainFrame_g.grid(row=0, column=2, sticky="nsew")
        mainFrame_g.grid_rowconfigure(0, weight=1)
        mainFrame_g.grid_columnconfigure(0, weight=1)

        # ttk.Notebook for the tabs
        ttk.Style(mainFrame).configure('leftalign.TNotebook', tabposition='nw')
        notebook = ttk.Notebook(mainFrame, style="leftalign.TNotebook")
        notebook.grid(row=0, column=0, sticky="nsew")

        # Initialize
        self.param = Parameters()
        self.Files = Files(master, notebook, mainFrame_g, self.param.paramlist)

        # Match canvas size to its frame
        master.update()
        #file.grph.canvas.config(width=file.Frames["quickoptions"].winfo_width(), height=file.Frames["quickoptions"].winfo_height())
