import tkinter as tk
from tkinter import ttk

import Gvars
from FileEdit import Tabs, Files
from Parameters import Parameters


class InitDialog:
    """Handles the initial menu upon start up.

    Attributes:
    -----------

    Methods:
    --------
    init_new() -> None
        Create an instance of MainApplication.
    close_tl(tl_id) -> None
        Remove top-level window object id from tl_windows. Destroy root if tl_windows is empty.
    """

    def __init__(self, root):
        """Create a root window as the start up menu.

        Parameters:
        -----------
        root : <class 'tkinter.Tk'>
            tkinter root
        """
        InitFrame = tk.Frame(root)
        InitFrame.grid(row=0, column=0, sticky="nsew")

        Gvars.tlwind_init()

        newButton = ttk.Button(InitFrame, text="New", takefocus=0, command=self.init_new)
        newButton.grid(row=0, column=0, sticky="nsew")

    def init_new(self):
        """Create an instance of MainApplication.

        Returns:
        --------
        None
        """
        Gvars.tl_windows.append(tk.Toplevel(Gvars.root))
        Gvars.tl_windows[-1].protocol("WM_DELETE_WINDOW", lambda tl_id=Gvars.tl_windows[-1]: self.close_tl(tl_id))
        MainApplication(Gvars.tl_windows[-1], "Untitled.obs", first=True)
        Gvars.root.withdraw()

    def close_tl(self, tl_id):
        """Remove top-level window object id from tl_windows. Destroy root if tl_windows is empty.

        Parameters:
        -----------
        tl_id : (<class 'tkinter.Toplevel'>)
            id of the top-level window instance

        Returns:
        --------
        None
        """
        Gvars.tl_windows.remove(tl_id)
        tl_id.destroy()
        if len(Gvars.tl_windows) == 0:
            Gvars.root.destroy()


class MainApplication:
    """Handles the initialization of the GUI.

    Attributes:
    -----------

    Methods:
    --------
    """
    def __init__(self, master, title, **kwargs):
        """Create the skeleton frame widgets for population and create instances of the necessary classes.

        Parameters:
        -----------
        master : <class 'tkinter.Tk'>
            Tkinter root instance
        title : str
            .obs filename on initialization

        Kwargs:
        -------
        "first" : bool
            True - indicates very first instance since run
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

        # Frame for the right side (Canvas, ...)
        mainFrame2 = tk.Frame(master, bg="#eaeaea")
        mainFrame2.grid(row=0, column=2, sticky="nsew")
        mainFrame2.grid_rowconfigure(0, weight=1)
        mainFrame2.grid_columnconfigure(0, weight=1)

        # ttk.Notebook for the tabs
        ttk.Style(mainFrame).configure('leftalign.TNotebook', tabposition='nw')
        notebook = ttk.Notebook(mainFrame, style="leftalign.TNotebook")
        notebook.grid(row=0, column=0, sticky="nsew")

        param = Parameters()

        # Fill in the tabs
        ## This list implementation is an artifact of a previous version
        ## tabsList is not further appended to after this
        tabsList = []
        tabsList.append(0)
        tabsList.append(Tabs(master, notebook, mainFrame2, param.paramlist, "Untitled"))

        # Include keyword if first instance
        try:
            Files(master, notebook, mainFrame2, tabsList, param.paramlist,
                               from_new=kwargs["first"])
        except KeyError:
            Files(master, notebook, mainFrame2, tabsList, param.paramlist)

        # Match canvas size to its frame
        master.update()
        tabsList[-1].grph.canvas.config(width=tabsList[-1].Frames[4].winfo_width(),
                                        height=tabsList[-1].Frames[4].winfo_height())
