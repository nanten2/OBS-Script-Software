import os, math, cmath
from itertools import chain
import numpy as np

import tkinter as tk
from tkinter import ttk, filedialog, simpledialog
from PIL import Image

import pyregion
from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord, Angle, FK4, Galactic, FK5
from astropy.wcs import WCS
from astropy.wcs._wcs import InvalidTransformError

from . import Gvars
from . import MainApp
from .Graphic import Graphic


class Files:
    """Create the application dropdown bar and link the Tabs objects to the Graphic objects.

    Attributes
    ----------
    *_index : tuple
        Attributes of this naming convention are 2-tuples of their necessary indices in paramlist
    Frames : dict
        Dictionary of frames to be populated in the OBS tab, FITS tab, and graphical component
    entry_list : list
        List of all entry boxes in the OBS tab
    fits_CurSize : tuple
    fits_OriSize : tuple
    fits_WCS : astropy.wcs.WCS
    fits_opened : bool
    fits_sys : any
        Coordinate system of the FITS file
    frame_angle : float
        Degree offset of the FITS image to the east of north of the rectangular pixel axes
    grph : Graphic
    mainFrame_g : tkinter.Frame
        Frame for the graphical components
    master : tkinter.Toplevel
    notebook : tkinter.ttk.Notebook
    obs_sys : list
        List of size 2 containing the previous and current COORDSYS entries
    paramlist : list
    relative : list
        List of size 2 containing the previous and current RELATIVE entries
    scale_var : list
        List of variables controlled by the sliders
    scan_direction : list
        List of size 2 containing the previous and current SCAN_DIRECTION entries
    tracers : list
    """

    def __init__(
        self, master: tk.Toplevel, notebook: ttk.Notebook, mainFrame_g: tk.Frame, paramlist: list
    ) -> None:
        """Create the necessary tk.Frame objects, populate the OBS tab, and initialize the graphical component.

        Parameters
        ----------
        master : tkinter.Toplevel
        notebook : tkinter.ttk.Notebook
        mainFrame_g : tkinter.Frame
        paramlist : list
        """
        self.master = master
        self.notebook = notebook
        self.mainFrame_g = mainFrame_g
        self.paramlist = paramlist

        # Set up frames
        self.frames_setup()
        self.fill_obstab(self.Frames["obstab"])
        self.grph = Graphic(self.master, self.Frames["canvas"], self.Frames["quickoptions"])

        # Create drop-down menu
        file_menu = tk.Menu(master)
        master.config(menu=file_menu)
        file_menu_dd = tk.Menu(file_menu)
        file_menu.add_cascade(label="File", menu=file_menu_dd)
        file_menu_dd.add_command(label="New...", command=self.new)
        file_menu_dd.add_command(label="Open...", command=self.open)
        file_menu_dd.add_command(label="Open FITS...", command=self.openFITS)
        file_menu_dd.add_command(label="Load ds9 Region File...", command=self.openREG)
        file_menu_dd.add_command(label="Save As...", command=self.save)

        self.fits_opened = False
        self.setup_cCu()

    def frames_setup(self) -> None:
        """Create the necessary tk.Frame widgets."""
        self.mainFrame_g.grid_rowconfigure([0], weight=1)
        self.mainFrame_g.grid_columnconfigure([0], weight=1)

        # Double up frames in the OBS tab to enable better padding for the notebook widget
        obs_tab_base = tk.Frame(self.notebook, padx=9, pady=0)
        obs_tab = tk.Frame(obs_tab_base)
        obs_tab.grid(row=0, column=0, pady=(0, 9))

        # FITS tab frame
        fits_tab = tk.Frame(self.notebook, padx=0, pady=0)
        fits_tab.grid_rowconfigure(0, weight=1)
        fits_tab.grid_columnconfigure(0, weight=1)

        # Create the notebook tabs
        self.notebook.add(obs_tab_base, text="OBS")
        self.notebook.add(fits_tab, text="FITS", state="hidden")

        # Create relevant frames for the OBS tab, FITS tab
        self.Frames = {}
        self.Frames["obstab"] = []
        for i in range(len(self.paramlist)):
            self.Frames["obstab"].append(tk.LabelFrame(
                obs_tab,
                labelwidget=tk.Label(obs_tab, text=self.paramlist[i][0], font="Calibri 9"),
                padx=3,
                pady=3))
            self.Frames["obstab"][-1].grid(row=i, column=0, columnspan=3, sticky="nsew")
        self.Frames["fitstab"] = tk.Frame(fits_tab, padx=3, pady=3)
        self.Frames["fitstab"].grid(row=0, column=0, sticky="nsew")

        # Create relevant frames for the canvas, quickotions
        self.Frames["canvas"] = tk.Frame(self.mainFrame_g, padx=3, pady=3)
        self.Frames["canvas"].grid(row=0, rowspan=1, column=0, sticky="nsew")
        ttk.Separator(self.mainFrame_g, orient="vertical").grid(row=0, rowspan=1, column=1, sticky="nsew")
        self.Frames["quickoptions"] = tk.Frame(self.mainFrame_g, padx=3, pady=3)
        self.Frames["quickoptions"].grid(row=2, column=0, rowspan=3, columnspan=3, sticky="nsew")
        ttk.Separator(self.mainFrame_g, orient="horizontal").grid(row=1, column=0, columnspan=3, sticky="nsew")

    def fill_obstab(self, obstabframe: list[tk.LabelFrame]) -> None:
        """Populate the OBS tab with the parameters from paramlist and keep a record of the relevant variables.

        Parameters
        ----------
        obstabframe : list
            List of tkinter.LabelFrame objects
        """
        # Font for the OBS tab
        font = "Calibri 9"

        # Create entry_list to store parameters in sets of [labelframe no.][name, stringvar, widget]
        self.entry_list: list[tuple[list[str], list[tk.StringVar], list[tk.Widget]]] = []
        for i in range(len(self.paramlist)):
            self.entry_list.append([[], [], []])

        label_list = []
        for f_num, frame in enumerate(obstabframe):
            # Unfocus entry boxes by clicking elsewhere or swapping tabs
            frame.bind("<Button-1>", lambda event: self.master.focus_set())
            frame.bind("<Visibility>", lambda event: self.master.focus_set())

            # Create the necessary widgets (tk.OptionMenu, tk.Entry) and record in entry_list
            for i in range(len(self.paramlist[f_num][1])):
                self.entry_list[f_num][0].append(self.paramlist[f_num][1][i]["name"])
                self.entry_list[f_num][1].append(tk.StringVar(frame, value=''))
                # Check for dropdown menu values and use tk.OptionMenu if present, otherwise tk.Entry
                if "values" in self.paramlist[f_num][1][i]:
                    temp_values = []
                    for option in self.paramlist[f_num][1][i]["values"]:
                        # This is to deal with some issue with tkinter where the text is out of bounds to the left
                        temp_values.append("   " + option)
                    self.entry_list[f_num][2].append(tk.OptionMenu(frame, self.entry_list[f_num][1][i], *temp_values))
                    self.entry_list[f_num][2][i].config(font=font, anchor="nw", bd=3)
                    self.entry_list[f_num][2][i].bind("<Button-1>", lambda event: self.master.focus_set())
                    self.entry_list[f_num][1][i].set(temp_values[0])
                    def_offset = "   "
                else:
                    self.entry_list[f_num][2].append(tk.Entry(frame, textvariable=self.entry_list[f_num][1][i],
                                                              font=font, disabledbackground="#d3cfd9", bd=3))
                    def_offset = 0
                # Check for default value
                try:
                    self.entry_list[f_num][1][i].set(def_offset + self.paramlist[f_num][1][i]["default"])
                except KeyError:
                    pass

                # Create tk.Label widgets and bind it to defocus entries
                label_list.append(tk.Label(frame, text=self.entry_list[f_num][0][i], font=font))
                label_list[-1].bind("<Button-1>", lambda event: self.master.focus_set())

                # Check for grid keyword in paramlist, if not present arrange in order
                try:
                    label_list[-1].grid(row=self.paramlist[f_num][1][i]["grid"][0],
                                        column=(2 * self.paramlist[f_num][1][i]["grid"][1]), sticky="e")
                    self.entry_list[f_num][2][i].grid(row=self.paramlist[f_num][1][i]["grid"][0],
                                                      column=(2 * self.paramlist[f_num][1][i]["grid"][1]) + 1,
                                                      sticky="ew")
                except KeyError:
                    label_list[-1].grid(row=i, column=0, sticky="e")
                    self.entry_list[f_num][2][i].grid(row=i, column=1, sticky="ew")

    def fill_fitstab(self, fitstabframe: tk.Frame):
        """Populate the FITS tab.

        Parameters
        ----------
        fitstabframe : tkinter.Frame
            tkinter.Frame object in the FITS tab
        """
        self.notebook.tab(1, state="normal")
        fitstabframe.bind("<Button-1>", lambda event: self.master.focus_set())
        fitstabframe.bind("<Visibility>", lambda event: self.master.focus_set())
        fitstabframe.grid_propagate(True)
        fitstabframe.grid_columnconfigure([1], weight=1)
        fitstabframe.grid_columnconfigure([0,2], weight=0)
        self.scale_var = [tk.StringVar(), tk.DoubleVar(), tk.DoubleVar(), tk.DoubleVar(), tk.DoubleVar(), tk.DoubleVar(),
                          tk.DoubleVar()]

        style = ttk.Style(self.master)
        style.configure('my.Horizontal.TScale', foreground='black')

        # Create scales and set default
        type_ddbox = tk.OptionMenu(fitstabframe, self.scale_var[0], "Regular", "Symmetric")
        type_ddbox.config(anchor="w", font="Calibri 9", width=6)
        type_ddbox.grid(row=0, column=0)
        self.scale_var[0].set("Regular")
        self.scale_var[0].trace("w", self.slider_callback)

        pixelmin = round(np.min(self.grph.fitsNp_ori), 1)
        pixelmax = round(np.max(self.grph.fitsNp_ori), 1)
        self.generate_sliders_fitstab(label="Gamma", variable=self.scale_var[1], from_=0, to=10, default=1, row=1)
        self.generate_sliders_fitstab(label="Gain", variable=self.scale_var[2], from_=0, to=10, default=1, row=2)
        self.generate_sliders_fitstab(label="X-Bias", variable=self.scale_var[3], from_=0, to=1, default=0.5, row=3)
        self.generate_sliders_fitstab(label="Y-Bias", variable=self.scale_var[4], from_=0, to=1, default=0.5, row=4)
        self.generate_sliders_fitstab(label="Lower Cutoff", variable=self.scale_var[5], from_=pixelmin, to=pixelmax,
                                      default=pixelmin, resolution=0.1, row=5)
        self.generate_sliders_fitstab(label="Upper Cutoff", variable=self.scale_var[6], from_=pixelmin, to=pixelmax,
                                      default=pixelmax, resolution=0.1, row=6)
        self.master.focus_set()

    def generate_sliders_fitstab(
        self, variable, row, default, from_, to, label: str = "", resolution: float = 0.0001
    ):
        """Create a FITS tab slider with label."""
        font = "Calibri 9"
        label = tk.Label(self.Frames["fitstab"], text=label, font=font)
        label.bind("<Button-1>", lambda event: self.master.focus_set())
        label.grid(row=row, column=0, sticky="ew")
        tk.Scale(self.Frames["fitstab"], from_=from_, to=to, orient="horizontal", variable=variable,
                 showvalue=0, resolution=resolution, repeatdelay=0,
                 command=self.slider_callback).grid(row=row, column=1, sticky="ew")
        entry = tk.Entry(self.Frames["fitstab"], textvariable=variable, font=font, bd=3, width=6)
        entry.grid(row=row, column=2, sticky="e")
        variable.set(default)
        variable.trace("w", lambda t1, t2, t3, entry=entry: self.slider_trace_callback(entry))

    def slider_trace_callback(self, entry: tk.Entry) -> None:
        entry.focus_set()
        self.slider_callback()

    def slider_callback(self, *args) -> None:
        """Pass the scaling function's variable values to Graphic.slider_master.

        Parameters
        ----------
        *args: iterable
            Dummy to accept arguments passed by the tkinter.StringVar.trace() method.
        """
        vartuple = (
            self.grph.fitsNp_ori, np.min(self.grph.fitsNp_ori), np.max(self.grph.fitsNp_ori), self.scale_var[1].get(),
            self.scale_var[2].get(), self.scale_var[3].get(), self.scale_var[4].get(), self.scale_var[5].get(),
            self.scale_var[6].get(), self.scale_var[0].get())
        self.grph.slider_master(vartuple)

    def new(self, **kwargs) -> None:
        """Create a new tkinter.Toplevel object.

        Parameters:
        -----------
        **kwargs : dict
            "title" : str
                Name of the new .obs file.
        """
        newtab_title = kwargs.get("title", simpledialog.askstring('Open...', 'New OBS file'))

        if not newtab_title.rsplit(".", 1)[-1] == "obs":
            newtab_title += ".obs"

        Gvars.tl_windows.append(tk.Toplevel(Gvars.root))
        Gvars.tl_windows[-1].protocol("WM_DELETE_WINDOW", lambda tl_id=Gvars.tl_windows[-1]: self.close_tl(tl_id))
        MainApp.MainApplication(Gvars.tl_windows[-1], newtab_title)

    def open(self) -> None:
        """Open a .obs file."""
        paths_list = filedialog.askopenfilenames(parent=self.master, initialdir=os.getcwd,
                                                 title='Please select a directory', filetypes=[("OBS files", "*.obs")])
        for newtab_path in paths_list:
            newtab_title = newtab_path.rsplit("/", 1)[-1]
            Gvars.tl_windows.append(tk.Toplevel(Gvars.root))
            Gvars.tl_windows[-1].protocol("WM_DELETE_WINDOW", lambda tl_id=Gvars.tl_windows[-1]: self.close_tl(tl_id))
            new_tl = MainApp.MainApplication(Gvars.tl_windows[-1], newtab_title)

            new_tl.Files.tracers_disable()
            param = new_tl.param.paramlist
            entry_list = new_tl.Files.entry_list

            with open(newtab_path, "r") as obs_file:
                temp_entry = []
                for script_line in obs_file:
                    temp_entry.append((script_line.rsplit("#", 1)[0]).split(" = ", 1))

            for frame in range(4):
                for i in range(len(param[frame][1])):
                    j = 0
                    if param[frame][1][i]["name"] == temp_entry[j][0]:
                        if param[frame][1][i]["string"] and not temp_entry[j][1].rstrip() == "{}":
                            unit = param[frame][1][i].get("unit", "")
                            temp_entry_final = temp_entry[j][1].rstrip()[1:-1 -  len(unit)]
                        else:
                            temp_entry_final = temp_entry[j][1].rstrip()
                        if "values" in param[frame][1][i].keys():
                            temp_entry_final = "   " + temp_entry_final
                        (entry_list[frame][1][i]).set(temp_entry_final)
                    else:
                        while param[frame][1][i]["name"] != temp_entry[j][0] and j < len(temp_entry) - 1:
                            j += 1
                            if param[frame][1][i]["name"] == temp_entry[j][0]:
                                if param[frame][1][i]["string"] and not temp_entry[j][1].rstrip() == "{}":
                                    try:
                                        temp_entry_final = temp_entry[j][1].rstrip()[1:-1-len(param[frame][1][i]["unit"])]
                                    except KeyError:
                                        temp_entry_final = temp_entry[j][1].rstrip()[1:-1]
                                else:
                                    temp_entry_final = temp_entry[j][1].rstrip()
                                if "values" in param[frame][1][i].keys():
                                    temp_entry_final = "   " + temp_entry_final
                                entry_list[frame][1][i].set(temp_entry_final)
                                break

            new_tl.Files.relative_trace_callback(forceconvert=True)
            new_tl.Files.tracers_init()

    def save(self) -> None:
        """Output the current entries as an .obs file."""
        filename = filedialog.asksaveasfilename(initialfile=self.master.title(), defaultextension=".obs",
                                                filetypes=(("OBS file", "*.obs"), ("All Files", "*.*")))
        with open(filename, "w") as obs_file:
            maxlength = len(max(chain.from_iterable(self.paramlist), key=len)) + len(
                max(chain.from_iterable(self.entry_list), key=len))

            for frame in range(len(self.paramlist)):
                obs_file.writelines("[" + self.paramlist[frame][0] + "]")
                obs_file.writelines('\n')
                for position in range(len(self.paramlist[frame][1])):
                    param_index = [index for index in range(len(self.paramlist[frame][1])) if
                                self.paramlist[frame][1][index]["no."] == position + 1][0]
                    unit = self.paramlist[frame][1][param_index].get("unit", None)
                    param_name = self.paramlist[frame][1][param_index]["name"]
                    value = self.entry_list[frame][1][param_index].get().lstrip()
                    if (unit is not None) and isinstance(value, str) and value:
                        temp = f"{param_name}[{unit}] = \"{value}\""
                    elif unit is not None:
                        temp = f"{param_name}[{unit}] = {value or '{}'}"
                    else:
                        temp = f"{param_name} = {value or '{}'}"

                    obs_file.writelines(temp.ljust(maxlength + 11))
                    obs_file.writelines("# ")
                    obs_file.writelines(self.paramlist[frame][1][param_index]["#"])
                    obs_file.writelines('\n')
                obs_file.writelines('\n')

    def openFITS(self) -> None:
        """Open a FITS file, convert the image bitmap to 8-bit,
        and pass the PIL.Image object to Graphic.fits_initialize.
        """
        temp_path = filedialog.askopenfilenames(
            parent=self.master, initialdir=os.getcwd, title='Open FITS...', filetypes=[("FITS files", "*.fits")])
        if not temp_path == "":
            self.grph.FITSpath = temp_path[0]

        print("pass")

        try:
            self.fits_WCS = WCS(fits.open(self.grph.FITSpath)[0].header, naxis=2)
        except InvalidTransformError:
            print("except")
            self.fits_WCS = WCS(fits.open(self.grph.FITSpath)[0].header, translate_units='s', naxis=2)
            ##test
            # print(self.fits_WCS)
            # print(fits.open(self.grph.FITSpath)[0].header["NAXIS"])
            ##
        print("pass1")

        kw_list = ["RADESYS", "CTYPE1", "CTYPE"]
        found_keyword = False
        for kw in kw_list:
            try:
                if not found_keyword:
                    self.fits_sys = fits.open(self.grph.FITSpath)[0].header[kw].lower()
                    found_keyword = True
            except KeyError:
                pass
        if self.fits_sys in ["ra--tan"]:
            self.fits_sys = "fk5"
        # try:
        #    self.fits_sys = fits.open(self.grph.FITSpath)[0].header["RADESYS"].lower()
        # except KeyError:
        #    self.fits_sys = fits.open(self.grph.FITSpath)[0].header["CTYPE"]

        if self.fits_sys == "fk4":
            self.fits_sys = FK4(equinox="B1950")
        img_array = fits.getdata(self.grph.FITSpath, ext=0)
        for _ in range(fits.open(self.grph.FITSpath)[0].header["NAXIS"] - 2):
            [img_array] = img_array.copy()
        self.grph.fitsNp_ori = img_array

        print(img_array.dtype.name)
        # Unnecessary?
        # if False:
        #     if img_array.dtype.name == "int16":
        #         lut = np.linspace(1 / (2 ** 16), 255, 2 ** 16).astype(np.uint16)
        #         img_array = lut[img_array].astype(np.uint8)
        #     elif img_array.dtype.name == "float32":
        #         img_array = (img_array / np.max(img_array) + 0.000001) * 255

        img_array = (img_array / np.max(img_array) + 0.000001) * 255


        self.fits_OriSize = img_array.shape[::-1]
        self.grph.fits_OriSize = self.fits_OriSize
        self.grph.fits_CurSize = self.fits_OriSize

        self.grph.fitsNp = img_array
        self.grph.fits_initialize(Image.fromarray(np.flip(img_array, 0)))
        self.grph.canvas.config(
            xscrollcommand=self.grph.hbar.set,
            yscrollcommand=self.grph.vbar.set,
            scrollregion=self.grph.canvas.bbox("all"))

        self.master.update()
        (cv_w, cv_h, im_w, im_h) = (
            self.grph.canvas.winfo_width() - 1,
            self.grph.canvas.winfo_height() - 1,
            self.grph.fitsPIL.size[0], self.grph.fitsPIL.size[1])
        self.grph.canvas.xview_moveto((im_w - cv_w) / im_w / 2)
        self.grph.canvas.yview_moveto((im_h - cv_h) / im_h / 2)
        self.grph.canvas.tag_raise("O")
        self.fill_fitstab(self.Frames["fitstab"])

        self.fits_opened = True
        self.currentCoords_update(None, forangle=True)
        self.lambet_trace_callback(None)

    def openREG(self) -> None:
        """Load a pyregion file and draw the boxes on the tkinter.Canvas instance."""
        paths_list = filedialog.askopenfilenames(
            parent=self.master,
            initialdir=os.getcwd,
            title='Load ds9 Region File...',
            filetypes=[("REG files", "*.reg")])
        reg = pyregion.open(paths_list[0])
        self.fits_WCS = WCS(fits.open(self.grph.FITSpath)[0].header)

        for self.sbox in reg:
            if self.sbox.coord_format == "fk4":
                self.sbox.coord_format = FK4(equinox="B1950")
            xmid, ymid = self.fits_WCS.world_to_pixel(
                SkyCoord(self.sbox.coord_list[0], self.sbox.coord_list[1], frame=self.sbox.coord_format, unit="deg"))
            x_nw, y_nw = self.fits_WCS.world_to_pixel(
                SkyCoord(
                    self.sbox.coord_list[0] + self.sbox.coord_list[2] / 2,
                    self.sbox.coord_list[1] + self.sbox.coord_list[3] / 2,
                    frame=self.sbox.coord_format, unit="deg"))
            x_se, y_se = self.fits_WCS.world_to_pixel(
                SkyCoord(
                    self.sbox.coord_list[0] - self.sbox.coord_list[2] / 2,
                    self.sbox.coord_list[1] - self.sbox.coord_list[3] / 2,
                    frame=self.sbox.coord_format, unit="deg"))

            (x_nw, x_se, xmid) = map(lambda x: x, (x_nw, x_se, xmid))
            (y_nw, y_se, ymid) = map(lambda y: self.grph.canvas.bbox("fits")[3] - y - 1, (y_nw, y_se, ymid))

            (x1, y1), (x2, y2), (x3, y3), (x4, y4) = map(
                lambda xy: self.REG_rotate(xy[0], xy[1], xmid, ymid, self.sbox.coord_list[4]),
                zip([x_nw, x_se, x_se, x_nw], [y_nw, y_nw, y_se, y_se]))

            self.grph.setBox(None, REG=(x1, y1, x2, y2, x3, y3, x4, y4), REGdeg=self.sbox.coord_list[4])
            self.grph.canvas.tag_raise("all")

    def REG_rotate(
        self, x: float, y: float, xmid: float, ymid: float, deg: float
    ) -> tuple[float, float]:
        """Rotate the pyregion box.

        Parameters
        ----------
        x : float
            X-coordinate of the corner of the box
        y : float
            Y-coordinate of the corner of the box
        xmid : float
            X-coordinate of the center of the box
        ymid : float
            Y-coordinate of the center of the box
        deg :
            Degrees to be rotated

        Returns
        -------
        tuple
            Rotated (x,y) coordinates
        """
        reg_offset = (0, 0)
        xy = complex(x - xmid, y - ymid) * cmath.exp(complex(0, -deg * math.pi / 180))
        return xy.real + xmid + reg_offset[0], xy.imag + ymid + reg_offset[1]

    def currentCoords_update(self, _, forangle: bool = False) -> None:
        """Update the parameters as the box object of Graphic changes.

        Parameters
        ----------
        _ : any
            Used to accept event argument from bind
        forangle : bool, optional
            True if only modyfing angle
        """
        self.tracers_disable()

        # Set parameters
        self.fits_CurSize = self.grph.fits_CurSize
        scale_ratio = (self.fits_OriSize[0] / self.fits_CurSize[0], self.fits_OriSize[1] / self.fits_CurSize[1])

        coordsys = self.entry_list[self.coordsys_index[0]][1][self.coordsys_index[1]].get()
        if self.grph.box_selected or forangle:
            try:
                if self.grph.box_selected:
                    # Calculate the coordinates of the on position and start position based on the drawn box
                    boxPos = self.grph.canvas.coords(self.grph.box_id)
                    box_degnow = self.grph.Box[self.grph.box_index][-1][1]
                    box_midx = (boxPos[0] + boxPos[2] + boxPos[4] + boxPos[6]) * scale_ratio[0] / 4
                    box_midy = (self.fits_CurSize[1] - 1 - (boxPos[1] + boxPos[3] + boxPos[5] + boxPos[7]) / 4) * scale_ratio[1]
                    lambet_on = self.fits_WCS.pixel_to_world(box_midx, box_midy)
                    startposx_1 = self.fits_WCS.pixel_to_world(
                        boxPos[self.grph.Box[self.grph.box_index][-1][2]] * scale_ratio[0],
                        (self.fits_CurSize[1] - 1 - boxPos[self.grph.Box[self.grph.box_index][-1][2] + 1]) * scale_ratio[1])
                elif forangle:
                    # Get the current values of the on position and start position and
                    # calculate the pixel coordiantes of the box center
                    lam_new = str(self.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].get())
                    bet_new = str(self.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].get())
                    startx_new = float(self.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].get())
                    starty_new = float(self.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].get())
                    lambet_on = SkyCoord(lam_new, bet_new, frame=self.obs_sys[-1])
                    startposx_1 = SkyCoord(
                        Angle(lam_new).to(u.arcsec) + startx_new * u.arcsec,
                        Angle(bet_new).to(u.arcsec) + starty_new * u.arcsec, frame=self.obs_sys[-1])
                    box_midx, box_midy = self.fits_WCS.world_to_pixel(lambet_on)
                    ###
                    box_midx, box_midy = box_midx / scale_ratio[0], box_midy / scale_ratio[1]
                    ###
            except ValueError:
                # Empty or invalid entry
                return

            # Set the on position and start poisiton values
            if coordsys.lstrip() == "Galactic":
                lambet_on = lambet_on.transform_to("galactic")
                self.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].set(lambet_on.l.to_string(u.hour))
                self.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].set(lambet_on.b.to_string(u.degree))
                startposx_1 = startposx_1.transform_to("galactic")
                self.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].set(
                    round(startposx_1.l.arcsec - lambet_on.l.arcsec, 1))
                self.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].set(
                    round(startposx_1.b.arcsec - lambet_on.b.arcsec, 1))
                compr_x, compr_y = self.fits_WCS.world_to_pixel(
                    SkyCoord(lambet_on.l.degree * u.deg, (lambet_on.b.degree + 1) * u.deg, frame=Galactic()))
            elif coordsys.lstrip() == "J2000":
                lambet_on = lambet_on.transform_to(FK5(equinox="J2000"))
                self.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].set(lambet_on.ra.to_string(u.hour))
                self.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].set(
                    lambet_on.dec.to_string(u.degree))
                startposx_1 = startposx_1.transform_to(FK5(equinox="J2000"))
                self.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].set(
                    round(startposx_1.ra.arcsec - lambet_on.ra.arcsec, 1))
                self.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].set(
                    round(startposx_1.dec.arcsec - lambet_on.dec.arcsec, 1))
                compr_x, compr_y = self.fits_WCS.world_to_pixel(
                    SkyCoord(lambet_on.ra.degree * u.deg, (lambet_on.dec.degree + 1) * u.deg,
                             frame=FK5(equinox="J2000")))
            elif coordsys.lstrip() == "B1950":
                lambet_on = lambet_on.transform_to(FK4(equinox="B1950"))
                self.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].set(lambet_on.ra.to_string(u.hour))
                self.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].set(
                    lambet_on.dec.to_string(u.degree))
                startposx_1 = startposx_1.transform_to(FK4(equinox="B1950"))
                self.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].set(
                    round(startposx_1.ra.arcsec - lambet_on.ra.arcsec, 1))
                self.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].set(
                    round(startposx_1.dec.arcsec - lambet_on.dec.arcsec, 1))
                compr_x, compr_y = self.fits_WCS.world_to_pixel(
                    SkyCoord(lambet_on.ra.degree * u.deg, (lambet_on.dec.degree + 1) * u.deg,
                             frame=FK4(equinox="B1950")))

            # Determine angle of the coordinate frame
            compr_x, compr_y = compr_x - box_midx, compr_y - box_midy
            self.frame_angle = np.arctan(compr_x / compr_y)

            # Adjust angle by 180 degrees when the box is inverted
            if not forangle:
                cur_angle_equator = box_degnow
                cross_comp = 0
                if self.grph.box_resize:
                    top_midx, top_midy = (boxPos[0] + boxPos[2]) / 2, (boxPos[1] + boxPos[3]) / 2
                    bottom_midx, bottom_midy = (boxPos[4] + boxPos[6]) / 2, (boxPos[5] + boxPos[7]) / 2
                    if (top_midx > bottom_midx and 0 < cur_angle_equator < math.pi) or (
                            top_midx < bottom_midx and cur_angle_equator > math.pi) or \
                            (round(cur_angle_equator, 2) == 0 and top_midy > bottom_midy) or (
                            round(cur_angle_equator, 2) == round(math.pi, 2) and top_midy < bottom_midy):
                        cross_comp = math.pi
                self.entry_list[self.angle_index[0]][1][self.angle_index[1]].set(
                    round(
                        round((cur_angle_equator + self.frame_angle + cross_comp + self.grph.degChange) * 180 / math.pi,
                              1) % 360, 1))

            # Determine otf_len
            try:
                corner = self.grph.Box[self.grph.box_index][-1][2]
                dirct = self.grph.Box[self.grph.box_index][-1][3]
                if ((corner == 0 or corner == 4) and dirct == 2) or ((corner == 2 or corner == 6) and dirct == -2):
                    self.entry_list[self.scanDirection_index[0]][1][self.scanDirection_index[1]].set("   X")
                else:
                    self.entry_list[self.scanDirection_index[0]][1][self.scanDirection_index[1]].set("   Y")
                sd_corner = (corner + dirct) % 8
                if self.grph.box_selected:
                    start_sc, turn_sc = self.fits_WCS.pixel_to_world(boxPos[corner] * scale_ratio[0],
                                                                     boxPos[corner + 1] * scale_ratio[1]), \
                                        self.fits_WCS.pixel_to_world(boxPos[sd_corner] * scale_ratio[0],
                                                                     boxPos[sd_corner + 1] * scale_ratio[1])
                    self.entry_list[self.otf_index[0]][1][self.otf_index[1]].set(
                        round(start_sc.separation(turn_sc).arcsec, 1))
                    self.Nspacing_trace_callback()
            except AttributeError:
                pass
        else:
            try:
                if self.fits_opened:
                    self.entry_list[self.angle_index[0]][1][self.angle_index[1]].set("")
                    self.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].set("")
                    self.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].set("")
                    self.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].set("")
                    self.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].set("")
                    self.entry_list[self.scanDirection_index[0]][1][self.scanDirection_index[1]].set(" ")
                    self.entry_list[self.otf_index[0]][1][self.otf_index[1]].set("")
                else:
                    pass
            except AttributeError:
                pass
        self.tracers_init()

    def setup_cCu(self) -> None:
        """Search for necessary indices of the paramaters using setup_cCu_iterator, set their initial states,
        and bind currentCoords_update to some mouse motions.
        """
        self.angle_index = self.setup_cCu_iterator("position_angle")
        self.lambet_index = self.setup_cCu_iterator("lambda_on")
        self.startpos_index = self.setup_cCu_iterator("start_position_x")
        self.coordsys_index = self.setup_cCu_iterator("coord_sys")
        self.offabs_index = self.setup_cCu_iterator("lambda_off")
        self.offrel_index = self.setup_cCu_iterator("delta_lambda")
        self.relative_index = self.setup_cCu_iterator("relative")
        self.scanDirection_index = self.setup_cCu_iterator("scan_direction")
        self.Nspacing_index = self.setup_cCu_iterator("n")
        self.otf_index = self.setup_cCu_iterator("scan_length")

        self.obs_sys = [FK5(equinox="J2000")]
        self.relative = ["true"]
        self.scan_direction = ["X"]
        self.entry_list[self.offabs_index[0]][1][self.offabs_index[1]].set("{}")
        self.entry_list[self.offabs_index[0]][1][self.offabs_index[1] + 1].set("{}")
        self.entry_list[self.offabs_index[0]][2][self.offabs_index[1]].config(state="disabled")
        self.entry_list[self.offabs_index[0]][2][self.offabs_index[1] + 1].config(state="disabled")
        self.entry_list[self.otf_index[0]][2][self.otf_index[1]].config(state="disabled")
        self.entry_list[self.Nspacing_index[0]][2][self.Nspacing_index[1]].config(state="disabled")

        self.grph.canvas.bind("<B1-Motion>", lambda event: self.currentCoords_update(event), add="+")
        self.grph.canvas.bind("<B2-Motion>", lambda event: self.currentCoords_update(event), add="+")
        self.grph.canvas.bind("<Shift-B4-Motion>", lambda event: self.currentCoords_update(event), add="+")

        # Not currently in use
        self.grph.canvas.bind(
            "<B5-Motion>",
            lambda event, disable=False: self.lambet_trace_callback(event, disabletracers=disable),
            add="+")

        #self.coordsys_trace_callback()
        self.tracers_init()

    def setup_cCu_iterator(self, name: str) -> None:
        """Search for the indices of a parameter in Parameters.paramlist.

        Parameters
        ----------
        name : str
            Name of the parameter

        Returns
        -------
        tuple
           2-tuple of the necessary paramlist indices
        """
        return [
            (frame, position)
            for frame, framelist in enumerate(self.entry_list)
            for position in range(len(framelist[0]))
            if self.entry_list[frame][0][position] == name][0]

    def tracers_init(self) -> None:
        """Bind functions to the changes made to the variables."""
        self.tracers = [
            [self.entry_list[self.coordsys_index[0]][1][self.coordsys_index[1]], lambda t1, t2, t3 : self.coordsys_trace_callback()],
            [self.entry_list[self.lambet_index[0]][1][self.lambet_index[1]], lambda t1, t2, t3 : self.lambet_trace_callback()],
            [self.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1], lambda t1, t2, t3 : self.lambet_trace_callback()],
            [self.entry_list[self.startpos_index[0]][1][self.startpos_index[1]], lambda t1, t2, t3 : self.lambet_trace_callback()],
            [self.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1], lambda t1, t2, t3 : self.lambet_trace_callback()],
            [self.entry_list[self.angle_index[0]][1][self.angle_index[1]], lambda t1, t2, t3 : self.lambet_trace_callback()],
            [self.entry_list[self.relative_index[0]][1][self.relative_index[1]], lambda t1, t2, t3 : self.relative_trace_callback()],
            [self.entry_list[self.scanDirection_index[0]][1][self.scanDirection_index[1]], lambda t1, t2, t3 : self.scanDirection_trace_callback()],
            [self.entry_list[self.Nspacing_index[0]][1][self.Nspacing_index[1] + 1], lambda t1, t2, t3 : self.Nspacing_trace_callback()]
        ]
        for i, var in enumerate(self.tracers):
            self.tracers[i].append(var[0].trace("w", var[1]))
            #variable.trace("w", lambda t1, t2, t3, entry=entry: self.slider_trace_callback(entry))
        #self.entry_list[self.coordsys_index[0]][1][self.coordsys_index[1]].trace("w", lambda t1, t2, t3:self.coordsys_trace_callback())

    def tracers_disable(self) -> None:
        """Delete all trace binds."""
        try:
            for i, (var, callback, tracer_id) in enumerate(self.tracers):
                var.trace_vdelete("w", tracer_id)
                self.tracers[i].remove(tracer_id)
        except ValueError:
            pass

    def relative_trace_callback(self, forceconvert: bool = False, *args) -> None:
        """Convert between relative and absolute coordiantes for the off position.

        Parameters
        ----------
        forceconvert : bool
            True to go through the process regardless if RELATIVE has actually been changed
        """
        relative = str(self.entry_list[self.relative_index[0]][1][self.relative_index[1]].get().lstrip())
        self.relative.append(relative)
        self.relative = self.relative[-2:]
        if self.relative[0] == self.relative[-1] and not forceconvert:
            return

        # Mark the entries to be disabled and enabled
        if relative == "false":
            active_index, disable_index = self.offabs_index, self.offrel_index
            offset1, offset2 = 0, 1
        elif relative == "true":
            active_index, disable_index = self.offrel_index, self.offabs_index
            offset1, offset2 = 1, 0

        try:
            # Get the current values
            lam_new = str(self.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].get())
            bet_new = str(self.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].get())
            lambet_new = SkyCoord(Angle(lam_new), Angle(bet_new).to(u.hourangle), frame=self.obs_sys[-1])
            if not forceconvert:
                off_point_lam = str(self.entry_list[disable_index[0]][1][disable_index[1]].get())
                off_point_bet = str(self.entry_list[disable_index[0]][1][disable_index[1] + 1].get())
                sys_index = -1
            else:
                off_point_lam = str(self.entry_list[active_index[0]][1][active_index[1]].get())
                off_point_bet = str(self.entry_list[active_index[0]][1][active_index[1] + 1].get())
                offset2 = offset1
                sys_index = 0
            lambet_new = lambet_new.transform_to(self.obs_sys[sys_index])

            # Convert between absolute and relative
            if type(self.obs_sys[sys_index]) is type(Galactic()):
                off_point_abs = SkyCoord(
                    Angle(off_point_lam) + offset2 * Angle(lambet_new.l.to_string(u.hour)),
                    Angle(off_point_bet) + offset2 * Angle(lambet_new.b.to_string(u.degree)),
                    frame=self.obs_sys[sys_index])
            elif type(self.obs_sys[sys_index]) is type(FK5(equinox="J2000")) or type(self.obs_sys[sys_index]) is type(FK4(equinox="B1950")):
                off_point_abs = SkyCoord(
                    Angle(off_point_lam) + offset2 * Angle(lambet_new.ra.to_string(u.hour)),
                    Angle(off_point_bet) + offset2 * Angle(lambet_new.dec.to_string(u.degree)),
                    frame=self.obs_sys[sys_index])
            off_point_abs = off_point_abs.transform_to(self.obs_sys[-1])
            lambet_new = lambet_new.transform_to(self.obs_sys[-1])

            if type(self.obs_sys[-1]) is type(Galactic()):
                self.entry_list[active_index[0]][1][active_index[1]].set(
                    (off_point_abs.l - offset1 * lambet_new.l).to_string(u.hour))
                self.entry_list[active_index[0]][1][active_index[1] + 1].set(
                    (off_point_abs.b - offset1 * lambet_new.b).to_string(u.degree))
            elif type(self.obs_sys[-1]) is type(FK5(equinox="J2000")) or type(self.obs_sys[-1]) is type(FK4(equinox="B1950")):
                self.entry_list[active_index[0]][1][active_index[1]].set(
                    (off_point_abs.ra - offset1 * lambet_new.ra).to_string(u.hour))
                self.entry_list[active_index[0]][1][active_index[1] + 1].set(
                    (off_point_abs.dec - offset1 * lambet_new.dec).to_string(u.degree))

        #except (ValueError, AttributeError):
        except KeyError:
            self.entry_list[active_index[0]][1][active_index[1]].set("")
            self.entry_list[active_index[0]][1][active_index[1] + 1].set("")

        # Configure the entry boxes
        self.entry_list[disable_index[0]][1][disable_index[1]].set("{}")
        self.entry_list[disable_index[0]][1][disable_index[1] + 1].set("{}")
        self.entry_list[disable_index[0]][2][disable_index[1]].config(state="disabled")
        self.entry_list[disable_index[0]][2][disable_index[1] + 1].config(state="disabled")
        self.entry_list[active_index[0]][2][active_index[1]].config(state="normal")
        self.entry_list[active_index[0]][2][active_index[1] + 1].config(state="normal")

    def coordsys_trace_callback(self, *args) -> None:
        """Convert coordinate values to the chosen frame."""
        # Record change
        coordsys = self.entry_list[self.coordsys_index[0]][1][self.coordsys_index[1]].get()
        if coordsys.lstrip() == "Galactic":
            self.obs_sys.append(Galactic())
        elif coordsys.lstrip() == "J2000":
            self.obs_sys.append(FK5(equinox="J2000"))
        elif coordsys.lstrip() == "B1950":
            self.obs_sys.append(FK4(equinox="B1950"))
        self.obs_sys = self.obs_sys[-2:]

        # Get and convert on position and start position
        if self.obs_sys[0] is not None and not self.obs_sys[0] == self.obs_sys[-1]:
            lam_new = str(self.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].get())
            bet_new = str(self.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].get())
            startx_new = float(self.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].get() or 0)
            starty_new = float(self.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].get() or 0)

            lambet_new = SkyCoord(Angle(lam_new), Angle(bet_new).to(u.hourangle), frame=self.obs_sys[0])
            startpos_new = SkyCoord(
                Angle(lam_new) + startx_new * u.arcsec,
                Angle(bet_new).to(u.hourangle) + starty_new * u.arcsec, frame=self.obs_sys[0])

            lambet_new = lambet_new.transform_to(self.obs_sys[-1])
            startpos_new = startpos_new.transform_to(self.obs_sys[-1])

            if coordsys.lstrip() == "Galactic":
                self.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].set(lambet_new.l.to_string(u.hour))
                self.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].set(lambet_new.b.to_string(u.degree))
                self.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].set(
                    round(startpos_new.l.arcsec - lambet_new.l.arcsec, 1))
                self.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].set(
                    round(startpos_new.b.arcsec - lambet_new.b.arcsec, 1))
            elif coordsys.lstrip() == "J2000" or coordsys.lstrip() == "B1950":
                self.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].set(lambet_new.ra.to_string(u.hour))
                self.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].set(lambet_new.dec.to_string(u.degree))
                self.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].set(
                    round(startpos_new.ra.arcsec - lambet_new.ra.arcsec, 1))
                self.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].set(
                    round(startpos_new.dec.arcsec - lambet_new.dec.arcsec, 1))

            # Convert off position
            self.relative_trace_callback(forceconvert=True)

    def lambet_trace_callback(
        self, event: tk.Event | None = None, disabletracers: bool = False, *args
    ) -> None:
        """Modify box when the relevant quantities are changed.

        Parameters
        ----------
        event : Tkinter.event
        disabletracers : bool
            True to disable tracers while making the changes

        Notes
        -----
        Deleting the traces does not seem to stop the bound methods from running.
        disabletracers not currently in use.
        """
        if not self.fits_opened:
            return
        try:
            #if disabletracers:
            #    self.tracers_disable()

            # Get the current values
            lam_new = str(self.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].get())
            bet_new = str(self.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].get())
            startx_new = float(self.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].get())
            starty_new = float(self.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].get())
            angle = float(self.entry_list[self.angle_index[0]][1][self.angle_index[1]].get()) * math.pi / 180 - self.frame_angle

            # Create the SkyCoord the transform to the FITS coordinate frame
            lambet_new = SkyCoord(Angle(lam_new), Angle(bet_new).to(u.hourangle), frame=self.obs_sys[-1])
            start_new = SkyCoord(
                Angle(lam_new) + startx_new * u.arcsec,
                Angle(bet_new).to(u.hourangle) + starty_new * u.arcsec, frame=self.obs_sys[-1])
            lambet_new = lambet_new.transform_to(self.fits_sys)
            start_new = start_new.transform_to(self.fits_sys)

            # Determine the relevant box pixel coordinates
            midPos_x, midPos_y = self.fits_WCS.world_to_pixel(lambet_new)
            start_x, start_y = self.fits_WCS.world_to_pixel(start_new)
            self.fits_CurSize = self.grph.fits_CurSize
            scale_ratio = (self.fits_OriSize[0] / self.fits_CurSize[0], self.fits_OriSize[1] / self.fits_CurSize[1])
            midPos_x, midPos_y = midPos_x / scale_ratio[0], midPos_y / scale_ratio[1]
            start_x, start_y = start_x / scale_ratio[0], start_y / scale_ratio[1]
            y_dim = round((self.grph.fitsNp_ori.shape[0] - 1) / scale_ratio[1])

            if not self.grph.boxDrawn:
                cur_angle = angle
            else:
                cur_angle = self.grph.Box[self.grph.box_index][-1][1]

            oristart = complex(start_x - midPos_x, start_y - midPos_y) * cmath.exp(complex(0, -cur_angle)) #startpos if position_angle=0
            width, height = -(oristart.real * 2), (oristart.imag * 2)
            rel_start_new = oristart * cmath.exp(complex(0, angle))
            start_x_new, start_y_new = rel_start_new.real + midPos_x, rel_start_new.imag + midPos_y
            height_x, height_y = height * math.sin(angle), height * math.cos(angle)
            width_x, width_y = width * math.cos(angle), width * math.sin(angle)

            finalNW = (start_x_new, start_y_new)
            finalNE = (start_x_new + width_x, start_y_new + width_y)
            finalSE = (start_x_new + width_x + height_x, start_y_new + width_y - height_y)
            finalSW = (start_x_new + height_x, start_y_new - height_y)

            # Match the corners of the box to the quadrant it belongs
            # This maintains the box's angle orientation as the ordering of the corners would be untouched
            oristart_angle = (np.angle(oristart) - math.pi / 2) % (2 * math.pi)
            if oristart_angle < math.pi / 2:
                self.grph.startxy_c = 0
                final = (*finalNW, *finalNE, *finalSE, *finalSW)
            elif oristart_angle < math.pi:
                self.grph.startxy_c = 6
                final = (*finalSW, *finalSE, *finalNE, *finalNW)
            elif oristart_angle < 1.5 * math.pi:
                self.grph.startxy_c = 4
                final = (*finalSE, *finalSW, *finalNW, *finalNE)
            elif oristart_angle < 2 * math.pi:
                self.grph.startxy_c = 2
                final = (*finalNE, *finalNW, *finalSW, *finalSE)

            # Modify existing box or draw if none exists
            if not self.grph.boxDrawn and not self.grph.box_selected:
                self.grph.canvas.create_polygon(
                    final[0], y_dim - final[1], final[2], y_dim - final[3],
                    final[4], y_dim - final[5], final[6], y_dim - final[7],
                    fill="", width=1, outline=self.grph.boxColor, tag="tempbox")
                self.grph.setBox(None, startxy_c=self.grph.startxy_c)
                self.grph.selectBox(None, simulate=True)
                self.grph.set_onpos(None, int(self.grph.startxy_c / 2))
                self.grph.setBox(None, startxy_c=self.grph.startxy_c)
            elif not self.grph.box_manip and self.grph.box_selected and not self.grph.clicked:
                self.grph.canvas.coords(
                    self.grph.Box[self.grph.box_index][0],
                    final[0], y_dim - final[1], final[2], y_dim - final[3],
                    final[4], y_dim - final[5], final[6], y_dim - final[7])
                self.grph.set_onpos(None, int(self.grph.startxy_c / 2))
                self.grph.setBox(None)

            #if disabletracers:
            #    self.tracers_init()
        except ValueError:
            pass

    def scanDirection_trace_callback(self, *args) -> None:
        """Update scan direction record and change the  position of the small circle."""
        try:
            self.scan_direction.append(
                str(self.entry_list[self.scanDirection_index[0]][1][self.scanDirection_index[1]].get()).lstrip())
            self.scan_direction = self.scan_direction[-2:]
            self.grph.scan_direction = self.scan_direction[-1]
            self.grph.set_onpos(None, corner=int(self.grph.Box[self.grph.box_index][-1][2] / 2))
            self.grph.setBox(None)
        except AttributeError:
            pass

    def Nspacing_trace_callback(self, *args) -> None:
        """Determine the scan_spacing value.

        Parameters
        ----------
        *args : iterable
            Dummy to accept arguments passed by bind
        """
        boxPos = self.grph.canvas.coords(self.grph.box_id)
        scan_spacing = float(self.entry_list[self.Nspacing_index[0]][1][self.Nspacing_index[1] + 1].get())
        corner = self.grph.Box[self.grph.box_index][-1][2]
        sd_corner = (corner - self.grph.Box[self.grph.box_index][-1][3]) % 8
        start_sc = self.fits_WCS.pixel_to_world(boxPos[corner]/self.grph.zoom, boxPos[corner + 1]/self.grph.zoom)
        end_sc = self.fits_WCS.pixel_to_world(boxPos[sd_corner]/self.grph.zoom, boxPos[sd_corner + 1]/self.grph.zoom)
        self.entry_list[self.Nspacing_index[0]][1][self.Nspacing_index[1]].set(
            math.ceil(start_sc.separation(end_sc).arcsec / scan_spacing))

    def close_tl(self, tl_id: tk.Toplevel) -> None:
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
