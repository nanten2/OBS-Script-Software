import os, math, cmath
import numpy as np
import tkinter as tk
from tkinter import ttk, filedialog, simpledialog
from itertools import chain
from PIL import Image, ImageTk
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
import pyregion
from astropy.coordinates import SkyCoord, Angle, FK4, Galactic, FK5
import ttkwidgets as ttkw

from Graphic import Graphic

class Tabs:
    def __init__(self, master, notebook, gframe_main, param, *arg):
        self.master = master
        for title in arg:
            self.tabtitle = title
            self.frames_setup(notebook, gframe_main, param)
            self.fill_obs_tab(self.Frames, param)
            self.grph = Graphic(master, self.Frames[5], self.Frames[4])

    def frames_setup(self, notebook, gframe_main, param):
        #self.obs_tab = tk.Frame(notebook, width=465, padx=9, pady=0, bg='#fafafa')
        self.obs_tab = tk.Frame(notebook, padx=9, pady=0, bg='#fafafa')
        #self.obs_tab.grid_propagate(False)
        notebook.add(self.obs_tab, text="OBS")

        self.fits_tab = tk.Frame(notebook, padx=0, pady=0)
        #self.fits_tab.grid_propagate(False)
        notebook.add(self.fits_tab, text="FITS", state="hidden")

        self.Frames = []

        for i in range(len(param)):
            self.Frames.append(tk.LabelFrame(self.obs_tab, labelwidget=tk.Label(self.obs_tab, text=param[i][0], font="Calibri 9", bg='#fafafa'), bg='#fafafa', padx=3, pady=3))
            self.Frames[-1].grid(row=i, column=0, columnspan=3, sticky="nsew")

        self.Frames.append(tk.Frame(gframe_main, bg='#fafafa', padx=3, pady=3))
        self.Frames.append(tk.Frame(gframe_main, bg='#ececec', padx=3, pady=3))
        self.Frames.append(tk.Frame(self.fits_tab, bg='#fafafa', padx=3, pady=3))

        self.fits_tab.grid_rowconfigure(0, weight=1)
        self.fits_tab.grid_columnconfigure(0, weight=1)
        gframe_main.grid_rowconfigure([0], weight=1)
        gframe_main.grid_columnconfigure([0], weight=1)

        self.Frames[4].grid(row=2, column=0, rowspan=3, columnspan=3, sticky="nsew")
        ttk.Separator(gframe_main, orient="horizontal").grid(row=1, column=0, columnspan=3, sticky="nsew")
        self.Frames[5].grid(row=0, rowspan=1, column=0, sticky="nsew")
        ttk.Separator(gframe_main, orient="vertical").grid(row=0, rowspan=1, column=1, sticky="nsew")
        self.Frames[6].grid(row=0, column=0, sticky="nsew")

    def fill_obs_tab(self, frames, param):
        font="Calibri 9"
        self.entry_list = []
        for i in range(len(param)):
            self.entry_list.append([[],[],[]])
        self.label_list = []
        for frame, frame_act in enumerate(frames):
            frame_act.bind("<Button-1>", lambda event: self.master.focus_set())
            frame_act.bind("<Visibility>", lambda event: self.master.focus_set())
            if frame < 4:
                for i in range(len(param[frame][1])):
                    self.entry_list[frame][0].append(param[frame][1][i]["name"])
                    self.entry_list[frame][1].append(tk.StringVar(frame_act, value=''))
                    try:
                        temp_values = []
                        for option in param[frame][1][i]["values"]:
                            temp_values.append("   " + option)
                        self.entry_list[frame][2].append(tk.OptionMenu(frame_act, self.entry_list[frame][1][i], *temp_values))
                        self.entry_list[frame][2][i].bind("<Button-1>", lambda event: self.master.focus_set())
                        self.entry_list[frame][2][i].config(font=font, anchor="nw", bg='#fafafa')
                        self.entry_list[frame][1][i].set(temp_values[0])
                    except KeyError:
                        self.entry_list[frame][2].append(tk.Entry(frame_act, textvariable=self.entry_list[frame][1][i], font=font, highlightbackground='#fafafa'))
                        try:
                            self.entry_list[frame][1][i].set(param[frame][1][i]["default"])
                        except KeyError:
                            pass
                    self.label_list.append(tk.Label(frame_act, text=self.entry_list[frame][0][i], font=font, bg='#fafafa'))
                    self.label_list[-1].bind("<Button-1>", lambda event: self.master.focus_set())
                    try:
                        self.label_list[-1].grid(row=param[frame][1][i]["grid"][0], column=(2*param[frame][1][i]["grid"][1]), sticky="e")
                        self.entry_list[frame][2][i].grid(row=param[frame][1][i]["grid"][0], column=(2*param[frame][1][i]["grid"][1])+1, sticky="ew")
                    except KeyError:
                        self.label_list[-1].grid(row=i, column=0, sticky="e")
                        self.entry_list[frame][2][i].grid(row=i, column=1, sticky="ew")

    def fill_fits_tab(self, notebook, Tframe):
        notebook.tab(1, state="normal")
        Tframe.grid_propagate(True)
        Tframe.grid_columnconfigure([1], weight=5)
        Tframe.grid_columnconfigure([2], weight=0)
        self.fits_gain, self.fits_gamma, self.fits_lowerb, self.fits_upperb, self.fits_biasx, self.fits_biasy = \
            tk.DoubleVar(), tk.DoubleVar(), tk.DoubleVar(), tk.DoubleVar(), tk.DoubleVar(), tk.DoubleVar()

        style = ttk.Style(self.master)
        #style.theme_use('clam')
        style.configure('my.Horizontal.TScale', background='#fafafa', foreground='black')

        self.gamma_scale = ttkw.TickScale(Tframe, from_=0, to=10, orient="horizontal", variable=self.fits_gamma, style="my.Horizontal.TScale",
                                          showvalue=False, tickinterval=0, digits=0, command=self.slider_callback)
        self.fits_gamma.set(1)
        self.gain_scale = ttkw.TickScale(Tframe, from_=0, to=10, orient="horizontal", variable=self.fits_gain, style="my.Horizontal.TScale",
                                         showvalue=False, tickinterval=0, digits=0, command=self.slider_callback)
        self.fits_gain.set(1)
        self.biasx_scale = ttkw.TickScale(Tframe, from_=0, to=1, orient="horizontal", variable=self.fits_biasx, style="my.Horizontal.TScale",
                                          showvalue=False, tickinterval=0, digits=0, command=self.slider_callback)
        self.fits_biasx.set(0.5)
        self.biasy_scale = ttkw.TickScale(Tframe, from_=0, to=1, orient="horizontal", variable=self.fits_biasy, style="my.Horizontal.TScale",
                                          showvalue=False, tickinterval=0, digits=0, command=self.slider_callback)
        self.fits_biasy.set(0.5)
        self.lowerb_scale = ttkw.TickScale(Tframe, from_=np.min(self.grph.fitsNp_ori), to=np.max(self.grph.fitsNp_ori), orient="horizontal", variable=self.fits_lowerb, style="my.Horizontal.TScale",
                                           showvalue=False, tickinterval=0, digits=0, command=self.slider_callback)
        self.fits_lowerb.set(np.min(self.grph.fitsNp_ori))
        self.upperb_scale = ttkw.TickScale(Tframe, from_=np.min(self.grph.fitsNp_ori), to=np.max(self.grph.fitsNp_ori), orient="horizontal", variable=self.fits_upperb, style="my.Horizontal.TScale",
                                           showvalue=False, tickinterval=0, digits=0, command=self.slider_callback)
        self.fits_upperb.set(np.max(self.grph.fitsNp_ori))

        self.templabel1 = tk.Label(Tframe, text="Gamma", bg="#fafafa")
        self.templabel1.grid(row=0, column=0, sticky="ew")
        self.gamma_scale.grid(row=0, column=1, sticky="ew")
        self.cmode = tk.StringVar()
        self.templbox1 = tk.OptionMenu(Tframe, self.cmode, "　Regular　", "　Symmetric")
        self.templbox1.config(anchor="nw")
        self.cmode.set("　Regular　")
        self.cmode.trace("w", self.slider_callback)
        self.templbox1.grid(row=0, column=2, sticky="ew")

        self.templabel2 = tk.Label(Tframe, text="Gain", bg="#fafafa")
        self.templabel2.grid(row=1, column=0, sticky="ew")
        self.gain_scale.grid(row=1, column=1, columnspan=2, sticky="ew")

        self.templabel3 = tk.Label(Tframe, text="X-Bias", bg="#fafafa")
        self.templabel3.grid(row=2, column=0, sticky="ew")
        self.biasx_scale.grid(row=2, column=1, columnspan=2, sticky="ew")

        self.templabel4 = tk.Label(Tframe, text="Y-Bias", bg="#fafafa")
        self.templabel4.grid(row=3, column=0, sticky="ew")
        self.biasy_scale.grid(row=3, column=1, columnspan=2, sticky="ew")

        self.templabel5 = tk.Label(Tframe, text="Lower Cutoff", bg="#fafafa")
        self.templabel5.grid(row=4, column=0, sticky="ew")
        self.lowerb_scale.grid(row=4, column=1, columnspan=2, sticky="ew")

        self.templabel6 = tk.Label(Tframe, text="Upper Cutoff", bg="#fafafa")
        self.templabel6.grid(row=5, column=0, sticky="ew")
        self.upperb_scale.grid(row=5, column=1, columnspan=2, sticky="ew")

    def slider_callback(self, *args):
        self.vartuple = (self.cmode.get(), self.grph.fitsNp_ori, self.fits_lowerb.get(), self.fits_upperb.get(), np.min(self.grph.fitsNp_ori), np.max(self.grph.fitsNp_ori),
                         self.fits_gamma.get(), self.fits_gain.get(), self.fits_biasx.get(), self.fits_biasy.get())
        self.grph.slider_master(self.vartuple)

class Files:
    def __init__(self, master, notebook, gframe_main, tabslist, param, **kwargs):
        self.master = master
        self.file_menu = tk.Menu(master)
        master.config(menu=self.file_menu)

        self.file_menu_dd = tk.Menu(self.file_menu)
        self.file_menu.add_cascade(label="File", menu=self.file_menu_dd)
        self.file_menu_dd.add_command(label="New...", command=lambda: self.new(notebook, gframe_main, tabslist, param))
        self.file_menu_dd.add_command(label="Open...", command= self.open)
        self.file_menu_dd.add_command(label="Open FITS...", command=lambda: self.openFITS(notebook, tabslist))
        self.file_menu_dd.add_command(label="Load ds9 Region File...", command=lambda: self.openREG(notebook, tabslist))
        self.file_menu_dd.add_command(label="Save As...", command=lambda: self.save(master, tabslist, param))

        try:
            self.window_from_new = kwargs["from_new"]
        except KeyError:
            self.window_from_new = False
        self.fits_opened = False
        self.setup_cCu(tabslist, notebook)


    def new(self, notebook, gframe_main, tabslist, param, **kwargs):
        try:
            newtab_title = kwargs["title"]
        except KeyError:
            newtab_title = simpledialog.askstring('Open...', 'New OBS file')

        if not newtab_title.rsplit(".", 1)[-1] == "obs":
            newtab_title += ".obs"

        tl_windows.append(tk.Toplevel(root))
        tl_windows[-1].protocol("WM_DELETE_WINDOW", lambda tl_id=tl_windows[-1]: self.close_tl(tl_id))
        MainApplication(tl_windows[-1], newtab_title, first=True)

    def open(self):
        paths_list = filedialog.askopenfilenames(parent=self.master, initialdir=os.getcwd, title='Please select a directory',
                                                 filetypes=[("OBS files", "*.obs")])
        for newtab_path in paths_list:
            newtab_title = newtab_path.rsplit("/", 1)[-1]
            tl_windows.append(tk.Toplevel(root))
            tl_windows[-1].protocol("WM_DELETE_WINDOW", lambda tl_id=tl_windows[-1]: self.close_tl(tl_id))
            self.new_tl = MainApplication(tl_windows[-1], newtab_title)

            tabslist = self.new_tl.tabsList
            param = self.new_tl.param.paramlist

            obs_file = open(newtab_path, "r")
            temp_entry = []
            for script_line in obs_file:
                temp_entry.append((script_line.rsplit("#", 1)[0]).split(" = ", 1))

            for frame in range(4):
                for i in range(len(param[frame][1])):
                    j = 0
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
                        (tabslist[-1].entry_list[frame][1][i]).set(temp_entry_final)
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
                                (tabslist[-1].entry_list[frame][1][i]).set(temp_entry_final)
                                break

            obs_file.close()

    def save(self, master, tabslist, param):
        self.entry_list = tabslist[1].entry_list
        filename = filedialog.asksaveasfilename(initialfile=master.title(), defaultextension=".obs", filetypes=(("OBS file", "*.obs"),("All Files", "*.*")))
        obs_file = open(filename, "w")
        maxlength = len(max(chain.from_iterable(param), key=len)) + len(max(chain.from_iterable(self.entry_list), key=len))

        for frame in range(len(param)):
            obs_file.writelines("["+param[frame][0]+"]")
            obs_file.writelines('\n')
            for position in range(len(param[frame][1])):
                self.param_index = [index for index in range(len(param[frame][1])) if param[frame][1][index]["no."] == position + 1][0]
                if param[frame][1][self.param_index]["string"] and not self.entry_list[frame][1][self.param_index].get().lstrip() == "{}":
                    try:
                        temp = param[frame][1][self.param_index]["name"] + " = " + '"' + self.entry_list[frame][1][self.param_index].get().lstrip() + param[frame][1][self.param_index]["unit"] + '"'
                    except KeyError:
                        temp = param[frame][1][self.param_index]["name"] + " = " + '"' + self.entry_list[frame][1][self.param_index].get().lstrip() + '"'
                else:
                    temp = param[frame][1][self.param_index]["name"] + " = " + self.entry_list[frame][1][self.param_index].get().lstrip()

                obs_file.writelines(temp.ljust(maxlength + 11))
                obs_file.writelines("# ")
                obs_file.writelines(param[frame][1][self.param_index]["#"])
                obs_file.writelines('\n')
            obs_file.writelines('\n')

        obs_file.close()

    def openFITS(self, notebook, tabslist):
        self.notebook_num = 1
        self.temp_path = filedialog.askopenfilenames(
            parent=self.master, initialdir=os.getcwd, title='Open FITS...', filetypes=[("FITS files", "*.fits")])
        if not self.temp_path == "":
            tabslist[self.notebook_num].grph.FITSpath = self.temp_path[0]
        self.fits_WSC = WCS(fits.open(tabslist[self.notebook_num].grph.FITSpath)[0].header)
        try:
            self.fits_sys = fits.open(tabslist[self.notebook_num].grph.FITSpath)[0].header["RADESYS"].lower()
        except KeyError:
            self.fits_sys = fits.open(tabslist[self.notebook_num].grph.FITSpath)[0].header["CTYPE"]
        if self.fits_sys == "fk4":
            self.fits_sys = FK4(equinox="B1950")
        self.img_array = fits.getdata(tabslist[self.notebook_num].grph.FITSpath, ext=0)
        tabslist[self.notebook_num].grph.fitsNp_ori = self.img_array

        if self.img_array.dtype.name == "int16":
            self.lut = np.linspace(1/(2**16), 255, 2**16).astype(np.uint16)
            self.img_array = self.lut[self.img_array].astype(np.uint8)
        elif self.img_array.dtype.name == "float32":
            self.img_array = (self.img_array / np.max(self.img_array) + 0.000001) * 255

        tabslist[self.notebook_num].grph.fitsNp = self.img_array
        tabslist[self.notebook_num].grph.fitsPIL = Image.fromarray(np.flip(self.img_array,0))
        tabslist[self.notebook_num].grph.fitsTk = ImageTk.PhotoImage(tabslist[self.notebook_num].grph.fitsPIL)
        tabslist[self.notebook_num].grph.fitsCanvas = tabslist[self.notebook_num].grph.canvas.create_image(
            tabslist[self.notebook_num].grph.fits_offset, image=tabslist[self.notebook_num].grph.fitsTk, anchor="nw", tag="fits")

        tabslist[self.notebook_num].grph.canvas.config(
            xscrollcommand=tabslist[self.notebook_num].grph.hbar.set,
            yscrollcommand=tabslist[self.notebook_num].grph.vbar.set,
            scrollregion=tabslist[self.notebook_num].grph.canvas.bbox("all"))

        self.master.update()
        (self.cv_w, self.cv_h, self.im_w, self.im_h) = (
            tabslist[self.notebook_num].grph.canvas.winfo_width()-1, tabslist[self.notebook_num].grph.canvas.winfo_height()-1,
            tabslist[self.notebook_num].grph.fitsPIL.size[0], tabslist[self.notebook_num].grph.fitsPIL.size[1])
        tabslist[self.notebook_num].grph.canvas.xview_moveto(-(self.cv_w-self.im_w)/self.im_w/2)
        tabslist[self.notebook_num].grph.canvas.yview_moveto(-(self.cv_h-self.im_h)/self.im_h/2)

        tabslist[self.notebook_num].grph.canvas.tag_raise("O")
        tabslist[self.notebook_num].fill_fits_tab(notebook, tabslist[self.notebook_num].Frames[6])

        self.fits_opened = True
        self.currentCoords_update(None, tabslist, notebook, forangle=True)
        self.lambet_trace_callback()

    def openREG(self, notebook, tabslist):
        self.notebook_num = notebook.index(notebook.select()) + 1
        self.paths_list = filedialog.askopenfilenames(parent=self.master, initialdir=os.getcwd, title='Load ds9 Region File...',
                                                      filetypes=[("REG files", "*.reg")])
        self.reg = pyregion.open(self.paths_list[0])
        self.fits_WSC = WCS(fits.open(tabslist[self.notebook_num].grph.FITSpath)[0].header)

        for self.sbox in self.reg:
            print(self.sbox.coord_format)
            print(self.sbox.coord_list[0], self.sbox.coord_list[1])
            print(SkyCoord(self.sbox.coord_list[0] + self.sbox.coord_list[2] / 2, self.sbox.coord_list[1] + self.sbox.coord_list[3] / 2,
                         frame=self.sbox.coord_format, unit="deg"))
            print(self.sbox.coord_list[0] - self.sbox.coord_list[2] / 2)
            if self.sbox.coord_format == "fk4":
                self.sbox.coord_format = FK4(equinox="B1950")
            xmid, ymid = self.fits_WSC.world_to_pixel(
                SkyCoord(self.sbox.coord_list[0], self.sbox.coord_list[1], frame=self.sbox.coord_format, unit="deg"))
            print("1:", xmid, ymid)
            print(self.fits_WSC.pixel_to_world(191.5,228.5))
            x_nw, y_nw = self.fits_WSC.world_to_pixel(
                SkyCoord(self.sbox.coord_list[0] + self.sbox.coord_list[2] / 2, self.sbox.coord_list[1] + self.sbox.coord_list[3] / 2,
                         frame=self.sbox.coord_format, unit="deg"))
            print("2:",x_nw, y_nw)
            x_se, y_se = self.fits_WSC.world_to_pixel(
                SkyCoord(self.sbox.coord_list[0] - self.sbox.coord_list[2] / 2, self.sbox.coord_list[1] - self.sbox.coord_list[3] / 2,
                         frame=self.sbox.coord_format, unit="deg"))
            print("3:", x_se, y_se)

            (x_nw, x_se, xmid) = map(lambda x: x, (x_nw, x_se, xmid))
            (y_nw, y_se, ymid) = map(lambda y: tabslist[self.notebook_num].grph.canvas.bbox("fits")[3] - y - 1, (y_nw, y_se, ymid))

            (x1, y1), (x2, y2), (x3, y3), (x4, y4) = map(lambda xy: self.REG_rotate(xy[0], xy[1], xmid, ymid, self.sbox.coord_list[4]),
                                                    zip([x_nw, x_se, x_se, x_nw],[y_nw, y_nw, y_se, y_se]))

            tabslist[self.notebook_num].grph.setBox(None, REG=(x1, y1, x2, y2, x3, y3, x4, y4), REGdeg=self.sbox.coord_list[4])
            tabslist[self.notebook_num].grph.canvas.tag_raise("all")

    def REG_rotate(self, x, y, xmid, ymid, deg):
        self.reg_offset = (0, 0)
        xy = complex(x-xmid, y-ymid) * cmath.exp(complex(0, -deg*math.pi/180))
        return xy.real + xmid + self.reg_offset[0], xy.imag + ymid + self.reg_offset[1]

    def offpos_callback(self):
        self.offabs_index = [(frame, position)
                            for frame, framelist in enumerate(self.current_tab.entry_list)
                            for position in range(len(framelist[0]))
                            if self.current_tab.entry_list[frame][0][position] == "lambda_off"][0]
        self.offrel_index = [(frame, position)
                            for frame, framelist in enumerate(self.current_tab.entry_list)
                            for position in range(len(framelist[0]))
                            if self.current_tab.entry_list[frame][0][position] == "lambdel_off"][0]

    def currentCoords_update(self, _, tabslist, notebook, forangle=False, SD=False):
        self.tracers_disable()

        self.current_tab = tabslist[1]
        coordsys = self.current_tab.entry_list[self.coordsys_index[0]][1][self.coordsys_index[1]].get()

        if self.current_tab.grph.box_selected or forangle:
            if self.current_tab.grph.box_selected:
                self.boxPos = self.current_tab.grph.canvas.coords(self.current_tab.grph.box_id)
                box_degnow = self.current_tab.grph.Box[self.current_tab.grph.box_index][-1][1]
                box_midx, box_midy = (self.boxPos[0]+self.boxPos[2]+self.boxPos[4]+self.boxPos[6])/4, tabslist[self.notebook_num].grph.fitsNp_ori.shape[0] - 1 - (self.boxPos[1]+self.boxPos[3]+self.boxPos[5]+self.boxPos[7])/4
                lambet_on = self.fits_WSC.pixel_to_world(box_midx, box_midy)
                startposx_1 = self.fits_WSC.pixel_to_world(self.boxPos[0], tabslist[self.notebook_num].grph.fitsNp_ori.shape[0] - 1 - self.boxPos[1])
            elif forangle:
                lam_new = str(self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].get())
                bet_new = str(self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].get())
                startx_new = float(self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].get())
                starty_new = float(self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].get())
                lambet_on = SkyCoord(lam_new, bet_new, frame=self.obs_sys[-1])
                startposx_1 = SkyCoord(Angle(lam_new).to(u.arcsec) + startx_new * u.arcsec, Angle(bet_new).to(u.arcsec) + starty_new * u.arcsec, frame=self.obs_sys[-1])
                box_midx, box_midy = self.fits_WSC.world_to_pixel(lambet_on)

            if coordsys.lstrip() == "Galactic":
                lambet_on = lambet_on.transform_to("galactic")
                self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].set(lambet_on.l.to_string(u.hour))
                self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].set(lambet_on.b.to_string(u.degree))
                startposx_1 = startposx_1.transform_to("galactic")
                self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].set(round(startposx_1.l.arcsec - lambet_on.l.arcsec,1))
                self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].set(round(startposx_1.b.arcsec - lambet_on.b.arcsec,1))
                compr_x, compr_y = self.fits_WSC.world_to_pixel(SkyCoord(lambet_on.l.degree * u.deg, (lambet_on.b.degree + 1) * u.deg, frame=Galactic()))
            elif coordsys.lstrip() == "J2000":
                lambet_on = lambet_on.transform_to(FK5(equinox="J2000"))
                self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].set(lambet_on.ra.to_string(u.hour))
                self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].set(lambet_on.dec.to_string(u.degree))
                startposx_1 = startposx_1.transform_to(FK5(equinox="J2000"))
                self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].set(round(startposx_1.ra.arcsec - lambet_on.ra.arcsec,1))
                self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].set(round(startposx_1.dec.arcsec - lambet_on.dec.arcsec,1))
                compr_x, compr_y = self.fits_WSC.world_to_pixel(SkyCoord(lambet_on.ra.degree * u.deg, (lambet_on.dec.degree + 1) * u.deg, frame=FK5(equinox="J2000")))
            elif coordsys.lstrip() == "B1950":
                lambet_on = lambet_on.transform_to(FK4(equinox="B1950"))
                self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].set(lambet_on.ra.to_string(u.hour))
                self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].set(lambet_on.dec.to_string(u.degree))
                startposx_1 = startposx_1.transform_to(FK4(equinox="B1950"))
                self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].set(round(startposx_1.ra.arcsec - lambet_on.ra.arcsec,1))
                self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].set(round(startposx_1.dec.arcsec - lambet_on.dec.arcsec,1))
                compr_x, compr_y = self.fits_WSC.world_to_pixel(SkyCoord(lambet_on.ra.degree * u.deg, (lambet_on.dec.degree + 1) * u.deg, frame=FK4(equinox="B1950")))

            compr_x, compr_y = compr_x - box_midx, compr_y - box_midy
            self.frame_angle = np.arctan(compr_x/compr_y)

            if not forangle:
                cur_angle_equator = box_degnow + self.current_tab.grph.degChange
                cross_comp = 0
                if self.current_tab.grph.box_resize:
                    top_midx, top_midy = (self.boxPos[0] + self.boxPos[2])/2, (self.boxPos[1] + self.boxPos[3])/2
                    cent_midx, cent_midy = top_midx, (self.boxPos[1] + self.boxPos[7])/2
                    if (top_midx > cent_midx and 0 < cur_angle_equator < math.pi) or (top_midx < cent_midx and cur_angle_equator > math.pi) or \
                            (round(cur_angle_equator,2) == 0 and top_midy > cent_midy) or (round(cur_angle_equator,2) == math.pi and top_midy < cent_midy):
                        cross_comp = math.pi
                self.current_tab.entry_list[self.angle_index[0]][1][self.angle_index[1]].set(round(round((cur_angle_equator + self.frame_angle + cross_comp)*180/math.pi, 1) % 360, 1))

            if SD or self.current_tab.grph.box_resize:
                start_corner = complex(self.boxPos[0] - box_midx, box_midy - self.boxPos[1]) \
                               * cmath.exp(complex(0, -(self.current_tab.grph.Box[self.current_tab.grph.box_index][-1][1] - self.frame_angle)))
                turn_corner = complex(self.boxPos[2] - box_midx, box_midy - self.boxPos[3]) \
                               * cmath.exp(complex(0, -(self.current_tab.grph.Box[self.current_tab.grph.box_index][-1][1] - self.frame_angle)))
                if start_corner.real < turn_corner.real:
                    if cross_comp == 0:
                        self.scan_direction.append("right")
                    else:
                        self.scan_direction.append("left")
                else:
                    if cross_comp == 0:
                        self.scan_direction.append("left")
                    else:
                        self.scan_direction.append("right")
                self.scan_direction = self.scan_direction[-2:]
                self.current_tab.entry_list[self.scanDirection_index[0]][1][self.scanDirection_index[1]].set("   "+self.scan_direction[-1])

            scan_d = self.current_tab.entry_list[self.scanDirection_index[0]][1][self.scanDirection_index[1]].get()
            if True:
                start_sc, turn_sc = self.fits_WSC.pixel_to_world(self.boxPos[0], self.boxPos[1]), self.fits_WSC.pixel_to_world(self.boxPos[2], self.boxPos[3])
                self.current_tab.entry_list[self.otf_index[0]][1][self.otf_index[1]].set(round(start_sc.separation(turn_sc).arcsec, 1))
                self.Nspacing_trace_callback()

        else:
            try:
                if self.fits_opened:
                    self.current_tab.entry_list[self.angle_index[0]][1][self.angle_index[1]].set("")
                    self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].set("")
                    self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].set("")
                    self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].set("")
                    self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].set("")
                    self.current_tab.entry_list[self.scanDirection_index[0]][1][self.scanDirection_index[1]].set(" ")
                    self.current_tab.entry_list[self.otf_index[0]][1][self.otf_index[1]].set("")
                else:
                    pass
            except AttributeError:
                pass
        self.tracers_init()

    def setup_cCu(self, tabslist, notebook):
        self.notebook_num =  1
        self.current_tab = tabslist[1]
        self.angle_index = [(frame, position)
                             for frame, framelist in enumerate(self.current_tab.entry_list)
                             for position in range(len(framelist[0]))
                             if self.current_tab.entry_list[frame][0][position] == "position_angle"][0]
        self.lambet_index = [(frame, position)
                            for frame, framelist in enumerate(self.current_tab.entry_list)
                            for position in range(len(framelist[0]))
                            if self.current_tab.entry_list[frame][0][position] == "lambda_on"][0]
        self.startpos_index = [(frame, position)
                             for frame, framelist in enumerate(self.current_tab.entry_list)
                             for position in range(len(framelist[0]))
                             if self.current_tab.entry_list[frame][0][position] == "start_pos_x"][0]
        self.coordsys_index = [(frame, position)
                               for frame, framelist in enumerate(self.current_tab.entry_list)
                               for position in range(len(framelist[0]))
                               if self.current_tab.entry_list[frame][0][position] == "coordsys"][0]
        self.offabs_index = [(frame, position)
                             for frame, framelist in enumerate(self.current_tab.entry_list)
                             for position in range(len(framelist[0]))
                             if self.current_tab.entry_list[frame][0][position] == "lambda_off"][0]
        self.offrel_index = [(frame, position)
                             for frame, framelist in enumerate(self.current_tab.entry_list)
                             for position in range(len(framelist[0]))
                             if self.current_tab.entry_list[frame][0][position] == "lamdel_off"][0]
        self.relative_index = [(frame, position)
                             for frame, framelist in enumerate(self.current_tab.entry_list)
                             for position in range(len(framelist[0]))
                             if self.current_tab.entry_list[frame][0][position] == "relative"][0]
        self.scanDirection_index = [(frame, position)
                             for frame, framelist in enumerate(self.current_tab.entry_list)
                             for position in range(len(framelist[0]))
                             if self.current_tab.entry_list[frame][0][position] == "scan_direction"][0]
        self.Nspacing_index = [(frame, position)
                             for frame, framelist in enumerate(self.current_tab.entry_list)
                             for position in range(len(framelist[0]))
                             if self.current_tab.entry_list[frame][0][position] == "N"][0]
        self.otf_index = [(frame, position)
                             for frame, framelist in enumerate(self.current_tab.entry_list)
                             for position in range(len(framelist[0]))
                             if self.current_tab.entry_list[frame][0][position] == "otflen"][0]

        self.obs_sys = [None]
        self.relative = [0]
        self.scan_direction = ["right"]
        self.current_tab.entry_list[self.offabs_index[0]][1][self.offabs_index[1]].set("{}")
        self.current_tab.entry_list[self.offabs_index[0]][1][self.offabs_index[1] + 1].set("{}")
        self.current_tab.entry_list[self.offabs_index[0]][2][self.offabs_index[1]].config(state="disabled")
        self.current_tab.entry_list[self.offabs_index[0]][2][self.offabs_index[1] + 1].config(state="disabled")
        self.current_tab.entry_list[self.otf_index[0]][2][self.otf_index[1]].config(state="disabled")
        self.current_tab.entry_list[self.Nspacing_index[0]][2][self.Nspacing_index[1]].config(state="disabled")
        self.current_tab.entry_list[self.scanDirection_index[0]][1][self.scanDirection_index[1]].set(" ")

        tabslist[self.notebook_num].grph.canvas.bind("<B1-Motion>",
                    lambda event, tabslist=tabslist, notebook=notebook: self.currentCoords_update(event, tabslist, notebook), add="+")
        tabslist[self.notebook_num].grph.canvas.bind("<B2-Motion>",
                    lambda event, tabslist=tabslist, notebook=notebook: self.currentCoords_update(event, tabslist, notebook), add="+")
        tabslist[self.notebook_num].grph.canvas.bind("<B4-Motion>",
                    lambda event, tabslist=tabslist, notebook=notebook: self.currentCoords_update(event, tabslist, notebook, SD=True), add="+")

        self.tracers_init()
        if self.window_from_new:
            self.coordsys_trace_callback()

    def tracers_init(self):
        self.tracers = []
        self.tracers = [
                        [self.current_tab.entry_list[self.coordsys_index[0]][1][self.coordsys_index[1]], self.coordsys_trace_callback],
                        [self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1]], self.lambet_trace_callback],
                        [self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1], self.lambet_trace_callback],
                        [self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1]], self.lambet_trace_callback],
                        [self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1], self.lambet_trace_callback],
                        [self.current_tab.entry_list[self.angle_index[0]][1][self.angle_index[1]], self.lambet_trace_callback],
                        [self.current_tab.entry_list[self.relative_index[0]][1][self.relative_index[1]], self.relative_trace_callback],
                        [self.current_tab.entry_list[self.scanDirection_index[0]][1][self.scanDirection_index[1]], self.scanDirection_trace_callback],
                        [self.current_tab.entry_list[self.Nspacing_index[0]][1][self.Nspacing_index[1]+1], self.Nspacing_trace_callback]
                       ]
        for i,var in enumerate(self.tracers):
            self.tracers[i].append(var[0].trace("w", var[1]))

    def tracers_disable(self):
        try:
            for i,(var,callback,tracer_id) in enumerate(self.tracers):
                var.trace_vdelete("w", tracer_id)
                self.tracers[i].remove(tracer_id)
        except ValueError:
            pass

    def relative_trace_callback(self, *args, **kwargs):
        relative = str(self.current_tab.entry_list[self.relative_index[0]][1][self.relative_index[1]].get().lstrip())
        self.relative.append(relative)
        self.relative = self.relative[-2:]
        if self.relative[0] == self.relative[-1] and not "convert" in kwargs:
            return None
        if relative == "false":
            active_index, disable_index = self.offabs_index, self.offrel_index
            offset1, offset2 = 0, 1
        elif relative == "true":
            active_index, disable_index = self.offrel_index, self.offabs_index
            offset1, offset2 = 1, 0

        try:
            lam_new = str(self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].get())
            bet_new = str(self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].get())
            lambet_new = SkyCoord(Angle(lam_new), Angle(bet_new).to(u.hourangle), frame=self.obs_sys[-1])
            if not "convert" in kwargs:
                off_point_lam = str(self.current_tab.entry_list[disable_index[0]][1][disable_index[1]].get())
                off_point_bet = str(self.current_tab.entry_list[disable_index[0]][1][disable_index[1] + 1].get())
                sys_index = -1
            else:
                off_point_lam = str(self.current_tab.entry_list[active_index[0]][1][active_index[1]].get())
                off_point_bet = str(self.current_tab.entry_list[active_index[0]][1][active_index[1] + 1].get())
                offset2 = offset1
                sys_index = 0
            lambet_new = lambet_new.transform_to(self.obs_sys[sys_index])

            if type(self.obs_sys[sys_index]) is type(Galactic()):
                off_point_abs = SkyCoord(Angle(off_point_lam) + offset2 * Angle(lambet_new.l.to_string(u.hour)), Angle(off_point_bet) + offset2 * Angle(lambet_new.b.to_string(u.degree)), frame=self.obs_sys[sys_index])
            elif type(self.obs_sys[sys_index]) is type(FK5(equinox="J2000")) or type(self.obs_sys[sys_index]) is type(FK4(equinox="B1950")):
                off_point_abs = SkyCoord(Angle(off_point_lam) + offset2 * Angle(lambet_new.ra.to_string(u.hour)), Angle(off_point_bet) + offset2 * Angle(lambet_new.dec.to_string(u.degree)), frame=self.obs_sys[sys_index])
            off_point_abs = off_point_abs.transform_to(self.obs_sys[-1])
            lambet_new = lambet_new.transform_to(self.obs_sys[-1])

            if type(self.obs_sys[-1]) is type(Galactic()):
                self.current_tab.entry_list[active_index[0]][1][active_index[1]].set((off_point_abs.l - offset1*lambet_new.l).to_string(u.hour))
                self.current_tab.entry_list[active_index[0]][1][active_index[1] + 1].set((off_point_abs.b - offset1*lambet_new.b).to_string(u.degree))
            elif type(self.obs_sys[-1]) is type(FK5(equinox="J2000")) or type(self.obs_sys[-1]) is type(FK4(equinox="B1950")):
                self.current_tab.entry_list[active_index[0]][1][active_index[1]].set((off_point_abs.ra - offset1*lambet_new.ra).to_string(u.hour))
                self.current_tab.entry_list[active_index[0]][1][active_index[1] + 1].set((off_point_abs.dec - offset1*lambet_new.dec).to_string(u.degree))

        except (ValueError, AttributeError):
            self.current_tab.entry_list[active_index[0]][1][active_index[1]].set("")
            self.current_tab.entry_list[active_index[0]][1][active_index[1] + 1].set("")

        self.current_tab.entry_list[disable_index[0]][1][disable_index[1]].set("{}")
        self.current_tab.entry_list[disable_index[0]][1][disable_index[1]+1].set("{}")
        self.current_tab.entry_list[disable_index[0]][2][disable_index[1]].config(state="disabled")
        self.current_tab.entry_list[disable_index[0]][2][disable_index[1]+1].config(state="disabled")
        self.current_tab.entry_list[active_index[0]][2][active_index[1]].config(state="normal")
        self.current_tab.entry_list[active_index[0]][2][active_index[1] + 1].config(state="normal")

    def coordsys_trace_callback(self, *args):
        coordsys = self.current_tab.entry_list[self.coordsys_index[0]][1][self.coordsys_index[1]].get()
        if coordsys.lstrip() == "Galactic":
            self.obs_sys.append(Galactic())
        elif coordsys.lstrip() == "J2000":
            self.obs_sys.append(FK5(equinox="J2000"))
        elif coordsys.lstrip() == "B1950":
            self.obs_sys.append(FK4(equinox="B1950"))
        self.obs_sys = self.obs_sys[-2:]

        if not type(self.obs_sys[0]) is type(None):
            lam_new = str(self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].get())
            bet_new = str(self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].get())
            startx_new = float(self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].get())
            starty_new = float(self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].get())

            lambet_new = SkyCoord(Angle(lam_new), Angle(bet_new).to(u.hourangle), frame=self.obs_sys[0])
            startpos_new = SkyCoord(Angle(lam_new) + startx_new * u.arcsec, Angle(bet_new).to(u.hourangle) + starty_new * u.arcsec, frame=self.obs_sys[0])

            lambet_new = lambet_new.transform_to(self.obs_sys[-1])
            startpos_new = startpos_new.transform_to(self.obs_sys[-1])

            if coordsys.lstrip() == "Galactic":
                self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].set(lambet_new.l.to_string(u.hour))
                self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1]+1].set(lambet_new.b.to_string(u.degree))
                self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].set(startpos_new.l.arcsec - lambet_new.l.arcsec)
                self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1]+1].set(startpos_new.b.arcsec - lambet_new.b.arcsec)
            elif coordsys.lstrip() == "J2000" or coordsys.lstrip() == "B1950":
                self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].set(lambet_new.ra.to_string(u.hour))
                self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].set(lambet_new.dec.to_string(u.degree))
                self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].set(startpos_new.ra.arcsec - lambet_new.ra.arcsec)
                self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1]+1].set(startpos_new.dec.arcsec - lambet_new.dec.arcsec)

            self.relative_trace_callback(convert=True)

    def lambet_trace_callback(self, *args, **kwargs):
        if not self.fits_opened:
            return None
        try:
            lam_new = str(self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1]].get())
            bet_new = str(self.current_tab.entry_list[self.lambet_index[0]][1][self.lambet_index[1] + 1].get())
            startx_new = float(self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1]].get())
            starty_new = float(self.current_tab.entry_list[self.startpos_index[0]][1][self.startpos_index[1] + 1].get())
            angle = float(self.current_tab.entry_list[self.angle_index[0]][1][self.angle_index[1]].get())*math.pi/180 - self.frame_angle

            lambet_new = SkyCoord(Angle(lam_new), Angle(bet_new).to(u.hourangle), frame=self.obs_sys[-1])
            start_new = SkyCoord(Angle(lam_new) + startx_new * u.arcsec, Angle(bet_new).to(u.hourangle) + starty_new * u.arcsec, frame=self.obs_sys[-1])
            lambet_new = lambet_new.transform_to(self.fits_sys)
            start_new = start_new.transform_to(self.fits_sys)

            midPos_x, midPos_y = self.fits_WSC.world_to_pixel(lambet_new)
            start_x, start_y = self.fits_WSC.world_to_pixel(start_new)
            y_dim = self.current_tab.grph.fitsNp_ori.shape[0] - 1

            if not self.current_tab.grph.boxDrawn:
                cur_angle = angle
            else:
                cur_angle = self.current_tab.grph.Box[self.current_tab.grph.box_index][-1][1]

            oristart = complex(start_x - midPos_x, start_y - midPos_y) * cmath.exp(complex(0, -cur_angle))
            width, height = -(oristart.real * 2), (oristart.imag * 2)
            rel_start_new = oristart * cmath.exp(complex(0, angle))
            start_x_new, start_y_new = rel_start_new.real + midPos_x, rel_start_new.imag + midPos_y
            height_x, height_y = height*math.sin(angle), height*math.cos(angle)
            width_x, width_y = width*math.cos(angle), width*math.sin(angle)

            finalNW = (start_x_new, start_y_new)
            finalNE = (start_x_new + width_x, start_y_new + width_y)
            finalSE = (start_x_new + width_x + height_x, start_y_new + width_y - height_y)
            finalSW = (start_x_new + height_x, start_y_new - height_y)

            if not self.current_tab.grph.boxDrawn and not self.current_tab.grph.box_selected:
                self.current_tab.grph.canvas.create_polygon(finalNW[0], y_dim - finalNW[1],
                                                            finalNE[0], y_dim - finalNE[1],
                                                            finalSE[0], y_dim - finalSE[1],
                                                            finalSW[0], y_dim - finalSW[1], fill="", width=1, outline=self.current_tab.grph.boxColor, tag="tempbox")
                self.current_tab.grph.setBox(None)
                self.current_tab.grph.selectBox(None, virtual=True)
            elif not self.current_tab.grph.box_manip and self.current_tab.grph.box_selected and not self.current_tab.grph.clicked:
                self.current_tab.grph.canvas.coords(self.current_tab.grph.Box[self.current_tab.grph.box_index][0],
                                                    finalNW[0], y_dim - finalNW[1],
                                                    finalNE[0], y_dim - finalNE[1],
                                                    finalSE[0], y_dim - finalSE[1],
                                                    finalSW[0], y_dim - finalSW[1])
                self.current_tab.grph.setBox(None)
        except ValueError:
            pass

    def scanDirection_trace_callback(self, *args):
        self.scan_direction.append(str(self.current_tab.entry_list[self.scanDirection_index[0]][1][self.scanDirection_index[1]].get()).lstrip())
        self.scan_direction = self.scan_direction[-2:]
        if self.current_tab.grph.box_selected:
            (NWPos, NEPos, SEPos, SWPos) = \
                tuple(self.current_tab.grph.canvas.coords(self.current_tab.grph.Box[self.current_tab.grph.box_index][0])[i:i + 2] for i in range(0, 8, 2))
            if self.scan_direction[-1] == "right":
                pass
            elif self.scan_direction[-1] == "left":
                pass
            if not self.scan_direction[0] == self.scan_direction[-1]:
                self.current_tab.grph.canvas.coords(self.current_tab.grph.Box[self.current_tab.grph.box_index][0],
                                                        NEPos[0], NEPos[1], NWPos[0], NWPos[1], SWPos[0], SWPos[1], SEPos[0], SEPos[1])
                self.current_tab.grph.setBox(None)

    def Nspacing_trace_callback(self, *args):
        scan_spacing = float(self.current_tab.entry_list[self.Nspacing_index[0]][1][self.Nspacing_index[1]+1].get())
        start_sc, end_sc = self.fits_WSC.pixel_to_world(self.boxPos[0], self.boxPos[1]), self.fits_WSC.pixel_to_world(self.boxPos[6], self.boxPos[7])
        self.current_tab.entry_list[self.Nspacing_index[0]][1][self.Nspacing_index[1]].set(math.ceil(start_sc.separation(end_sc).arcsec/scan_spacing))

    def close_tl(self, tl_id):
        tl_windows.remove(tl_id)
        tl_id.destroy()
        if len(tl_windows) == 0:
            root.destroy()
