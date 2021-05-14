import os, math, cmath
import numpy as np
import tkinter as tk
from tkinter import ttk, filedialog, simpledialog
from itertools import chain
from colour import Color
from PIL import Image, ImageDraw, ImageFont, ImageTk, ImageEnhance, ImageFilter
from astropy.io import fits
from astropy import units as u
from astropy.utils.data import get_pkg_data_filename
from astropy.wcs import WCS
from astropy import wcs
import pyregion
from astropy.coordinates import SkyCoord, Angle, FK4, Galactic, FK5
import ttkwidgets as ttkw


class InitDialog:
    def __init__(self):
        InitFrame = tk.Frame(root)
        InitFrame.grid(row=0, column=0, sticky="nsew")

        global tl_windows
        tl_windows = []

        self.newButton = ttk.Button(InitFrame, text="New", takefocus=0, command=self.init_new)
        self.newButton.grid(row=0, column=0, sticky="nsew")

    def init_new(self):
        tl_windows.append(tk.Toplevel(root))
        tl_windows[-1].protocol("WM_DELETE_WINDOW", lambda tl_id=tl_windows[-1]: self.close_tl(tl_id))
        MainApplication(tl_windows[-1], "Untitled.obs", first=True)
        root.withdraw()

    def close_tl(self, tl_id):
        tl_windows.remove(tl_id)
        tl_id.destroy()
        if len(tl_windows) == 0:
            root.destroy()


class MainApplication:
    def __init__(self, master, title, **kwargs):
        self.master = master
        master.title(title)
        master.geometry("1005x610")
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
        self.tabsList = []
        self.tabsList.append(0)
        self.tabsList.append(Tabs(master, self.notebook, self.mainFrame2, self.param.paramlist, "Untitled"))
        try:
            self.Files = Files(master, self.notebook, self.mainFrame2, self.tabsList, self.param.paramlist, from_new=kwargs["first"])
        except KeyError:
            self.Files = Files(master, self.notebook, self.mainFrame2, self.tabsList, self.param.paramlist)
        master.update()
        self.tabsList[-1].grph.canvas.config(width=self.tabsList[-1].Frames[4].winfo_width(),
                                                       height=self.tabsList[-1].Frames[4].winfo_height())


class Tabs:
    def __init__(self, master, notebook, gframe_main, param, *arg):
        self.master = master
        for title in arg:
            self.tabtitle = title
            self.frames_setup(notebook, gframe_main, param)
            self.fill_obs_tab(self.Frames, param)
            self.grph = Graphic(master, self.Frames[5], self.Frames[4], self.Frames[6])

    def frames_setup(self, notebook, gframe_main, param):
        self.obs_tab = tk.Frame(notebook, width=465, padx=9, pady=0, bg='#fafafa')
        self.obs_tab.grid_propagate(False)
        notebook.add(self.obs_tab, text="OBS")

        self.fits_tab = tk.Frame(notebook, padx=0, pady=0)
        self.fits_tab.grid_propagate(False)
        notebook.add(self.fits_tab, text="FITS", state="hidden")

        self.Frames = []

        for i in range(len(param)):
            self.Frames.append(tk.LabelFrame(self.obs_tab, labelwidget=tk.Label(self.obs_tab, text=param[i][0], font="TimesNewRoman 9", bg='#fafafa'), bg='#fafafa', padx=3, pady=3))
            self.Frames[-1].grid(row=i, column=0, columnspan=3, sticky="nsew")
            #ttk.Separator(self.obs_tab, orient="horizontal").grid(row=(2*i)+1, column=0, columnspan=2, sticky="nsew")

        #self.Frames.append(tk.Frame(self.obs_tab, bg='white', padx=3, pady=3))
        #self.Frames.append(tk.Frame(self.obs_tab, bg='white', padx=3, pady=3))
        self.Frames.append(tk.Frame(gframe_main, bg='#fafafa', height=100, padx=3, pady=3))
        self.Frames.append(tk.Frame(gframe_main, bg='#ececec', padx=3, pady=3))
        #self.Frames.append(tk.Frame(gframe_main, bg='#ececec', padx=3, pady=3))
        self.Frames.append(tk.Frame(self.fits_tab, bg='#fafafa', padx=3, pady=3))

        #self.obs_tab.grid_rowconfigure([1,2,3,4], weight=1)
        #self.obs_tab.grid_columnconfigure([0,1,2,3,4], weight=1)
        self.fits_tab.grid_rowconfigure(0, weight=1)
        self.fits_tab.grid_columnconfigure(0, weight=1)
        gframe_main.grid_rowconfigure([0,2], weight=1)
        gframe_main.grid_columnconfigure([0], weight=1)

        #self.Frames[0].grid(row=0, column=0, columnspan=3, sticky="nsew")
        #ttk.Separator(self.obs_tab, orient="horizontal").grid(row=1, column=0, columnspan=2, sticky="nsew")
        #self.Frames[1].grid(row=2, rowspan=3, column=0, sticky="nsew")
        #ttk.Separator(self.obs_tab, orient="vertical").grid(row=0, rowspan=5, column=1, sticky="nsew")
        #self.Frames[2].grid(row=0, rowspan=5, column=2, sticky="nsew")

        self.Frames[4].grid(row=2, column=0, rowspan=3, columnspan=3, sticky="nsew")
        ttk.Separator(gframe_main, orient="horizontal").grid(row=1, column=0, columnspan=3, sticky="nsew")
        self.Frames[5].grid(row=0, rowspan=1, column=0, sticky="nsew")
        ttk.Separator(gframe_main, orient="vertical").grid(row=0, rowspan=1, column=1, sticky="nsew")
        #self.Frames[5].grid(row=0, rowspan=1, column=2, sticky="nsew")

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


class Graphic:
    def __init__(self, master, Gframe, Sframe, Tframe):
        self.master=master
        self.boxColor, self.regColor = "red", "red"
        self.degChange = 0.
        self.Box = []
        self.box_selected,  self.over_object, self.over_selected, self.clicked, self.box_manip, self.box_resize, self.boxDrawn = \
            False, False, False, False, False, False, False
        self.fits_offset = (0, 0)

        Gframe.grid_rowconfigure(0, weight=1)
        Gframe.grid_columnconfigure(0, weight=1)
        self.canvas = tk.Canvas(Gframe, bg="gray", cursor="cross-hair", highlightthickness=0, yscrollincrement = 1, xscrollincrement = 1)
        self.canvas.grid(row=0, column=0, sticky="nsew")
        self.vbar = ttk.Scrollbar(Gframe, orient="vertical", command=self.canvas.yview)
        self.vbar.grid(row=0, column=1, rowspan=1, sticky="ns")
        self.hbar = ttk.Scrollbar(Gframe, orient="horizontal", command=self.canvas.xview)
        self.hbar.grid(row=1, column=0, columnspan=2, sticky="ew")

        Sframe.grid_rowconfigure([0,1], weight=0)
        Sframe.grid_columnconfigure([0], weight=1)
        self.clearButton = ttk.Button(Sframe, text="Clear", state=tk.DISABLED, takefocus=0, command=lambda _=None, mode="select": self.resetBox(_, mode))
        self.clearButton.grid(row=0, column=1, rowspan=2, sticky="nsew")
        self.canvas.bind("<BackSpace>", lambda event, mode="select": self.resetBox(event, mode))
        self.slidercanvas_1 = tk.Canvas(Sframe, bg="gray", width=200, height=10, highlightthickness=0, highlightcolor="black",
                                        borderwidth=0, relief=tk.GROOVE)
        self.slidercanvas_2 = tk.Canvas(Sframe, bg="gray", width=200, height=10, highlightthickness=0, highlightcolor="black",
                                        borderwidth=0, relief=tk.GROOVE)
        self.slidercanvas_1.grid(row=0, column=0, rowspan=1, columnspan=1, sticky="ew")
        self.slidercanvas_2.grid(row=1, column=0, rowspan=1, columnspan=1, sticky="ew")

        self.scBg = []
        for num, s_set in enumerate(((self.slidercanvas_1, ["box","marker"]), (self.slidercanvas_2, ["reg","regmarker"]))):
            self.colorSlider(s_set[0], "red", "pink")
            s_set[0].bind("<Button-1>",
                          lambda event, canvas=s_set[0], target=[s_set[1],num], pad=1: self.sliderNob(event, canvas, target, pad))
            s_set[0].bind("<B1-Motion>",
                          lambda event, canvas=s_set[0], target=[s_set[1],num], pad=1: self.sliderNob(event, canvas, target, pad))

        self.canvas.bind("<Enter>", lambda event, mode="default": self.cursors(event, mode))
        self.canvas.bind("<Button-1>", self.B12_callback)
        self.canvas.bind("<B1-Motion>", self.B1M_callback)
        self.canvas.bind("<ButtonRelease-1>", self.B1R_callback)
        self.canvas.bind("<ButtonRelease-2>", self.B2R_callback)
        self.canvas.bind('<MouseWheel>', self.v_scroll)
        self.canvas.bind('<Shift-MouseWheel>', self.h_scroll)

    def B12_callback(self, event):
        self.canvas.focus_set()
        self.bad_start = False
        self.clicked = True
        if not False:
            if self.box_selected or self.boxDrawn:
                self.deselectBox(event, "B12")
            else:
                self.setStart_draw(event)
        else:
            self.bad_start = True

    def B1M_callback(self, event):
        if not self.box_selected and not self.boxDrawn and self.clicked:
            self.drawBox(event)

    def B1R_callback(self, event):
        if self.box_selected and not self.over_object:
            self.deselectBox(event, "B1R")
        elif self.box_manip:
            self.setBox(event)
        self.clicked = False

    def B2R_callback(self, event):
        if self.box_selected and self.box_manip:
            self.setBox(event)

    def slider_master(self, vartuple):
        self.fits_offset = (0, 0)
        self.canvas.focus_set()
        try:
            self.temp = self.color_func(*vartuple)
            self.temp2 = self.temp.copy()
            self.temp2 = np.clip(self.temp2, np.min(self.fitsNp_ori), np.max(self.fitsNp_ori))
            self.temp2 = (self.temp2 / np.max(self.fitsNp_ori)) * 255
            self.temp = Image.fromarray(np.flip(self.temp2, 0))
            self.temp = ImageTk.PhotoImage(self.temp)
            self.canvas.delete("fits")
            self.canvas.create_image(self.fits_offset, image=self.temp, anchor="nw", tag="fits")
            self.canvas.tag_raise("O")
        except AttributeError:
            pass

    def color_func(self, mode, pixel, lowerb, upperb, pixel_min, pixel_max, gamma, gain, bias_x, bias_y):
        f1 = pixel_min
        f4 = pixel_max
        if mode == "　Symmetric":
            bias_sep_x = lowerb + (bias_x * (upperb - lowerb))
            bias_sep_y = pixel_min + (bias_y * (pixel_max - pixel_min))
            f2 = lambda pixel: -((-(pixel-bias_sep_x)/(0.5*(upperb-lowerb)))**gamma)*(0.5*(pixel_max-pixel_min))*gain + bias_sep_y
            f3 = lambda pixel: (((pixel-bias_sep_x)/(0.5*(upperb-lowerb)))**gamma)*(0.5*(pixel_max-pixel_min))*gain + bias_sep_y
            return np.piecewise(pixel, [(pixel<lowerb), (lowerb<=pixel)*(pixel<bias_sep_x),
                                        (bias_sep_x<=pixel)*(pixel<=upperb),(upperb<pixel)], [f1, f2, f3, f4])
        elif mode == "　Regular　":
            bias_sep_x = lowerb + ((bias_x-0.5) * (upperb - lowerb))
            bias_sep_y = pixel_min + ((bias_y-0.5) * (pixel_max - pixel_min))
            f23 = lambda pixel: (((pixel-bias_sep_x)/(upperb-lowerb))**gamma)*(pixel_max-pixel_min)*gain + bias_sep_y
            return np.piecewise(pixel, [(pixel<lowerb), (lowerb <= pixel) * (pixel <= upperb), (upperb < pixel)], [f1, f23, f4])

    def v_scroll(self, event):
        self.canvas.yview_scroll(-3 * event.delta, 'units')

    def h_scroll(self, event):
        self.canvas.xview_scroll(-3 * event.delta, 'units')

    def setStart_draw(self, event):
        self.startPos = (self.canvas.canvasx(event.x), self.canvas.canvasy(event.y))

    def drawBox(self, event):
        self.box_manip = True
        self.canvas.delete("tempbox")
        self.endPos = (self.canvas.canvasx(event.x), self.canvas.canvasy(event.y))
        self.NEPos = (self.canvas.canvasx(event.x), self.startPos[1])
        self.SWPos = (self.startPos[0], self.canvas.canvasy(event.y))
        self.tempBox = self.canvas.create_polygon(self.startPos, self.NEPos, self.endPos, self.SWPos,
                                                  fill="", width=1, outline=self.boxColor, tag="tempbox")

    def setBox(self, _, **kwargs):
        try:
            self.Box.append([self.canvas.create_polygon(kwargs["REG"], fill="", width=1, outline=self.regColor, tag="newbox")])
            self.box_id = self.canvas.find_withtag("newbox")[0]
            self.box_index = -1
            self.canvas.itemconfig("newbox", tag=("O", "reg"))
        except KeyError:
            if not self.box_selected and not self.boxDrawn:
                tempbox_coords = self.canvas.coords("tempbox")
                self.canvas.delete("tempbox")
                self.Box.append([self.canvas.create_polygon(tempbox_coords, fill="", width=1, outline=self.boxColor, tag="newbox")])
                self.box_id = self.canvas.find_withtag("newbox")[0]
                self.box_index = -1
                self.canvas.itemconfig("newbox", tag=("O", "box"))
                self.over_selected = False
                self.boxDrawn = True
        if True:
            try:
                self.id_str = str(self.box_id)
                self.canvas.delete("mid"+self.id_str, "resize"+self.id_str,"marker"+self.id_str,"regmarker"+self.id_str)

                (self.NWPos, self.NEPos, self.SEPos, self.SWPos) = tuple(self.canvas.coords(self.box_id)[i:i + 2] for i in range(0, 8, 2))
                self.UPos = ((self.NWPos[0] + self.NEPos[0]) / 2, (self.NWPos[1] + self.NEPos[1]) / 2)
                self.DPos = ((self.SEPos[0] + self.SWPos[0]) / 2, (self.SEPos[1] + self.SWPos[1]) / 2)
                self.LPos = ((self.NWPos[0] + self.SWPos[0]) / 2, (self.NWPos[1] + self.SWPos[1]) / 2)
                self.RPos = ((self.SEPos[0] + self.NEPos[0]) / 2, (self.SEPos[1] + self.NEPos[1]) / 2)
                self.midPos = ((self.LPos[0] + self.RPos[0]) / 2, (self.UPos[1] + self.DPos[1]) / 2)
                try:
                    self.degNow = kwargs["REGdeg"]*math.pi/180
                except KeyError:
                    self.degNow = -(((0.5 * math.pi) - cmath.phase(complex(self.UPos[0] - self.midPos[0], self.midPos[1] - self.UPos[1])) - (2*math.pi)) % (-2*math.pi))

                self.resizeFill = "white"

                self.Box[self.box_index] = [self.box_id]

                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.NWPos[0] - 2, self.NWPos[1] - 2, self.NWPos[0] + 2, self.NWPos[1] + 2,
                                            width=1, fill=self.boxColor, outline=self.boxColor, tag=("O", "box", "marker", "marker"+self.id_str)))
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.NEPos[0] - 1, self.NEPos[1] - 1, self.NEPos[0] + 1, self.NEPos[1] + 1,
                                            width=1, fill=self.boxColor, outline=self.boxColor, tag=("O", "box", "marker", "marker"+self.id_str)))
                if "reg" in self.canvas.gettags(self.box_id):
                    self.canvas.itemconfig(self.Box[self.box_index][1], fill=self.regColor, outline=self.regColor, tag=("O", "reg", "regmarker", "regmarker"+self.id_str))
                    self.canvas.itemconfig(self.Box[self.box_index][2], fill=self.regColor, outline=self.regColor, tag=("O", "reg", "regmarker", "regmarker"+self.id_str))
                if self.box_selected:
                    self.canvas.itemconfig(self.Box[self.box_index][1], fill="")
                    self.canvas.itemconfig(self.Box[self.box_index][2], fill="")

                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.NWPos[0] - 6, self.NWPos[1] + 5, self.NWPos[0] + 5, self.NWPos[1] - 6,
                                        width=0, fill=self.resizeFill, tag=("O", "resize", "resize" + self.id_str, "C", "NW")))
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.NEPos[0] - 6, self.NEPos[1] + 5, self.NEPos[0] + 5, self.NEPos[1] - 6,
                                        width=0, fill=self.resizeFill, tag=("O", "resize", "resize" + self.id_str, "C", "NE")))
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.SEPos[0] - 6, self.SEPos[1] + 5, self.SEPos[0] + 5, self.SEPos[1] - 6,
                                        width=0, fill=self.resizeFill, tag=("O", "resize", "resize" + self.id_str, "C", "SE")))
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.SWPos[0] - 6, self.SWPos[1] + 5, self.SWPos[0] + 5, self.SWPos[1] - 6,
                                        width=0, fill=self.resizeFill, tag=("O", "resize", "resize" + self.id_str, "C", "SW")))
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.UPos[0] - 6, self.UPos[1] + 5, self.UPos[0] + 5, self.UPos[1] - 6,
                                        width=0, fill=self.resizeFill, tag=("O", "resize", "resize" + self.id_str, "UD", "U")))
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.DPos[0] - 6, self.DPos[1] + 5, self.DPos[0] + 5, self.DPos[1] - 6,
                                        width=0, fill=self.resizeFill, tag=("O", "resize", "resize" + self.id_str, "UD", "D")))
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.LPos[0] - 6, self.LPos[1] + 5, self.LPos[0] + 5, self.LPos[1] - 6,
                                        width=0, fill=self.resizeFill, tag=("O", "resize", "resize" + self.id_str, "LR", "L")))
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.RPos[0] - 6, self.RPos[1] + 5, self.RPos[0] + 5, self.RPos[1] - 6,
                                        width=0, fill=self.resizeFill, tag=("O", "resize", "resize" + self.id_str, "LR", "R")))
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.midPos[0] - 6, self.midPos[1] + 5, self.midPos[0] + 5, self.midPos[1] - 6,
                                        width=0, fill=self.resizeFill, tag=("O", "mid", "mid" + self.id_str)))
                self.Box[self.box_index].append((self.midPos, self.degNow))

                self.canvas.tag_bind("all", "<Enter>", lambda event, mode="Enter": self.hover_detect(event, mode))
                self.canvas.tag_bind("all", "<Leave>", lambda event, mode="Leave": self.hover_detect(event, mode))
                self.canvas.tag_bind("all", "<Button-1>", lambda event: self.setStart_manip(event))
                self.canvas.tag_bind("C", "<Button-2>", lambda event: self.setStart_manip(event))

                self.canvas.tag_bind("C", "<B2-Motion>", lambda event, mode=("rotate", None): self.manipulateBox(event, mode))
                self.canvas.tag_bind("C", "<Leave>", lambda event, mode="default": self.cursors(event, mode))
                self.canvas.tag_bind("all", "<Leave>", lambda event, mode="default": self.cursors(event, mode))

                self.canvas.tag_bind("box", "<Enter>", lambda event, mode="move": self.cursors(event, mode))
                self.canvas.tag_bind("box", "<B1-Motion>", lambda event, mode=("move", None): self.manipulateBox(event, mode))
                self.canvas.tag_bind("reg", "<Enter>", lambda event, mode="move": self.cursors(event, mode))
                self.canvas.tag_bind("reg", "<B1-Motion>", lambda event, mode=("move", None): self.manipulateBox(event, mode))
                self.canvas.tag_bind("mid", "<Enter>", lambda event, mode="move": self.cursors(event, mode))
                self.canvas.tag_bind("mid", "<B1-Motion>", lambda event, mode=("move", None): self.manipulateBox(event, mode))

                self.canvas.tag_bind("UD", "<Enter>", lambda event, mode="hand": self.cursors(event, mode))
                self.canvas.tag_bind("U", "<B1-Motion>", lambda event, mode=("U", "stretch"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("D", "<B1-Motion>", lambda event, mode=("D", "stretch"): self.manipulateBox(event, mode))

                self.canvas.tag_bind("LR", "<Enter>", lambda event, mode="hand": self.cursors(event, mode))
                self.canvas.tag_bind("L", "<B1-Motion>", lambda event, mode=("L", "stretch"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("R", "<B1-Motion>", lambda event, mode=("R", "stretch"): self.manipulateBox(event, mode))

                self.canvas.tag_bind("C", "<Enter>", lambda event, mode="hand": self.cursors(event, mode))
                self.canvas.tag_bind("C", "<B2-Motion>", lambda event, mode=("rotate", None): self.manipulateBox(event, mode))
                self.canvas.tag_bind("NW", "<B1-Motion>", lambda event, mode=("NW", "free"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("NW", "<Shift-B1-Motion>", lambda event, mode=("NW", "ratio"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("NE", "<B1-Motion>", lambda event, mode=("NE", "free"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("NE", "<Shift-B1-Motion>", lambda event, mode=("NE", "ratio"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("SW", "<B1-Motion>", lambda event, mode=("SW", "free"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("SW", "<Shift-B1-Motion>", lambda event, mode=("SW", "ratio"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("SE", "<B1-Motion>", lambda event, mode=("SE", "free"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("SE", "<Shift-B1-Motion>", lambda event, mode=("SE", "ratio"): self.manipulateBox(event, mode))

                self.canvas.tag_bind("all", "<ButtonRelease-1>", lambda event: self.selectBox(event))

                self.box_manip, self.box_resize = False, False
                self.box_index = self.box_index_selected
                self.box_id = self.Box[self.box_index][0]
                self.degChange = 0.
                self.canvas.event_generate("<B4-Motion>")
                self.canvas.event_generate("<Enter>")
            except (AttributeError, IndexError):
                pass

    def B12_leave(self, event):
        try:
            if self.boxDrawn:
                self.endPos_temp = (self.canvas.canvasx(event.x), self.canvas.canvasy(event.y))
                self.boxPos_temp = self.canvas.coords(self.Box[self.box_index_hover][0])
                if max(self.boxPos_temp) in self.endPos_temp or min(self.boxPos_temp) in self.endPos_temp:
                    self.canvas.event_generate("<Leave>")
        except AttributeError:
            pass

    def hover_detect(self, event, mode):
        self.over_selected = False
        if mode == "Enter":
            self.over_object = True
            try:
                if self.canvas.find_withtag("current")[0] in self.Box[self.box_index] and self.box_selected:
                    self.over_selected = True
                else:
                    self.box_index_hover = [index for index, box in enumerate(self.Box) if self.canvas.find_withtag("current")[0] in box][0]
            except IndexError:
                pass
        elif mode == "Leave":
            self.over_object = False
            self.cursors(event, "default")

    def selectBox(self, event, **kwargs):
        try:
            if not self.over_selected:
                try:
                    if kwargs["virtual"]:
                        self.box_index = len(self.Box)-1
                        physical = False
                        self.over_object = True
                except KeyError:
                    self.box_index = [index for index, box in enumerate(self.Box) if self.canvas.find_withtag("current")[0] in box][0]
                    physical = True
                self.box_index_selected = self.box_index
                self.box_id = self.Box[self.box_index][0]
                if self.over_object and not self.box_manip:
                    self.canvas.tag_raise(self.box_id)
                    self.canvas.tag_raise("marker"+str(self.box_id))
                    self.canvas.tag_raise("resize"+str(self.box_id))
                    self.canvas.itemconfig("box", dash=())
                    self.canvas.itemconfig("reg", dash=())
                    self.canvas.itemconfig("marker", fill=self.boxColor)
                    self.canvas.itemconfig("regmarker", fill=self.regColor)
                    self.canvas.itemconfig(self.box_id, dash=(4, 5))
                    self.canvas.itemconfig(self.Box[self.box_index][1], fill="")
                    self.canvas.itemconfig(self.Box[self.box_index][2], fill="")
                    self.box_selected = True
                    self.over_object = physical
                    self.over_selected = physical
                    self.canvas.event_generate("<B4-Motion>")
                    self.cursors(event, "move")
                    self.clearButton.config(state=tk.ACTIVE)
        except IndexError:
            self.deselectBox(event, "B12")
            if self.box_selected:
                self.deselectBox(event, "B1R")

    def deselectBox(self, event, mode):
        if mode == "B12":
            self.clicked = True
        elif mode == "B1R":
            if self.clicked:
                self.canvas.itemconfig("box", dash=())
                self.canvas.itemconfig("reg", dash=())
                self.canvas.itemconfig("marker", fill=self.boxColor)
                self.canvas.itemconfig("regmarker", fill=self.regColor)
                self.box_selected = False
                self.over_selected = False
                self.clicked = False
                self.clearButton.config(state=tk.DISABLED)
                self.canvas.event_generate("<B4-Motion>")

    def resetBox(self, _, mode):
        self.canvas.focus_set()
        if mode == "select":
            if self.box_selected:
                self.box_selected = False
                self.canvas.event_generate("<B4-Motion>")
                if not "reg" in self.canvas.gettags(self.box_id):
                    self.boxDrawn = False
                self.over_selected = False
                self.clicked = False
                self.bad_start = False
                self.clearButton.config(state=tk.DISABLED)
                for i in range(len(self.Box[self.box_index])):
                    self.canvas.delete(self.Box[self.box_index][i])

    def setStart_manip(self, event):
        if self.over_selected:
            try:
                self.startPos = (self.canvas.canvasx(event.x), self.canvas.canvasy(event.y))
                (self.NWPos, self.NEPos, self.SEPos, self.SWPos) = tuple(self.canvas.coords(self.box_id)[i:i + 2] for i in range(0, 8, 2))
                self.UPos = ((self.NWPos[0] + self.NEPos[0]) / 2, (self.NWPos[1] + self.NEPos[1]) / 2)
                self.DPos = ((self.SEPos[0] + self.SWPos[0]) / 2, (self.SEPos[1] + self.SWPos[1]) / 2)
                self.LPos = ((self.NWPos[0] + self.SWPos[0]) / 2, (self.NWPos[1] + self.SWPos[1]) / 2)
                self.RPos = ((self.SEPos[0] + self.NEPos[0]) / 2, (self.SEPos[1] + self.NEPos[1]) / 2)
                self.boxPos = self.canvas.coords(self.box_id)
            except AttributeError:
                pass

    def manipulateBox(self, event, mode):
        try:
            self.canvas.tag_unbind("all", "<Leave>")
            self.finalPos = self.boxPos
            self.midPos = self.Box[self.box_index][-1][0]
            self.degNow = self.Box[self.box_index][-1][1]
            if self.over_selected:
                self.box_manip = True
                if mode[0] == "move":
                    self.posChange = [self.canvas.canvasx(event.x) - self.startPos[0],
                                      self.canvas.canvasy(event.y) - self.startPos[1]]
                    self.finalPos = (self.boxPos[0] + self.posChange[0], self.boxPos[1] + self.posChange[1],
                                     self.boxPos[2] + self.posChange[0], self.boxPos[3] + self.posChange[1],
                                     self.boxPos[4] + self.posChange[0], self.boxPos[5] + self.posChange[1],
                                     self.boxPos[6] + self.posChange[0], self.boxPos[7] + self.posChange[1])
                elif mode[0] == "rotate":
                    self.degChange = -(
                        ((cmath.phase(complex(self.startPos[0] - self.midPos[0], self.midPos[1] - self.startPos[1]))
                          - cmath.phase(complex(self.canvas.canvasx(event.x) - self.midPos[0],
                                                self.midPos[1] - self.canvas.canvasy(event.y))))
                         % (2 * math.pi)) / (math.pi / 900)) * (math.pi / 900) + (2*math.pi)
                    self.finalPos = []
                    for i in range(0, 8, 2):
                        self.dummyPos = complex(self.boxPos[i] - self.midPos[0], self.midPos[1] - self.boxPos[i + 1])
                        self.dummyPos = self.dummyPos * cmath.exp(complex(0, self.degChange))
                        self.finalPos.append(self.dummyPos.real + self.midPos[0])
                        self.finalPos.append(self.midPos[1] - self.dummyPos.imag)
                    self.finalPos = tuple(self.finalPos)
                elif mode[0] == "U":
                    self.posChange = self.manipulateVar(event, mode[1], self.UPos, self.DPos, 0.5, 1,-1,0,0, 0,1,0,0,0,0,0,0, 1)
                    self.finalPos = (self.boxPos[0] + self.posChange[0], self.boxPos[1] + self.posChange[1],
                                     self.boxPos[2] + self.posChange[0], self.boxPos[3] + self.posChange[1],
                                     self.boxPos[4], self.boxPos[5], self.boxPos[6], self.boxPos[7])
                elif mode[0] == "D":
                    self.posChange = self.manipulateVar(event, mode[1], self.DPos, self.UPos, 0.5, 1,-1,0,0, 0,1,0,0,0,0,0,0, 1)
                    self.finalPos = (self.boxPos[0], self.boxPos[1], self.boxPos[2], self.boxPos[3],
                                     self.boxPos[4] + self.posChange[0], self.boxPos[5] + self.posChange[1],
                                     self.boxPos[6] + self.posChange[0], self.boxPos[7] + self.posChange[1])
                elif mode[0] == "L":
                    self.posChange = self.manipulateVar(event, mode[1], self.LPos, self.RPos, 1, -1,-1,0,0, 0,0,0,1,0,0,0,0, 1)
                    self.finalPos = (self.boxPos[0] + self.posChange[0], self.boxPos[1] + self.posChange[1], self.boxPos[2], self.boxPos[3],
                                     self.boxPos[4], self.boxPos[5], self.boxPos[6] + self.posChange[0], self.boxPos[7] + self.posChange[1])
                elif mode[0] == "R":
                    self.posChange = self.manipulateVar(event, mode[1], self.RPos, self.LPos, 1, -1,-1,0,0, 0,0,0,1,0,0,0,0, 1)
                    self.finalPos = (self.boxPos[0], self.boxPos[1], self.boxPos[2] + self.posChange[0], self.boxPos[3] + self.posChange[1],
                                     self.boxPos[4] + self.posChange[0], self.boxPos[5] + self.posChange[1], self.boxPos[6], self.boxPos[7])
                elif mode[0] == "NW":
                    self.posChange = self.manipulateVar(event, mode[1], self.SEPos, self.NWPos, 0.5, 1,-1,1,1, 0,1,0,0,1,0,1,1, 0)
                    self.finalPos = (self.endPos[0], self.endPos[1], self.boxPos[2] + self.posChange[0], self.boxPos[3] + self.posChange[1],
                                     self.boxPos[4], self.boxPos[5], self.boxPos[6] + self.posChange[2], self.boxPos[7] + self.posChange[3])
                elif mode[0] == "NE":
                    self.posChange = self.manipulateVar(event, mode[1], self.SWPos, self.NEPos, 0.5, 1,-1,1,1, 0,1,0,0,1,0,1,1, 0)
                    self.finalPos = (self.boxPos[0] + self.posChange[0], self.boxPos[1] + self.posChange[1], self.endPos[0], self.endPos[1],
                                     self.boxPos[4] + self.posChange[2], self.boxPos[5] + self.posChange[3], self.boxPos[6], self.boxPos[7])
                elif mode[0] == "SW":
                    self.posChange = self.manipulateVar(event, mode[1], self.NEPos, self.SWPos, 0, 1,1,-1,1, 0,0,0,1,1,1,1,0, 0)
                    self.finalPos = (self.boxPos[0] + self.posChange[0], self.boxPos[1] + self.posChange[1], self.boxPos[2], self.boxPos[3],
                                     self.boxPos[4] + self.posChange[2], self.boxPos[5] + self.posChange[3], self.endPos[0], self.endPos[1])
                elif mode[0] == "SE":
                    self.posChange = self.manipulateVar(event, mode[1], self.NWPos, self.SEPos, 0, 1,1,-1,1, 0,0,0,1,1,1,1,0, 0)
                    self.finalPos = (self.boxPos[0], self.boxPos[1], self.boxPos[2] + self.posChange[0], self.boxPos[3] + self.posChange[1],
                        self.endPos[0], self.endPos[1], self.boxPos[6] + self.posChange[2], self.boxPos[7] + self.posChange[3])

                self.canvas.coords(self.box_id, self.finalPos)
                self.canvas.coords(self.Box[self.box_index][1], self.finalPos[0]-2, self.finalPos[1]-2, self.finalPos[0]+2, self.finalPos[1]+2)
                self.canvas.coords(self.Box[self.box_index][2], self.finalPos[2]-1, self.finalPos[3]-1, self.finalPos[2]+1, self.finalPos[3]+1)
        except AttributeError:
            pass

    def manipulateVar(self, event, mode1, ref1, ref2, o, a1, a2, a3, a4, phi11, phi12, phi21, phi22, phi31, phi32, phi41, phi42, r):
        self.box_resize = True
        self.endPos = (self.canvas.canvasx(event.x), self.canvas.canvasy(event.y))
        if mode1 == "stretch":
            self.refPos = [ref1[0] - ref2[0], ref1[1] - ref2[1]]
        elif mode1 == "free":
            self.refPos = [0, 0]
        elif mode1 == "ratio":
            self.refPos = [0, 0]
            try:
                self.diagVector = [ref2[0] - ref1[0], ref2[1] - ref1[1]]
                self.posChange = [self.endPos[0] - ref2[0], self.endPos[1] - ref2[1]]
                self.diagDeg = cmath.phase(complex(self.diagVector[0], self.diagVector[1]))
                self.projDeg = self.diagDeg - cmath.phase(complex(self.posChange[0], self.posChange[1]))
                self.posChangeNorm = np.linalg.norm(self.posChange)
                self.endPos = (ref2[0] + self.posChangeNorm * math.cos(self.projDeg) * math.cos(self.diagDeg),
                               ref2[1] + self.posChangeNorm * math.cos(self.projDeg) * math.sin(self.diagDeg))
            except ZeroDivisionError:
                pass

        self.posChange = [self.endPos[0] - ref2[0], ref2[1] - self.endPos[1]]
        self.projDeg = (o * math.pi) - (cmath.phase(complex(self.posChange[0], self.posChange[1]))) + self.degNow
        self.posChangeNorm = np.linalg.norm(self.posChange)
        self.posChange = \
            [a1 * (self.posChangeNorm * math.cos(self.projDeg - phi11*0.5*math.pi) * math.cos(-self.degNow - phi12*0.5*math.pi)) - (r * self.refPos[0]),
             a2 * (self.posChangeNorm * math.cos(self.projDeg - phi21*0.5*math.pi) * math.cos(-self.degNow - phi22*0.5*math.pi)) - (r * self.refPos[1]),
             a3 * (self.posChangeNorm * math.cos(self.projDeg - phi31*0.5*math.pi) * math.cos(-self.degNow - phi32*0.5*math.pi)),
             a4 * (self.posChangeNorm * math.cos(self.projDeg - phi41*0.5*math.pi) * math.cos(-self.degNow - phi42*0.5*math.pi))]

        return self.posChange

    def colorSlider(self, canvas, color1, color2):
        canvas.update()
        self.sliderwidth, self.sliderheight = canvas.winfo_width(), canvas.winfo_height()
        self.scBg.append([])
        self.colorGradient = list(Color(color1).range_to(Color(color2), self.sliderwidth))
        self.gradBg = Image.new("RGB", (self.sliderwidth, self.sliderheight), "#FFFFFF")
        self.gradBgDraw = ImageDraw.Draw(self.gradBg)
        for x, color in enumerate(self.colorGradient):
            self.gradBgDraw.line((x, 0, x, self.sliderheight), fill=str(color), width=1)
        self.scBg[-1].append(ImageTk.PhotoImage(self.gradBg))
        self.scBg[-1].append(canvas.create_image(0, 0, image=self.scBg[-1][0], anchor=tk.NW))
        self.scBg[-1].append(canvas.create_line(1, 0, 1, self.sliderheight, width=2, fill="#444444"))

    def sliderNob(self, event, canvas, target, pad):
        self.canvas.focus_set()
        width, height = self.sliderwidth, self.sliderheight
        try:
            if width > event.x > pad:
                canvas.coords(self.scBg[target[1]][2], event.x, pad - 1, event.x, height)
            elif event.x <= pad:
                canvas.coords(self.scBg[target[1]][2], pad, pad - 1, pad, height)
            elif width <= event.x:
                canvas.coords(self.scBg[target[1]][2], width - pad, pad - 1, width - pad, height)
            if target[0][0] == "box":
                self.boxColor = str(self.colorGradient[int(canvas.coords(self.scBg[target[1]][2])[0] + pad - 1)])
                temp_color = self.boxColor
            elif target[0][0] == "reg":
                self.regColor = str(self.colorGradient[int(canvas.coords(self.scBg[target[1]][2])[0] + pad - 1)])
                temp_color = self.regColor
            self.canvas.itemconfig(target[0][0], outline=temp_color)
            self.canvas.itemconfig(target[0][1], fill=temp_color)
            if self.box_selected:
                self.canvas.itemconfig(self.Box[self.box_index][1], fill="")
                self.canvas.itemconfig(self.Box[self.box_index][2], fill="")
            if not self.resizeFill == "":
                self.canvas.itemconfig("resize", fill=temp_color)
        except AttributeError:
            pass

    def cursors(self, _, mode):
        if self.over_selected:
            if mode == "default":
                self.canvas.config(cursor="cross-hair")
            elif mode == "hand":
                self.canvas.config(cursor="openhand")
            elif mode == "move":
                self.canvas.config(cursor="fleur")
            elif mode == "rotate":
                self.canvas.config(cursor="exchange")
        else:
            self.canvas.config(cursor="cross-hair")


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


class Parameters:
    def __init__(self):
        self.paramlist = []
        self.paramlist.append(["observation_property",[
            {"name": "observer", "#": "name of observer", "string": True, "no.": 1},
            {"name": "object", "#": "name of target object", "string": True, "no.": 2},
            {"name": "molecule_1", "#": "line identifier of target molecule", "string": True, "no.": 3}
        ]])
        self.paramlist.append(["coordinate",[
            {"name": "lambda_on", "#": "x-coordinate of map center (degree)", "string": True, "grid": (1,0), "no.": 1},
            {"name": "beta_on", "#": "y-coordinate of map center (degree)", "string": True, "grid": (2,0), "no.": 2},
            {"name": "relative", "#": "OFF position is given by relative (0) or absolute (1)", "string": False, "grid": (4, 0), "values": ["true", "false"], "default":"false", "no.": 3},
            {"name": "lambda_off", "#": "off point x-coordinate (= 0 if relative == 1) (degree)", "string": True, "grid": (5,0), "no.": 4},
            {"name": "beta_off", "#": "off point y-coordinate (= 0 if relative == 1) (degree)", "string": True, "grid": (6,0), "no.": 5},
            {"name": "lamdel_off", "#": "offset_off (= 0 if relative == 0)", "string": True, "grid": (5,1), "no.": 6},
            {"name": "betdel_off", "#": "offset_off (= 0 if relative == 0)", "string": True, "grid": (6,1), "no.": 7},
            {"name": "position_angle", "#": "position angle of the scan direction (default = 0) (degree)", "string": True, "unit": "deg", "grid": (3,0), "no.": 8},
            {"name": "otadel", "#": "dcos (Y/N)", "string": False, "values": ["true", "false"], "grid": (4,1), "no.": 9},
            {"name": "start_pos_x", "#": "start position X relative to lambda_on (arcsec) (sw)", "string": True, "unit": "arcsec", "grid": (1,1), "no.": 10},
            {"name": "start_pos_y", "#": "start position Y relative to lambda_on (arcsec) (sw)", "string": True, "unit": "arcsec", "grid": (2,1), "no.": 11},
            {"name": "coordsys", "#": "(J2000/b1950/galactic/horizontal)", "string": True, "values": ["J2000", "B1950", "Galactic"], "grid": (0,0), "no.": 12},
        ]])
        self.paramlist.append(["scan_property",[
            {"name": "scan_direction", "#": "X (=0) or Y (=1) scan", "values": ["right", "left"], "string": True, "grid": (0,0), "no.": 1},
            {"name": "N", "#": "Number of scanlines (integer)", "string": False, "grid": (1,0), "no.": 2},
            {"name": "scan_spacing", "#": "scan spacing (arcsec) (default = 60)", "string": True, "unit": "arcsec", "default": 60, "grid": (2,0), "no.": 3},
            {"name": "otflen", "#": "time duration for 1 scan  (s) ", "string": True, "unit": "s", "grid": (1,1), "no.": 4},
            {"name": "otfvel", "#": "scan velocity (arcsec/s) (default = 600) (<=600)", "string": True, "unit": "arcsec/s", "default": 600, "grid": (2,1), "no.": 5},
            {"name": "ramp_pixel", "#": "approach length (pixel) (default = 4)", "string": False, "grid": (3,0), "no.": 6},
            {"name": "integ_on", "#": "ON point integratoin time (s)", "string": True, "unit": "s", "grid": (4,0), "no.": 7}
        ]])
        self.paramlist.append(["calibration",[
            {"name": "integ_off", "#": "OFF point integration time (s)", "string": True, "unit": "s", "no.": 1},
            {"name": "integ_hot", "#": "hot-load integration time (s)", "string": True, "unit": "s", "no.": 2},
            {"name": "off_interval", "#": "off interval (number of scan lines: default = 1)", "string": False, "default": 1, "no.": 3},
            {"name": "load_interval", "#": "time interval of hot-load measurements (min)", "string": True, "unit": "min", "no.": 4}
        ]])


if __name__ == '__main__':
    root = tk.Tk()
    initdialog = InitDialog()
    root.mainloop()


