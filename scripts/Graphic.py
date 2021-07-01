import math, cmath
import numpy as np
from PIL import Image, ImageDraw, ImageTk
from colour import Color
import tkinter as tk
from tkinter import ttk

class Graphic:
    def __init__(self, master, Gframe, Sframe):
        self.master=master
        self.boxColor, self.regColor = "red", "red"
        self.degChange = 0.
        self.scan_direction = "X"
        self.dirct = 2
        self.startxy_c_de = 0
        self.startxy_c = self.startxy_c_de
        self.fits_offset = (0, 0)
        self.Box = []
        self.box_selected,  self.over_object, self.over_selected, self.clicked, self.box_manip, self.box_resize, self.boxDrawn = \
            False, False, False, False, False, False, False

        Gframe.grid_rowconfigure(0, weight=1)
        Gframe.grid_columnconfigure(0, weight=1)
        self.canvas = tk.Canvas(Gframe, bg="gray", cursor="cross-hair", highlightthickness=0, yscrollincrement=1, xscrollincrement=1)
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
        self.slidercanvas_1 = tk.Canvas(Sframe, bg="gray", height=10, highlightthickness=0, highlightcolor="black", borderwidth=0, relief=tk.GROOVE)
        self.slidercanvas_2 = tk.Canvas(Sframe, bg="gray", width=200, height=10, highlightthickness=0, highlightcolor="black", borderwidth=0, relief=tk.GROOVE)
        self.slidercanvas_1.grid(row=0, column=0, rowspan=1, columnspan=1, sticky="ew")
        self.slidercanvas_2.grid(row=1, column=0, rowspan=1, columnspan=1, sticky="ew")
        self.slider_temp()
        Sframe.bind("<Configure>", self.slider_temp)

        self.canvas.bind("<Enter>", lambda event, mode="default": self.cursors(event, mode))
        self.canvas.bind("<Button-1>", self.B12_callback)
        self.canvas.bind("<B1-Motion>", self.B1M_callback)
        self.canvas.bind("<ButtonRelease-1>", self.B1R_callback)
        self.canvas.bind("<ButtonRelease-2>", self.B2R_callback)
        self.canvas.bind('<MouseWheel>', self.v_scroll)
        self.canvas.bind('<Shift-MouseWheel>', self.h_scroll)

        self.canvas.bind("<Command-MouseWheel>", lambda event: self.testzoom(event))
        self.canvas.bind("<KeyRelease-Meta_L>", lambda event: self.endzoom(event))
        self.canvas.bind("a", lambda event: self.testpos(event))
        self.canvas.bind("d", lambda event: self.deletetestpos(event))
        self.zoom = 1

    def h_scroll(self, event):
        self.canvas.xview_scroll(-3 * event.delta, 'units')

    def v_scroll(self, event):
        self.canvas.yview_scroll(-3 * event.delta, 'units')

    def testzoom(self, event):
        self.canvas.focus_set()
        px,py = event.widget.winfo_pointerxy()
        rx,ry = (event.widget.winfo_rootx(), event.widget.winfo_rooty())
        cx,cy = (px-rx, py-ry)
        self.eventx, self.eventy = self.canvas.canvasx(cx),self.canvas.canvasy(cy)
        print(self.canvas.canvasx(cx),self.canvas.canvasy(cy))

        scale = 1 - np.sign(event.delta)*0.02
        self.zoom *= scale

        imgx, imgy = self.fitsPIL.size
        print(imgx, imgy)
        finalx, finaly = imgx*self.zoom, imgy*self.zoom
        print(finalx, finaly)
        fitsPILtemp = self.fitsPIL.resize((int(round(finalx)), int(round(finaly))), Image.NEAREST)
        self.fitsTk = ImageTk.PhotoImage(fitsPILtemp)

        #fitspos = (self.canvas.canvasx(cx)*(1-scale), self.canvas.canvasy(cy)*(1-scale))
        self.canvas.delete("fits")
        self.fitsCanvas = self.canvas.create_image((0,0), image=self.fitsTk, anchor="nw", tag="fits")
        self.fits_CurSize = (self.fitsTk.width(), self.fitsTk.height())
        print(self.fits_CurSize)

        self.canvas.config(scrollregion=self.canvas.bbox("all"))

        #self.canvas.xview_moveto(-(self.cv_w-self.im_w)/self.im_w/2-(self.canvas.canvasx(cx)*self.zoom))
        #self.canvas.yview_moveto(-(self.cv_h-self.im_h)/self.im_h/2-(self.canvas.canvasy(cy)*self.zoom))

        self.canvas.xview_scroll(int(round(self.eventx*(scale-1))), 'units')
        self.canvas.yview_scroll(int(round(self.eventy*(scale-1))), 'units')

        self.canvas.tag_lower("fits")

        self.canvas.scale("box",0,0,scale,scale)

    def endzoom(self, *args):
        if self.box_selected:
            ori_selected_state = True
            ori_selected_box_id = self.box_id
        else:
            ori_selected_state = False
        self.box_selected = True
        for i in range(len(self.Box)):
            self.box_id = self.Box[i][0]
            self.setBox(None, generate_B4=False)
        if ori_selected_state:
            self.box_id = ori_selected_box_id
            self.canvas.event_generate("<B5-Motion>")
        else:
            self.box_selected = False


    def testpos(self, event):
        px,py = event.widget.winfo_pointerxy()
        rx,ry = (event.widget.winfo_rootx(), event.widget.winfo_rooty())
        cx,cy = (px-rx, py-ry)
        self.eventx, self.eventy = self.canvas.canvasx(cx),self.canvas.canvasy(cy)
        print(self.eventx, self.eventy)
        self.canvas.create_line(5,0,50,0, fill="orange", tag=("testline","box"))
        a = self.canvas.create_line(5,0.4,100,0.4, fill="pink", tag=("testline","box"))
        self.canvas.create_line(5,0.6,150,0.6, fill="cyan", tag=("testline","box"))
        self.canvas.create_line(5,1,200,1, fill="white", tag=("testline","box"))
        print(self.canvas.coords(a))


    def deletetestpos(self, *args):
        self.canvas.delete("testline")
        print(self.eventx*self.zoom, self.eventy*self.zoom)
        self.canvas.xview_moveto(-(self.cv_w-self.im_w)/self.im_w/2-(self.eventx*(self.zoom-1)))
        self.canvas.yview_moveto(-(self.cv_h-self.im_h)/self.im_h/2-(self.eventy*(self.zoom-1)))



    def slider_temp(self, *args):
        """Initialize the box color slider."""
        self.scBg = []
        self.colorSlider([self.slidercanvas_1, self.slidercanvas_2], "red", "pink")
        for num, s_set in enumerate(((self.slidercanvas_1, ["box","marker"]), (self.slidercanvas_2, ["reg","regmarker"]))):
            s_set[0].bind("<Button-1>",
                          lambda event, canvas=s_set[0], target=[s_set[1],num], pad=1: self.sliderNob(event, canvas, target, pad))
            s_set[0].bind("<B1-Motion>",
                          lambda event, canvas=s_set[0], target=[s_set[1],num], pad=1: self.sliderNob(event, canvas, target, pad))

    def colorSlider(self, canvas_list, color1, color2):
        """Draw the color gradient and slider nob."""
        canvas_list[0].update()
        self.sliderwidth, self.sliderheight = canvas_list[0].winfo_width(), canvas_list[0].winfo_height()
        self.colorGradient = list(Color(color1).range_to(Color(color2), self.sliderwidth))
        self.gradBg = Image.new("RGB", (self.sliderwidth, self.sliderheight), "#FFFFFF")
        self.gradBgDraw = ImageDraw.Draw(self.gradBg)
        for x, color in enumerate(self.colorGradient):
            self.gradBgDraw.line((x, 0, x, self.sliderheight), fill=str(color), width=1)
        for canvas in canvas_list:
            self.scBg.append([])
            self.scBg[-1].append(ImageTk.PhotoImage(self.gradBg))
            self.scBg[-1].append(canvas.create_image(0, 0, image=self.scBg[-1][0], anchor=tk.NW))
            self.scBg[-1].append(canvas.create_line(1, 0, 1, self.sliderheight, width=2, fill="#444444"))

    def sliderNob(self, event, canvas, target, pad):
        """Determine the color chosen based on the position of the nob."""
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
        """Scale and update the FITS image.

        Args:
            vartuple (tuple) : Contains the variable values for the scaling function to be use in color_func.
                               c.f. FileEdit.Tabs.slider_callback
        """
        self.fits_offset = (0, 0)
        self.canvas.focus_set()
        try:
            self.temp = self.color_func(*vartuple)
            self.temp2 = self.temp.copy()
            self.temp2 = np.clip(self.temp2, np.min(self.fitsNp_ori), np.max(self.fitsNp_ori))
            self.temp2 = (self.temp2 / np.max(self.fitsNp_ori)) * 255
            self.temp = Image.fromarray(np.flip(self.temp2, 0))
            self.fitsPIL = self.temp
            self.temp = ImageTk.PhotoImage(self.temp.resize(self.fits_CurSize,Image.NEAREST))
            self.canvas.delete("fits")
            self.canvas.create_image(self.fits_offset, image=self.temp, anchor="nw", tag="fits")
            self.canvas.tag_raise("O")
        except AttributeError:
            pass

    def color_func(self, mode, pixel, lowerb, upperb, pixel_min, pixel_max, gamma, gain, bias_x, bias_y):
        """Scale the input pixel."""
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

    def setBox(self, _, generate_B4=True, **kwargs):
        """Initializes box manipulation functions.

        Kwargs:
            "REG" : Separate initialization for pyregion boxes
        """
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
                self.startxy_c = 0
        if True:
            try:
                self.id_str = str(self.box_id)
                self.canvas.delete("mid"+self.id_str, "resize"+self.id_str,"marker"+self.id_str,"regmarker"+self.id_str)

                box_coords = self.canvas.coords(self.box_id)
                (self.NWPos, self.NEPos, self.SEPos, self.SWPos) = tuple(box_coords[i:i + 2] for i in range(0, 8, 2))
                self.UPos = ((self.NWPos[0] + self.NEPos[0]) / 2, (self.NWPos[1] + self.NEPos[1]) / 2)
                self.DPos = ((self.SEPos[0] + self.SWPos[0]) / 2, (self.SEPos[1] + self.SWPos[1]) / 2)
                self.LPos = ((self.NWPos[0] + self.SWPos[0]) / 2, (self.NWPos[1] + self.SWPos[1]) / 2)
                self.RPos = ((self.SEPos[0] + self.NEPos[0]) / 2, (self.SEPos[1] + self.NEPos[1]) / 2)
                self.midPos = ((self.LPos[0] + self.RPos[0]) / 2, (self.UPos[1] + self.DPos[1]) / 2)
                try:
                    self.degNow = kwargs["REGdeg"]*math.pi/180
                except KeyError:
                    self.degNow = -(((0.5 * math.pi) - cmath.phase(complex(self.UPos[0] - self.midPos[0], self.midPos[1] - self.UPos[1])) - (2*math.pi)) % (-2*math.pi))

                self.Box[self.box_index] = [self.box_id]
                self.Box[self.box_index].append(
                    self.canvas.create_oval(box_coords[self.startxy_c] - 2, box_coords[self.startxy_c+1] - 2,
                                            box_coords[self.startxy_c] + 2, box_coords[self.startxy_c+1] + 2,
                                            width=1, fill=self.boxColor, outline=self.boxColor, tag=("O", "box", "marker", "marker"+self.id_str)))
                small_circle_index_temp = (self.startxy_c+self.dirct)%8
                self.Box[self.box_index].append(
                    self.canvas.create_oval(box_coords[small_circle_index_temp] - 1, box_coords[small_circle_index_temp+1] - 1,
                                            box_coords[small_circle_index_temp] + 1, box_coords[small_circle_index_temp+1] + 1,
                                            width=1, fill=self.boxColor, outline=self.boxColor, tag=("O", "box", "marker", "marker"+self.id_str)))

                if "reg" in self.canvas.gettags(self.box_id):
                    self.canvas.itemconfig(self.Box[self.box_index][1], fill=self.regColor, outline=self.regColor, tag=("O", "reg", "regmarker", "regmarker"+self.id_str))
                    self.canvas.itemconfig(self.Box[self.box_index][2], fill=self.regColor, outline=self.regColor, tag=("O", "reg", "regmarker", "regmarker"+self.id_str))
                if self.box_selected:
                    self.canvas.itemconfig(self.Box[self.box_index][1], fill="")
                    self.canvas.itemconfig(self.Box[self.box_index][2], fill="")

                self.resizeFill = ""
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.NWPos[0] - 6, self.NWPos[1] + 5, self.NWPos[0] + 5, self.NWPos[1] - 6,
                                        width=0, fill="", tag=("O", "resize", "resize" + self.id_str, "C", "NW")))
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.NEPos[0] - 6, self.NEPos[1] + 5, self.NEPos[0] + 5, self.NEPos[1] - 6,
                                        width=0, fill="", tag=("O", "resize", "resize" + self.id_str, "C", "NE")))
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.SEPos[0] - 6, self.SEPos[1] + 5, self.SEPos[0] + 5, self.SEPos[1] - 6,
                                        width=0, fill="", tag=("O", "resize", "resize" + self.id_str, "C", "SE")))
                self.Box[self.box_index].append(
                    self.canvas.create_oval(self.SWPos[0] - 6, self.SWPos[1] + 5, self.SWPos[0] + 5, self.SWPos[1] - 6,
                                        width=0, fill="", tag=("O", "resize", "resize" + self.id_str, "C", "SW")))
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
                self.Box[self.box_index].append([self.midPos, self.degNow, self.startxy_c, self.dirct])

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
                self.canvas.tag_bind("NW", "<Double-Button-1>", lambda event, corner=0, manual=True: self.set_onpos(event, corner, manual))
                self.canvas.tag_bind("NE", "<B1-Motion>", lambda event, mode=("NE", "free"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("NE", "<Shift-B1-Motion>", lambda event, mode=("NE", "ratio"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("NE", "<Double-Button-1>", lambda event, corner=1, manual=True: self.set_onpos(event, corner, manual))
                self.canvas.tag_bind("SE", "<B1-Motion>", lambda event, mode=("SE", "free"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("SE", "<Shift-B1-Motion>", lambda event, mode=("SE", "ratio"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("SE", "<Double-Button-1>", lambda event, corner=2, manual=True: self.set_onpos(event, corner, manual))
                self.canvas.tag_bind("SW", "<B1-Motion>", lambda event, mode=("SW", "free"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("SW", "<Shift-B1-Motion>", lambda event, mode=("SW", "ratio"): self.manipulateBox(event, mode))
                self.canvas.tag_bind("SW", "<Double-Button-1>", lambda event, corner=3, manual=True: self.set_onpos(event, corner, manual))

                self.canvas.tag_bind("all", "<ButtonRelease-1>", lambda event: self.selectBox(event))

                #
                #self.canvas.tag_raise(self.Box[self.box_index][1])
                #
                self.box_manip, self.box_resize = False, False
                self.box_index = self.box_index_selected
                self.box_id = self.Box[self.box_index][0]
                self.degChange = 0.
                if generate_B4:
                    self.canvas.event_generate("<B4-Motion>")
                self.canvas.event_generate("<Enter>")
            except (AttributeError, IndexError):
                pass

    def set_onpos(self, _, corner, manual=False):
        if self.box_selected:
            box_coord = tuple(self.canvas.coords(self.box_id)[i:i + 2] for i in range(0, 8, 2))
            if self.scan_direction == "X":
                if corner == 0 or corner == 2:
                    self.dirct = 2
                else:
                    self.dirct = -2
            else:
                if corner == 1 or corner == 3:
                    self.dirct = 2
                else:
                    self.dirct = -2
            (temp1x, temp1y) = box_coord[corner]
            self.canvas.coords(self.Box[self.box_index][1], temp1x - 2, temp1y - 4, temp1x + 20, temp1y + 4)
            self.Box[self.box_index][-1][2] = corner*2
            self.startxy_c = corner*2
        if manual:
            self.setBox(None)

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
        """Detect cursor status."""
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
        """Assign special 'selected' state to a clicked box."""
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
                self.startxy_c = self.Box[self.box_index][-1][2]
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
        """Delete 'selected' box."""
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
                self.startxy_c, self.scan_direction, self.dirct = 0, "X", 2
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
        """Calculate and redraw box when interactively manipulated."""
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
                small_circle_index_temp = (self.startxy_c+self.dirct)%8
                self.canvas.coords(self.Box[self.box_index][1], self.finalPos[self.startxy_c]-2, self.finalPos[self.startxy_c+1]-2,
                                                                self.finalPos[self.startxy_c]+2, self.finalPos[self.startxy_c+1]+2)
                self.canvas.coords(self.Box[self.box_index][2], self.finalPos[small_circle_index_temp]-1, self.finalPos[small_circle_index_temp+1]-1,
                                                                self.finalPos[small_circle_index_temp]+1, self.finalPos[small_circle_index_temp+1]+1)
        except AttributeError:
            pass

    def manipulateVar(self, event, mode1, ref1, ref2, o, a1, a2, a3, a4, phi11, phi12, phi21, phi22, phi31, phi32, phi41, phi42, r):
        """Calcaulte final position of box vertices."""
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
