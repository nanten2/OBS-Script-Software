import math, cmath
import numpy as np

import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageDraw, ImageTk
from colour import Color

import Gvars


class Graphic:
    """Handles the displaying and scaling of the FITS image, as well as the drawing, selection, and manipulation of the
    boxes.

    Attributes
    ----------
    Box : list
        list containing all Tkinter.
    boxColor : str
    boxDrawn : bool
    box_id : int
    box_index : int
    box_index_hover : int
    box_index_selected : int
    box_manip : bool
    box_resize : bool
    box_selected : bool
    canvas : tkinter.Canvas
        Main canvas for the FITS image boxes
    clearButton : tkinter.ttk.Button
        Button to delete selected box
    clicked : bool
    colorGradient : list
    degChange : float
        Degrees changed while rotating, before setBox is called
    degNow : float
    delta_coeff : int
        Modifier for event.delta
    dirct : int
        Index difference from the start position corner to the end corner of the scan line
    endPos : tuple
    fitsCanvas : int
        Canvas image object
    fitsPIL : PIL.Image.Image
    fitsTk : PIL.ImageTk.PhotoImage
    fits_CurSize : tuple
    fits_region : int
        Bounding box Canvas rectangle object for the FITS image
    hbar : tkinter.ttk.Scrollbar
    manipBox_initvars : tuple
    master : tkinter.Toplevel
    over_object : bool
    over_selected : bool
    regColor : str
    resizeFill : str
    scBg : list
        List of the slider objects and the gradient image
    scan_direction : str
    slidercanvas_1 : tkinter.Canvas
    slidercanvas_2 : tkinter.Canvas
    sliderheight : int
    sliderwidth : int
    startPos_draw : tuple
    startxy_c : int
        Index of the start position corner
    vbar : tkinter.ttk.Scrollbar
    zoom : float

    Notes
    -----
    Events bound to FileEdit.Files.currentCoords_update:
    <B1-Motion> (left click drag) :
        To update parameter values while moving or resizing box
    <B2-Motion> (right click drag) :
        To update parameter values while rotating box
    <Shift-B4-Motion> :
        To execute FileEdit.Files.currentCoords_update
        This might cause issues on operating systems that use B4 for the mousewheel.
    """
    def __init__(self, master, Gframe, Sframe):
        self.master = master
        self.degChange = 0.
        self.boxColor, self.regColor = "red", "red"
        self.scan_direction = "X"
        self.dirct = 2
        self.startxy_c = 0
        self.Box = []
        self.box_selected,  self.over_object, self.over_selected, self.clicked, self.box_manip, self.box_resize, self.boxDrawn = \
            False, False, False, False, False, False, False
        self.zoom = 1
        self.delta_coeff_init()

        Gframe.grid_rowconfigure(0, weight=1)
        Gframe.grid_columnconfigure(0, weight=1)
        self.canvas = tk.Canvas(Gframe, bg="gray", highlightthickness=0, yscrollincrement=1, xscrollincrement=1)
        self.cursors(None, "default")
        self.canvas.grid(row=0, column=0, sticky="nsew")
        self.vbar = ttk.Scrollbar(Gframe, orient="vertical", command=self.canvas.yview)
        self.vbar.grid(row=0, column=1, rowspan=1, sticky="ns")
        self.hbar = ttk.Scrollbar(Gframe, orient="horizontal", command=self.canvas.xview)
        self.hbar.grid(row=1, column=0, columnspan=2, sticky="ew")

        Sframe.grid_rowconfigure([0,1], weight=0)
        Sframe.grid_columnconfigure([0], weight=1)
        self.clearButton = ttk.Button(Sframe, text="Clear", state=tk.DISABLED, takefocus=0, command=lambda _=None, mode="select": self.resetBox(_, mode))
        self.clearButton.grid(row=0, column=1, rowspan=2, sticky="nsew")
        self.slidercanvas_1 = tk.Canvas(Sframe, bg="gray", height=10, highlightthickness=0, highlightcolor="black", borderwidth=0, relief=tk.GROOVE)
        self.slidercanvas_1.grid(row=0, column=0, rowspan=1, columnspan=1, sticky="ew")
        self.slidercanvas_2 = tk.Canvas(Sframe, bg="gray", width=200, height=10, highlightthickness=0, highlightcolor="black", borderwidth=0, relief=tk.GROOVE)
        self.slidercanvas_2.grid(row=1, column=0, rowspan=1, columnspan=1, sticky="ew")
        self.slider_temp(None)
        Sframe.bind("<Configure>", self.slider_temp)

        self.draw_keybind()
        self.zoom_keybind()

    def draw_keybind(self):
        self.canvas.bind("<BackSpace>", lambda event, mode="select": self.resetBox(event, mode))
        self.canvas.bind("<Enter>", lambda event, mode="default": self.cursors(event, mode))
        self.canvas.bind("<Button-1>", self.B12_callback)
        self.canvas.bind("<B1-Motion>", self.B1M_callback)
        self.canvas.bind("<ButtonRelease-1>", self.B1R_callback)
        self.canvas.bind("<ButtonRelease-2>", self.B2R_callback)
        if Gvars.curOS == "Linux":
            self.canvas.bind('<Button-4>', lambda event: self.v_scroll(event, delta=120))
            self.canvas.bind('<Button-5>', lambda event: self.v_scroll(event, delta=-120))
            self.canvas.bind('<Shift-Button-4>', lambda event: self.h_scroll(event, delta=120))
            self.canvas.bind('<Shift-Button-5>', lambda event: self.h_scroll(event, delta=-120))
        else:
            self.canvas.bind('<MouseWheel>', self.v_scroll)
            self.canvas.bind('<Shift-MouseWheel>', self.h_scroll)

    def zoom_keybind(self):
        self.canvas.bind("<Configure>", lambda event, target=self.canvas: self.update_proxy(event, self.canvas))
        self.canvas.bind("<KeyRelease-Meta_L>", lambda event: self.endzoom(event))
        if Gvars.curOS == "Linux":
            self.canvas.bind("<Control-Button-4>", lambda event: self.fits_zoom(event, delta=120))
            self.canvas.bind("<Control-Button-5>", lambda event: self.fits_zoom(event, delta=-120))
        else:
            self.canvas.bind("<Command-MouseWheel>", lambda event: self.fits_zoom(event))

    def update_proxy(self, _, target):
        target.update()

    def B12_callback(self, event):
        self.canvas.focus_set()
        self.clicked = True
        if self.box_selected or self.boxDrawn:
            self.deselectBox(event, "B12")
        else:
            self.setStart_draw(event)

    def B12_leave(self, event):
        try:
            if self.boxDrawn:
                endPos_temp = (self.canvas.canvasx(event.x), self.canvas.canvasy(event.y))
                boxPos_temp = self.canvas.coords(self.Box[self.box_index_hover][0])
                if max(boxPos_temp) in endPos_temp or min(boxPos_temp) in endPos_temp:
                    self.canvas.event_generate("<Leave>")
        except AttributeError:
            pass

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

    def delta_coeff_init(self):
        if Gvars.curOS == "Darwin":
            self.delta_coeff = 1
        elif Gvars.curOS == "Windows" or Gvars.curOS == "Linux":
            self.delta_coeff = 120

    def h_scroll(self, event, delta=None):
        if delta is not None:
            event.delta = delta
        self.canvas.xview_scroll(int(-3 * event.delta / self.delta_coeff), 'units')
        self.fits_zoom(event, zoom=False)

    def v_scroll(self, event, delta=None):
        if delta is not None:
            event.delta = delta
        self.canvas.yview_scroll(int(-3 * event.delta / self.delta_coeff), 'units')
        self.fits_zoom(event, zoom=False)

    def fits_initialize(self, PILimage):
        """Create the necessary image attributes when a FITS file is loaded."""
        # Keep reference to the image
        self.fitsPIL = PILimage
        self.fitsTk = ImageTk.PhotoImage(self.fitsPIL)
        self.fitsCanvas = self.canvas.create_image(0, 0, image=self.fitsTk, anchor="nw", tag="fits")

        # Update current size of image
        self.fits_CurSize = (self.fitsTk.width(), self.fitsTk.height())

        # Create bounding box for the image used in determining its scaled size
        self.fits_region = self.canvas.create_rectangle(0, 0, self.fits_CurSize[0], self.fits_CurSize[1],
                                                        width=0, outline="yellow", fill="", tag="fregion")

        # Add binding for window resize
        self.canvas.bind("<Configure>", lambda event: self.fits_zoom(event), add="+")

    def fits_zoom(self, event, zoom=True, delta=None):
        """Handles zooming and scrolling of the zoomed image."""
        if zoom:
            # Obtain cursor position on the canvas
            # Due to a bug in MacOS the more direct self.canvas.canvasx(event.x) is not used
            self.canvas.focus_set()
            px, py = event.widget.winfo_pointerxy()
            rx, ry = (event.widget.winfo_rootx(), event.widget.winfo_rooty())
            cx, cy = (px-rx, py-ry)
            eventx, eventy = self.canvas.canvasx(cx), self.canvas.canvasy(cy)

            if delta is not None:
                event.delta = delta

            # Set current and overall scaling factor
            scale = 1 - event.delta * 0.01 / self.delta_coeff
            self.zoom *= scale
            if self.zoom < 0.1:
                self.zoom /= scale
                return

            self.canvas.delete("marker")
        else:
            scale = 1
            eventx, eventy = 0, 0

        # Determine bounding region of the zoomed image
        # Because of the integer rounding required for some methods, the bounding box is not entirely accurate
        # which causes some of the edge to be trimmed off when zoomed beyond the canvas size
        self.canvas.scale("fregion", eventx, eventy, scale, scale)
        fits_region_bbox = self.canvas.coords(self.fits_region)
        self.canvas.config(scrollregion=fits_region_bbox)

        # Determine display region of the zoomed tile
        display_region_bbox = fits_region_bbox.copy()
        if display_region_bbox[0] < self.canvas.canvasx(0):
            display_region_bbox[0] = self.canvas.canvasx(0)
        if display_region_bbox[1] < self.canvas.canvasy(0):
            display_region_bbox[1] = self.canvas.canvasy(0)
        if display_region_bbox[2] > self.canvas.canvasx(self.canvas.winfo_width()):
            display_region_bbox[2] = self.canvas.canvasx(self.canvas.winfo_width())
        if display_region_bbox[3] > self.canvas.canvasy(self.canvas.winfo_height()):
            display_region_bbox[3] = self.canvas.canvasy(self.canvas.winfo_height())

        # Determine cropping area of original image and execute crop
        crop_area = [max(int(round((display_region_bbox[0]-fits_region_bbox[0])/self.zoom)), 0),
                     max(int(round((display_region_bbox[1]-fits_region_bbox[1])/self.zoom)), 0),
                     min(int(round(self.fits_OriSize[0]-0-(fits_region_bbox[2]-display_region_bbox[2])/self.zoom)), self.fits_OriSize[0]-0),
                     min(int(round(self.fits_OriSize[1]-0-(fits_region_bbox[3]-display_region_bbox[3])/self.zoom)), self.fits_OriSize[1]-0)]
        fitsPIL_cropped = self.fitsPIL.crop(crop_area)

        # Resize cropped tile and redraw image
        final_size = (int(round((crop_area[2]-crop_area[0]+1)*self.zoom))-1,
                      int(round((crop_area[3]-crop_area[1]+1)*self.zoom))-1)
        fitsPIL_cropped_zoomed = fitsPIL_cropped.resize(final_size, Image.NEAREST)
        self.fitsTk = ImageTk.PhotoImage(fitsPIL_cropped_zoomed)
        self.canvas.delete("fits")
        self.fitsCanvas = self.canvas.create_image(-fits_region_bbox[0]+display_region_bbox[0],
                                                    -fits_region_bbox[1]+display_region_bbox[1],
                                                    image=self.fitsTk, anchor="nw", tag="fits")

        if zoom:
            # Adjust position to keep cursor on target while zooming
            self.canvas.xview_scroll(int(round(eventx*(scale-1))), 'units')
            self.canvas.yview_scroll(int(round(eventy*(scale-1))), 'units')

            # Match the bounding box's position to the image
            self.canvas.moveto("fregion", 0, 0)

            # Update image information
            self.fits_CurSize = (int(round(fits_region_bbox[2]-fits_region_bbox[0])),
                                 int(round(fits_region_bbox[3]-fits_region_bbox[1])))

            # Zoom canvas objects
            self.canvas.scale("box", 0, 0, scale, scale)

        self.canvas.tag_lower("fits")

    def endzoom(self, event):
        """Update fits tab after zooming to account for pixel rounding errors."""
        if self.box_selected:
            ori_selected_state = True
            ori_selected_box_id = self.box_id
        else:
            ori_selected_state = False
        self.box_selected = True
        for i in range(len(self.Box)):
            self.box_id = self.Box[i][0]
            self.setBox(None)
        if ori_selected_state:
            self.box_id = ori_selected_box_id
            self.canvas.event_generate("<Shift-B4-Motion>")
        else:
            self.box_selected = False

    def slider_temp(self, _):
        """Initialize the box color slider."""
        # Create list to hold references
        self.scBg = []

        # Create sliders
        self.colorSlider([self.slidercanvas_1, self.slidercanvas_2], "red", "pink")

        # Bind slider
        for num, s_set in enumerate(((self.slidercanvas_1, ["box","marker"]), (self.slidercanvas_2, ["reg","regmarker"]))):
            s_set[0].bind("<Button-1>",
                          lambda event, canvas=s_set[0], target=[s_set[1],num], pad=1: self.sliderNob(event, canvas, target, pad))
            s_set[0].bind("<B1-Motion>",
                          lambda event, canvas=s_set[0], target=[s_set[1],num], pad=1: self.sliderNob(event, canvas, target, pad))

    def colorSlider(self, canvas_list, color1, color2):
        """Draw the color gradient and slider nob."""
        # Get current slider dimensions
        canvas_list[0].update()
        self.sliderwidth, self.sliderheight = canvas_list[0].winfo_width(), canvas_list[0].winfo_height()

        # Create color gradient image
        self.colorGradient = list(Color(color1).range_to(Color(color2), self.sliderwidth))
        gradBg = Image.new("RGB", (self.sliderwidth, self.sliderheight), "#FFFFFF")
        gradBgDraw = ImageDraw.Draw(gradBg)
        for x, color in enumerate(self.colorGradient):
            gradBgDraw.line((x, 0, x, self.sliderheight), fill=str(color), width=1)

        # Create canvas images and keep references
        for canvas in canvas_list:
            self.scBg.append([])
            self.scBg[-1].append(ImageTk.PhotoImage(gradBg))
            self.scBg[-1].append(canvas.create_image(0, 0, image=self.scBg[-1][0], anchor=tk.NW))
            self.scBg[-1].append(canvas.create_line(1, 0, 1, self.sliderheight, width=2, fill="#444444"))

    def sliderNob(self, event, canvas, target, pad):
        """Determine the color chosen based on the position of the nob."""
        width, height = self.sliderwidth, self.sliderheight
        try:
            # Determine corresponding position of the nob on the color gradient
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

            # Apply color change
            self.canvas.itemconfig(target[0][0], outline=temp_color)
            self.canvas.itemconfig(target[0][1], fill=temp_color)

            # Hollow out markers if the box is in the "selected" state
            if self.box_selected:
                self.canvas.itemconfig(self.Box[self.box_index][1], fill="")
                self.canvas.itemconfig(self.Box[self.box_index][2], fill="")

            # For debugging with self.resizeFill
            if not self.resizeFill == "":
                self.canvas.itemconfig("resize", fill=temp_color)
        except AttributeError:
            pass

    def slider_master(self, vartuple):
        """Scale and update the FITS image.

        Parameters
        ----------
        vartuple : tuple
            Contains the variable values for the scaling function to be use in color_func.
            c.f. FileEdit.Tabs.slider_callback
        """
        self.canvas.focus_set()
        try:
            temp_bitmap = self.color_func(*vartuple)
            temp_bitmap_2 = temp_bitmap.copy()

            # Trim the values larger or smaller than the original's maximum and minimum
            temp_bitmap_2 = np.clip(temp_bitmap_2, np.min(self.fitsNp_ori), np.max(self.fitsNp_ori))

            # Convert to 8-bit and invert to have Y-axis pointing up
            temp_bitmap_2 = (temp_bitmap_2 / np.max(self.fitsNp_ori)) * 255
            temp_bitmap = Image.fromarray(np.flip(temp_bitmap_2, 0))

            self.fitsPIL = temp_bitmap
            self.fits_zoom(None, zoom=False)
        except ValueError:
            pass

    def color_func(self, pixel, pixel_min, pixel_max, gamma, gain, bias_x, bias_y, lowerb, upperb, mode):
        """Scale the input pixel by evaluating a piecewise scaling function.

        Returns
        -------
        numpy.ndarray
            Scaled pixel value
        """
        f1 = pixel_min
        f4 = pixel_max
        if mode == "　Symmetric":
            bias_sep_x = lowerb + (bias_x * (upperb - lowerb))
            bias_sep_y = pixel_min + (bias_y * (pixel_max - pixel_min))
            bound_diff = upperb - lowerb
            global_diff = pixel_max - pixel_min
            f2 = lambda pixel: -((-(pixel-bias_sep_x)/(0.5*bound_diff))**gamma)*(0.5*global_diff)*gain + bias_sep_y
            f3 = lambda pixel: (((pixel-bias_sep_x)/(0.5*bound_diff))**gamma)*(0.5*global_diff)*gain + bias_sep_y
            return np.piecewise(pixel, [(pixel<lowerb), (lowerb<=pixel)*(pixel<bias_sep_x),
                                        (bias_sep_x<=pixel)*(pixel<=upperb),(upperb<pixel)], [f1, f2, f3, f4])
        elif mode == "　Regular　":
            bias_sep_x = lowerb + ((bias_x-0.5) * (upperb - lowerb))
            bias_sep_y = pixel_min + ((bias_y-0.5) * (pixel_max - pixel_min))
            bound_diff = upperb - lowerb
            global_diff = pixel_max - pixel_min
            print(gamma)
            f23 = lambda pixel: (((pixel-bias_sep_x)/bound_diff)**gamma)*global_diff*gain + bias_sep_y
            return np.piecewise(pixel, [(pixel<lowerb), (lowerb <= pixel) * (pixel <= upperb), (upperb < pixel)], [f1, f23, f4])

    def setStart_draw(self, event):
        """Set the starting click position."""
        self.startPos_draw = (self.canvas.canvasx(event.x), self.canvas.canvasy(event.y))

    def drawBox(self, event):
        """Draw the box specified by the starting click position and the cursors current position."""
        self.box_manip = True
        self.canvas.delete("tempbox")
        endPos = (self.canvas.canvasx(event.x), self.canvas.canvasy(event.y))
        NEPos = (self.canvas.canvasx(event.x), self.startPos_draw[1])
        SWPos = (self.startPos_draw[0], self.canvas.canvasy(event.y))
        self.canvas.create_polygon(self.startPos_draw, NEPos, endPos, SWPos,
                                   fill="", width=1, outline=self.boxColor, tag="tempbox")

    def setBox(self, _, **kwargs):
        """Records references to box objects and initializes box manipulation functions.

        Parameters
        ----------
        _ : any
            Used to accept the Tkinter.Event argument passed by bind
        **kwargs: dict
            "REG" : tuple
                8-tuple of the pyregion box's polygon's corners
            "REGdeg" : float
                Position angle of the pyregion box

        Notes
        -----
        This method is generally called every time the box is physically changed.
        """
        # Create box with the appropriate tags
        try:
            self.Box.append([self.canvas.create_polygon(kwargs["REG"],
                                                        fill="", width=1, outline=self.regColor, tag="newbox")])
            self.box_id = self.canvas.find_withtag("newbox")[0]
            self.box_index = -1
            self.canvas.itemconfig("newbox", tag=("O", "reg"))
        except KeyError:
            if not self.box_selected and not self.boxDrawn:
                # Redraw tempbox from drawBox
                tempbox_coords = self.canvas.coords("tempbox")
                self.canvas.delete("tempbox")
                self.Box.append([self.canvas.create_polygon(tempbox_coords,
                                                            fill="", width=1, outline=self.boxColor, tag="newbox")])
                self.box_id = self.canvas.find_withtag("newbox")[0]
                self.box_index = -1
                self.canvas.itemconfig("newbox", tag=("O", "box"))
                self.over_selected = False
                self.boxDrawn = True
                try :
                    self.startxy_c = kwargs["startxy_c"]
                except KeyError:
                    self.startxy_c = 0
        try:
            id_str = str(self.box_id)

            # Delete peripheral items from the old box
            self.canvas.delete("mid"+id_str, "resize"+id_str, "marker"+id_str, "regmarker"+id_str)

            # Calculate relevant quantities in advance
            box_coords = self.canvas.coords(self.box_id)
            (NWPos, NEPos, SEPos, SWPos) = tuple(box_coords[i:i + 2] for i in range(0, 8, 2))
            UPos = ((NWPos[0] + NEPos[0]) / 2, (NWPos[1] + NEPos[1]) / 2)
            DPos = ((SEPos[0] + SWPos[0]) / 2, (SEPos[1] + SWPos[1]) / 2)
            LPos = ((NWPos[0] + SWPos[0]) / 2, (NWPos[1] + SWPos[1]) / 2)
            RPos = ((SEPos[0] + NEPos[0]) / 2, (SEPos[1] + NEPos[1]) / 2)
            midPos = ((LPos[0] + RPos[0]) / 2, (UPos[1] + DPos[1]) / 2)
            try:
                self.degNow = kwargs["REGdeg"]*math.pi/180
            except KeyError:
                self.degNow = -(((0.5 * math.pi) - cmath.phase(complex(UPos[0] - midPos[0], midPos[1] - UPos[1])) - (2*math.pi)) % (-2*math.pi))

            # Save references to self.Box
            self.Box[self.box_index] = [self.box_id]

            # Scan direction markers
            self.Box[self.box_index].append(
                self.canvas.create_oval(box_coords[self.startxy_c] - 2, box_coords[self.startxy_c+1] - 2,
                                        box_coords[self.startxy_c] + 2, box_coords[self.startxy_c+1] + 2,
                                        width=1, fill=self.boxColor, outline=self.boxColor, tag=("O", "box", "marker", "marker"+id_str)))
            small_circle_index_temp = (self.startxy_c+self.dirct) % 8
            self.Box[self.box_index].append(
                self.canvas.create_oval(box_coords[small_circle_index_temp] - 1, box_coords[small_circle_index_temp+1] - 1,
                                        box_coords[small_circle_index_temp] + 1, box_coords[small_circle_index_temp+1] + 1,
                                        width=1, fill=self.boxColor, outline=self.boxColor, tag=("O", "box", "marker", "marker"+id_str)))

            # Maintain appropriate tags to differentiate pyregion boxes
            if "reg" in self.canvas.gettags(self.box_id):
                self.canvas.itemconfig(self.Box[self.box_index][1], fill=self.regColor, outline=self.regColor, tag=("O", "reg", "regmarker", "regmarker"+id_str))
                self.canvas.itemconfig(self.Box[self.box_index][2], fill=self.regColor, outline=self.regColor, tag=("O", "reg", "regmarker", "regmarker"+id_str))
            if self.box_selected:
                self.canvas.itemconfig(self.Box[self.box_index][1], fill="")
                self.canvas.itemconfig(self.Box[self.box_index][2], fill="")

            # Items on the perimeter for box manipulation
            self.resizeFill = ""
            self.Box[self.box_index].append(
                self.canvas.create_oval(NWPos[0] - 6, NWPos[1] + 5, NWPos[0] + 5, NWPos[1] - 6,
                                    width=0, fill="", tag=("O", "resize", "resize" + id_str, "C", "NW")))
            self.Box[self.box_index].append(
                self.canvas.create_oval(NEPos[0] - 6, NEPos[1] + 5, NEPos[0] + 5, NEPos[1] - 6,
                                    width=0, fill="", tag=("O", "resize", "resize" + id_str, "C", "NE")))
            self.Box[self.box_index].append(
                self.canvas.create_oval(SEPos[0] - 6, SEPos[1] + 5, SEPos[0] + 5, SEPos[1] - 6,
                                    width=0, fill="", tag=("O", "resize", "resize" + id_str, "C", "SE")))
            self.Box[self.box_index].append(
                self.canvas.create_oval(SWPos[0] - 6, SWPos[1] + 5, SWPos[0] + 5, SWPos[1] - 6,
                                    width=0, fill="", tag=("O", "resize", "resize" + id_str, "C", "SW")))
            self.Box[self.box_index].append(
                self.canvas.create_oval(UPos[0] - 6, UPos[1] + 5, UPos[0] + 5, UPos[1] - 6,
                                    width=0, fill=self.resizeFill, tag=("O", "resize", "resize" + id_str, "UD", "U")))
            self.Box[self.box_index].append(
                self.canvas.create_oval(DPos[0] - 6, DPos[1] + 5, DPos[0] + 5, DPos[1] - 6,
                                    width=0, fill=self.resizeFill, tag=("O", "resize", "resize" + id_str, "UD", "D")))
            self.Box[self.box_index].append(
                self.canvas.create_oval(LPos[0] - 6, LPos[1] + 5, LPos[0] + 5, LPos[1] - 6,
                                    width=0, fill=self.resizeFill, tag=("O", "resize", "resize" + id_str, "LR", "L")))
            self.Box[self.box_index].append(
                self.canvas.create_oval(RPos[0] - 6, RPos[1] + 5, RPos[0] + 5, RPos[1] - 6,
                                    width=0, fill=self.resizeFill, tag=("O", "resize", "resize" + id_str, "LR", "R")))
            self.Box[self.box_index].append(
                self.canvas.create_oval(midPos[0] - 6, midPos[1] + 5, midPos[0] + 5, midPos[1] - 6,
                                    width=0, fill=self.resizeFill, tag=("O", "mid", "mid" + id_str)))
            self.Box[self.box_index].append([midPos, self.degNow, self.startxy_c, self.dirct])

            # Bind the canvas items to their respective functions
            self.canvas.tag_bind("all", "<Enter>", lambda event, mode="Enter": self.hover_detect(event, mode))
            self.canvas.tag_bind("all", "<Leave>", lambda event, mode="Leave": self.hover_detect(event, mode))
            self.canvas.tag_bind("all", "<Button-1>", lambda event: self.manipBox_callback(event))
            self.canvas.tag_bind("C", "<Button-2>", lambda event: self.manipBox_callback(event))

            # self.canvas.tag_bind("C", "<B2-Motion>", lambda event, mode=("rotate", None): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("C", "<Leave>", lambda event, mode="default": self.cursors(event, mode))
            self.canvas.tag_bind("all", "<Leave>", lambda event, mode="default": self.cursors(event, mode))

            self.canvas.tag_bind("box", "<Enter>", lambda event, mode="move": self.cursors(event, mode))
            self.canvas.tag_bind("box", "<B1-Motion>", lambda event, mode=("move", None): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("reg", "<Enter>", lambda event, mode="move": self.cursors(event, mode))
            self.canvas.tag_bind("reg", "<B1-Motion>", lambda event, mode=("move", None): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("mid", "<Enter>", lambda event, mode="move": self.cursors(event, mode))
            self.canvas.tag_bind("mid", "<B1-Motion>", lambda event, mode=("move", None): self.manipBox_callback(event, mode=mode))

            self.canvas.tag_bind("UD", "<Enter>", lambda event, mode="hand": self.cursors(event, mode))
            self.canvas.tag_bind("U", "<B1-Motion>", lambda event, mode=("U", "stretch"): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("D", "<B1-Motion>", lambda event, mode=("D", "stretch"): self.manipBox_callback(event, mode=mode))

            self.canvas.tag_bind("LR", "<Enter>", lambda event, mode="hand": self.cursors(event, mode))
            self.canvas.tag_bind("L", "<B1-Motion>", lambda event, mode=("L", "stretch"): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("R", "<B1-Motion>", lambda event, mode=("R", "stretch"): self.manipBox_callback(event, mode=mode))

            self.canvas.tag_bind("C", "<Enter>", lambda event, mode="hand": self.cursors(event, mode))
            self.canvas.tag_bind("C", "<B2-Motion>", lambda event, mode=("rotate", None): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("C", "<B3-Motion>", lambda event, mode=("rotate", None): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("NW", "<B1-Motion>", lambda event, mode=("NW", "free"): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("NW", "<Shift-B1-Motion>", lambda event, mode=("NW", "ratio"): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("NW", "<Double-Button-1>", lambda event, corner=0, manual=True: self.set_onpos(event, corner, manual))
            self.canvas.tag_bind("NE", "<B1-Motion>", lambda event, mode=("NE", "free"): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("NE", "<Shift-B1-Motion>", lambda event, mode=("NE", "ratio"): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("NE", "<Double-Button-1>", lambda event, corner=1, manual=True: self.set_onpos(event, corner, manual))
            self.canvas.tag_bind("SE", "<B1-Motion>", lambda event, mode=("SE", "free"): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("SE", "<Shift-B1-Motion>", lambda event, mode=("SE", "ratio"): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("SE", "<Double-Button-1>", lambda event, corner=2, manual=True: self.set_onpos(event, corner, manual))
            self.canvas.tag_bind("SW", "<B1-Motion>", lambda event, mode=("SW", "free"): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("SW", "<Shift-B1-Motion>", lambda event, mode=("SW", "ratio"): self.manipBox_callback(event, mode=mode))
            self.canvas.tag_bind("SW", "<Double-Button-1>", lambda event, corner=3, manual=True: self.set_onpos(event, corner, manual))

            self.canvas.tag_bind("all", "<ButtonRelease-1>", lambda event: self.selectBox(event))

            # Final miscellaneous updates
            self.degChange = 0.
            self.box_manip, self.box_resize = False, False
            self.box_index = self.box_index_selected
            self.box_id = self.Box[self.box_index][0]
            self.canvas.event_generate("<Shift-B4-Motion>")
            self.canvas.event_generate("<Enter>")
        except (AttributeError, IndexError):
            pass

    def set_onpos(self, event, corner, manual=False):
        """Insert markers at the on position and the ending corner of the first scan.

        Parameters
        ----------
        event : Tkinter.Event
        corner : int
            Takes values 1, 2, 3, or 4 for each corner in order
        manual : bool
            True if changing interactively
        """
        if self.box_selected:
            box_coord = tuple(self.canvas.coords(self.box_id)[i:i + 2] for i in range(0, 8, 2))

            # Determine relative position of the turning corner
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

            # Relocate start_pos marker
            (temp1x, temp1y) = box_coord[corner]
            self.canvas.coords(self.Box[self.box_index][1], temp1x - 2, temp1y - 4, temp1x + 2, temp1y + 4)

            # Update attributes
            self.Box[self.box_index][-1][2] = corner*2
            self.startxy_c = corner*2

        if manual:
            self.setBox(None)

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

    def selectBox(self, event, simulate=False):
        """Assign special 'selected' state to a clicked box.

        Parameters
        ----------
        event : Tkinter.Event
        simulate : bool, optional
            True if simulating box selection

        Notes
        -----
        A 'selected' box would have a dotted perimeter and hollow markers
        """
        try:
            if not self.over_selected:
                # Check conditions and determine the box index
                if simulate:
                    self.box_index = len(self.Box) - 1
                    physical = False
                    self.over_object = True
                else:
                    self.box_index = [index for index, box in enumerate(self.Box) if self.canvas.find_withtag("current")[0] in box][0]
                    physical = True

                # Set selected box details
                self.box_index_selected = self.box_index
                self.box_id = self.Box[self.box_index][0]
                self.startxy_c = self.Box[self.box_index][-1][2]

                # Configure to be 'selected' state
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
                    self.canvas.event_generate("<Shift-B4-Motion>")
                    self.cursors(event, "move")
                    self.clearButton.config(state=tk.ACTIVE)
        except IndexError:
            self.deselectBox(event, "B12")
            if self.box_selected:
                self.deselectBox(event, "B1R")

    def deselectBox(self, event, mode):
        """Deselect the selected box on release of <Button-1>."""
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
                self.canvas.event_generate("<Shift-B4-Motion>")

    def resetBox(self, _, mode):
        """Delete the selected box."""
        self.canvas.focus_set()
        if mode == "select":
            if self.box_selected:
                self.box_selected = False
                self.canvas.event_generate("<Shift-B4-Motion>")
                if "reg" not in self.canvas.gettags(self.box_id):
                    self.boxDrawn = False
                self.over_selected = False
                self.clicked = False
                self.startxy_c, self.scan_direction, self.dirct = 0, "X", 2
                self.clearButton.config(state=tk.DISABLED)
                for i in range(len(self.Box[self.box_index])):
                    self.canvas.delete(self.Box[self.box_index][i])

    def manipBox_callback(self, event, **kwargs):
        """Callback for box manipulation.

        Parameters
        ----------
        event : Tkinter.event
        **kwargs : dict
            "mode" : tuple
                2-tuple of the activation point and manipulation type to be passed to manipulateBox
        """
        self.box_resize = True
        try:
            # Manipulate box
            self.manipulateBox(event, kwargs["mode"], *self.manipBox_initvars)
        except (KeyError, AttributeError):
            # Initiliaze manipulation
            if self.over_selected:
                try:
                    # Record cursor's starting position
                    startPos = (self.canvas.canvasx(event.x), self.canvas.canvasy(event.y))

                    # Calculate revelant quantities in advance
                    (NWPos, NEPos, SEPos, SWPos) = tuple(self.canvas.coords(self.box_id)[i:i + 2] for i in range(0, 8, 2))
                    UPos = ((NWPos[0] + NEPos[0]) / 2, (NWPos[1] + NEPos[1]) / 2)
                    DPos = ((SEPos[0] + SWPos[0]) / 2, (SEPos[1] + SWPos[1]) / 2)
                    LPos = ((NWPos[0] + SWPos[0]) / 2, (NWPos[1] + SWPos[1]) / 2)
                    RPos = ((SEPos[0] + NEPos[0]) / 2, (SEPos[1] + NEPos[1]) / 2)
                    boxPos = self.canvas.coords(self.box_id)
                    self.manipBox_initvars = (startPos, NWPos, NEPos, SEPos, SWPos, UPos, DPos, LPos, RPos, boxPos)
                except AttributeError:
                    pass

    def manipulateBox(self, event, mode, startpos, nwpos, nepos, sepos, swpos, upos, dpos, lpos, rpos, boxpos):
        """Calculate and redraw box when interactively manipulated."""
        try:
            self.canvas.tag_unbind("all", "<Leave>")
            finalPos = boxpos
            midPos = self.Box[self.box_index][-1][0]
            self.degNow = self.Box[self.box_index][-1][1]

            if self.over_selected:
                self.box_manip = True
                if mode[0] == "move":
                    posChange = [self.canvas.canvasx(event.x) - startpos[0],
                                 self.canvas.canvasy(event.y) - startpos[1]]
                    finalPos = (boxpos[0] + posChange[0], boxpos[1] + posChange[1],
                                boxpos[2] + posChange[0], boxpos[3] + posChange[1],
                                boxpos[4] + posChange[0], boxpos[5] + posChange[1],
                                boxpos[6] + posChange[0], boxpos[7] + posChange[1])
                elif mode[0] == "rotate":
                    self.degChange = -(
                        ((cmath.phase(complex(startpos[0] - midPos[0], midPos[1] - startpos[1]))
                          - cmath.phase(complex(self.canvas.canvasx(event.x) - midPos[0],
                                                midPos[1] - self.canvas.canvasy(event.y))))
                         % (2 * math.pi)) / (math.pi / 900)) * (math.pi / 900) + (2*math.pi)
                    finalPos = []
                    for i in range(0, 8, 2):
                        dummyPos = complex(boxpos[i] - midPos[0], midPos[1] - boxpos[i + 1])
                        dummyPos = dummyPos * cmath.exp(complex(0, self.degChange))
                        finalPos.append(dummyPos.real + midPos[0])
                        finalPos.append(midPos[1] - dummyPos.imag)
                    finalPos = tuple(finalPos)
                elif mode[0] == "U":
                    posChange = self.manipulateVar(event, mode[1], upos, dpos, 0.5, 1,-1,0,0, 0,1,0,0,0,0,0,0, 1)
                    finalPos = (boxpos[0] + posChange[0], boxpos[1] + posChange[1],
                                boxpos[2] + posChange[0], boxpos[3] + posChange[1],
                                boxpos[4], boxpos[5], boxpos[6], boxpos[7])
                elif mode[0] == "D":
                    posChange = self.manipulateVar(event, mode[1], dpos, upos, 0.5, 1,-1,0,0, 0,1,0,0,0,0,0,0, 1)
                    finalPos = (boxpos[0], boxpos[1], boxpos[2], boxpos[3],
                                boxpos[4] + posChange[0], boxpos[5] + posChange[1],
                                boxpos[6] + posChange[0], boxpos[7] + posChange[1])
                elif mode[0] == "L":
                    posChange = self.manipulateVar(event, mode[1], lpos, rpos, 1, -1,-1,0,0, 0,0,0,1,0,0,0,0, 1)
                    finalPos = (boxpos[0] + posChange[0], boxpos[1] + posChange[1], boxpos[2], boxpos[3],
                                boxpos[4], boxpos[5], boxpos[6] + posChange[0], boxpos[7] + posChange[1])
                elif mode[0] == "R":
                    posChange = self.manipulateVar(event, mode[1], rpos, lpos, 1, -1,-1,0,0, 0,0,0,1,0,0,0,0, 1)
                    finalPos = (boxpos[0], boxpos[1], boxpos[2] + posChange[0], boxpos[3] + posChange[1],
                                boxpos[4] + posChange[0], boxpos[5] + posChange[1], boxpos[6], boxpos[7])
                elif mode[0] == "NW":
                    posChange = self.manipulateVar(event, mode[1], sepos, nwpos, 0.5, 1,-1,1,1, 0,1,0,0,1,0,1,1, 0)
                    finalPos = (self.endPos[0], self.endPos[1], boxpos[2] + posChange[0], boxpos[3] + posChange[1],
                                boxpos[4], boxpos[5], boxpos[6] + posChange[2], boxpos[7] + posChange[3])
                elif mode[0] == "NE":
                    posChange = self.manipulateVar(event, mode[1], swpos, nepos, 0.5, 1,-1,1,1, 0,1,0,0,1,0,1,1, 0)
                    finalPos = (boxpos[0] + posChange[0], boxpos[1] + posChange[1], self.endPos[0], self.endPos[1],
                                boxpos[4] + posChange[2], boxpos[5] + posChange[3], boxpos[6], boxpos[7])
                elif mode[0] == "SW":
                    posChange = self.manipulateVar(event, mode[1], nepos, swpos, 0, 1,1,-1,1, 0,0,0,1,1,1,1,0, 0)
                    finalPos = (boxpos[0] + posChange[0], boxpos[1] + posChange[1], boxpos[2], boxpos[3],
                                boxpos[4] + posChange[2], boxpos[5] + posChange[3], self.endPos[0], self.endPos[1])
                elif mode[0] == "SE":
                    posChange = self.manipulateVar(event, mode[1], nwpos, sepos, 0, 1,1,-1,1, 0,0,0,1,1,1,1,0, 0)
                    finalPos = (boxpos[0], boxpos[1], boxpos[2] + posChange[0], boxpos[3] + posChange[1],
                                self.endPos[0], self.endPos[1], boxpos[6] + posChange[2], boxpos[7] + posChange[3])

                self.canvas.coords(self.box_id, finalPos)
                small_circle_index_temp = (self.startxy_c+self.dirct) % 8
                self.canvas.coords(self.Box[self.box_index][1],
                                   finalPos[self.startxy_c]-2, finalPos[self.startxy_c+1]-2,
                                   finalPos[self.startxy_c]+2, finalPos[self.startxy_c+1]+2)
                self.canvas.coords(self.Box[self.box_index][2],
                                   finalPos[small_circle_index_temp]-1, finalPos[small_circle_index_temp+1]-1,
                                   finalPos[small_circle_index_temp]+1, finalPos[small_circle_index_temp+1]+1)
        except AttributeError:
            pass

    def manipulateVar(self, event, mode1, ref1, ref2, o, a1, a2, a3, a4, phi11, phi12, phi21, phi22, phi31, phi32, phi41, phi42, r):
        """Calcaulte final position of box vertices.

        Returns
        -------
        list
            Coordiante change of the box's corners, in the order of the first position of the polygon, excluding the
            corner(s) that do not change and the activation corner if applicable.
        """
        # Offset final position if streatch box sides
        self.endPos = (self.canvas.canvasx(event.x), self.canvas.canvasy(event.y))
        if mode1 == "stretch":
            pos_offset = [ref1[0] - ref2[0], ref1[1] - ref2[1]]
        elif mode1 == "free":
            pos_offset = [0, 0]
        elif mode1 == "ratio":
            pos_offset = [0, 0]
            # Different end position to cursor if maintaining box dimensions
            try:
                diagVector = [ref2[0] - ref1[0], ref2[1] - ref1[1]]
                refEnd = [self.endPos[0] - ref2[0], self.endPos[1] - ref2[1]]
                diagDeg = cmath.phase(complex(diagVector[0], diagVector[1]))
                projDeg = diagDeg - cmath.phase(complex(refEnd[0], refEnd[1]))
                refEndNorm = np.linalg.norm(refEnd)
                self.endPos = (ref2[0] + refEndNorm * math.cos(projDeg) * math.cos(diagDeg),
                               ref2[1] + refEndNorm * math.cos(projDeg) * math.sin(diagDeg))
            except ZeroDivisionError:
                pass

        # Calculate the projection degree and norm of the end position relative to the reference points
        refEnd = [self.endPos[0] - ref2[0], ref2[1] - self.endPos[1]]
        refEndNorm = np.linalg.norm(refEnd)
        projDeg = (o * math.pi) - (cmath.phase(complex(refEnd[0], refEnd[1]))) + self.degNow

        # Calculate and return final change in coordinates of the relevant points
        pi_half = 0.5*math.pi
        posChange = \
            [a1 * (refEndNorm * math.cos(projDeg - phi11*pi_half) * math.cos(-self.degNow - phi12*pi_half)) - (r * pos_offset[0]),
             a2 * (refEndNorm * math.cos(projDeg - phi21*pi_half) * math.cos(-self.degNow - phi22*pi_half)) - (r * pos_offset[1]),
             a3 * (refEndNorm * math.cos(projDeg - phi31*pi_half) * math.cos(-self.degNow - phi32*pi_half)),
             a4 * (refEndNorm * math.cos(projDeg - phi41*pi_half) * math.cos(-self.degNow - phi42*pi_half))]

        return posChange

    def cursors(self, _, mode):
        if self.over_selected:
            if mode == "default":
                if Gvars.curOS == "Darwin":
                    cursor = "cross-hair"
                else:
                    cursor = "crosshair"
            elif mode == "hand":
                if Gvars.curOS == "Darwin":
                    cursor = "openhand"
                else:
                    cursor = "hand1"
            elif mode == "move":
                if Gvars.curOS == "Darwin":
                    cursor = "fleur"
                else:
                    cursor = "fleur"
            elif mode == "rotate":
                if Gvars.curOS == "Darwin":
                    cursor = "exchange"
                else:
                    cursor = "exchange"
            self.canvas.config(cursor=cursor)
        else:
            if Gvars.curOS == "Darwin":
                cursor = "cross-hair"
            else:
                cursor = "crosshair"
        self.canvas.config(cursor=cursor)