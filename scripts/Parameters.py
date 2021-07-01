class Parameters:
    """List of parameters.

    Each top-level element is a tk.LabelFrame widget, for which the first element is to be it's Label (str) and the
    second element is to be a single list of dictionaries for the observational parameters and their respective
    properties.

    Dict Keys:
        "name" (str) : Name of the parameter
        "#" (str) : Description of the parameter
        "string" (bool) : Status as string literal in the output OBS
        "no." (int) : Position within each section in the output OBS
        "grid" (tuple) : .grid coordinates used when the parameters are inserted into the GUI interface
                         Recommended to be used for either all or none of the parameters in a section
        "default" : Default entry
        "values" (list) : Options for tk.OptionMenu dropdown menu widget
                          First option is default
        "unit" (str) : Unit suffix in the output OBS
    """

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
            {"name": "scan_direction", "#": "X (=0) or Y (=1) scan", "values": ["X", "Y"], "string": True, "grid": (0,0), "no.": 1},
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

