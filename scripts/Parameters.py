from n_const import REST_FREQ

class Parameters:
    """List of parameters to be parsed and have widgets created for in FileEdit.Files.setup_obs_frame().

    Notes
    -----
    Each top-level element is a tk.LabelFrame widget, for which the first element is to be it's Label (str) and the
    second element is to be a single list of dictionaries for the observational parameters and their respective
    properties.

    Dict Keys:
        "name" : str
            Name of the parameter
        "#" : str
            Description of the parameter
        "string" : bool
            Status as string literal in the output OBS
        "no." : int
            Position within each section in the output OBS
        "grid" : tuple
            .grid coordinates used when the parameters are inserted into the GUI interface.
            Reccomended to use for either all or none of the parameters in a section.
        "default" :
            Default entry
        "values" : list
            Options for tk.OptionMenu dropdown menu widget. First option is default.
        "unit" : str
            Unit suffix in the output OBS
    """

    def __init__(self):
        self.paramlist = []
        self.paramlist.append(["observation_property", [
            {"name": "observer", "#": "Name of observer", "string": True, "no.": 1},
            {"name": "target", "#": "Name of target object", "string": True, "no.": 2},
        ]])
        self.paramlist.append(["coordinate", [
            {"name": "lambda_on", "#": "X-coordinate of map center", "string": True, "grid": (1,0), "no.": 1},
            {"name": "beta_on", "#": "Y-coordinate of map center", "string": True, "grid": (2,0), "no.": 2},
            {"name": "relative", "#": "If true, OFF position is given in offset from ON point", "string": False, "grid": (4, 0), "values": ["true", "false"], "default":"true", "no.": 3},
            {"name": "lambda_off", "#": "OFF point x-coordinate (degree)", "string": True, "grid": (5,0), "no.": 4},
            {"name": "beta_off", "#": "OFF point y-coordinate (degree)", "string": True, "grid": (6,0), "no.": 5},
            {"name": "delta_lambda", "#": "offset_off", "string": True, "grid": (5,1), "no.": 6},
            {"name": "delta_beta", "#": "offset_off", "string": True, "grid": (6,1), "no.": 7},
            {"name": "position_angle", "#": "Position angle of the scan direction (default = 0deg)", "string": True, "unit": "deg", "grid": (3,0), "no.": 8},
            {"name": "start_position_x", "#": "X start position relative to lambda_on (arcsec) (sw)", "string": True, "unit": "arcsec", "grid": (1,1), "no.": 9},
            {"name": "start_position_y", "#": "Y start position relative to lambda_on (arcsec) (sw)", "string": True, "unit": "arcsec", "grid": (2,1), "no.": 10},
            {"name": "coord_sys", "#": "(J2000/b1950/galactic/horizontal)", "string": True, "values": ["J2000", "B1950", "Galactic"], "grid": (0,0), "no.": 11},
        ]])
        self.paramlist.append(["scan_property", [
            {"name": "scan_direction", "#": "X (=0) or Y (=1) scan", "values": ["X", "Y"], "string": True, "grid": (0,0), "no.": 1},
            {"name": "n", "#": "Number of scanlines (integer)", "string": False, "grid": (1,0), "no.": 2},
            {"name": "scan_spacing", "#": "Scan spacing (default = 60arcsec)", "string": True, "unit": "arcsec", "default": 60, "grid": (2,0), "no.": 3},
            {"name": "scan_length", "#": "time duration for 1 scan  (s) ", "string": True, "unit": "s", "grid": (1,1), "no.": 4},
            {"name": "scan_velocity", "#": "scan velocity (default = 600arcsec/s) (<=600)", "string": True, "unit": "arcsec/s", "default": 600, "grid": (2,1), "no.": 5},
            {"name": "ramp_pixel", "#": "approach length (pixel) (default = 4)", "string": False, "grid": (3,0), "no.": 6},
            {"name": "integ_on", "#": "ON point integratoin time (s)", "string": True, "unit": "s", "grid": (4,0), "no.": 7}
        ]])
        self.paramlist.append(["calibration", [
            {"name": "integ_off", "#": "OFF point integration time (s)", "string": True, "unit": "s", "no.": 1},
            {"name": "integ_hot", "#": "hot-load integration time (s)", "string": True, "unit": "s", "no.": 2},
            {"name": "off_interval", "#": "off interval (number of scan lines: default = 1)", "string": False, "default": 1, "no.": 3},
            {"name": "load_interval", "#": "time interval of hot-load measurements (min)", "string": True, "unit": "min", "no.": 4}
        ]])

