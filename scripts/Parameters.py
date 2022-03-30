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
            {"name": "OBSERVER", "#": "name of observer", "string": True, "no.": 1},
            {"name": "OBJECT", "#": "name of target object", "string": True, "no.": 2},
            {"name": "MOLECULE_1", "#": "line identifier of target molecule", "string": True, "no.": 3, "values": REST_FREQ.keys()}
        ]])
        self.paramlist.append(["coordinate", [
            {"name": "LambdaOn", "#": "x-coordinate of map center (degree)", "string": True, "grid": (1,0), "no.": 1},
            {"name": "BetaOn", "#": "y-coordinate of map center (degree)", "string": True, "grid": (2,0), "no.": 2},
            {"name": "RELATIVE", "#": "OFF position is given by relative (0) or absolute (1)", "string": False, "grid": (4, 0), "values": ["true", "false"], "default":"true", "no.": 3},
            {"name": "LambdaOff", "#": "off point x-coordinate (= 0 if relative == 1) (degree)", "string": True, "grid": (5,0), "no.": 4},
            {"name": "BetaOff", "#": "off point y-coordinate (= 0 if relative == 1) (degree)", "string": True, "grid": (6,0), "no.": 5},
            {"name": "deltaLambda", "#": "offset_off (= 0 if relative == 0)", "string": True, "grid": (5,1), "no.": 6},
            {"name": "deltaBeta", "#": "offset_off (= 0 if relative == 0)", "string": True, "grid": (6,1), "no.": 7},
            {"name": "position_angle", "#": "position angle of the scan direction (default = 0) (degree)", "string": True, "unit": "deg", "grid": (3,0), "no.": 8},
            {"name": "OTADEL", "#": "dcos (Y/N)", "string": False, "values": ["true", "false"], "grid": (4,1), "no.": 9},
            {"name": "StartPositionX", "#": "start position X relative to lambda_on (arcsec) (sw)", "string": True, "unit": "arcsec", "grid": (1,1), "no.": 10},
            {"name": "StartPositionY", "#": "start position Y relative to lambda_on (arcsec) (sw)", "string": True, "unit": "arcsec", "grid": (2,1), "no.": 11},
            {"name": "COORD_SYS", "#": "(J2000/b1950/galactic/horizontal)", "string": True, "values": ["J2000", "B1950", "Galactic"], "grid": (0,0), "no.": 12},
        ]])
        self.paramlist.append(["scan_property", [
            {"name": "SCAN_DIRECTION", "#": "X (=0) or Y (=1) scan", "values": ["X", "Y"], "string": True, "grid": (0,0), "no.": 1},
            {"name": "n", "#": "Number of scanlines (integer)", "string": False, "grid": (1,0), "no.": 2},
            {"name": "scan_spacing", "#": "scan spacing (arcsec) (default = 60)", "string": True, "unit": "arcsec", "default": 60, "grid": (2,0), "no.": 3},
            {"name": "scan_length", "#": "time duration for 1 scan  (s) ", "string": True, "unit": "s", "grid": (1,1), "no.": 4},
            {"name": "scan_velocity", "#": "scan velocity (arcsec/s) (default = 600) (<=600)", "string": True, "unit": "arcsec/s", "default": 600, "grid": (2,1), "no.": 5},
            {"name": "ramp_pixel", "#": "approach length (pixel) (default = 4)", "string": False, "grid": (3,0), "no.": 6},
            {"name": "integ_on", "#": "ON point integratoin time (s)", "string": True, "unit": "s", "grid": (4,0), "no.": 7}
        ]])
        self.paramlist.append(["calibration", [
            {"name": "integ_off", "#": "OFF point integration time (s)", "string": True, "unit": "s", "no.": 1},
            {"name": "integ_hot", "#": "hot-load integration time (s)", "string": True, "unit": "s", "no.": 2},
            {"name": "off_interval", "#": "off interval (number of scan lines: default = 1)", "string": False, "default": 1, "no.": 3},
            {"name": "load_interval", "#": "time interval of hot-load measurements (min)", "string": True, "unit": "min", "no.": 4}
        ]])

