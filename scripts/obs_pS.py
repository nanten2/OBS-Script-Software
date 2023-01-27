from . import Gvars
from . import MainApp

def main():
    Gvars.gvar_init()
    initdialog = MainApp.InitDialog(Gvars.root)
    Gvars.root.mainloop()


if __name__ == '__main__':
    main()
