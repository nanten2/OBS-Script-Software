import Gvars
import MainApp


if __name__ == '__main__':
    Gvars.gvar_init()
    initdialog = MainApp.InitDialog(Gvars.root)
    Gvars.root.mainloop()
