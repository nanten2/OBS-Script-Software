import tkinter as tk
from guppy import hpy
import MainApp

if __name__ == '__main__':
    root = tk.Tk()
    h = hpy()
    initdialog = MainApp.InitDialog(root)
    root.mainloop()

print(h.heap())