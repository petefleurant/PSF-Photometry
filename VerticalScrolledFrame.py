# Based on
#   https://web.archive.org/web/20170514022131id_/http://tkinter.unpythonic.net/wiki/VerticalScrolledFrame

try:  # Python 2
    import tkinter as tk
    import tkinter.ttk as ttk
    from tkinter.constants import *
except ImportError:  # Python 2
    import Tkinter as tk
    import ttk
    from tkinter.constants import *


class VerticalScrolledFrame(ttk.Frame):
    """A pure Tkinter scrollable frame that actually works!
    * Use the 'interior' attribute to place widgets inside the scrollable frame.
    * Construct and pack/place/grid normally.
    * This frame only allows vertical scrolling.
    """
    def __init__(self, parent, height, *args, **kw):
        ttk.Frame.__init__(self, parent, *args, **kw)

        # Create a canvas object and a vertical scrollbar for scrolling it.
        vscrollbar = ttk.Scrollbar(self, orient=VERTICAL)
        vscrollbar.pack(fill=Y, side=RIGHT, expand=FALSE)
        #canvas = tk.Canvas(self, bd=0, highlightthickness=0,yscrollcommand=vscrollbar.set)
        canvas = tk.Canvas(self, bd=0, height=height, highlightthickness=0,yscrollcommand=vscrollbar.set)
        #canvas.pack(side=LEFT, fill=BOTH, expand=TRUE)
        canvas.pack(side=LEFT, fill=NONE, expand=FALSE)
        canvas.configure(background='black')
        vscrollbar.config(command=canvas.yview)

        # Reset the view
        canvas.xview_moveto(0)
        canvas.yview_moveto(0)

        # Create a frame inside the canvas which will be scrolled with it.
        self.interior = interior = ttk.Frame(canvas, height=height)
        interior_id = canvas.create_window(0, 0, window=interior, anchor=NW)

        # Track changes to the canvas and frame width and sync them,
        # also updating the scrollbar.
        def _configure_interior(event):
            # Update the scrollbars to match the size of the inner frame.
            size = (interior.winfo_reqwidth(), interior.winfo_reqheight())
            canvas.config(scrollregion="0 0 %s %s" % size)
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # Update the canvas's width to fit the inner frame.
                canvas.config(width=interior.winfo_reqwidth())
            #if interior.winfo_reqheight() != canvas.winfo_height():
            #    # Update the canvas's width to fit the inner frame.
            #    canvas.config(height=interior.winfo_reqheight())
        interior.bind('<Configure>', _configure_interior)

        def _configure_canvas(event):
            if interior.winfo_reqwidth() != canvas.winfo_width():
                # Update the inner frame's width to fill the canvas.
                canvas.itemconfigure(interior_id, width=canvas.winfo_width())
        canvas.bind('<Configure>', _configure_canvas)


if __name__ == "__main__":

    class SampleApp(tk.Tk):
        def __init__(self, *args, **kwargs):
            root = tk.Tk.__init__(self, *args, **kwargs)

            self.frame = VerticalScrolledFrame(root)
            self.frame.pack()
            self.label = ttk.Label(self, text="Shrink the window to activate the scrollbar.")
            self.label.pack()
            buttons = []
            for i in range(100):
                buttons.append(ttk.Button(self.frame.interior, text="Button " + str(i)))
                buttons[-1].pack()

    app = SampleApp()
    app.mainloop()