#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_tk.py                                                        #
#                                                                             #
# PURPOSE:  Functions and classes for TK graphical elements.                  #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  ScrolledTreeTab     ... scrolled Treeview widget mimicking a table         #
#  ScrolledTreeView    ... standard scrolled TreeView widget                  #
#  ScrolledCanvasFrame ... widget frame embedded in a scrolling canvas        #
#  ScatterPlot         ... simple canvas-based scatter plot widget            #
#  ScrolledListBox     ... listbox with scrollbars and convenience functions  #
#  SingleFigFrame      ... canvas widget to display a matplotlib figure       #
#  DoubleScale         ... double-handled slider to define a value range      #
#  ImageBrowser        ... show images in a directory with a preview          #
#                                                                             #
#=============================================================================#
#                                                                             #
# The MIT License (MIT)                                                       #
#                                                                             #
# Copyright (c) 2017 Cormac R. Purcell                                        #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the "Software"),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
#=============================================================================#
import os
import glob
try:
    # Python 2.7x
    import Tkinter as tk
    import ttk
    import tkFont
except Exception:
    # Python 3.x
    import tkinter as tk
    from tkinter import ttk
    import tkinter.font as tkFont
import numpy as np
from PIL import Image, ImageTk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


#-----------------------------------------------------------------------------#
class ScrolledTreeTab(ttk.Frame):
    """Use a ttk.Treeview as a multicolumn ListBox with scrollbars."""
    
    def __init__(self, parent, virtEvent="<<tab_row_selected>>", strPad=10,
                 *args, **kw):
        ttk.Frame.__init__(self, parent, *args, **kw)
        self.parent = parent
        self.rowSelected = None
        self.textSelected = None
        self.virtEvent = virtEvent
        self.strPad = strPad
        
        # Create the treeview and the scrollbars
        self.tree = ttk.Treeview(self, show="headings")
        vsb = ttk.Scrollbar(self, orient="vertical",
                            command=self.tree.yview)
        hsb = ttk.Scrollbar(self, orient="horizontal",
                            command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        # Grid the tree and scrollbars within the container frame
        self.tree.grid(column=0, row=0, sticky="NSWE")
        vsb.grid(column=1, row=0, sticky="NS")
        hsb.grid(column=0, row=1, sticky="WE")
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.rowconfigure(0, weight=1)
        
        # Limit to single row selections
        self.tree.configure(selectmode="browse")
        self.tree.bind("<<TreeviewSelect>>", self._on_row_select)
        
    def _on_row_select(self, event=None):
        """Store the index of the selected row and generate an event.
        The original index of the column is stored in the 'text' property
        of the TreeView widget"""
        
        item =  event.widget.selection()
        if not item=="":
            indx = event.widget.item(item, "text")
            self.rowSelected = int(indx)
            self.textSelected = event.widget.item(item, "value")
            self.event_generate(self.virtEvent)

    def _sortby(self, tree, col, descending):
        """Sort tree contents when a column header is clicked."""

        # Fetch column IDs and data values to sort on
        data = [(tree.set(child, col), child)
                for child in tree.get_children('')]

        # If the data to be sorted is numeric change to float
        data = self._change_numeric_onestep(data)
        
        # now sort the data in place
        data.sort(reverse=descending)
        for i, item in enumerate(data):
            tree.move(item[1], '', i)
            
        # Toggle the sort function
        tree.heading(col, command=lambda col=col: \
                     self._sortby(tree, col, int(not descending)))
        
    def _change_numeric_onestep(self, data):
        """If the data to be sorted is numeric change to float."""
        
        newData = []
        try:
            for child, col in data:
                if child=="None":
                    child = "-inf"   # Regard text "None" as -infinity
                newData.append((float(child), col))
            #print("Sorting column as numeric data.") # DEBUG
            return newData
        except Exception:
            #print("Sorting column as ascii data.") # DEBUG
            return data

    def name_columns(self, colNames):
        """Insert the column headings"""
        
        self.tree['columns'] = colNames
        for col in colNames:
            self.tree.heading(col, text=col, command=lambda c=col:
                              self._sortby(self.tree, c, 0))
            
            # Set the column width to the width of the header string
            strWidth = tkFont.Font().measure(col.title())
            self.tree.column(col, width=strWidth + self.strPad)
            self.tree.column(col, minwidth=strWidth + self.strPad)
        
    def insert_rows(self, rows, colNames=None):
        """Insert rows from a 2D iterable object."""

        # Populate the headers and set the sort function
        if colNames is None:
            colNames = ["Row "+ str(x+1) for x in range(len(rows[0]))]
        if len(self.tree['columns'])==0:
            self.tree['columns'] = colNames
            for col in colNames:
                self.tree.heading(col, text=col, command=lambda c=col:
                                  self._sortby(self.tree, c, 0))
        
            # Set the column width to the width of the header string
            strWidth = tkFont.Font().measure(col.title())
            self.tree.column(col, width=strWidth + self.strPad)
            self.tree.column(col, minwidth=strWidth + self.strPad)
            
        # Populate the rows
        rowIndx = 0
        for row in rows:
            row = [str(x) for x in row]   # Convert to plain strings
            self.tree.insert('', 'end', values=row, text=str(rowIndx))
            rowIndx += 1
            
            # Adjust the column width (& minwidth) to fit each value
            for i, val in enumerate(row):
                strWidth = tkFont.Font().measure(val.title())
                if self.tree.column(colNames[i], width=None)<\
                   (strWidth + self.strPad):
                    self.tree.column(colNames[i], width=strWidth +
                                     self.strPad)
                    self.tree.column(colNames[i], minwidth=strWidth +
                                     self.strPad)
                    
    def get_indx_selected(self):
        """Return the index of the last row selected."""

        if self.rowSelected is None:
            return None
        else:
            return int(self.rowSelected)

    def get_text_selected(self):
        """Return the text of the selected row"""

        if self.textSelected is None:
            return None
        else:
            return self.textSelected
            
    def get_all_text(self):
#        """Return a list of all text entries in the listbox"""
        try:
            itemLst = self.tree.get_children()
            valueLst = []
            for item in itemLst:
                valueLst.append(self.tree.item(item, "value"))
            return valueLst
        except Exception:
            return None
        
    def insert_recarray(self, arr):
        """Insert a numpy.recarray into the treeview widget"""
        
        colNames = arr.dtype.names
        self.insert_rows(arr, colNames)

    def clear_selected(self):
        """Clear all the entries from the table."""        
        try:
            idx = self.tree.selection()
            self.tree.delete(idx)
            self.rowSelected = None
            self.textSelected = None
            return idx
        except Exception:
            return None
        
    def clear_entries(self):
        """Clear all the entries from the table."""
        
        try:
            x = self.tree.get_children() 
            for entry in x:
                self.tree.delete(entry)
            self.rowSelected = None
        except Exception:
            pass

        
#-----------------------------------------------------------------------------#
class ScrolledTreeView(ttk.Frame):
    """A ttk.Treeview with scrollbars. Selecting a row generates a virtual
    event so that the parent object can query the selection."""
    
    def __init__(self, parent, virtEvent="<<tree_selected>>", *args, **kw):
        ttk.Frame.__init__(self, parent, *args, **kw)
        self.parent = parent
        self.textSelected = None
        self.textRootSelected = None
        self.virtEvent = virtEvent
    
        # Create the treeview and the scrollbars
        self.tree = ttk.Treeview(self)
        vsb = ttk.Scrollbar(self, orient="vertical",
                            command=self.tree.yview)
        hsb = ttk.Scrollbar(self, orient="horizontal",
                            command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        
        # Grid the tree and scrollbars within the container frame
        self.tree.grid(column=0, row=0, sticky="NSEW")
        vsb.grid(column=1, row=0, sticky="NS")
        hsb.grid(column=0, row=1, sticky="WE")
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.rowconfigure(0, weight=1)

        # Binding to handle row selections
        self.tree.bind("<<TreeviewSelect>>", self._on_row_select)
        
    def _on_row_select(self, event=None):
        """Triggered when a row is selected: fills the class variables with
        the text selected."""
        
        item =  event.widget.selection()
        if not item=="":
            self.textSelected = event.widget.item(item, "text")
            self.textRootSelected = event.widget.item(item, "text")
            while True:                
                parentItem =  event.widget.parent(item)
                if parentItem=="":
                    break
                self.textRootSelected = event.widget.item(parentItem, "text")
                item = parentItem
            self.event_generate(self.virtEvent)

    def get_text_selected(self):
        return self.textSelected, self.textRootSelected


#-----------------------------------------------------------------------------#
class ScrolledCanvasFrame(ttk.Frame):
    """Canvas with embedded frame, window & scrollbar. Used to provide a
    scrolled pane in which to layout widgets. EXPERIMENTAL."""
    
    def __init__(self, parent, *args, **kw):
        ttk.Frame.__init__(self, parent, *args, **kw)
        self.parent = parent

        # Create the canvas and the scrollbars
        self.canvas = tk.Canvas(self, border=0, highlightthickness=0)
        vsb = ttk.Scrollbar(self, orient="vertical",
                            command=self.canvas.yview)
        #hsb = ttk.Scrollbar(self, orient="horizontal",
        #                    command=self.canvas.xview)
        self.canvas.configure(yscrollcommand=vsb.set)#, xscrollcommand=hsb.set)
        
        self.canvas.xview_moveto(0)
        self.canvas.yview_moveto(0)
        
        # Grid the canvas and scrollbars within the container frame
        self.canvas.grid(column=0, row=0, sticky="NSEW")
        vsb.grid(column=1, row=0, sticky="NS")
        #hsb.grid(column=0, row=1, sticky="WE")
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.rowconfigure(0, weight=1)
        
        # Now create a frame within the canvas
        self.interior = ttk.Frame(self.canvas)
        self.winID = self.canvas.create_window((0,0), window=self.interior,
                                         anchor="nw", tags="self.interior")
        
        self.interior.bind('<Configure>', self._configure_interior)
        self.canvas.bind('<Configure>', self._configure_canvas)

    def _configure_interior(self, event):
        """Set the interior frame properties"""
        
        size = (self.interior.winfo_reqwidth(),
                self.interior.winfo_reqheight())
        self.canvas.config(scrollregion="0 0 %s %s" % size)
        #self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        if self.interior.winfo_reqwidth() != self.canvas.winfo_width():
            self.canvas.config(width=self.interior.winfo_reqwidth())

    def _configure_canvas(self, event):
        """Set the canvas window properties."""
        
        if self.interior.winfo_reqwidth() != self.canvas.winfo_width():
            self.canvas.itemconfigure(self.winID,
                                      width=self.canvas.winfo_width())
        

#-----------------------------------------------------------------------------#
class ScatterPlot(ttk.Frame):
    """Canvas configured as a simple scatterplot widget."""
    
    #                                                left, right, bottom, top
    def __init__(self, parent, width=500, height=500, axPad=(100,25,70,25),
                 tickLen=10, nXticks=3, nYticks=3, padF=0.05, aspect="free",
                 pntSize=2.0, interactive=False, *args, **kw):
        ttk.Frame.__init__(self, parent, *args, **kw)
        self.parent = parent

        # Create the canvas and grid
        self.canvas = tk.Canvas(self, background="white", width=width,
                                height=height, highlightthickness=0)
        self.canvas.grid(column=0, row=0, padx=0, pady=0)
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        
        # Data and plot parameters in world coordinates
        self.xArr = None 
        self.yArr = None
        self.xPlotMax = None
        self.xPlotMin = None
        self.yPlotMax = None
        self.yPlotMin = None
        self.padF = padF
        
        # Data and plot parameters in canvas coordinates
        self.xCanArr = None 
        self.yCanArr = None
        self.xCanMax = None
        self.xCanMin = None
        self.yCanMax = None
        self.yCanMin = None

        # Canvas and axis layout parameters
        self.width = width
        self.height = height
        self.axPad = axPad
        self.aspect = aspect
        self.pntSize = pntSize

        # Tickmarks and labels
        self.xTicks = None
        self.yTicks = None
        self.xTickVals = None
        self.xTickVals = None
        self.xPwr = None
        self.yPwr = None
        self.nXticks = nXticks
        self.nYticks = nYticks
        self.tickLen = tickLen
        
        # Bindings & location variables.
        self.oldx = 0
        self.oldy = 0
        if interactive:
            self._create_bindings()

    def _init_axes(self):
        """Draw the axis lines"""
        
        self.canvas.delete("all")
        self.xCanMin = self.axPad[0]
        self.xCanMax = self.width - self.axPad[1]
        self.yCanMin = self.height - self.axPad[2]
        self.yCanMax = self.axPad[3]
        self.canvas.create_line(self.xCanMin, self.yCanMin,
                                self.xCanMin, self.yCanMax, width=1)
        self.canvas.create_line(self.xCanMin, self.yCanMin,
                                self.xCanMax, self.yCanMin, width=1)
        self.canvas.create_line(self.xCanMin, self.yCanMax,
                                self.xCanMax, self.yCanMax, width=1)
        self.canvas.create_line(self.xCanMax, self.yCanMin,
                                self.xCanMax, self.yCanMax, width=1)
        
    def _scale_data_axes(self):
        """Autoscale the axes and convert the data to canvas coordinates.
        Leave a border of range*padF around the data."""
        
        xMin, xMax = (np.min(self.xArr), np.max(self.xArr))
        yMin, yMax = (np.min(self.yArr), np.max(self.yArr))
        xRng = (xMax - xMin)
        yRng = (yMax - yMin)
        if self.aspect=="equal":
            rng = max(xRng, yRng)
            xRng = rng
            yRng = rng
        else:
            xRng = 1 if xRng==0 else xRng
            yRng = 1 if yRng==0 else yRng
        if xMin==xMax:
            xMin = xMin - xRng/2.0
            xMax = xMax + xRng/2.0
        if yMin==yMax:
            yMin = yMin - yRng/2.0
            yMax = yMax + yRng/2.0
        self.xPlotMin = xMin - xRng * self.padF
        self.xPlotMax = xMax + xRng * self.padF
        self.yPlotMin = yMin - yRng * self.padF
        self.yPlotMax = yMax + yRng * self.padF
        
        self.xCanArr, self.yCanArr = self._world2canvas(self.xArr, self.yArr)

    def _world2canvas(self, x=None, y=None):
        """Convert x-y arrays of world coordinates to canvas coordinates"""
        
        l = None
        m = None
        xCanRng = self.xCanMax - self.xCanMin
        yCanRng = self.yCanMax - self.yCanMin
        xPlotRng = self.xPlotMax - self.xPlotMin
        yPlotRng = self.yPlotMax - self.yPlotMin
        if x is not None:
            l = (x-self.xPlotMin)*xCanRng/xPlotRng + self.xCanMin
        if y is not None:
            m = (y-self.yPlotMin)*yCanRng/yPlotRng + self.yCanMin
        return l, m

    def _canvas2world(self, l=None, m=None):
        """Convert l-m arrays of canvas coordinates to world coordinates"""
        
        x = None
        y = None
        xCanRng = self.xCanMax-self.xCanMin
        yCanRng = self.yCanMax-self.yCanMin
        xPlotRng = self.xPlotMax-self.xPlotMin
        yPlotRng = self.yPlotMax-self.yPlotMin
        if l is not None:
            x = (l-self.xCanMin)*xPlotRng/xCanRng + self.xPlotMin
        if m is not None:
            y = (m-self.yCanMin)*yPlotRng/yCanRng + self.yPlotMin
        return x, y
        
    def _layout_ticks(self):
        """Calculate positions of ticks with rounded values."""

        # Calculate the tick values based on the number of requested ticks
        # and the 10-power of the axis range. The number of ticks will be
        # likely different than requested due to the rounding, but the
        # labels will be numbers with finite significant figures.
        def calc_tick_vals(xMin, xMax, nTicks):
            rng = xMax-xMin
            pwr = int(np.floor(np.log10(rng))-1)
            d = np.round(rng/(float(nTicks)*10.0**pwr)) * 10.0**pwr
            start = np.round(xMin/10.0**pwr) * 10.0**pwr
            start = np.round(start/d)*d  # Force ticks to shift to zero
            ticks = []
            i = 0
            while True:
                tick = start + i*d
                i += 1
                if tick>=xMin and tick<=xMax:
                    ticks.append(tick)
                if tick>xMax:
                    break
            return np.array(ticks), pwr

        # Calculate the tick values and convert to canvas coordinates
        self.xTickVals, self.xPwr = calc_tick_vals(self.xPlotMin,
                                                   self.xPlotMax,
                                                   self.nXticks)
        self.yTickVals, self.yPwr = calc_tick_vals(self.yPlotMin,
                                                   self.yPlotMax,
                                                   self.nYticks)
        self.xTicks, dummy = self._world2canvas(x=self.xTickVals)
        dummy, self.yTicks = self._world2canvas(y=self.yTickVals)
        
    def _draw_ticks(self):
        """Draw the tick-marks on the X and Y axis"""

        # Set the string formatting codes
        if self.xPwr<=0:
            xFmt = "{:" +".{:s}f".format(str(abs(self.xPwr))) + "}"
        else:
            xFmt = "{}"
        if self.yPwr<=0:
            yFmt = "{:" +".{:s}f".format(str(abs(self.yPwr))) + "}"
        else:
            yFmt = "{}"
            
        # Draw X ticks
        for i in range(len(self.xTicks)):
            self.canvas.create_line(self.xTicks[i],
                                    self.yCanMin,
                                    self.xTicks[i],
                                    self.yCanMin-self.tickLen,
                                    width=1)
            self.canvas.create_text(self.xTicks[i],
                                    self.yCanMin+2*self.tickLen,
                                    text=xFmt.format(self.xTickVals[i]))
            
        # Draw Y ticks
        for i in range(len(self.yTicks)):
            self.canvas.create_line(self.xCanMin,
                                    self.yTicks[i],
                                    self.xCanMin+self.tickLen,
                                    self.yTicks[i],
                                    width=1)
            self.canvas.create_text(self.xCanMin-2*self.tickLen,
                                    self.yTicks[i],
                                    text=yFmt.format(self.yTickVals[i]),
                                    anchor=tk.E,)

    def _draw_points(self, pntSize=3, colour="blue"):
        """Draw the points as circles on the canvas."""
            
        for i in range(len(self.xCanArr)):
            x = self.xCanArr[i]
            y = self.yCanArr[i]
            item = self.canvas.create_oval(x-pntSize, y-pntSize,
                                           x+pntSize, y+pntSize,
                                           width=1, outline='black',
                                           fill=colour)
            self.canvas.addtag_withtag('point', item)

    # Event handling code ----------------------------------------------------#

    def _point_enter(self, evt):   
        self.canvas.itemconfigure('current', fill='red')
        self.oldx = evt.x
        self.oldy = evt.y
        
    def _point_leave(self, evt):   
        self.canvas.itemconfigure('current', fill='SkyBlue2')
        
    def _point_drag(self, evt):
        x, y = self.canvas.canvasx(evt.x), self.canvas.canvasy(evt.y)
        self.canvas.move('current', x-self.oldx, y-self.oldy)
        self.oldx, self.oldy = x, y

    def _create_bindings(self):
        self.canvas.tag_bind('point', '<Any-Enter>', self._point_enter)
        self.canvas.tag_bind('point', '<Any-Leave>', self._point_leave)
        self.canvas.tag_bind('point', '<B1-Motion>', self._point_drag)
        
    #-------------------------------------------------------------------------#
            
    def draw_zerolines(self):
        """Draw Lines at (0,0) in world coordinates"""
        
        zero = np.array([0.0, 0.0])
        x, y = self._world2canvas(zero[0], zero[1])
        if x>self.xCanMin and x<self.xCanMax:
            item = self.canvas.create_line(x, self.yCanMin, x, self.yCanMax,
                                           width=1, fill="LightGreen")
            self.canvas.tag_lower(item)
        if y>self.yCanMax and y<self.yCanMin:
            item = self.canvas.create_line(self.xCanMin, y, self.xCanMax, y,
                                           width=1, fill="LightGreen")
            self.canvas.tag_lower(item)
        
    def set_xlabel(self, label=""):
        """Label the X-axis"""
        
        self.canvas.create_text(self.axPad[0]+(self.xCanMax-self.xCanMin)/2,
                                self.yCanMin+self.axPad[2]/1.5,
                                text=label, fill='red')
        
    def set_ylabel(self, label=""):
        """Label the Y-axis"""
        try:
            self.canvas.create_text(self.axPad[2]/4,
                            self.axPad[3] + (self.yCanMin-self.yCanMax)/2,
                            text=label, fill='red', anchor=tk.N, angle=90.0)
        except Exception:
            pass
    
    def display_message(self, text=""):
        """Display text in the centre of the canvas."""
        
        self.canvas.create_text(self.width/2, self.height/2, justify="center",
                                text=text, fill='black')
        
    def load_data(self, xArr, yArr):
        """Load the xy data arrays and draw the data points"""

        self.xArr = xArr
        self.yArr = yArr

        # (Re)initialise the canvas and draw the axis lines
        self._init_axes()
        
        # Scale the data and axis limits to canvas coordinates
        self._scale_data_axes()

        # Layout the tick-marks and plot
        self._layout_ticks()
        self._draw_ticks()
        
        # Draw the data
        self._draw_points(pntSize=self.pntSize, colour="SkyBlue2")
            

#-----------------------------------------------------------------------------#
class ScrolledListBox(ttk.Frame):

    def __init__(self, parent, selectmode="single",
                 virtEvent="<<list_row_selected>>", *args, **kw):
        ttk.Frame.__init__(self, parent, *args, **kw)
        self.parent = parent
        self.virtEvent = virtEvent
        self.rowSelected = None
        self.textSelected = None
        
        # Create the listbox and the scrollbars
        self.listBox = tk.Listbox(self, selectmode=selectmode)
        vsb = ttk.Scrollbar(self, orient="vertical",
                            command=self.listBox.yview)
        hsb = ttk.Scrollbar(self, orient="horizontal",
                            command=self.listBox.xview)
        self.listBox.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        # Grid the listbox and scrollbars within the container frame
        self.listBox.grid(column=0, row=0, sticky="NSWE")
        vsb.grid(column=1, row=0, sticky="NS")
        hsb.grid(column=0, row=1, sticky="WE")
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=0)
        self.rowconfigure(0, weight=1)
        
        # Binding to handle row selections: generates a virtual event
        self.listBox.bind("<ButtonRelease-1>", self._on_row_select)

    def insert(self, item):
        """Insert a single items"""
        
        self.listBox.insert(tk.END, str(item))
    
    def insert_list(self, items):
        """Insert items from a 1D iterable object"""

        for item in items:
            self.listBox.insert(tk.END, str(item))
            
    def clear(self):
        """Clear all items from the listbox"""
        
        self.listBox.delete(0, tk.END)  
        
    def clear_selected(self):
        """Clear the current selection"""
        
        try:            
            index = self.listBox.curselection()[0]
            self.listBox.delete(index)
        except Exception:
            pass

    def get_text_selected(self):
        """Convenience function to get selected text"""
        try:
            index = self.listBox.curselection()[0]
            return self.listBox.get(index)
        except Exception:
            return None

    def get_row_selected(self):
        """Convenience function to get selected index"""
        try:
            return self.listBox.curselection()[0]
        except Exception:
            return None

    def get_all_text(self):
        """Return a list of all text entries in the listbox"""
        try:
            return self.listBox.get(0, tk.END)
        except Exception:
            return None
        
    def _on_row_select(self, event=None):
        """Triggered when a row is selected. Stores the last row index
        and text selected, and generates a virtual event."""
        
        try:
            self.rowSelected = self.listBox.curselection()[0]
            self.textSelected = self.listBox.get(self.rowSelected)
            self.event_generate(self.virtEvent)
        except Exception:
            pass
        

#-----------------------------------------------------------------------------#
class SingleFigFrame(ttk.Frame):
    """Use a canvas widget to display a matplotlib figure. EXPERIMENTAL"""

    def __init__(self, parent, *args, **kw):
        ttk.Frame.__init__(self, parent, *args, **kw)
        self.parent = parent

        # Create the blank figure canvas and grid its tk canvas
        self.fig = Figure(figsize=(4.5, 4.0))
        self.figCanvas = FigureCanvasTkAgg(self.fig, master=self)
        self.figCanvas.show()
        self.canvas = self.figCanvas.get_tk_widget()
        self.canvas.grid(column=0, row=0, padx=0, pady=0, sticky="NSEW")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.ax = None

    def add_axis(self):
        if self.ax is None:
            self.ax = self.fig.add_subplot(111)
        return self.ax

    def show(self):
        self.figCanvas.show()

    
#-----------------------------------------------------------------------------#
class DoubleScale(ttk.Frame):
    """Create a custom double slider widget in a canvas"""

    def __init__(self, parent, width=400, handlesize=7, from_=0, to=100,
                 initLeft=None, initRight=None, tickIntMajor=20,
                 tickIntMinor=None, linewidth=1, ticklen=20, yPad=30, xPad=30,
                 labFmt="{}"):
        ttk.Frame.__init__(self, parent)
        self.parent = parent        
        bgColour = ttk.Style().lookup("TFrame", "background")

        # Canvas & ruler layout parameters
        self.width = width
        self.yZero = yPad + ticklen
        self.height = self.yZero + handlesize*3 + linewidth
        xPad = max(xPad, handlesize/2.0)
        self.canMin = xPad
        self.canMax = width-xPad
        self.lineWidth = linewidth
        self.tickLen = ticklen
        self.handleSize = handlesize
        
        # Data range and intervals
        self.xMin = from_
        self.xMax = to
        self.dX = tickIntMajor
        self.labFmt = labFmt
        if tickIntMinor is None:
            tickIntMinor = tickIntMajor/5.
        self.dx = tickIntMinor
        self.dxCan = (float(self.canMax-self.canMin)*self.dx/
                      float(self.xMax-self.xMin))
        
        # Store the current handle values
        self.valueLeft = tk.DoubleVar()
        self.valueRight = tk.DoubleVar()
        self.limitLeft = None
        self.limitRight = None
        self.x = None
        
        # Insert the canvas
        self.canvas = tk.Canvas(self, background=bgColour, width=self.width,
                                height=self.height, highlightthickness=0)
        self.canvas.grid(row=0, column=0, columnspan=3,padx=0, pady=0)
        
        # Draw the axes
        self._draw_ruler()

        # Draw the handles & create the bindings
        if initLeft is None:
            initLeft = self.canMin
        if initRight is None:
            initRight = self.canMax
        hLeft = self._draw_handle(self._world2canvas(initLeft), 'left')
        hRight = self._draw_handle(self._world2canvas(initRight), 'right')
        self._create_bindings()

        # Draw the limits as entry boxes
        cWidth = tkFont.Font().measure("0")
        self.lowLab = ttk.Label(self, textvariable=self.valueLeft, width=7,
                                anchor="center")
        self.lowLab.grid(row=1, column=0, padx=0, pady=5, sticky="W")
        self.highLab = ttk.Label(self, textvariable=self.valueRight, width=7,
                                 anchor="center")
        self.highLab.grid(row=1, column=2, padx=0, pady=0, sticky="E")
        
    def _draw_ruler(self):
        """Draw the ruler line with tick marks"""

        # Draw the ruler line
        self.canvas.delete("all")
        self.canvas.create_line(self.canMin, self.yZero,
                                self.canMax, self.yZero,
                                width=self.lineWidth)
        
        # Lay out the major ticks and draw labels
        self.tickVals = np.arange(self.xMin, self.xMax+self.dX, self.dX)
        self.ticks = self._world2canvas(self.tickVals)
        self._draw_ticks(yLine=self.yZero,
                         xLst=self.ticks,
                         tickLen=self.tickLen,
                         lineWidth=self.lineWidth,
                         yLab=(self.yZero-self.tickLen)/2.,
                         valLst=self.tickVals)

        # Lay out the minor ticks
        mintickVals = np.arange(self.xMin, self.xMax+self.dx, self.dx)
        self.minticks = self._world2canvas(mintickVals)
        self._draw_ticks(yLine=self.yZero,
                         xLst=self.minticks,
                         tickLen=self.tickLen/2,
                         lineWidth=self.lineWidth)

    def _draw_ticks(self, yLine, xLst, tickLen=20, lineWidth=1, yLab=None,
                    valLst=None):
        """Draw the tick marks on the ruler"""
        
        for i in range(len(xLst)):
            self.canvas.create_line(xLst[i],
                                    yLine + np.round(lineWidth/2.),
                                    xLst[i],
                                    yLine + np.round(lineWidth/2.)-tickLen,
                                    width=lineWidth)
            if yLab is not None and valLst is not None:
                self.canvas.create_text(xLst[i],
                                        yLab,
                                        text='{}'.format(valLst[i]))
    
    def _world2canvas(self, x):
        """Convert an array of world coordinates to canvas coordinates"""

        canRng = float(self.canMax - self.canMin)
        xRng = float(self.xMax - self.xMin)
        l = (x-self.xMin)*canRng/xRng + self.canMin
        
        return l

    def _canvas2world(self, l):
        """Convert an array of canvas coordinates to world coordinates"""

        canRng = float(self.canMax - self.canMin)
        xRng = float(self.xMax - self.xMin)
        x = (l-self.canMin)*xRng/canRng + self.xMin
        
        return x
    
    def _draw_handle(self, x, tag):
        """Draw a handle as a polygon."""
        
        y = self.yZero-self.tickLen/3.0
        size = self.handleSize
        if tag=="left":
            polygon = (x, y,
                       x, y+size*2.,
                       x, y+size*3.,
                       x-size*2, y+size*3.,
                       x-size*2, y+size*2.,
                       x, y)
        else:
            polygon = (x, y,
                       x, y+size*2.,
                       x, y+size*3.,
                       x+size*2, y+size*3.,
                       x+size*2, y+size*2.,
                       x, y)
        item = self.canvas.create_polygon(polygon, fill="lightblue",
                                          outline="black")
        
        self.canvas.itemconfigure(item, tag=('handle', tag))
        self._set_handle_value(item, x, self.labFmt)

        return item
    
    #-------------------------------------------------------------------------#
    def _create_bindings(self):
        self.canvas.tag_bind('handle', '<1>', self._sel_handle)
        self.canvas.bind('<B1-Motion>', self._move_handle)
        self.canvas.bind('<Any-ButtonRelease-1>', self._release_handle)
            
    def _sel_handle(self, evt):
        """Triggered when the user clicks on a handle"""
        
        self.x = self.canvas.canvasx(evt.x)
        self.canvas.addtag_withtag('active', 'current')
        self.canvas.itemconfigure('active',{'fill': 'red', 'stipple': ''})

    def _move_handle(self, evt):
        """Triggered when the user clicks and drags handle"""

        if not self.canvas.find_withtag('active'):
            return

        # Current cursor X
        cx = self.canvas.canvasx(evt.x)         

        # Find the difference between the clicked position & edge
        item = self.canvas.find_withtag('active')
        edgeX = self.canvas.coords(item)[0]
        dX = self.x - edgeX
        
        # Prevent collisions
        if "left" in self.canvas.gettags(item):
            limLeft = self.canMin
            limRight = self.limitRight
        if "right" in self.canvas.gettags(item):
            limLeft = self.limitLeft
            limRight = self.canMax      
        if cx <= limLeft +dX:
            cx = limLeft +dX
        if cx >= limRight +dX:
            cx = limRight +dX           
        self.canvas.move('active', cx - self.x, 0)
        self.x = cx

        # Set the handle variable
        self._set_handle_value(item, edgeX)

    def _release_handle(self, evt):
        """Triggered when the user releases the mouse"""

        # Don't trigger outside the handle polygon
        if not self.canvas.find_withtag('active'):
            return

        # Change the colour and remove 'active' tag
        self.canvas.itemconfigure('active',
                                  {'fill': 'lightblue', 'stipple': ''})
        self.canvas.dtag('active')

    def _set_handle_value(self, item, value, fmt="{}"):
        
        valFmt = "{0:.2f}".format(self._canvas2world(value))
        if "left" in self.canvas.gettags(item):
            self.limitLeft = value
            self.valueLeft.set(valFmt)
        if "right" in self.canvas.gettags(item):
            self.limitRight = value
            self.valueRight.set(valFmt)


#-----------------------------------------------------------------------------#
class ImageBrowser(ttk.Frame):    
    
    def __init__(self, parent, virtEvent="<<tab_row_selected>>", path=None,
                 thumbSize=(300, 300), *args, **kw):
        ttk.Frame.__init__(self, parent, *args, **kw)
        self.parent = parent
        if path is None:
            path = os.getcwd()
        self.path = path
        self.thumbSize = thumbSize
        self.filePath = None        
        self.virtEvent = virtEvent

        # Allowed filetypes (can be read by PIL)
        self.FILETYPES = ['.bmp', '.dib', '.dcx', '.gif', '.im', '.jpg',
                          '.jpe', '.jpeg', '.pcd', '.pcx', '.png', '.pbm',
                          '.pgm', '.ppm', '.psd', '.tif', '.tiff', '.xbm', 
                          '.xpm']
        self.FILETYPES += [x.upper() for x in self.FILETYPES]
        
        # Create the directory box & browse button
        self.dirFrm = ttk.Frame(self)
        self.dirFrm.grid(column=0, row=0, columnspan=2, padx=0, pady=0,
                         sticky="EW")
        self.dirFrm.columnconfigure(1, weight=1)
        self.dirLab = ttk.Label(self.dirFrm, text="Image Directory: ")
        self.dirLab.grid(row=0, column=0, padx=5, pady=5, sticky="W")
        self.dirName = tk.StringVar()
        self.dirEnt = ttk.Entry(self.dirFrm, width=50,
                                textvariable=self.dirName)
        self.dirEnt.grid(row=0, column=1, padx=5, pady=5, sticky="EW")
        self.dirEnt.bind('<Return>', self._load_dir)
        self.browseBtn = ttk.Button(self.dirFrm, text="Browse",
                                    command=self._browse_dir)
        self.browseBtn.grid(row=0, column=2, padx=5, pady=5, sticky="E")
        
        # Create the listbox and the scrollbars
        self.files = ScrolledListBox(self)
        self.files.grid(column=0, row=1, sticky="NSW")
        self.files.listBox.bind('<ButtonRelease-1>', self._load_preview)

        # Create the image preview pane
        self.imgFrm = ttk.Labelframe(self, text=" Preview ")
        self.imgFrm.grid(column=1, row=1, padx=(10,5), pady=0, sticky="NSWE")
        self.imgThumb = None       # Image handle to prevent garbage collection
        self.imgLab = ttk.Label(self.imgFrm, image=self.imgThumb, border=2)
        self.imgLab.grid(row=0, column=0, padx=5, pady=5),# sticky="NSEW")
        self.imgFrm.columnconfigure(0, weight=1)
        self.imgFrm.rowconfigure(0, weight=1)

        # Set the expansion properties
        self.rowconfigure(1, weight=1)
        self.columnconfigure(1, weight=1)
        
        # Populate the file list
        if os.path.exists(self.path):
            self.dirName.set(self.path)
        self._load_dir()
                
    def _browse_dir(self):
        imgPath = tkFileDialog.askdirectory(initialdir=self.dirName.get(),
                                        parent=self,
                                        title='Select a new image directory',
                                        mustexist=True)
        if imgPath:
            self.dirName.set(imgPath)
            self._load_dir()
    
    def _load_dir(self, *args):
        dn = os.path.join(self.dirName.get(), '*')
        self.files.clear()
        
        # Find the files in the given directory, insert into listbox
        for f in glob.iglob(dn):
            if os.path.splitext(f)[1] in self.FILETYPES:
                self.files.insert(os.path.basename(f))
        
    def _load_preview(self, event):
        item = self.files.listBox.curselection()
        fn = self.files.listBox.get(item)
        imagePath = os.path.join(self.dirName.get(), fn)

        # Get the current size of the preview frame
        self.thumbSize = (self.imgFrm.winfo_width(), self.imgFrm.winfo_height())
        
        #try:
        if True:
            imgPIL = Image.open(imagePath)
            imgPIL.thumbnail(self.thumbSize, Image.ANTIALIAS)
            if imgPIL.mode == '1':
                self.imgThumb = ImageTk.BitmapImage(imgPIL)
            else:
                self.imgThumb = ImageTk.PhotoImage(imgPIL)
            self.imgLab.configure(image=self.imgThumb)
            self.imagePath = imagePath
            self.event_generate(self.virtEvent)
            
        #except Exception:
        #    
        #    # Colour failed files in red
        #    self.files.listBox.itemconfig(item, background='red',
        #                                  selectbackground='red')

    def display_image(self, imagePath):
            imgPIL = Image.open(imagePath)
            imgPIL.thumbnail(self.thumbSize, Image.ANTIALIAS)
            if imgPIL.mode == '1':
                self.imgThumb = ImageTk.BitmapImage(imgPIL)
            else:
                self.imgThumb = ImageTk.PhotoImage(imgPIL)
            self.imgLab.configure(image=self.imgThumb)
            self.imagePath = imagePath
            self.event_generate(self.virtEvent)

