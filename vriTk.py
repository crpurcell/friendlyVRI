#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     vriTK.py                                                          #
#                                                                             #
# PURPOSE:  A virtual interferometer application written in Tkinter.          #
#                                                                             #
# REQUIRED: Requires numpy, tkinter, matplotlib                               #
#                                                                             #
# MODIFIED: 03-Mar-2017 by cpurcell                                           #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# App                ... class to create the root and array config windows    #
# ArrayChooserFrame  ... class defining the array chooser interface           #
# ModelImageFrame    ... class defining the model image interface             #
# ModelFFTframe      ...                                                      #
# uvCoverageFrame    ...                                                      #
# ObservedFFTframe   ...                                                      #
# SynthBeamFrame     ...                                                      #
# ObservedImageFrame ...                                                      #
# StatusFrame        ...                                                      #
#                                                                             #
#=============================================================================#
#                                                                             #
# The MIT License (MIT)                                                       #
#                                                                             #
# Copyright (c) 2015 Cormac R. Purcell                                        #
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

# Window geometry
geometryMainWin = "1024x800"
geometryAConfWin = "1360x450"

#-----------------------------------------------------------------------------#

import os
import sys
try:               # Python 2.7x
    import Tkinter as tk
    import ttk
    import tkFont
    import tkMessageBox
    import tkFileDialog
    import tkSimpleDialog
except Exception:  # Python 3.x
    import tkinter as tk
    from tkinter import ttk
    import tkinter.font as tkFont
    import tkinter.messagebox as tkMessageBox
    import tkinter.filedialog as tkFileDialog
    import tkinter.simpledialog as tkSimpleDialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from Imports.util_tk import *
from Imports.util_plot import *
from Imports.vriCalc import *


#-----------------------------------------------------------------------------#
class App:
    """
    Class defining the Virtual Interferometer application.

    This class creates a root window to show the inputs and outputs of the 
    virtual interferometer and a secondary TopLevel window, used to choose and
    accumulate array configurations, and set observing parameters.
    
    """

    def __init__(self, root):
        self.root = root
        self.root.title("The Friendly Virtual Interferometer")
        self.root.geometry(geometryMainWin)
        self.root.resizable(True, True)
        self.root.protocol("WM_DELETE_WINDOW", self.applicationExit)
        self.obsManager = None
        
        # Create the model image panel
        self.modelImgFrm = ModelImageFrame(self.root)
        self.modelImgFrm.grid(row=0, column=0, rowspan=2, padx=5, pady=5,
                              sticky="NS")
        self.modelImgFrm.rowconfigure(0, weight=1)
        
        # Create the model FFT panel
        self.modelFFTfrm = ModelFFTframe(self.root)
        self.modelFFTfrm.grid(row=2, column=0, padx=5, pady=5, sticky="NSEW")
        
        # Create the uv-coverage panel
        self.uvCovFrm = uvCoverageFrame(self.root)
        self.uvCovFrm.grid(row=1, column=1, padx=5, pady=5, sticky="SEW")
        
        # Create the observed FFT panel
        self.obsFFTfrm = ObservedFFTframe(self.root)
        self.obsFFTfrm.grid(row=2, column=1, padx=5, pady=5, sticky="NSEW")

        # Create the synthesised beam panel
        self.beamFrm = SynthBeamFrame(self.root)
        self.beamFrm.grid(row=1, column=2, padx=5, pady=5, sticky="SEW")

        # Create the observed image panel
        self.obsImgFrm = ObservedImageFrame(self.root)
        self.obsImgFrm.grid(row=2, column=2, padx=5, pady=5, sticky="NSEW")

        # Create the status panel
        self.statusFrm = StatusFrame(self.root)
        self.statusFrm.grid(row=0, column=1, columnspan=2, padx=5, pady=5,
                          sticky="NSEW")
        self.statusFrm.columnconfigure(0, weight=1)
        
        # Create the array chooser window and set the focus back to root
        self.aConfWin = tk.Toplevel(self.root)
        self.aConfWin.title("Array Configuration Chooser")
        self.aConfWin.geometry(geometryAConfWin)
        self.aConfWin.resizable(True, True)
        self.aConfWin.protocol("WM_DELETE_WINDOW", self.applicationExit)

        # Create the array chooser frame
        self.chooserFrame = ArrayChooserFrame(self.aConfWin)
        self.chooserFrame.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")

        # Load the back-end and populate the array configuration list
        self.obsManager = observationManager()
        vals = self.obsManager.arrsAvailable.values()
        configLst = zip([x["telescope"] for x in vals],
                        [x["config"] for x in vals])
        
        self.chooserFrame.configInTab.name_columns(
            ("Telescope", "Array"))
        self.chooserFrame.configInTab.insert_rows(configLst,
                                                   ("Telescope", "Array"))

        # Populate the observing parameters widgets
        self.chooserFrame.freq_MHz.set(1420.0)
        self.chooserFrame.sampRt_s.set(300)
        self.chooserFrame.dec_deg.set(5.5)
        self.modelImgFrm.pixScale_asec.set(0.5)
        self.modelImgFrm.modelFile.set("models/radio_galaxy.png")
        
        # Bind virtual events generated by sub-widgets
        self.aConfWin.bind("<<config_in_selected>>",
                           lambda event : self._on_select_config(event))
        self.aConfWin.bind("<<selection_changed>>",
                           lambda event : self._on_sel_change(event))
        self.aConfWin.bind("<<plot_uvcoverage>>",
                           lambda event : self._on_plot_uvcov(event))
        self.aConfWin.bind("<<plot_elevation>>",
                           lambda event : self._on_plot_elevation(event))
        self.root.bind("<<do_observation>>",
                           lambda event : self._on_do_observation(event))
        self.root.bind("<<load_model_image>>",
                       lambda event : self._on_load_model(event))
        
        # Force a minimum size on the windows
        self.root.update()
        self.root.minsize(self.root.winfo_width(), 
                          self.root.winfo_height())
        self.aConfWin.update()
        self.aConfWin.minsize(self.aConfWin.winfo_width(),
                            self.aConfWin.winfo_height())
        
        # Raise the config window and start the event loop
        self.aConfWin.after(1, lambda: self.aConfWin.focus_force())
        self.aConfWin.after(1, lambda: self.aConfWin.lift())
        self.root.mainloop()
        
    def applicationExit(self):
        """Exit the application cleanly."""
        
        self.root.destroy()

    # Event handlers bound to virtual events ---------------------------------#
    #
    # _on_select_config 
    # _on_sel_change
    # _on_plot_uvcov
    # _on_plot_elevation
    # _on_load_model
    #
    #-------------------------------------------------------------------------#
    
    def _on_select_config(self, event=None):
        """Plot the antenna layout for the selected array configuration"""
        
        # Fetch the antenna layout of the selected configuration
        row = event.widget.get_indx_selected()
        xArr, yArr = self.obsManager.get_ant_coordinates(row)
        
        # Plot the antenna positions 
        self.chooserFrame.antPosPlot.load_data(xArr/1000.0, yArr/1000.0)
        self.chooserFrame.antPosPlot.draw_zerolines()
        self.chooserFrame.antPosPlot.set_xlabel("East-West (km)")
        self.chooserFrame.antPosPlot.set_ylabel("North-South (km)")

    def _on_sel_change(self, event=None):
        self.obsManager.set_obs_parms(self.chooserFrame.freq_MHz.get(),
                                      self.chooserFrame.sampRt_s.get(),
                                      self.chooserFrame.dec_deg.get())
        self.obsManager.clear_selection()
        
    def _on_plot_uvcov(self, event=None):
        """Plot the uv-coverage for all selected array configurations"""

        # Clear any previous selection in the observationManager
        self.obsManager.clear_selection()
        
        # Re-do the selection based on the arrays and HAs in the GUI
        for selection in self.chooserFrame.configOutTab.get_all_text():
            key = "_".join(selection[:2])
            haStart = float(selection[2])
            haEnd = float(selection[3])
            self.obsManager.select_array(key, haStart, haEnd)
            
        # Calculate the uv-coverage for the selected observation parameters
        self.obsManager.calc_selected_uvcoverage()

        # Plot the uv-coverage as an external MPL figure
        #fig = plt.figure(figsize=(10,10))
        #ax = fig.add_subplot(111)
        #plot_uvcov_ax(ax, self.obsManager.arrsSelected)    
        #fig.show()

        # ALT: plot in the main window frame
        ax = self.uvCovFrm.pltFrame.add_axis()
        plot_uvcov_ax(ax, self.obsManager.arrsSelected)
        self.uvCovFrm.pltFrame.show()
        
    def _on_plot_elevation(self, event=None):
        # DEBUGING 1
        #for k, v in self.obsManager.arrsAvailable.items():
        #    print
        #    print k
        #    for k1, v1 in v.items():
        #        print "\t", k1, v1

        # DEBUGGING 2
        #print "\nSELECTED:\n"
        #for entry in self.obsManager.arrsSelected:
        #    for key, val in entry.items():
        #        print "\t", key, "=", val
        #    print
        #print

        # DEBUG 3
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        for selection in self.chooserFrame.configOutTab.get_all_text():
            key = "_".join(selection[:2])
            ar =  self.obsManager.arrsAvailable[key]["antArray"]
            ax.scatter(x=ar.Bx_m/1000, y=ar.By_m/1000, marker="o",
                       edgecolor='k', s=15)
            ax.set_xlabel(u"x baseline (m)")
            ax.set_ylabel(u"y baseline (m)")
            ax.set_aspect('equal', 'datalim')
            ax.margins(0.02)
        fig.show()

    def _on_load_model(self, event=None):

        # Load the model into the observationManager
        modelFile = self.modelImgFrm.modelFile.get()
        pixScale_asec = self.modelImgFrm.pixScale_asec.get()
        self.obsManager.load_model_image(modelFile, pixScale_asec)

        # Calculate the FFT of the model
        self.obsManager.invert_model()        
        
        # Plot the model image
        ax = self.modelImgFrm.pltFrame.add_axis()
        plot_image_ax(ax, self.obsManager.modelImgArr)
        self.modelImgFrm.pltFrame.show()
        
        # Plot the model FFT
        ax = self.modelFFTfrm.pltFrame.add_axis()
        plot_fft_ax(ax, self.obsManager.modelFFTarr)
        self.modelFFTfrm.pltFrame.show()

    def _on_do_observation(self, event=None):

        # Calculate the uv-coverage if not cached.
        
        # Grid the uv-coverage to make a mask
        self.obsManager.grid_uvcoverage()
        
        # Show the observed uv-coverage
        ax = self.obsFFTfrm.pltFrame.add_axis()
        plot_fft_ax(ax, self.obsManager.obsFFTarr)
        self.obsFFTfrm.pltFrame.show()
        
        # Calculate the PSF
        self.obsManager.calc_beam()
        
        # Show the observed beam
        ax = self.beamFrm.pltFrame.add_axis()
        plot_image_ax(ax, self.obsManager.beamArr)
        self.beamFrm.pltFrame.show()

        # Invert the masked FFT
        self.obsManager.invert_observation()
        
        # Show the observed image
        ax = self.obsImgFrm.pltFrame.add_axis()
        plot_image_ax(ax, self.obsManager.obsImgArr)
        self.obsImgFrm.pltFrame.show()
        
        
#-----------------------------------------------------------------------------#
class ArrayChooserFrame(tk.Frame):
    """Frame presenting an interface to choose one or more array contigurations
    and set the observing parameters."""
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent

        # Array configuration chooser frame
        self.chooseFrm = ttk.LabelFrame(self,
                                    text="  Array Configurations  ")
        self.chooseFrm.grid(column=0, row=0, padx=5, pady=5, rowspan=4,
                            sticky="NSEW")
        
        # Canvas illustrating the antenna positions of selected array
        self.antPosPlot = ScatterPlot(self.chooseFrm, width=400, height=400,
                                      axPad=(75,25))
        self.antPosPlot.grid(column=0, row=0, padx=5, pady=5)
        
        # Listbox showing the available array configurations
        self.configInTab = ScrolledTreeTab(self.chooseFrm,
                                            virtEvent="<<config_in_selected>>")
        self.configInTab.grid(column=1, row=0, padx=5, pady=5, sticky="NSEW")
        
        # Create hour angle slider
        self.haLab = ttk.Label(self,
                                 text="Hour Angle Range:")
        self.haLab.grid(column=1, row=0, padx=5, pady=15, sticky="NEW")
        self.haScale = DoubleScale(self, from_=-12.0, to=12.0,
                                   tickIntMajor=6, tickIntMinor=2, width=270)
        self.haScale.grid(column=1, row=1, padx=0, pady=0, sticky="NEW")
        
        # Fancy add button with strike-through arrow
        self.addButFrm = tk.Frame(self)
        self.addButFrm.grid(column=1, row=2, padx=5, pady=15,sticky="NS")
        self.rowconfigure(2, weight=1)
        self.addButFrm.columnconfigure(0, weight=1)
        self.addButFrm.columnconfigure(2, weight=1)
        self.arrowL = tk.Canvas(self.addButFrm, width=70, height=20)
        self.arrowL.create_line(0,10,70,10, width=2)
        self.arrowL.grid(column=0, row=0, padx=0, pady=0, sticky="E")
        self.arrowR = tk.Canvas(self.addButFrm, width=70, height=20)
        self.arrowR.create_line(0,10,70,10, width=2,arrow=tk.LAST,
                                arrowshape=(15,20,7))
        self.arrowR.grid(column=2, row=0, padx=0, pady=0, sticky="W")
        self.addBtn = ttk.Button(self.addButFrm, text="\nADD\n", width=10,
                                 command=self._handler_add_button)
        self.addBtn.grid(column=1, row=0, padx=0, pady=0, sticky="S" )
        
        # Create the observing parameter frame
        self.obsParmFrm = ttk.LabelFrame(self,
                                    text="  Common Observing Parameters  ")
        self.obsParmFrm.grid(column=1, row=3, padx=5, pady=5, sticky="SEW")
        #self.obsParmFrm.rowconfigure(3, weight=1)

        # Frequency entry box
        self.freqLab = ttk.Label(self.obsParmFrm, text="Frequency (MHz):")
        self.freqLab.grid(column=0, row=0, padx=5, pady=5, sticky="W")
        self.freq_MHz = tk.StringVar()
        self.freqEnt = ttk.Entry(self.obsParmFrm, width=10,
                                 textvariable=self.freq_MHz)
        self.freqEnt.grid(column=1, row=0, padx=5, pady=5, sticky="E")
        
        # Sampling rate entry box
        self.sampRtLab = ttk.Label(self.obsParmFrm, text="Sampling Rate (s):")
        self.sampRtLab.grid(column=0, row=1, padx=5, pady=5, sticky="W")
        self.sampRt_s = tk.StringVar()
        self.sampRtEnt = ttk.Entry(self.obsParmFrm, width=10,
                                   textvariable=self.sampRt_s)
        self.sampRtEnt.grid(column=1, row=1, padx=5, pady=5, sticky="E")
        
        # Source declination entry box
        self.decSrcLab = ttk.Label(self.obsParmFrm,
                                   text="Source Declination (deg):")
        self.decSrcLab.grid(column=0, row=3, padx=5, pady=5, sticky="W")
        self.dec_deg = tk.StringVar()
        self.decSrcEnt = ttk.Entry(self.obsParmFrm, width=10,
                                   textvariable=self.dec_deg)
        self.decSrcEnt.grid(column=1, row=3, padx=5, pady=5, sticky="E")

        # Accumulated array configuration frame
        self.selectedFrm = ttk.LabelFrame(self,
                                    text="  Selected Configurations  ")
        self.selectedFrm.grid(column=2, row=0, padx=5, pady=5, rowspan=4,
                              sticky="NSEW")
        self.selectedFrm.rowconfigure(0, weight=1)
        self.selectedFrm.columnconfigure(0, weight=1)
        
        # Listbox showing the selected array configurations
        self.configOutTab = ScrolledTreeTab(self.selectedFrm,
                                        virtEvent="<<config_out_selected>>")
        self.configOutTab.name_columns(("Telescope", "Array", "Start", "End"))
        self.configOutTab.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW")

        # Delete buttons
        self.delBtn = ttk.Button(self.selectedFrm, text="Clear Selected",
                                 command=self._handler_clear_button)
        self.delBtn.grid(column=0, row=1, padx=5, pady=5, sticky="SEW" )
        self.delAllBtn = ttk.Button(self.selectedFrm, text="Clear All",
                                    command=self._handler_clear_all_button)
        self.delAllBtn.grid(column=0, row=2, padx=5, pady=5, sticky="SEW" )

        # Plot uv-coverage button
        self.pltuvBtn = ttk.Button(self.selectedFrm, text="Update uv-Coverage",
                                   command=self._handler_plt_uvcov)
        self.pltuvBtn.grid(column=0, row=3, padx=5, pady=5, sticky="SEW" )
        
        # Plot elevation button
        self.pltElBtn = ttk.Button(self.selectedFrm, text="Plot Elevation",
                                   command=self._handler_plt_elevation)
        self.pltElBtn.grid(column=0, row=4, padx=5, pady=5, sticky="SEW" )
        
    def _handler_add_button(self):
        """Add the selected configuration to the list box"""
        selConf = self.configInTab.get_text_selected()
        if selConf is not None:
            selConf = list(selConf)
            selConf.append(str(self.haScale.valueLeft))
            selConf.append(str(self.haScale.valueRight))
            self.configOutTab.insert_rows([selConf],
                                    ("Telescope", "Array", "Start", "End"))

        prevAddedLst = self.configOutTab.get_all_text()
        self.event_generate("<<selection_changed>>")
    
    def _handler_clear_button(self):
        """Delete the selected configuration from the list box"""

        self.configOutTab.clear_selected()   
        self.event_generate("<<selection_changed>>")
        
    def _handler_clear_all_button(self):
        """Delete all configurations from the list box"""
        
        self.configOutTab.clear_entries()   
        self.event_generate("<<selection_changed>>")

    def _handler_plt_uvcov(self):
        """Generate a virtual event to plot the uv-coverage"""
        
        self.event_generate("<<plot_uvcoverage>>")

    def _handler_plt_elevation(self):
        """Generate a virtual event to plot the elevation for each array"""
        
        self.event_generate("<<plot_elevation>>")
        
#-----------------------------------------------------------------------------#
class ModelImageFrame(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent

        # Model chooser frame
        self.modelFrm = ttk.LabelFrame(self, text="  Model Image  ")
        self.modelFrm.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW")
        self.modelFrm.rowconfigure(2, weight=1)

        # Input file
        self.fileLab = ttk.Label(self.modelFrm, text="File Name:")
        self.fileLab.grid(column=0, row=0, padx=5, pady=5, sticky="W")
        self.modelFile = tk.StringVar()
        self.fileEnt = ttk.Entry(self.modelFrm, width=20,
                                 textvariable=self.modelFile)
        self.fileEnt.grid(column=1, row=0, columnspan=2, padx=5, pady=5,
                          sticky="EW")
        self.browseBtn = ttk.Button(self.modelFrm, text="Browse",
                                    command=self._handler_browse_button)
        self.browseBtn.grid(column=0, row=1, padx=5, pady=2, sticky="EW")
        self.loadBtn = ttk.Button(self.modelFrm, text="Load",
                                  command=self._handler_load_button)
        self.loadBtn.grid(column=2, row=1, padx=5, pady=5, sticky="E")

        # Pixel Scale
        self.pixScaLab = ttk.Label(self.modelFrm,
                                   text="Pixel Scale (arcseconds):")
        self.pixScaLab.grid(column=0, row=2, columnspan=2, padx=5, pady=5,
                            sticky="NW")
        self.pixScale_asec = tk.DoubleVar()
        self.pixScaEnt = ttk.Entry(self.modelFrm, width=10,
                                   textvariable=self.pixScale_asec)
        self.pixScaEnt.grid(column=2, row=2, padx=5, pady=5, sticky="NE")

        # Image plot
        self.pltFrame = SingleFigFrame(self.modelFrm, width=300, height=300)
        self.pltFrame.grid(column=0, row=3, columnspan=3,
                                  padx=5, pady=5, sticky="S")
        
    def _handler_browse_button(self):
        """Open the file selection dialog and set the session dir."""
        
        modelFile = tkFileDialog.askopenfilename(parent=self,
                                    title="Choose a model file")
        if not len(modelFile)==0:
            self.modelFile.set(modelFile)
            
    def _handler_load_button(self):      
        """Raise a <<load_session>> virtual event in the parent window."""
        
        self.event_generate("<<load_model_image>>")


#-----------------------------------------------------------------------------#
class ModelFFTframe(tk.Frame):
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent

        # Model FFT
        self.modelFrm = ttk.LabelFrame(self, text="  Model FFT  ")
        self.modelFrm.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW")
        self.modelFrm.rowconfigure(0, weight=1)
        self.modelFrm.columnconfigure(0, weight=1)
        
        # Image plot
        self.pltFrame = SingleFigFrame(self.modelFrm, width=300, height=300)
        self.pltFrame.grid(column=0, row=1, padx=5, pady=5)


#-----------------------------------------------------------------------------#
class uvCoverageFrame(tk.Frame):
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent

        # uv-Coverage
        self.uvcovFrm = ttk.LabelFrame(self, text="  uv-Coverage  ")
        self.uvcovFrm.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW")
        self.uvcovFrm.rowconfigure(0, weight=1)
        self.uvcovFrm.columnconfigure(0, weight=1)
        
        # uv-coverage plot
        self.pltFrame = SingleFigFrame(self.uvcovFrm, width=300, height=300)
        self.pltFrame.grid(column=0, row=1, padx=5, pady=5)

        
#-----------------------------------------------------------------------------#
class ObservedFFTframe(tk.Frame):
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent

        # Observed FFT
        self.obsFFTfrm = ttk.LabelFrame(self, text="  Observed FFT  ")
        self.obsFFTfrm.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW")
        self.obsFFTfrm.rowconfigure(0, weight=1)
        self.obsFFTfrm.columnconfigure(0, weight=1)
        
        # Observed FFT plot
        self.pltFrame = SingleFigFrame(self.obsFFTfrm, width=300, height=300)
        self.pltFrame.grid(column=0, row=1, padx=5, pady=5)


#-----------------------------------------------------------------------------#
class SynthBeamFrame(tk.Frame):
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent

        # Synthesised Beam
        self.synthBeamFrm = ttk.LabelFrame(self, text="  Synthesised Beam  ")
        self.synthBeamFrm.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW")
        self.synthBeamFrm.rowconfigure(0, weight=1)
        self.synthBeamFrm.columnconfigure(0, weight=1)
        
        # Synthesised Beam plot
        self.pltFrame = SingleFigFrame(self.synthBeamFrm,
                                        width=300, height=300)
        self.pltFrame.grid(column=0, row=1, padx=5, pady=5)

        
#-----------------------------------------------------------------------------#
class ObservedImageFrame(tk.Frame):
    
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent

        # Synthesised Beam
        self.obsImgFrm = ttk.LabelFrame(self, text="  Observed Image  ")
        self.obsImgFrm.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW")
        self.obsImgFrm.rowconfigure(0, weight=1)
        self.obsImgFrm.columnconfigure(0, weight=1)
        
        # Synthesised Beam plot
        self.pltFrame = SingleFigFrame(self.obsImgFrm,
                                        width=300, height=300)
        self.pltFrame.grid(column=0, row=1, padx=5, pady=5)


#-----------------------------------------------------------------------------#
class StatusFrame(tk.Frame):
        
    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent

        # Container
        self.statFrm = ttk.LabelFrame(self, text="  Status  ")
        self.statFrm.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW")
        self.statFrm.columnconfigure(1, weight=1)
        self.statFrm.columnconfigure(2, weight=1)

        # Model loaded
        self.modelLab = ttk.Label(self.statFrm, justify="center",
                                  foreground="black",
                                  background="orange",
                                  text="Not Loaded", relief="solid",
                                  anchor="center", padding=5, width=15)
        self.modelStat = ttk.Label(self.statFrm, justify="center",
                                   text="Model Image:",
                                   anchor="center", padding=0)
        self.modelLab.grid(row=1, column=0, padx=15, pady=5, sticky="NS")
        self.modelStat.grid(row=0, column=0, padx=15, pady=5, sticky="NS")

        # Array loaded
        self.arrayLab = ttk.Label(self.statFrm, justify="center",
                                  foreground="black",
                                  background="orange",
                                  text="Not Loaded", relief="solid",
                                  anchor="center", padding=5, width=15)
        self.arrayStat = ttk.Label(self.statFrm, justify="center",
                                   text="Array Configurations:",
                                   anchor="center", padding=0)
        self.arrayLab.grid(row=1, column=1, padx=15, pady=5, sticky="NS")
        self.arrayStat.grid(row=0, column=1, padx=15, pady=5, sticky="NS")

        # Observe Button
        self.observeBtn = ttk.Button(self.statFrm, text="\nDO OBSERVATION\n",
                                     width=20,
                                     command=self._handler_observe_button)
        self.observeBtn.grid(row=0, column=2, rowspan=2, padx=5, pady=5,
                             sticky="NSEW")
        
    def _handler_observe_button(self):
        """Run the observation"""
        
        self.event_generate("<<do_observation>>")

        
#-----------------------------------------------------------------------------#
if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    
