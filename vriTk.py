#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     vriTK.py                                                          #
#                                                                             #
# PURPOSE:  A virtual interferometer application written in Tkinter.          #
#                                                                             #
# REQUIRED: Requires numpy, tkinter, matplotlib                               #
#                                                                             #
# MODIFIED: 13-Mar-2017 by cpurcell                                           #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# App                ... class containing the main application logic          #
# PlotFrame          ... class defining the plotting window                   #
# ObsControlFrame    ... class defining the observation controller interface  #
# StatusFrame        ... class defining status indicators and action buttons  #
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
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from Imports.util_tk import *
from Imports.util_plot import *
from Imports.vriCalc import *


#-----------------------------------------------------------------------------#
class App(ttk.Frame):
    """Class defining the Virtual Interferometer application.

    This class creates a root window to show the inputs and outputs of the 
    virtual interferometer and a secondary TopLevel window, used to choose and
    accumulate array configurations, and set observing parameters.
    
    """

    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.parent.title("Friendly VRI: Control Window")
        self.obsManager = None
        
        # Create the display window and set the focus back to root
        self.displayWin = tk.Toplevel(self)
        self.displayWin.title("Friendly VRI: Display Window")
        self.displayWin.protocol("WM_DELETE_WINDOW", self.applicationExit)
        self.displayWin.columnconfigure(0, weight=1)
        self.displayWin.rowconfigure(0, weight=1)
        
        # Create the display interface
        self.pltFrm = PlotFrame(self.displayWin)
        self.pltFrm.grid(row=0, column=0, padx=0, pady=0, sticky="NSEW")
        
        # Create the observation controller interface
        self.ctrlFrm = ObsControlFrame(self)
        self.ctrlFrm.grid(row=0, column=0, padx=0, pady=0, sticky="NSEW")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        # DEBUG
        #self.testWin = tk.Toplevel(self)
        #self.testWin.title(" TEST WINDOW ")
        #self.testWin.protocol("WM_DELETE_WINDOW", self.applicationExit)
        #self.testWin.resizable(True, True)
        #self.testWin.columnconfigure(0, weight=1)
        #self.testWin.rowconfigure(0, weight=1)
        #self.a =  ScaleInformation(self.testWin)
        #self.a.grid(row=0, column=0, padx=0, pady=0)
        
        # Load the back-end and populate the array configuration list
        self.obsManager = observationManager(verbose=True, debug=True)
        vals = self.obsManager.arrsAvailable.values()
        configLst = zip([x["telescope"] for x in vals],
                        [x["config"] for x in vals])        
        #self.ctrlFrm.configInTab.name_columns(("Telescope", "Array"))
        self.ctrlFrm.selector.configInTab.insert_rows(
            configLst,("Telescope", "Array"))

        # Populate the observing parameters widgets with defaults
        self.ctrlFrm.obsParms.freq_MHz.set(1420.0)
        self.ctrlFrm.obsParms.sampRt_s.set(300)
        self.ctrlFrm.obsParms.dec_deg.set(20.0)
        self.ctrlFrm.modelParms.pixScale_asec.set(1.0)
        self.ctrlFrm.modelParms.modelPath = "models/radio_galaxy.png"
        self.ctrlFrm.modelParms.modelFile.set("radio_galaxy.png")
        
        
        # Bind virtual events generated by the control widgets
        self.parent.bind("<<config_in_selected>>",
                  lambda event : self._on_select_config(event))
        self.parent.bind("<<selection_changed>>",
                  lambda event : self._on_sel_change(event))
        self.parent.bind("<<plot_modFFT>>",
                  lambda event : self._on_plot_modFFT(event))
        self.parent.bind("<<plot_uvcoverage>>",
                  lambda event : self._on_plot_uvcov(event))
        self.parent.bind("<<plot_elevation>>",
                  lambda event : self._on_plot_elevation(event))
        self.parent.bind("<<do_observation>>",
                  lambda event : self._on_do_observation(event))
        self.parent.bind("<<load_model_image>>",
                       lambda event : self._on_load_model(event))
        
        # Force a minimum size on the windows
        self.parent.update()
        self.parent.minsize(self.parent.winfo_width(),
                            self.parent.winfo_height())
        self.parent.resizable(True, True)
        self.displayWin.update()
        self.displayWin.minsize(self.displayWin.winfo_width(),
                                self.displayWin.winfo_height())
        self.displayWin.resizable(True, True)
        
    def applicationExit(self):
        """Exit the application cleanly if the display window is closed."""
        
        self.parent.destroy()
        
    # Event handlers bound to virtual events ---------------------------------#
    #
    # _on_select_config 
    # _on_sel_change
    # _on_plot_modFFT
    # _on_plot_uvcov
    # _on_plot_elevation
    # _on_load_model
    # _on_do_observation
    #
    #-------------------------------------------------------------------------#
    
    def _on_select_config(self, event=None):
        """Plot the antenna layout for the selected array configuration"""
        
        # Fetch the antenna layout of the selected configuration
        row = event.widget.get_indx_selected()
        xArr, yArr = self.obsManager.get_ant_coordinates(row=row)
        
        # Plot the antenna positions 
        self.ctrlFrm.selector.antPosPlot.load_data(xArr/1000.0, yArr/1000.0)
        self.ctrlFrm.selector.antPosPlot.draw_zerolines()
        self.ctrlFrm.selector.antPosPlot.set_xlabel("East-West (km)")
        self.ctrlFrm.selector.antPosPlot.set_ylabel("North-South (km)")

        # Update the antenna diameter and array latitude
        d = self.obsManager.get_telescope_diameter(row=row)
        self.ctrlFrm.selector.antD_m.set(str(d) + " m" )
        lat = self.obsManager.get_telescope_latitude(row=row)
        self.ctrlFrm.selector.antL_deg.set(str(lat) + u"\u00B0")
        
    def _on_sel_change(self, event=None):
        """When the selection changes in the GUI, refresh the observation
        parameters in the observation manager."""

        # Reset the common parameters
        self.obsManager.set_obs_parms(self.ctrlFrm.obsParms.freq_MHz.get(),
                                      self.ctrlFrm.obsParms.sampRt_s.get(),
                                      self.ctrlFrm.obsParms.dec_deg.get())

        # Reset the array an hour-angle selection
        self.obsManager.clear_all_selections()
        for selection in self.ctrlFrm.selector.configOutTab.get_all_text():
            key = "_".join(selection[:2])
            haStart = float(selection[2])
            haEnd = float(selection[3])
            self.obsManager.select_array(key, haStart, haEnd)
        
        # Update the status indicators
        stateDict = self.obsManager.get_status()
        self.ctrlFrm.statFrm.set_state_dict(stateDict)
        
    def _on_plot_modFFT(self, event=None):
        """Show the FFT of the model image."""

        # Invert the model image
        self.obsManager.invert_model()
        
        # Plot the model FFT
        ax = self.pltFrm.modelFFTfrm.add_axis()
        plot_fft_ax(ax, self.obsManager.modelFFTarr)
        self.pltFrm.modelFFTfrm.show()
        
        # Update the status indicators
        stateDict = self.obsManager.get_status()
        self.ctrlFrm.statFrm.set_state_dict(stateDict)
        
    def _on_plot_uvcov(self, event=None):
        """Plot the uv-coverage for all selected array configurations"""
        
        # Calculate the uv-coverage for the selected observation parameters
        self.obsManager.calc_uvcoverage()

        # Plot the uv-coverage in the display window
        ax = self.pltFrm.uvCovFrm.add_axis()
        plot_uvcov_ax(ax, self.obsManager.arrsSelected)
        self.pltFrm.uvCovFrm.show()
        
        # ALT: Plot the uv-coverage as an external MPL figure
        #fig = plt.figure(figsize=(10,10))
        #ax = fig.add_subplot(111)
        #plot_uvcov_ax(ax, self.obsManager.arrsSelected)    
        #fig.show()
        
        # Update the status indicators
        stateDict = self.obsManager.get_status()
        self.ctrlFrm.statFrm.set_state_dict(stateDict)
        
    def _on_plot_elevation(self, event=None):
        """Create a plot showing the elevation of the source as seen from 
        the currently selected telescopes."""
        
        # DEBUG
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111)
        for selection in self.ctrlFrm.selector.configOutTab.get_all_text():
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
        """Load a model image into the observation manager"""
        
        # Load the model image
        modelFile = self.ctrlFrm.modelParms.modelPath
        pixScale_asec = self.ctrlFrm.modelParms.pixScale_asec.get()
        self.obsManager.load_model_image(modelFile, pixScale_asec)
        
        # Plot the model image
        ax = self.pltFrm.modelImgFrm.add_axis()
        plot_image_ax(ax, self.obsManager.modelImgArr)
        self.pltFrm.modelImgFrm.show()
        
        # Update the status indicators
        stateDict = self.obsManager.get_status()
        self.ctrlFrm.statFrm.set_state_dict(stateDict)

        # Change the icon of the load model button
        self.ctrlFrm.modelParms._change_icon()
        
    def _on_do_observation(self, event=None):
        """Perform the bulk of the observing steps"""
        
        # Calculate the Fourier transform of the model if not cached
        stateDict = self.obsManager.get_status()
        if not stateDict['statusModelFFT']:
            self.obsManager.invert_model()
        
            # Plot the model FFT
            ax = self.pltFrm.modelFFTfrm.add_axis()
            plot_fft_ax(ax, self.obsManager.modelFFTarr)
            self.pltFrm.modelFFTfrm.show()
        
            # Update the status indicators
            stateDict = self.obsManager.get_status()
            self.ctrlFrm.statFrm.set_state_dict(stateDict)
        
        # Calculate the uv-coverage if not cached
        stateDict = self.obsManager.get_status()
        if not stateDict['statusuvCalc']:
            self.obsManager.calc_uvcoverage()

            # Plot the uv-coverage
            ax = self.pltFrm.uvCovFrm.add_axis()
            plot_uvcov_ax(ax, self.obsManager.arrsSelected)
            self.pltFrm.uvCovFrm.show()

            # Update the status indicators
            stateDict = self.obsManager.get_status()
            self.ctrlFrm.statFrm.set_state_dict(stateDict)
            
        # Grid the uv-coverage to make a mask
        self.obsManager.grid_uvcoverage()
        
        # Update the status indicators
        stateDict = self.obsManager.get_status()
        self.ctrlFrm.statFrm.set_state_dict(stateDict)
        
        # Calculate the PSF
        self.obsManager.calc_beam()
        
        # Show the synthesised beam
        ax = self.pltFrm.beamFrm.add_axis()
        plot_image_ax(ax, self.obsManager.beamArr)
        self.pltFrm.beamFrm.show()
        
        # Update the status indicators
        stateDict = self.obsManager.get_status()
        self.ctrlFrm.statFrm.set_state_dict(stateDict)
        
        # Apply the gridded uv-coverage and invert
        self.obsManager.invert_observation()
        
        # Show the observed FFT
        ax = self.pltFrm.obsFFTfrm.add_axis()
        plot_fft_ax(ax, self.obsManager.obsFFTarr)
        self.pltFrm.obsFFTfrm.show()
        
        # Show the observed image
        ax = self.pltFrm.obsImgFrm.add_axis()
        plot_image_ax(ax, self.obsManager.obsImgArr)
        self.pltFrm.obsImgFrm.show()
        
        # Update the status indicators
        stateDict = self.obsManager.get_status()
        self.ctrlFrm.statFrm.set_state_dict(stateDict)
        

#-----------------------------------------------------------------------------#
class PlotFrame(ttk.Frame):
    """Frame showing the plots produced by the virtual interferometer."""
    
    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)

        # Set the expansion
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)
        
        # Create the model image panel
        self.modelImgFrm = SingleFigFrame(self)
        self.modelImgFrm.grid(column=0, row=0, columnspan=1, padx=5, pady=5,
                              sticky="NSEW")

        # Create the model FFT panel
        self.modelFFTfrm = SingleFigFrame(self)
        self.modelFFTfrm.grid(column=0, row=1, padx=5, pady=5, sticky="NSEW")
        
        # Create the uv-coverage panel
        self.uvCovFrm = SingleFigFrame(self)
        self.uvCovFrm.grid(column=1, row=0, padx=5, pady=5, sticky="NSEW")
        
        # Create the observed FFT panel
        self.obsFFTfrm = SingleFigFrame(self)
        self.obsFFTfrm.grid(column=1, row=1, padx=5, pady=5, sticky="NSEW")

        # Create the synthesised beam panel
        self.beamFrm = SingleFigFrame(self)
        self.beamFrm.grid(column=2, row=0,  padx=5, pady=5, sticky="NSEW")

        # Create the observed image panel
        self.obsImgFrm = SingleFigFrame(self)
        self.obsImgFrm.grid(column=2, row=1,  padx=5, pady=5, sticky="NSEW")

        
#-----------------------------------------------------------------------------#
class ObsControlFrame(ttk.Frame):
    """Frame presenting an interface to control the virtual radio 
    interferometer."""
    
    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        # Array selector interface
        self.selector =  ArraySelector(self)
        self.selector.grid(column=0, row=0, columnspan=3, padx=10, pady=0,
                           sticky="NW")

        sep = ttk.Separator(self, orient="horizontal")
        sep.grid(column=0, row=1, columnspan=3, padx=0, pady=5, sticky="EW")
        
        # Common Observing parameters
        self.obsParms = ObsSettings(self)
        self.obsParms.grid(column=0, row=2, padx=10, pady=0, sticky="NW")
        
        sep = ttk.Separator(self, orient="vertical")
        sep.grid(column=1, row=2, columnspan=1, padx=5, pady=0, sticky="NS")

        # Model chooser frame
        self.modelParms = ModelChooser(self)
        self.modelParms.grid(row=2, column=2, padx=5, pady=0, sticky="NW")

        sep = ttk.Separator(self, orient="horizontal")
        sep.grid(column=0, row=3, columnspan=3, padx=5, pady=10, sticky="EW")

        # uv-coverage and model FFT information
        self.scaleInfo = ScaleInformation(self)
        self.scaleInfo.grid(column=0, row=4, columnspan=3, padx=5, pady=10,
                            sticky="EW")

        sep = ttk.Separator(self, orient="horizontal")
        sep.grid(column=0, row=5, columnspan=3, padx=5, pady=10, sticky="EW")
        
        # New Status frame
        self.statFrm = StatusFrame(self, boxWidth=20, gapWidth=150)
        self.statFrm.grid(column=0, row=6, columnspan=4, padx=5, pady=5)
        
    def _handler_observe_button(self):
        """Run the observation"""
        
        self.event_generate("<<do_observation>>")
    
    def _handler_plt_uvcov(self):
        """Generate a virtual event to plot the uv-coverage"""
        
        self.event_generate("<<plot_uvcoverage>>")

    def _handler_plt_elevation(self):
        """Generate a virtual event to plot the elevation for each array"""
        
        self.event_generate("<<plot_elevation>>")

        
#-----------------------------------------------------------------------------#
class ArraySelector(ttk.Frame):
    """Two multi-column listboxes an hour-angle slider and a 'Select' button
    that make up the array selection interface."""
    
    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        # Set the expansion properties
        self.columnconfigure(2, weight=1)
        self.columnconfigure(6, weight=2)
        self.columnconfigure(7, weight=2)
        self.rowconfigure(3, weight=2)
        
        # Scatter plot of antenna locations
        self.antPosPlot = ScatterPlot(self, width=340, height=300,
                                      axPad=(100,25,70,25), aspect="equal",
                                      pntSize=2.5)
        self.antPosPlot.grid(column=0, row=0, rowspan=6,
                             padx=5, pady=5)

        # Listbox showing the available array configurations
        self.configInTab = ScrolledTreeTab(self,
                                           virtEvent="<<config_in_selected>>")
        self.configInTab.name_columns(("Telescope", "Array"))
        self.configInTab.grid(column=1, row=0, columnspan=4, rowspan=4,
                              padx=5, pady=5, sticky="NSEW")
        
        # Antenna Diameter and array latitude labels
        self.antDlab1 = ttk.Label(self, text=u"Antenna \u0398:   ")
        self.antDlab1.grid(column=1, row=4, padx=0, pady=0, sticky="E")
        self.antD_m = tk.StringVar()
        self.antDlab2 = ttk.Label(self, textvariable=self.antD_m)
        self.antDlab2.grid(column=2, row=4, padx=0, pady=0, sticky="W")
        self.antLlab1 = ttk.Label(self, text=u"Telescope \u03c6:   ")
        self.antLlab1.grid(column=1, row=5, padx=0, pady=0, sticky="E")
        self.antL_deg = tk.StringVar()
        self.antLlab2 = ttk.Label(self, textvariable=self.antL_deg)
        self.antLlab2.grid(column=2, row=5, padx=0, pady=0, sticky="W")
        
        # Hour angle slider
        self.haLab = ttk.Label(self, text="Hour Angle Range (hours):")
        self.haLab.grid(column=5, row=0, padx=0, pady=5, sticky="NW")
        self.haScale = DoubleScale(self, from_=-12.0, to=12.0,
                                   initLeft=-1.0, initRight=+1.0,
                                   tickIntMajor=6, tickIntMinor=1, width=270)
        self.haScale.grid(column=5, row=1, padx=0, pady=5, sticky="N")

        # Fancy add button with strike-through arrow
        bgColour = ttk.Style().lookup("TFrame", "background")
        self.canvas = tk.Canvas(self, width=270, height=30,
                               background=bgColour, highlightthickness=0)
        self.canvas.create_line(10,15,260,15, width=2, arrow=tk.LAST,
                                arrowshape=(10,15,5), fill="black")
        self.addBtn = ttk.Button(self, text="Add", width=10,
                                 command=self._handler_add_button)
        self.canvas.create_window(135, 15, window=self.addBtn)
        self.canvas.grid(column=5, row=2, padx=0, pady=5, sticky="NS")

        # Listbox showing the selected array configurations
        self.configOutTab = ScrolledTreeTab(self,
                                        virtEvent="<<config_out_selected>>")
        self.configOutTab.name_columns(("Telescope", "  Array  ",
                                        "HA-Start", " HA-End "))
        self.configOutTab.grid(column=6, row=0, columnspan=2, rowspan=4,
                               padx=5, pady=5, sticky="NSEW")
        
        # Delete buttons
        self.delBtn = ttk.Button(self, text="Clear Selected", width=20,
                                 command=self._handler_clear_button)
        self.delBtn.grid(column=6, row=5, padx=5, pady=5, sticky="EW" )
        self.delAllBtn = ttk.Button(self, text="Clear All", width=20,
                                    command=self._handler_clear_all_button)
        self.delAllBtn.grid(column=7, row=5, padx=5, pady=5, sticky="EW" )
                
    def _handler_add_button(self):
        """Add the selected configuration to the list box"""
        
        selConf = self.configInTab.get_text_selected()
        if selConf is not None:
            selConf = list(selConf)
            selConf.append(self.haScale.valueLeft.get())
            selConf.append(self.haScale.valueRight.get())
            self.configOutTab.insert_rows(
                [selConf], ("Telescope", "  Array  ", "HA-Start", " HA-End "))
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


#-----------------------------------------------------------------------------#
class ObsSettings(ttk.Frame):
    """Common settings for the observation (observing frequency, data sampling
    cadence, Declination of source, Robust weighting factor"""
    
    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        self.title = ttk.Label(self, text="Observation Settings:")
        self.title.grid(row=0, column=0, padx=0, pady=0, columnspan=4,
                        sticky="NW")

        # Observing frequency
        self.freqLab = ttk.Label(self,
                                 text="Observing Frequency (MHz):")
        self.freqLab.grid(column=0, row=1, padx=(15,5), pady=5, sticky="E")
        self.freq_MHz = tk.StringVar()
        self.freqEnt = ttk.Entry(self, width=10, textvariable=self.freq_MHz)
        self.freqEnt.grid(column=1, row=1, padx=5, pady=5, sticky="EW")

        # Source Declination
        self.decSrcLab = ttk.Label(self,
                                   text="Source Declination (deg):")
        self.decSrcLab.grid(column=0, row=2, padx=(15,5), pady=5, sticky="E")
        self.dec_deg = tk.StringVar()
        self.decSrcEnt = ttk.Entry(self, width=10, textvariable=self.dec_deg)
        self.decSrcEnt.grid(column=1, row=2, padx=5, pady=5, sticky="EW")

        # Sampling cadence
        self.sampRtLab = ttk.Label(self,
                                   text="Sampling Rate (s):")
        self.sampRtLab.grid(column=2, row=1, padx=(15,5), pady=5, sticky="E")
        sampRtLst_s = ["10", "60", "60", "100", "300", "600", "1200", "1800",
                       "3600"]
        self.sampRt_s = tk.StringVar()
        self.sampRtComb = ttk.Combobox(self, textvariable=self.sampRt_s,
                                       values=sampRtLst_s, width=7)
        self.sampRtComb.current(4)
        self.sampRtComb.grid(column=3, row=1, padx=5, pady=5, sticky="EW")
        
        # Weighting factor
        self.robustLab = ttk.Label(self,
                                   text="Robust Weighting Factor:")
        self.robustLab.grid(column=2, row=2, padx=(15,5), pady=5, sticky="E")
        #robustLst = ["None", "-2.0", "-1.5", "-1.0", "-0.5", "0.0",
        #             "0.5", "1.0", "1.5", "2.0"]
        robustLst = ["None"]
        self.robust = tk.StringVar()
        self.robustComb = ttk.Combobox(self, textvariable=self.robust,
                                       values=robustLst, width=15)
        self.robustComb.current(0)
        self.robustComb.grid(column=3, row=2, padx=5, pady=5, sticky="EW")

        
#-----------------------------------------------------------------------------#
class ModelChooser(ttk.Frame):
    """Interface to choose a model image file."""

    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.modelFile = None
        self.modelPath = None

        # Title
        self.title = ttk.Label(self, text="Model Image:")
        self.title.grid(row=0, column=0, padx=0, pady=0, columnspan=4,
                        sticky="W")
        
        # Input file name
        self.fileLab = ttk.Label(self, text="File:")
        self.fileLab.grid(column=0, row=1, padx=(15,5), pady=5, sticky="E")
        self.modelFile = tk.StringVar()
        self.fileEnt = ttk.Entry(self, width=30,
                                 textvariable=self.modelFile)
        self.fileEnt.configure(state="readonly")        
        self.fileEnt.grid(column=1, row=1, columnspan=1, padx=5, pady=5,
                          sticky="EW")

        # Browse and load buttons
        buttonImage = Image.open('Imports/folder_ico.png')
        browsePhoto = ImageTk.PhotoImage(buttonImage)
        self.browseBtn = ttk.Button(self, image=browsePhoto,
                                    command=self._handler_browse_button)
        self.browseBtn.grid(column=2, row=1, padx=5, pady=5, sticky="E")
        self.browseBtn.image = browsePhoto
        buttonImage = Image.open('Imports/load_ico.png')
        self.loadPhoto = ImageTk.PhotoImage(buttonImage)
        self.loadBtn = ttk.Button(self, image=self.loadPhoto,
                                  command=self._handler_load_button)
        self.loadBtn.grid(column=3, row=1, padx=5, pady=5, sticky="E")
        self.loadBtn.image = self.loadPhoto

        # Pixel scale
        self.pixScaLab = ttk.Label(self, text="Pixel Scale (arcsec):")
        self.pixScaLab.grid(column=1, row=3, columnspan=1, padx=5, pady=5,
                            sticky="E")
        self.pixScale_asec = tk.DoubleVar()
        self.pixScaEnt = ttk.Entry(self, width=3,
                                   textvariable=self.pixScale_asec)
        self.pixScaEnt.grid(column=2, row=3, columnspan=2,padx=5, pady=5,
                            sticky="EW")
        
        # Information
        #text = "Supported image formats: PNG, JPEG, TIFF."
        #self.formatInfoLab = ttk.Label(self,text=text)
        #self.formatInfoLab.grid(column=0, row=2, padx=5, pady=5,
        #                        columnspan=4, sticky="W")
        #self.columnconfigure(1, weight=1)

    def _change_icon(self):
        buttonImage = Image.open('Imports/reload_ico.png')
        self.reloadPhoto = ImageTk.PhotoImage(buttonImage)
        self.loadBtn.configure(image=self.reloadPhoto)
        self.loadBtn.image = self.reloadPhoto
        
    def _handler_browse_button(self):
        """Open the file selection dialog and set the session dir."""
        
        modelPath = tkFileDialog.askopenfilename(parent=self,
                            initialdir="models", title="Choose a model file")
        if not len(modelPath)==0:
            if os.path.exists(modelPath):
                s = os.path.split(modelPath)
                if len(s)==2:
                    self.modelFile.set(s[-1])
                    self.modelPath = modelPath
        self.event_generate("<<load_model_image>>")
        
    def _handler_load_button(self):      
        """Raise a <<load_session>> virtual event in the parent window."""

        self.event_generate("<<load_model_image>>")

        
#-----------------------------------------------------------------------------#
class ScaleInformation(ttk.Frame):
    """Frame presenting the parameters of the selected arrays and input model.
    """

    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        self.title = ttk.Label(self, text="Properties of uv-Coverage:")
        self.title.grid(row=0, column=0, padx=0, pady=0, columnspan=2,
                        sticky="NW")

        # Resolution, maximum scale and field of view
        self.resLab = ttk.Label(self, text=u"Resolution (arcsec):")
        self.resLab.grid(column=0, row=1, padx=(15,5), pady=5, sticky="E")
        self.resolution_asec = tk.StringVar()
        self.resVal = ttk.Label(self, width=20,
                                textvariable=self.resolution_asec)
        self.resVal.grid(column=1, row=1, padx=5, pady=5, sticky="W")
        self.maxScaleLab = ttk.Label(self, text=u"Maximum Scale (arcmin):")
        self.maxScaleLab.grid(column=0, row=2, padx=(15,5), pady=5, sticky="E")
        self.maxScale_arcmin = tk.StringVar()
        self.maxScaleVal = ttk.Label(self, width=20,
                                     textvariable=self.maxScale_arcmin)
        self.maxScaleVal.grid(column=1, row=2, padx=5, pady=5, sticky="W")
        self.fovLab = ttk.Label(self, text=u"Field of View: (arcmin):")
        self.fovLab.grid(column=0, row=3, padx=(15,5), pady=5, sticky="E")
        self.fov_arcmin = tk.StringVar()
        self.fovVal = ttk.Label(self, width=20,
                                textvariable=self.fov_arcmin)
        self.fovVal.grid(column=1, row=3, padx=5, pady=5, sticky="W")

        
        sep = ttk.Separator(self, orient="vertical")
        sep.grid(column=2, row=0, rowspan=4, padx=5, pady=0, sticky="NS")
        
        self.title = ttk.Label(self, text="Properties of Image and FFT:")
        self.title.grid(column=3, row=0, padx=0, pady=0, columnspan=2,
                        sticky="NW")
        
        # Image size in arcsec, pixel scale of FFT, Image limits
        self.imgSizeLab = ttk.Label(self, text=u"Image Extent [x, y] (arcmin):")
        self.imgSizeLab.grid(column=3, row=1, padx=(15,5), pady=5, sticky="E")
        self.imgSize_arcmin = tk.StringVar()
        self.imgSizeVal = ttk.Label(self, width=20,
                                    textvariable=self.imgSize_arcmin)
        self.imgSizeVal.grid(column=4, row=1, padx=5, pady=5, sticky="W")

        self.pixScaleLab = ttk.Label(self,
                                     text=u"FFT Pixel Scale [x, y] (k\u03bb):")
        self.pixScaleLab.grid(column=3, row=2, padx=(15,5), pady=5, sticky="E")
        self.pixScale_kl = tk.StringVar()
        self.pixScaleVal = ttk.Label(self, width=20,
                                     textvariable=self.pixScale_kl)
        self.pixScaleVal.grid(column=4, row=2, padx=5, pady=5, sticky="W")
        
        self.fftLimLab = ttk.Label(self,
                                     text=u"FFT Image Limits (k\u03bb):")
        self.fftLimLab.grid(column=3, row=3, padx=(15,5), pady=5, sticky="E")
        self.fftLim_kl = tk.StringVar()
        self.fftLimVal = ttk.Label(self, width=20,
                                     textvariable=self.fftLim_kl)
        self.fftLimVal.grid(column=4, row=3, padx=5, pady=5, sticky="W")
        
        

#-----------------------------------------------------------------------------#
class StatusFrame(ttk.Frame):
    """Frame presenting status indicators lights for the individual steps in
    the application and the two main control buttons."""
    
    def __init__(self, parent, boxWidth=20, gapWidth=30, yPad=25):
        ttk.Frame.__init__(self, parent)
        self.parent = parent
        bgColour = ttk.Style().lookup("TFrame", "background")
        #bgColour = "white"

        # Properties of the status boxes
        self.tagLst = ['statusModel', 'statusSelection', 'statusModelFFT',
                       'statusuvCalc', 'statusuvGrid', 'statusBeam',
                       'statusObsDone']
        labLst = ['Model', 'Array', 'Model FFT', 'uv-Coverage', 'uv-Grid',
                  'Beam', 'Observation']
        self.state = [0] * len(self.tagLst)
        
        # Calculate canvas dimensions and internal grid coordinates
        self.boxWidth = boxWidth
        self.gapWidth = gapWidth
        self.yPad = yPad
        self.nBox = len(self.tagLst)
        self.width = self.nBox * boxWidth + self.nBox * gapWidth
        self.height = boxWidth + yPad * 4.8
        yCentLab = yPad * 0.75
        self.yCentBox = boxWidth * 2.0
        yCentBrkt = self.yCentBox + boxWidth / 2.0 + yPad
        yCentBtn = yCentBrkt + yPad * 1.5
        self.xCentLst = []
        for i in range(self.nBox):
            x = gapWidth/2.0 + boxWidth/2.0+ i*(boxWidth + gapWidth)
            self.xCentLst.append(x)
            
        # Insert the canvas
        self.canvas = tk.Canvas(self, width=self.width, height=self.height,
                                background=bgColour, highlightthickness=0)
        self.canvas.grid(row=0, column=0, padx=0, pady=0)

        # Draw the boxes and labels in order
        for x, tag, text, in zip(self.xCentLst, self.tagLst, labLst):
            self._draw_checkbox(x, self.yCentBox, self.boxWidth, tag)
            self.canvas.create_text(x, yCentLab, text=text)

        # Draw the model line & plotting button
        x1 = self.xCentLst[0]
        y1 = self.yCentBox + boxWidth / 2.0
        y2 = yCentBrkt
        y3 = yCentBtn
        self.canvas.create_line(x1, y1, x1, y3, width=2, fill="black",
                                joinstyle=tk.MITER)
        self.pltModFFTbtn = ttk.Button(self, text = "Plot Model FFT",
                                       width=18,
                                       command=self._handler_plt_modFFT)
        self.pltModFFTbtn.configure(state="disabled")
        self.canvas.create_window(x1, y2, window=self.pltModFFTbtn)
        self.pltPwrSpecbtn = ttk.Button(self, text = "Plot Power Spectrum",
                                        width=18,
                                        command=self._handler_plt_pwrspec)
        self.pltPwrSpecbtn.configure(state="disabled")
        self.canvas.create_window(x1, y3, window=self.pltPwrSpecbtn)
        
        # Draw the array line & the plotting buttons
        x1 = self.xCentLst[1]
        y1 = self.yCentBox + boxWidth / 2.0
        y2 = yCentBrkt
        y3 = yCentBtn
        self.canvas.create_line(x1, y1, x1, y3, width=2, fill="black",
                                joinstyle=tk.MITER)
        self.pltuvCovBtn = ttk.Button(self, text = "Plot uv-Coverage",
                                      width=18,
                                      command=self._handler_plt_uvcov)
        self.pltuvCovBtn.configure(state="disabled")
        self.canvas.create_window(x1, y2, window=self.pltuvCovBtn)
        self.pltElBtn = ttk.Button(self, text = "Plot Elevation",
                                   width=18,
                                   command=self._handler_plt_elevation)
        self.pltElBtn.configure(state="disabled")
        self.elBtnW = self.canvas.create_window(x1, y3, window=self.pltElBtn)
        
        # Draw the line seperating inputs from outputs
        x1 = self.xCentLst[1] + (self.xCentLst[2] - self.xCentLst[1]) / 2.0
        x2 = self.xCentLst[0] + (self.xCentLst[1] - self.xCentLst[0]) / 2.0
        self.canvas.create_line(x1, 0, x1, self.height, width=2, dash=4)

        # Draw the observe button & bracket
        x1 = self.xCentLst[2] - gapWidth / 3.0
        x2 = self.xCentLst[-1] + gapWidth / 3.0
        x3 = self.xCentLst[4]
        y1 = self.yCentBox
        y2 = yCentBrkt
        y3 = yCentBtn
        self.canvas.create_line(x1, y1, x1, y2, x2, y2, x2, y1, width=2,
                                fill="black", joinstyle=tk.MITER)
        self.canvas.create_line(x3, y2, x3, y3, width=2, fill="black",
                                joinstyle=tk.MITER)
        self.obsBtn = ttk.Button(self, text = "Do Observation", width=30,
                                 command=self._handler_observe_button)
        self.obsBtn.configure(state="disabled")
        obsBtnW = self.canvas.create_window(x3, y3, window=self.obsBtn)
        
        
    def _draw_checkbox(self, xCent, yCent, size, tag, state=0, lw=3):
        """Draw a checkbox on the canvas."""

        r = float(size)/2.0
        
        # Draw box & ticks
        item = self.canvas.create_rectangle(xCent-r, yCent-r,
                                            xCent+r, yCent+r,
                                            outline="black", fill="white",
                                            width=1.0, tag="checkbox")
        self.canvas.addtag_withtag(tag, item)
        if state:
            item = self.canvas.create_line(xCent-r+lw, yCent+r/2-lw,
                                           xCent-r/2+lw, yCent+r-lw,
                                           xCent+r-lw, yCent-r+lw,
                                           fill="green", capstyle=tk.ROUND,
                                           width=lw, tag="tick")
            self.canvas.addtag_withtag(tag, item)
        else:
            item = self.canvas.create_line(xCent-r+lw, yCent-r+lw,
                                           xCent+r-lw, yCent+r-lw,
                                           fill="red", capstyle=tk.ROUND,
                                           width=lw, tag="tick")
            self.canvas.addtag_withtag(tag, item)
            item = self.canvas.create_line(xCent-r+lw, yCent+r-lw,
                                           xCent+r-lw, yCent-r+lw,
                                           fill="red", capstyle=tk.ROUND,
                                           width=lw, tag="tick")
            self.canvas.addtag_withtag(tag, item)

    def set_state_dict(self, stateDict):
        for k, v in stateDict.items():
            self.set_state(k, v)
            
    def set_state(self, tag, state=0):

        # Get the index of the box with the tag & set the new state
        indx = self.tagLst.index(tag)
        self.state[indx] = state

        # Find and delete all items with requested tag
        itemLst = self.canvas.find_withtag(tag)
        for item in itemLst:
            self.canvas.delete(item)
            
        # Re-draw the checkbox with the new state
        self._draw_checkbox(self.xCentLst[indx], self.yCentBox, self.boxWidth,
                            tag, state=self.state[indx])

        # Enable or dissable buttons based on state
        if self.state[1]:
            self.pltuvCovBtn.configure(state="enabled")
            self.pltElBtn.configure(state="enabled")
        else:
            self.pltuvCovBtn.configure(state="disabled")
            self.pltElBtn.configure(state="disabled")
        if self.state[0]:
            self.pltModFFTbtn.configure(state="enabled")
        else:
            self.pltModFFTbtn.configure(state="disabled")
        if self.state[0] and self.state[1]:
            self.obsBtn.configure(state="enabled")
        else:
            self.obsBtn.configure(state="disabled")
    
    # Event handlers ---------------------------------------------------------#

    def _handler_plt_modFFT(self):
        self.event_generate("<<plot_modFFT>>")
        
    def _handler_plt_pwrspec(self):
        self.event_generate("<<plot_powerspec>>")
        
    def _handler_plt_elevation(self):
        self.event_generate("<<plot_elevation>>")
        
    def _handler_plt_uvcov(self):
        self.event_generate("<<plot_uvcoverage>>")
        
    def _handler_observe_button(self):
        self.event_generate("<<do_observation>>")
        
        
#-----------------------------------------------------------------------------#
if __name__ == "__main__":
    root = tk.Tk()

    # Scaling tests
    #root.tk.call('tk', 'scaling', 4.0)
    #root.tk.call('tk', 'scaling', '-displayof', '.', 50)

    # Force colours 
    bgColour = ttk.Style().lookup("TFrame", "background")
    ttk.Style().configure("TFrame", background=bgColour)
    ttk.Style().configure("TLabelframe", background=bgColour)

    # Force font
    default_font = tkFont.nametofont("TkDefaultFont")
    default_font.configure(size=10)
    root.option_add("*Font", default_font)

    # Grid the main window and start mainloop
    app = App(root).pack(side="top", fill="both", expand=True)
    root.mainloop()
    
