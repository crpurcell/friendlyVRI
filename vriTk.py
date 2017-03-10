#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     vriTK.py                                                          #
#                                                                             #
# PURPOSE:  A virtual interferometer application written in Tkinter.          #
#                                                                             #
# REQUIRED: Requires numpy, tkinter, matplotlib                               #
#                                                                             #
# MODIFIED: 08-Mar-2017 by cpurcell                                           #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# App                ... class to create the root and array config windows    #
# ObsControlFrame    ... class defining the array chooser interface           #
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

# Window geometry
geometryMainWin = "1280x800"
geometryAConfWin = "1260x650"

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
import matplotlib as mpl
mpl.use("TkAgg")
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
        self.root.title("Friendly VRI: Display")
        #self.root.geometry(geometryMainWin)
        self.root.resizable(True, True)
        self.root.protocol("WM_DELETE_WINDOW", self.applicationExit)
        self.obsManager = None

        # Set the expansion
        self.root.columnconfigure(0, weight=1)
        self.root.columnconfigure(1, weight=1)
        self.root.columnconfigure(2, weight=1)
        self.root.rowconfigure(0, weight=1)
        self.root.rowconfigure(1, weight=1)
        self.root.rowconfigure(2, weight=1)
        
        # Create the model image panel
        self.modelImgFrm = SingleFigFrame(self.root)
        self.modelImgFrm.grid(column=0, row=0, columnspan=1, padx=5, pady=5,
                              sticky="NSEWS")

        # Create the model FFT panel
        self.modelFFTfrm = SingleFigFrame(self.root)
        self.modelFFTfrm.grid(column=0, row=1, padx=5, pady=5, sticky="NSEW")
        
        # Create the uv-coverage panel
        self.uvCovFrm = SingleFigFrame(self.root)
        self.uvCovFrm.grid(column=1, row=0, padx=5, pady=5, sticky="NSEW")
        
        # Create the observed FFT panel
        self.obsFFTfrm = SingleFigFrame(self.root)
        self.obsFFTfrm.grid(column=1, row=1, padx=5, pady=5, sticky="NSEW")

        # Create the synthesised beam panel
        self.beamFrm = SingleFigFrame(self.root)
        self.beamFrm.grid(column=2, row=0,  padx=5, pady=5, sticky="NSEW")

        # Create the observed image panel
        self.obsImgFrm = SingleFigFrame(self.root)
        self.obsImgFrm.grid(column=2, row=1,  padx=5, pady=5, sticky="NSEW")

        # Create the array chooser window and set the focus back to root
        self.aConfWin = tk.Toplevel(self.root)
        self.aConfWin.title("Friendly VRI: Observation Controller")
        #self.aConfWin.geometry(geometryAConfWin)
        self.aConfWin.resizable(True, True)
        self.aConfWin.protocol("WM_DELETE_WINDOW", self.applicationExit)
        
        # Create the array chooser frame
        self.controlFrm = ObsControlFrame(self.aConfWin)
        self.controlFrm.grid(row=0, column=0, padx=0, pady=0, sticky="NSEW")

        # Load the back-end and populate the array configuration list
        self.obsManager = observationManager(verbose=True, debug=True)
        vals = self.obsManager.arrsAvailable.values()
        configLst = zip([x["telescope"] for x in vals],
                        [x["config"] for x in vals])
        
        self.controlFrm.configInTab.name_columns(
            ("Telescope", "Array"))
        self.controlFrm.configInTab.insert_rows(configLst,
                                                   ("Telescope", "Array"))

        # Populate the observing parameters widgets
        self.controlFrm.freq_MHz.set(1420.0)
        self.controlFrm.sampRt_s.set(300)
        self.controlFrm.dec_deg.set(20.0)
        self.controlFrm.pixScale_asec.set(0.5)
        self.controlFrm.modelFile.set("models/radio_galaxy.png")
        
        # Bind virtual events generated by sub-widgets
        self.aConfWin.bind("<<config_in_selected>>",
                           lambda event : self._on_select_config(event))
        self.aConfWin.bind("<<selection_changed>>",
                           lambda event : self._on_sel_change(event))
        self.aConfWin.bind("<<plot_uvcoverage>>",
                           lambda event : self._on_plot_uvcov(event))
        self.aConfWin.bind("<<plot_elevation>>",
                           lambda event : self._on_plot_elevation(event))
        self.aConfWin.bind("<<do_observation>>",
                           lambda event : self._on_do_observation(event))
        self.aConfWin.bind("<<load_model_image>>",
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
        xArr, yArr = self.obsManager.get_ant_coordinates(row=row)
        
        # Plot the antenna positions 
        self.controlFrm.antPosPlot.load_data(xArr/1000.0, yArr/1000.0)
        self.controlFrm.antPosPlot.draw_zerolines()
        self.controlFrm.antPosPlot.set_xlabel("East-West (km)")
        self.controlFrm.antPosPlot.set_ylabel("North-South (km)")
       
    def _on_sel_change(self, event=None):
        self.obsManager.set_obs_parms(self.controlFrm.freq_MHz.get(),
                                      self.controlFrm.sampRt_s.get(),
                                      self.controlFrm.dec_deg.get())
        self.obsManager.clear_all_selections()

        # Re-do the selection based on the arrays and HAs in the GUI
        for selection in self.controlFrm.configOutTab.get_all_text():
            key = "_".join(selection[:2])
            haStart = float(selection[2])
            haEnd = float(selection[3])
            self.obsManager.select_array(key, haStart, haEnd)
        
        # Update the status
        stateDict = self.obsManager.get_status()
        self.controlFrm.nst.set_state_dict(stateDict)
        
    def _on_plot_uvcov(self, event=None):
        """Plot the uv-coverage for all selected array configurations"""

        # Clear any previous selection in the observationManager
        #self.obsManager.clear_all_selections()
        
#        # Re-do the selection based on the arrays and HAs in the GUI
#        for selection in self.controlFrm.configOutTab.get_all_text():
#            key = "_".join(selection[:2])
#            haStart = float(selection[2])
#            haEnd = float(selection[3])
#            self.obsManager.select_array(key, haStart, haEnd)
            
        # Calculate the uv-coverage for the selected observation parameters
        self.obsManager.calc_uvcoverage()

        # Plot in the display window
        ax = self.uvCovFrm.add_axis()
        plot_uvcov_ax(ax, self.obsManager.arrsSelected)
        self.uvCovFrm.show()
        
        # Update the status
        stateDict = self.obsManager.get_status()
        self.controlFrm.nst.set_state_dict(stateDict)
        
        # ALT: Plot the uv-coverage as an external MPL figure
        #fig = plt.figure(figsize=(10,10))
        #ax = fig.add_subplot(111)
        #plot_uvcov_ax(ax, self.obsManager.arrsSelected)    
        #fig.show()
        
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
        for selection in self.controlFrm.configOutTab.get_all_text():
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
        modelFile = self.controlFrm.modelFile.get()
        pixScale_asec = self.controlFrm.pixScale_asec.get()
        self.obsManager.load_model_image(modelFile, pixScale_asec)

        # Calculate the FFT of the model
        #self.obsManager.invert_model()        
        
        # Plot the model image
        ax = self.modelImgFrm.add_axis()
        plot_image_ax(ax, self.obsManager.modelImgArr)
        self.modelImgFrm.show()
        
        # Plot the model FFT
        ax = self.modelFFTfrm.add_axis()
        plot_fft_ax(ax, self.obsManager.modelFFTarr)
        self.modelFFTfrm.show()

        # Update the status
        stateDict = self.obsManager.get_status()
        self.controlFrm.nst.set_state_dict(stateDict)
        
    def _on_do_observation(self, event=None):

        # Calculate the uv-coverage if not cached.
        stateDict = self.obsManager.get_status()
        if not stateDict['statusuvCalc']:
            self.obsManager.calc_uvcoverage()            

            # Plot in the display window
            ax = self.uvCovFrm.add_axis()
            plot_uvcov_ax(ax, self.obsManager.arrsSelected)
            self.uvCovFrm.show()
            
            stateDict = self.obsManager.get_status()
            self.controlFrm.nst.set_state_dict(stateDict)
            
        # Grid the uv-coverage to make a mask
        self.obsManager.grid_uvcoverage()
        stateDict = self.obsManager.get_status()
        self.controlFrm.nst.set_state_dict(stateDict)
        
        # Calculate the PSF
        self.obsManager.calc_beam()
        stateDict = self.obsManager.get_status()
        self.controlFrm.nst.set_state_dict(stateDict)
        
        # Show the observed beam
        ax = self.beamFrm.add_axis()
        plot_image_ax(ax, self.obsManager.beamArr)
        self.beamFrm.show()
        
        # Apply the gridded uv-coverage and invert
        self.obsManager.invert_observation()
        
        # Show the observed uv-coverage
        ax = self.obsFFTfrm.add_axis()
        plot_fft_ax(ax, self.obsManager.obsFFTarr)
        self.obsFFTfrm.show()
        
        # Show the observed image
        ax = self.obsImgFrm.add_axis()
        plot_image_ax(ax, self.obsManager.obsImgArr)
        self.obsImgFrm.show()
        
        # Update the status
        stateDict = self.obsManager.get_status()
        self.controlFrm.nst.set_state_dict(stateDict)
        
#-----------------------------------------------------------------------------#
class ObsControlFrame(ttk.Frame):
    """Frame presenting an interface to choose one or more array contigurations
    and set the observing parameters."""
    
    def __init__(self, parent):
        ttk.Frame.__init__(self, parent)
        self.parent = parent

        # Row & column weightings        
        self.rowconfigure(1, weight=1)
        
        # Scatter plot of antenna locations
        self.antPosPlot = ScatterPlot(self, width=380, height=350,
                                      axPad=(100,25,70,25), aspect="equal")
        self.antPosPlot.grid(column=0, row=0, rowspan=3, padx=5, pady=5)
        
        # Listbox showing the available array configurations
        self.configInTab = ScrolledTreeTab(self,
                                           virtEvent="<<config_in_selected>>")
        self.configInTab.grid(column=1, row=0, rowspan=3, padx=5, pady=5,
                              sticky="NSEW")
        
        # Hour angle slider
        self.haFrm = ttk.Frame(self)
        self.haFrm.grid(column=2, row=0, padx=5, pady=5)
        self.haLab = ttk.Label(self.haFrm, text="Hour Angle Range (hours):")
        self.haLab.grid(column=0, row=0, padx=0, pady=0, sticky="NW")
        self.haScale = DoubleScale(self.haFrm, from_=-12.0, to=12.0,
                                   initLeft=-1.0, initRight=+1.0,
                                   tickIntMajor=6, tickIntMinor=1, width=270)
        self.haScale.grid(column=0, row=1, padx=0, pady=5, sticky="NEW")
        
        # Fancy add button with strike-through arrow
        bgColour = ttk.Style().lookup("TFrame", "background")
        self.addButFrm = ttk.Frame(self)
        self.addButFrm.grid(column=2, row=1, padx=5, pady=5, sticky="NSEW")
        self.arrow = tk.Canvas(self.addButFrm, width=270, height=30,
                               background=bgColour, highlightthickness=0)
        self.arrow.create_line(10,15,260,15, width=2, arrow=tk.LAST,
                                arrowshape=(10,15,5))
        self.arrow.grid(column=0, row=0, columnspan=3, padx=0, pady=0,
                        sticky="EW")
        self.addBtn = ttk.Button(self.addButFrm, text="ADD", width=10,
                                 command=self._handler_add_button)
        self.addBtn.grid(column=1, row=0, padx=0, pady=0)

        # Action button frame: delete selection, plot uv-coverage & elevation
        self.actFrm = ttk.Frame(self)
        self.actFrm.grid(column=2, row=2, padx=5, pady=5, sticky="SEW")
        self.actFrm.columnconfigure(0, weight=1)
        self.delAllBtn = ttk.Button(self.actFrm, text="Clear All",
                                    command=self._handler_clear_all_button)
        self.delAllBtn.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW" )
        self.delBtn = ttk.Button(self.actFrm, text="Clear Selected",
                                 command=self._handler_clear_button)
        self.delBtn.grid(column=0, row=1, padx=5, pady=5, sticky="EW" )
        
        # Listbox showing the selected array configurations
        self.configOutTab = ScrolledTreeTab(self,
                                        virtEvent="<<config_out_selected>>")
        self.configOutTab.name_columns(("Telescope", "  Array  ",
                                        "HA-Start", " HA-End "))
        self.configOutTab.grid(column=3, row=0, rowspan=3, padx=5, pady=5,
                               sticky="NSEW")
        
        # Common observing parameter frame
        self.obsParmFrm = ttk.LabelFrame(self,
                                    text="  Observation Settings  ")
        self.obsParmFrm.grid(column=0, row=3, padx=5, pady=5, sticky="NSEW")
        self.freqLab = ttk.Label(self.obsParmFrm,
                                 text="Observing Frequency (MHz):")
        self.freqLab.grid(column=0, row=0, padx=5, pady=5, sticky="E")
        self.freq_MHz = tk.StringVar()
        self.freqEnt = ttk.Entry(self.obsParmFrm, width=10,
                                 textvariable=self.freq_MHz)
        self.freqEnt.grid(column=1, row=0, padx=5, pady=5, sticky="EW")
        self.sampRtLab = ttk.Label(self.obsParmFrm,
                                   text="Data Sampling Rate (s):")
        self.sampRtLab.grid(column=0, row=1, padx=5, pady=5, sticky="E")
        sampRtLst_s = ["10", "60", "60", "100", "300", "600", "1200", "1800",
                       "3600"]
        self.sampRt_s = tk.StringVar()
        self.sampRtComb = ttk.Combobox(self.obsParmFrm,
                                       textvariable=self.sampRt_s,
                                       values=sampRtLst_s, width=15)
        self.sampRtComb.current(4)
        self.sampRtComb.grid(column=1, row=1, padx=5, pady=5, sticky="EW")
        self.decSrcLab = ttk.Label(self.obsParmFrm,
                                   text="Source Declination (degrees):")
        self.decSrcLab.grid(column=0, row=2, padx=5, pady=5, sticky="E")
        self.dec_deg = tk.StringVar()
        self.decSrcEnt = ttk.Entry(self.obsParmFrm, width=15,
                                   textvariable=self.dec_deg)
        self.decSrcEnt.grid(column=1, row=2, padx=5, pady=5, sticky="EW")
        
        self.robustLab = ttk.Label(self.obsParmFrm,
                                  text="Robust Weighting Factor:")
        self.robustLab.grid(column=0, row=3, padx=5, pady=5, sticky="E")
        robustLst = ["None", "-2.0", "-1.5", "-1.0", "-0.5", "0.0",
                     "0.5", "1.0", "1.5", "2.0"]
        self.robust = tk.StringVar()
        self.robustComb = ttk.Combobox(self.obsParmFrm, width=15,
                                       textvariable=self.robust,
                                       values=robustLst)
        self.robustComb.grid(column=1, row=3, padx=5, pady=5, sticky="EW")
        self.robustComb.config(state="disabled")      # DISABLED UNTIL BUILT
        self.obsParmFrm.columnconfigure(0, weight=1)
        self.obsParmFrm.columnconfigure(1, weight=1)

        # Model chooser frame
        self.modelLoadFrm = ttk.LabelFrame(self, text="  Model Image  ")
        self.modelLoadFrm.grid(column=1, row=3, padx=5, pady=5,
                               columnspan=2, sticky="NSEW")
        self.modelFile = tk.StringVar()
        self.fileLab = ttk.Label(self.modelLoadFrm, text="File:")
        self.fileLab.grid(column=0, row=0, padx=5, pady=5, sticky="E")
        self.browseBtn = ttk.Button(self.modelLoadFrm, text="Browse",
                                    width=10,
                                    command=self._handler_browse_button)
        self.browseBtn.grid(column=2, row=0, padx=5, pady=2, sticky="E")
        self.loadBtn = ttk.Button(self.modelLoadFrm, text="Load",
                                  width=10,
                                  command=self._handler_load_button)
        self.loadBtn.grid(column=3, row=0, padx=5, pady=5, sticky="E")
        self.modelFile = tk.StringVar()
        self.fileEnt = ttk.Entry(self.modelLoadFrm, width=30,
                                 textvariable=self.modelFile)
        self.fileEnt.grid(column=0, row=1, columnspan=4, padx=5, pady=5,
                          sticky="EW")
        self.modelLoadFrm.columnconfigure(1, weight=1)
        text = "Supported image formats: PNG, JPEG, TIFF."
        self.formatInfoLab = ttk.Label(self.modelLoadFrm,text=text)
        self.formatInfoLab.grid(column=0, row=2, padx=5, pady=5,
                                columnspan=4, sticky="E")
        self.pixScaLab = ttk.Label(self.modelLoadFrm,
                                   text="Pixel Scale (arcseconds):")
        self.pixScaLab.grid(column=1, row=3, columnspan=2, padx=5, pady=5,
                            sticky="E")
        self.pixScale_asec = tk.DoubleVar()
        self.pixScaEnt = ttk.Entry(self.modelLoadFrm, width=10,
                                   textvariable=self.pixScale_asec)
        self.pixScaEnt.grid(column=3, row=3, padx=5, pady=5, sticky="EW")
        


        # Information panel
        self.infoFrm = ttk.LabelFrame(self, text="  Information  ")
        self.infoFrm.grid(column=0, row=4, columnspan=4, padx=5, pady=5,
                          sticky="NSEW")
        self.uvcovLab = ttk.Label(self.infoFrm,
                                  text="Parameters from uv-coverage:")
        self.uvcovLab.grid(row=0, column=0, columnspan=2, padx=5, pady=5,
                           sticky="NSEW")
        
        sep = ttk.Separator(self.infoFrm, orient="vertical")
        sep.grid(row=0, column=2, rowspan=4, padx=5, pady=5, sticky="NS")
        self.imgLab = ttk.Label(self.infoFrm,
                                text="Parameters from model image:")
        self.imgLab.grid(row=0, column=3, columnspan=2, padx=5, pady=5,
                         sticky="NSEW")
        
        # New Status frame
        self.nst = StatusFrame(self, boxWidth=20, gapWidth=150)
        self.nst.grid(column=0, row=5, columnspan=4, padx=5, pady=5,
                      sticky="NSEW")

        
    def _handler_observe_button(self):
        """Run the observation"""
        
        self.event_generate("<<do_observation>>")
        
    def _handler_browse_button(self):
        """Open the file selection dialog and set the session dir."""
        
        modelFile = tkFileDialog.askopenfilename(parent=self,
                                    title="Choose a model file")
        if not len(modelFile)==0:
            self.modelFile.set(modelFile)
        
    def _handler_load_button(self):      
        """Raise a <<load_session>> virtual event in the parent window."""
        
        self.event_generate("<<load_model_image>>")

    def _handler_add_button(self):
        """Add the selected configuration to the list box"""
        selConf = self.configInTab.get_text_selected()
        if selConf is not None:
            selConf = list(selConf)
            selConf.append(self.haScale.valueLeft.get())
            selConf.append(self.haScale.valueRight.get())
            self.configOutTab.insert_rows([selConf],
                            ("Telescope", "  Array  ", "HA-Start", " HA-End "))
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
            
        # Draw the line & the plotting buttons
        x1 = self.xCentLst[1] - boxWidth / 2.0
        x2 = x1 - gapWidth * 0.5
        y1 = self.yCentBox
        y2 = yCentBrkt
        y3 = yCentBtn
        self.canvas.create_line(x1, y1, x2, y1, x2, y3, width=2, fill="black",
                                joinstyle=tk.MITER)
        
        self.pltuvCovBtn = ttk.Button(self, text = "Plot uv-Coverage",
                                      width=20,
                                      command=self._handler_plt_uvcov)
        self.pltuvCovBtn.configure(state="disabled")        
        self.canvas.create_window(x2, y2, window=self.pltuvCovBtn)
        
        self.pltElBtn = ttk.Button(self, text = "Plot Elevation", width=20,
                                   command=self._handler_plt_elevation)
        self.pltElBtn.configure(state="disabled")
        self.elBtnW = self.canvas.create_window(x2, y3, window=self.pltElBtn)
        
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
            
        if self.state[0] and self.state[1]:
            self.obsBtn.configure(state="enabled")
        else:
            self.obsBtn.configure(state="disabled")
            
        
                
            
    # Event handlers ---------------------------------------------------------#
    
    def _handler_plt_elevation(self):
        self.event_generate("<<plot_elevation>>")
        
    def _handler_plt_uvcov(self):
        self.event_generate("<<plot_uvcoverage>>")
        
    def _handler_observe_button(self):
        self.event_generate("<<do_observation>>")
        
        
#-----------------------------------------------------------------------------#
if __name__ == "__main__":
    root = tk.Tk()
    #root.tk.call('tk', 'scaling', 4.0)
    #root.tk.call('tk', 'scaling', '-displayof', '.', 50)
    bgColour = ttk.Style().lookup("TFrame", "background")
    ttk.Style().configure("TFrame", background=bgColour)
    ttk.Style().configure("TLabelframe", background=bgColour)
    
    default_font = tkFont.nametofont("TkDefaultFont")
    default_font.configure(size=10)
    root.option_add("*Font", default_font)
    app = App(root)
    
