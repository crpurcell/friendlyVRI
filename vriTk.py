#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     vriTK.py                                                          #
#                                                                             #
# PURPOSE:  A virtual interferometer application written in Tkinter.          #
#                                                                             #
# REQUIRED: Requires numpy, tkinter, matplotlib                               #
#                                                                             #
# MODIFIED: 29-Mar-2017 by cpurcell                                           #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# App (class)        ... class containing the main application logic          #
#     _applicationExit                                                        #
#     _show_about                                                             #
#     _show_help                                                              #
#     _update_status                                                          #
#     _on_select_config                                                       #
#     _on_sel_change                                                          #
#     _on_obsparm_change                                                      #
#     _on_pixscale_change                                                     #
#     _on_plot_modFFT                                                         #
#     _on_plot_uvcov                                                          #
#     _on_plot_elevation                                                      #
#     _on_load_model                                                          #
#     _on_do_observation                                                      #
#                                                                             #
# ArraySelector      ... class defining the array selection interface         #
#     _handler_add_button                                                     #
#     _handler_clear_button                                                   #
#     _handler_clear_all_button                                               #
#                                                                             #
# ObsInputs          ... class exposing the remaining observation inputs      #
#     _handler_browse_button                                                 #
#                                                                             #
# StatusFrame        ... class defining status indicators and action buttons  #
#     _draw_checkbox                                                          #
#     set_state_by_dict                                                       #
#     set_state                                                               #
#     _handler_plt_modFFT                                                     #
#     _handler_plt_pwrspec                                                    #
#     _handler_plt_elevation                                                  #
#     _handler_plt_uvcov                                                      #
#     _handler_observe_button                                                 #
#                                                                             #
# InformationPanel   ... class showing derived properties of observation      #
#     update                                                                  #
#                                                                             #
# PlotFrame          ... class defining the plotting window                   #
#     _plot_image                                                             #
#     _plot_fft                                                               #
#     _plot_uvcov                                                             #
#     plot_model_image                                                        #
#     clear_model_image                                                       #
#     plot_model_fft                                                          #
#     clear_model_fft                                                         #
#     plot_uvcov                                                              #
#     clear_uvcov                                                             #
#     plot_beam                                                               #
#     clear_beam                                                              #
#     plot_obs_fft                                                            #
#     clear_obs_fft                                                           #
#     plot_obs_image                                                          #
#     clear_obs_image                                                         #
#     show                                                                    #
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
    from ScrolledText import ScrolledText as tkScrolledText
except Exception:  # Python 3.x
    import tkinter as tk
    from tkinter import ttk
    import tkinter.font as tkFont
    import tkinter.messagebox as tkMessageBox
    import tkinter.filedialog as tkFileDialog
    import tkinter.simpledialog as tkSimpleDialog
    import tkinter.scrolledtext as tkScrolledText
import numpy as np
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from Imports.util_tk import *
from Imports.vriCalc import *


#-----------------------------------------------------------------------------#
class App(ttk.Frame):
    """Class defining the Virtual Radio Interferometer application.

    This class creates a root window used to set observation parameters, choose
    and accumulate array configurations and load a model image. A secondary 
    Toplevel window is used to display the input, output and intermediate
    Fourier transforms and images."""
    
    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.parent.title("Friendly VRI: Control Window")
        self.obsManager = None

        # Set the grid expansion properties
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)
        self.rowconfigure(2, weight=1)

        # Menu bar
        self.menuBar = tk.Menu(self)
        self.menuBar.add_command(label="About", command=self._show_about)
        self.menuBar.add_command(label="Help", command=self._show_help)
        self.menuBar.add_command(label="Quit", command=self._applicationExit)
        self.parent.config(menu=self.menuBar)
        
        # Array selector interface
        self.selector =  ArraySelector(self)
        self.selector.grid(column=0, row=0, padx=10, pady=5, sticky="EW")
        
        sep = ttk.Separator(self, orient="horizontal")
        sep.grid(column=0, row=1, padx=10, pady=5, sticky="EW")
        
        # Observation settings
        self.inputs = ObsInputs(self)
        self.inputs.grid(column=0, row=2, padx=10, pady=5, sticky="EW")
        
        sep = ttk.Separator(self, orient="horizontal")
        sep.grid(column=0, row=3, padx=10, pady=5, sticky="EW")
        
        # Status frame
        self.statFrm = StatusFrame(self, boxWidth=20, gapWidth=150)
        self.statFrm.grid(column=0, row=4, padx=10, pady=5)
        
        # Create the display window
        self.dispWin = tk.Toplevel(self)
        self.dispWin.title("Friendly VRI: Display Window")
        self.dispWin.protocol("WM_DELETE_WINDOW", self._applicationExit)
        self.dispWin.columnconfigure(0, weight=1)
        self.dispWin.rowconfigure(0, weight=1)

        # Draw the display interface
        self.pltFrm = PlotFrame(self.dispWin)
        self.pltFrm.grid(row=0, column=0, padx=0, pady=0, sticky="NSEW")
        
        # DEBUG
        if False:
            self.testWin = tk.Toplevel(self)
            self.testWin.title(" TEST WINDOW ")
            self.testWin.protocol("WM_DELETE_WINDOW", self._applicationExit)
            self.testWin.resizable(True, True)
            self.testWin.columnconfigure(0, weight=1)
            self.testWin.rowconfigure(0, weight=1)
            self.a =PlotFrame(self.testWin)
            self.a.grid(row=0, column=0, padx=0, pady=0, sticky="NSEW")
            
        # Load the back-end and populate the array configuration list
        # The observationManager class can be used stand-alone from an
        # iPython terminal and replicates the functionality of the GUI.
        # The tkinter GUI serves as a controller for this class.
        self.obsManager = observationManager(verbose=False, debug=False)
        vals = self.obsManager.arrsAvailable.values()
        configLst = zip([x["telescope"] for x in vals],
                        [x["config"] for x in vals])        
        self.selector.configInTab.insert_rows(configLst,
                                                      ("Telescope", "Array"))
        
        # Bind virtual events generated by the control widgets
        self.parent.bind("<<config_in_selected>>",
                  lambda event : self._on_select_config(event))
        self.parent.bind("<<selection_changed>>",
                  lambda event : self._on_sel_change(event))
        self.parent.bind("<<obsparm_changed>>",
                  lambda event : self._on_obsparm_change(event))    
        self.parent.bind("<<pixscale_changed>>",
                  lambda event : self._on_pixscale_change(event))        
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
        
        # Force a minimum size on the windows & set resize properties
        self.parent.update()
        self.parent.minsize(self.parent.winfo_width(),
                            self.parent.winfo_height())
        self.parent.resizable(False, False)
        self.dispWin.update()
        self.dispWin.minsize(self.dispWin.winfo_width(),
                                self.dispWin.winfo_height())
        self.dispWin.resizable(True, True)
        
    def _applicationExit(self):
        """Exit the application cleanly if the window is closed."""
        
        self.parent.destroy()

    def _show_about(self):
        """Show the about text in a new window."""
        
        self.aboutWin = tk.Toplevel(self)
        self.aboutWin.title("About Friendly VRI")
        self.aboutTxt = tkScrolledText(self.aboutWin)
        self.aboutTxt.config(state="normal")
        with open('README.txt','r') as f:
            text = f.read()
        self.aboutTxt.insert('1.0', text)
        self.aboutTxt.config(state="disabled")
        self.aboutTxt.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW")
        self.closeBtn = ttk.Button(self.aboutWin, text='Close',
                           command=self.aboutWin.destroy)
        self.closeBtn.grid(column=0, row=1, padx=5, pady=5, sticky="E")
        self.aboutWin.rowconfigure(0, weight=1)
        self.aboutWin.columnconfigure(0, weight=1)
        
    def _show_help(self):
        """Show the help text in a new window."""
        
        self.helpWin = tk.Toplevel(self)
        self.helpWin.title("Help File for Friendly VRI")
        self.helpTxt = tkScrolledText(self.helpWin)
        self.helpTxt.config(state="normal")
        with open('HELP.txt','r') as f:
            text = f.read()
        self.helpTxt.insert('1.0', text)
        self.helpTxt.config(state="disabled")
        self.helpTxt.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW")
        self.closeBtn = ttk.Button(self.helpWin, text='Close',
                           command=self.helpWin.destroy)
        self.closeBtn.grid(column=0, row=1, padx=5, pady=5, sticky="E")
        self.helpWin.rowconfigure(0, weight=1)
        self.helpWin.columnconfigure(0, weight=1)
        
    def _update_status(self):
        """Update the status of the user interface, including the checkbox
        indicators, plot axes and information panels. The status of each step
        is queried from the instance of the observationManager class."""

        # Query the status and set the indicators
        stateDict = self.obsManager.get_status()
        self.statFrm.set_state_by_dict(stateDict)
        
        # Clear the plots based on the status
        self.pltFrm.clear_by_dict(stateDict)
        
        # Query the scales and update the information panel
        parmDict = self.obsManager.get_scales()
        self.pltFrm.infoPanel.update(parmDict)
        
    # Event handlers bound to virtual events ---------------------------------#
    
    def _on_select_config(self, event=None):
        """Plot the antenna layout for the selected array configuration"""
        
        # Query the antenna parameters of the selected configuration
        row = event.widget.get_indx_selected()
        d = self.obsManager.get_array_params(row=row)
        
        # Plot the antenna positions 
        self.selector.antPosPlot.load_data(d["x"]/1000.0, d["y"]/1000.0)
        self.selector.antPosPlot.draw_zerolines()
        self.selector.antPosPlot.set_xlabel("East-West (km)")
        self.selector.antPosPlot.set_ylabel("North-South (km)")
        
        # Show the antenna diameter and array latitude
        text = "{:.1f} m".format(d["diameter_m"])
        self.selector.antD_m.set(text)
        text = u"{:.6f}\u00B0".format(d["latitude_deg"])
        self.selector.antL_deg.set(text)
        text = u"{:.3f} km".format(d["baseMin_m"]/1000.0)
        self.selector.minBase_km.set(text)
        text = u"{:.3f} km".format(d["baseMax_m"]/1000.0)
        self.selector.maxBase_km.set(text)
        
    def _on_sel_change(self, event=None):
        """When the arrays selected change in the GUI, refresh the observation
        parameters in the observation manager."""

        # Reset the common parameters
        self.obsManager.set_obs_parms(self.inputs.freq_MHz.get(),
                                      self.inputs.sampRt_s.get(),
                                      self.inputs.dec_deg.get())

        # Reset the array and hour-angle selections
        self.obsManager.clear_all_selections()
        for selection in self.selector.configOutTab.get_all_text():
            key = "_".join(selection[:2])
            haStart = float(selection[2])
            haEnd = float(selection[3])
            self.obsManager.select_array(key, haStart, haEnd)
        
        # Update the status
        self._update_status()

    def _on_obsparm_change(self, event=None):
        """When the observation parameters change in the GUI, reset the 
        calculations."""

        # Reset the common parameters
        self.obsManager.set_obs_parms(self.inputs.freq_MHz.get(),
                                      self.inputs.sampRt_s.get(),
                                      self.inputs.dec_deg.get())
        
        # Update the status
        self._update_status()
        
    def _on_pixscale_change(self, event=None):
        """When the pixel scale is changed in the GUI, clear the FFT plot."""

        # Only operate if a model is already loaded
        stateDict = self.obsManager.get_status()
        if not stateDict["statusModel"]:
            return

        # Reset the pixel scale
        # TODO: check for float
        pixScale_asec = self.inputs.pixScale_asec.get()
        self.obsManager.set_pixscale(pixScale_asec)
        
        # Set the extent label in the load frame
        pixScaleImg_deg = self.obsManager.pixScaleImg_asec / 3600.0
        nX = self.obsManager.nX
        nY = self.obsManager.nY
        text = " %d x %d pix  /  %s x %s" % (nX, nY,
                                             ang2str(nX * pixScaleImg_deg),
                                             ang2str(nY * pixScaleImg_deg))
        self.inputs.extent.set(text)
        
        # Update the status
        self._update_status()
        
    def _on_plot_modFFT(self, event=None):
        """Show the FFT of the model image."""

        # Invert the model image
        self.obsManager.invert_model()
        
        # Plot the model FFT
        parmDict = self.obsManager.get_scales()
        lim_kl = parmDict["fftScale_lam"]/1e3
        self.pltFrm.plot_model_fft(self.obsManager.modelFFTarr, limit=lim_kl)
        
        # Update the status
        self._update_status()
        
    def _on_plot_uvcov(self, event=None):
        """Plot the uv-coverage for all selected array configurations"""
        
        # Calculate the uv-coverage for the selected observation parameters
        stateDict = self.obsManager.get_status()
        if not stateDict['statusuvCalc']:        
            self.obsManager.calc_uvcoverage()

            # Plot the uv-coverage in the display window
            self.pltFrm.plot_uvcov(self.obsManager.arrsSelected)
        
            # Update the status
            self._update_status()
        
    def _on_plot_elevation(self, event=None):
        """Create a plot showing the elevation of the source as seen from 
        the currently selected telescopes."""
        pass 
        
        # DEBUG
        if False:
            fig = plt.figure(figsize=(10,10))
            ax = fig.add_subplot(111)
            for selection in self.selector.configOutTab.get_all_text():
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
        modelFile = self.inputs.modelPath
        pixScale_asec = self.inputs.pixScale_asec.get()
        self.obsManager.load_model_image(modelFile, pixScale_asec)

        # Set the extent label in the load frame
        pixScaleImg_deg = self.obsManager.pixScaleImg_asec / 3600.0
        nX = self.obsManager.nX
        nY = self.obsManager.nY
        text = " %d x %d pix  /  %s x %s" % (nX, nY,
                                             ang2str(nX * pixScaleImg_deg),
                                             ang2str(nY * pixScaleImg_deg))
        self.inputs.extent.set(text)
        
        # Plot the model image
        self.pltFrm.plot_model_image(self.obsManager.modelImgArr)
        
        # Update the status
        self._update_status()
        
    def _on_do_observation(self, event=None):
        """Perform the bulk of the observing steps"""
        
        # Calculate the Fourier transform of the model if not cached
        stateDict = self.obsManager.get_status()
        if not stateDict['statusModelFFT']:
            self.obsManager.invert_model()
        
            # Plot the model FFT
            parmDict = self.obsManager.get_scales()
            lim_kl = parmDict["fftScale_lam"]/1e3
            self.pltFrm.plot_model_fft(self.obsManager.modelFFTarr,
                                       limit=lim_kl)
            
            # Update the status
            self._update_status()
        
        # Calculate the uv-coverage if not cached
        stateDict = self.obsManager.get_status()
        if not stateDict['statusuvCalc']:
            self.obsManager.calc_uvcoverage()

            # Plot the uv-coverage
            self.pltFrm.plot_uvcov(self.obsManager.arrsSelected)
            
            # Update the status
            self._update_status()
            
        # Grid the uv-coverage to make a mask
        self.obsManager.grid_uvcoverage()
        
        # Apply the gridded uv-coverage and invert
        self.obsManager.invert_observation()
        
        # Update the status
        self._update_status()
        
        # Show the observed FFT
        parmDict = self.obsManager.get_scales()
        lim_kl = parmDict["fftScale_lam"]/1e3
        self.pltFrm.plot_obs_fft(self.obsManager.obsFFTarr, limit=lim_kl)
        
        # Calculate the PSF
        self.obsManager.calc_beam()
        
        # Show the synthesised beam
        self.pltFrm.plot_beam(np.abs(self.obsManager.beamArr))
        
        # Update the status
        self._update_status()
        
        # Show the observed image
        self.pltFrm.plot_obs_image(np.abs(self.obsManager.obsImgArr))
        
        # Update the status
        self._update_status()

        
#-----------------------------------------------------------------------------#
class ArraySelector(ttk.Frame):
    """Two multi-column listboxes an hour-angle slider and a 'Select' button
    that make up the array selection interface."""
    
    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        # Set the expansion properties
        self.columnconfigure(2, weight=10)
        self.columnconfigure(4, weight=1)
        self.columnconfigure(5, weight=20)
        self.columnconfigure(6, weight=20)
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
        self.configInTab.grid(column=1, row=0, columnspan=2, rowspan=4,
                              padx=5, pady=5, sticky="NSEW")
        
        # Antenna Diameter and array latitude labels
        self.antDlab = ttk.Label(self, text=u"Antenna \u0398:   ")
        self.antDlab.grid(column=1, row=4, padx=0, pady=0, sticky="E")
        self.antD_m = tk.StringVar()
        self.antDval = ttk.Label(self, textvariable=self.antD_m)
        self.antDval.grid(column=2, row=4, padx=0, pady=0, sticky="W")
        self.antLlab = ttk.Label(self, text=u"Telescope \u03c6:   ")
        self.antLlab.grid(column=1, row=5, padx=0, pady=0, sticky="E")
        self.antL_deg = tk.StringVar()
        self.antLval = ttk.Label(self, textvariable=self.antL_deg)
        self.antLval.grid(column=2, row=5, padx=0, pady=0, sticky="W")

        # Minimum and maximum baseline labels
        self.minBaseLab = ttk.Label(self, text=u"Min Baseline:   ")
        self.minBaseLab.grid(column=3, row=5, padx=0, pady=0, sticky="E")
        self.minBase_km = tk.StringVar()
        self.minBaseVal = ttk.Label(self, textvariable=self.minBase_km)
        self.minBaseVal.grid(column=4, row=5, padx=0, pady=0, sticky="W")
        self.maxBaseLab = ttk.Label(self, text=u"Max Baseline:   ")
        self.maxBaseLab.grid(column=3, row=4, padx=0, pady=0, sticky="E")
        self.maxBase_km = tk.StringVar()
        self.maxBaseVal = ttk.Label(self, textvariable=self.maxBase_km )
        self.maxBaseVal.grid(column=4, row=4, padx=0, pady=0, sticky="W")
        
        # Hour angle slider
        self.haLab = ttk.Label(self, text="Hour Angle Range (hours):")
        self.haLab.grid(column=3, row=0, columnspan=2, padx=0, pady=5,
                        sticky="NW")
        self.haScale = DoubleScale(self, from_=-12.0, to=12.0,
                                   initLeft=-1.0, initRight=+1.0,
                                   tickIntMajor=6, tickIntMinor=1, width=270)
        self.haScale.grid(column=3, row=1, columnspan=2, padx=0, pady=5,
                          sticky="N")

        # Fancy add button with strike-through arrow
        bgColour = ttk.Style().lookup("TFrame", "background")
        self.canvas = tk.Canvas(self, width=270, height=30,
                               background=bgColour, highlightthickness=0)
        self.canvas.create_line(10,15,260,15, width=2, arrow=tk.LAST,
                                arrowshape=(10,15,5), fill="black")
        self.addBtn = ttk.Button(self, text="Add", width=10,
                                 command=self._handler_add_button)
        self.canvas.create_window(135, 15, window=self.addBtn)
        self.canvas.grid(column=3, row=2, columnspan=2, padx=0, pady=5,
                         sticky="NS")

        # Listbox showing the selected array configurations
        self.configOutTab = ScrolledTreeTab(self,
                                        virtEvent="<<config_out_selected>>")
        self.configOutTab.name_columns(("Telescope", "  Array  ",
                                        "HA-Start", " HA-End "))
        self.configOutTab.grid(column=5, row=0, columnspan=2, rowspan=4,
                               padx=5, pady=5, sticky="NSEW")
        
        # Delete buttons
        self.delBtn = ttk.Button(self, text="Clear Selected", width=20,
                                 command=self._handler_clear_button)
        self.delBtn.grid(column=5, row=5, padx=5, pady=5, sticky="EW" )
        self.delAllBtn = ttk.Button(self, text="Clear All", width=20,
                                    command=self._handler_clear_all_button)
        self.delAllBtn.grid(column=6, row=5, padx=5, pady=5, sticky="EW" )
                
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

        if self.configOutTab.clear_selected():
            self.event_generate("<<selection_changed>>")
        
    def _handler_clear_all_button(self):
        """Delete all configurations from the list box"""
        
        self.configOutTab.clear_entries()   
        self.event_generate("<<selection_changed>>")


#-----------------------------------------------------------------------------#
class ObsInputs(ttk.Frame):
    """Input settings for the observation (observing frequency, data sampling
    rate, Declination of source, Robust weighting factor, model image."""

    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        # Set the expansion properties
        self.columnconfigure(3, weight=1)
        self.columnconfigure(6, weight=1)
        
        # Observing frequency
        self.freqLab = ttk.Label(self, text="Observing Frequency (MHz):")
        self.freqLab.grid(column=1, row=0, padx=5, pady=5, sticky="E")
        self.freq_MHz = tk.StringVar()
        self.freq_MHz.trace("w", lambda dummy1, dummy2, dummy3:
                            self.event_generate("<<obsparm_changed>>"))
        self.freq_MHz.set(1420.0)
        self.freqEnt = ttk.Entry(self, width=10, textvariable=self.freq_MHz)
        self.freqEnt.grid(column=2, row=0, padx=5, pady=5, sticky="EW")

        # Source Declination
        self.decSrcLab = ttk.Label(self,
                                   text="Source Declination (deg):")
        self.decSrcLab.grid(column=1, row=1, padx=5, pady=5, sticky="E")
        self.dec_deg = tk.StringVar()
        self.dec_deg.trace("w", lambda dummy1, dummy2, dummy3:
                           self.event_generate("<<obsparm_changed>>"))
        self.dec_deg.set(20.0)
        self.decSrcEnt = ttk.Entry(self, width=10, textvariable=self.dec_deg)
        self.decSrcEnt.grid(column=2, row=1, padx=5, pady=5, sticky="EW")
        
        # Sampling rate
        self.sampRtLab = ttk.Label(self, text="Sampling Rate (s):")
        self.sampRtLab.grid(column=4, row=0, padx=5, pady=5, sticky="E")
        sampRtLst_s = ["10", "30", "60", "100", "300", "600", "1200", "1800",
                       "3600"]
        self.sampRt_s = tk.StringVar()
        self.sampRt_s.trace("w", lambda dummy1, dummy2, dummy3:
                            self.event_generate("<<obsparm_changed>>"))
        self.sampRtComb = ttk.Combobox(self, state="readonly",
                                       textvariable=self.sampRt_s,
                                       values=sampRtLst_s, width=7)
        self.sampRtComb.current(4)
        self.sampRtComb.grid(column=5, row=0, padx=5, pady=5, sticky="EW")
        
        # Weighting factor
        self.robustLab = ttk.Label(self, state="readonly" ,
                                   text="Robust Weighting Factor:")
        self.robustLab.grid(column=4, row=1, padx=(15,5), pady=5, sticky="E")
        #robustLst = ["None", "-2.0", "-1.5", "-1.0", "-0.5", "0.0",
        #             "0.5", "1.0", "1.5", "2.0"]
        robustLst = ["None"]
        self.robust = tk.StringVar()
        self.robust.trace("w", lambda dummy1, dummy2, dummy3:
                          self.event_generate("<<obsparm_changed>>"))
        self.robustComb = ttk.Combobox(self, state="readonly",
                                       textvariable=self.robust,
                                       values=robustLst, width=15)
        self.robustComb.current(0)
        self.robustComb.grid(column=5, row=1, padx=5, pady=5, sticky="EW")
    
        # Model image
        self.fileLab = ttk.Label(self, text="Model Image:")
        self.fileLab.grid(column=7, row=0, padx=(15,5), pady=5, sticky="E")
        self.modelFile = tk.StringVar()
        self.fileEnt = ttk.Entry(self, width=30,
                                 textvariable=self.modelFile)
        self.fileEnt.configure(state="readonly")        
        self.fileEnt.grid(column=8, row=0, columnspan=4, padx=5, pady=5,
                          sticky="EW")

        # Browse and load buttons
        self.browsePhoto = tk.PhotoImage(file='Imports/folder.gif')
        self.browseBtn = ttk.Button(self, image=self.browsePhoto,
                                    command=self._handler_browse_button)
        self.browseBtn.grid(column=12, row=0, padx=5, pady=5, sticky="E")
        self.loadPhoto = tk.PhotoImage(file='Imports/reload.gif')
        self.loadBtn = ttk.Button(self, image=self.loadPhoto,
                    command=lambda: self.event_generate("<<load_model_image>>")
        self.loadBtn.grid(column=13, row=0, padx=5, pady=5, sticky="E")
        
        # Pixel scale
        self.pixScaLab = ttk.Label(self, text="Pixel Scale (arcsec):")
        self.pixScaLab.grid(column=7, row=1, columnspan=1, padx=5, pady=5,
                            sticky="E")
        self.pixScale_asec = tk.DoubleVar()
        self.pixScale_asec.trace("w", lambda dummy1, dummy2, dummy3:
                            self.event_generate("<<pixscale_changed>>"))
        self.pixScale_asec.set(1.0)
        self.pixScaEnt = ttk.Entry(self, width=5,
                                   textvariable=self.pixScale_asec)
        self.pixScaEnt.grid(column=8, row=1, columnspan=1, padx=5, pady=5,
                            sticky="W")
        
        # Scale information
        self.extent = tk.StringVar()
        self.extent.set("")
        self.extentLab = ttk.Label(self, textvariable=self.extent)
        self.extentLab.grid(column=9, row=1, columnspan=5, padx=5, pady=5,
                            sticky="E")
        
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
        
        
#-----------------------------------------------------------------------------#
class StatusFrame(ttk.Frame):
    """Frame presenting status indicators lights and action buttons for the
    steps in the observation process."""

    # TODO: Redefine the grid labels for neatness.
    
    def __init__(self, parent, boxWidth=20, gapWidth=30, yPad=25):
        ttk.Frame.__init__(self, parent)
        self.parent = parent
        bgColour = ttk.Style().lookup("TFrame", "background")

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
        self.pltModFFTbtn = ttk.Button(self, text = "Plot Model FFT", width=18,
                        command=lambda: self.event_generate("<<plot_modFFT>>"))
        self.pltModFFTbtn.configure(state="disabled")
        self.canvas.create_window(x1, y2, window=self.pltModFFTbtn)
        self.pltPwrSpecbtn = ttk.Button(self, text = "Plot Power Spectrum",
                                        width=18,
                    command=lambda: self.event_generate("<<plot_powerspec>>"))
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
                    command=lambda: self.event_generate("<<plot_uvcoverage>>"))
        self.pltuvCovBtn.configure(state="disabled")
        self.canvas.create_window(x1, y2, window=self.pltuvCovBtn)
        self.pltElBtn = ttk.Button(self, text = "Plot Elevation",
                                   width=18,
                    command=lambda: self.event_generate("<<plot_elevation>>"))
        self.pltElBtn.configure(state="disabled")
        self.elBtnW = self.canvas.create_window(x1, y3, window=self.pltElBtn)
        
        # Draw the line seperating inputs from outputs
        x1 = self.xCentLst[1] + (self.xCentLst[2] - self.xCentLst[1]) / 1.7
        x2 = self.xCentLst[0] + (self.xCentLst[1] - self.xCentLst[0]) / 1.7
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
                    command=lambda: self.event_generate("<<do_observation>>"))
        self.obsBtn.configure(state="disabled")
        obsBtnW = self.canvas.create_window(x3, y3, window=self.obsBtn)
        
    def _draw_checkbox(self, xCent, yCent, size, tag, state=0, lw=3):
        """Draw a large checkbox on the canvas with a specified state."""

        r = float(size)/2.0
        
        # Draw box
        item = self.canvas.create_rectangle(xCent-r, yCent-r,
                                            xCent+r, yCent+r,
                                            outline="black", fill="white",
                                            width=1.0, tag="checkbox")
        self.canvas.addtag_withtag(tag, item)
        
        if state:
            # Draw green tickmark
            item = self.canvas.create_line(xCent-r+lw, yCent+r/2-lw,
                                           xCent-r/2+lw, yCent+r-lw,
                                           xCent+r-lw, yCent-r+lw,
                                           fill="green", capstyle=tk.ROUND,
                                           width=lw, tag="tick")
            self.canvas.addtag_withtag(tag, item)
        else:
            # Draw red X
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

    def set_state_by_dict(self, stateDict):
        """Use a dictionary to set the status of the indicator checkboxes."""
        
        for k, v in stateDict.items():
            self.set_state(k, v)
            
    def set_state(self, tag, state=0):
        """Re-draw the status indicator checkboxes with either red crosses or
        green ticks, depending on the state. Update the self.state list, which
        keeps a record of the current status of each step. Enable or disable
        plotting buttons as appropriate."""
        
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
            
#-----------------------------------------------------------------------------#
class InformationPanel(ttk.Frame):
    """Frame presenting the information on the selected arrays, input model
    and results"""

    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent

        # Scales from uv-coverage
        self.uvTitle = ttk.Label(self, text="Properties of uv-Coverage")
        self.uvTitle.grid(row=0, column=0, padx=0, pady=0, columnspan=2,
                          sticky="EW")
        
        self.uvScaleLab = ttk.Label(self, text=u"Scale Range:")
        self.uvScaleLab.grid(column=0, row=1, padx=5, pady=5, sticky="E")
        self.uvScale = tk.StringVar()
        self.uvScaleVal = ttk.Label(self, width=20, textvariable=self.uvScale)
        self.uvScaleVal.grid(column=1, row=1, padx=5, pady=5, sticky="W")

        self.uvRangeLab = ttk.Label(self, text=u"uv-Range:")
        self.uvRangeLab.grid(column=0, row=2, padx=5, pady=5, sticky="E")
        self.uvRange = tk.StringVar()
        self.uvRangeVal = ttk.Label(self, width=20, textvariable=self.uvRange)
        self.uvRangeVal.grid(column=1, row=2, padx=5, pady=5, sticky="W")
        
        #self.fovLab = ttk.Label(self, text=u"Field of View:")
        #self.fovLab.grid(column=0, row=2, padx=5, pady=5, sticky="E")
        #self.fov = tk.StringVar()
        #self.fovVal = ttk.Label(self, width=20, textvariable=self.fov)
        #self.fovVal.grid(column=1, row=2, padx=5, pady=5, sticky="W")

        sep = ttk.Separator(self, orient="vertical")
        sep.grid(column=2, row=0, rowspan=4, padx=5, pady=0, sticky="NS")
        
        # Scales from model FFT
        self.title = ttk.Label(self, text="Properties of Model FFT:")
        self.title.grid(column=3, row=0, padx=0, pady=0, columnspan=2,
                        sticky="EW")
        
        self.pixScaleLab = ttk.Label(self,
                                     text="FFT Pixel Scale:")
        self.pixScaleLab.grid(column=3, row=1, padx=(15,5), pady=5, sticky="E")
        self.pixScale = tk.StringVar()
        self.pixScaleVal = ttk.Label(self, width=20,
                                     textvariable=self.pixScale)
        self.pixScaleVal.grid(column=4, row=1, padx=5, pady=5, sticky="W")
        
        self.fftLimLab = ttk.Label(self, text="FFT Image Limit:")
        self.fftLimLab.grid(column=3, row=2, padx=(15,5), pady=5, sticky="E")
        self.fftLim = tk.StringVar()
        self.fftLimVal = ttk.Label(self, width=20, textvariable=self.fftLim)
        self.fftLimVal.grid(column=4, row=2, padx=5, pady=5, sticky="W")

    def update(self, parmDict):
        """Update the values using a dictionary."""
        
        p = parmDict
        if p["scaleMin_deg"] is None or p["scaleMin_deg"] is None:
            self.uvScale.set("")
        else:
            self.uvScale.set(u"%s to %s" % (ang2str(p["scaleMin_deg"]),
                                            ang2str(p["scaleMax_deg"])))
            
        if p["uvRngMin_lam"] is None or p["uvRngMax_lam"] is None:
            self.uvRange.set("")
        else:
            self.uvRange.set(u"%.3f to %.3f k\u03bb" %
                              (p["uvRngMin_lam"]/1000.0,
                               p["uvRngMax_lam"]/1000.0))
            
        #if "priBeamMax_deg" is None:
        #    self.fov.set("")
        #else:
        #    self.fov.set(ang2str(p["priBeamMax_deg"]))
            
        if p["pixScaleFFTX_lam"] is None or p["pixScaleFFTY_lam"] is None:
            self.pixScale.set("")
        else:
            self.pixScale.set(u"%.3f x %.3f k\u03bb" %
                              (p["pixScaleFFTX_lam"]/1000.0,
                               p["pixScaleFFTY_lam"]/1000.0))
        if p["fftScale_lam"] is None:
            self.fftLim.set("")
        else:
            self.fftLim.set(u"%.3f k\u03bb" % (p["fftScale_lam"]/1e3))
            

#-----------------------------------------------------------------------------#
class PlotFrame(ttk.Frame):
    """Frame showing the plots produced by the virtual interferometer."""

    def __init__(self, parent, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)

        # Track which plots are active
        self.pltDict = {"modelImg": 0,
                        "modelFFT": 0,
                        "uvCov": 0,
                        "beam": 0,
                        "obsFFT": 0,
                        "obsImg": 0}
        
        # Create the blank figure canvas and grid its tk canvas
        self.fig = Figure(figsize=(15.0, 9.0))
        self.figCanvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas = self.figCanvas.get_tk_widget()
        self.canvas.grid(column=0, row=0, padx=0, pady=0, sticky="NSEW")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)
        
        # Add the axes for the plots
        self.modAx = None
        self.modAx = None
        self.uvAx = None
        self.beamAx = None
        self.fftAx = None
        self.gridAx = None
        self.obsAx = None
        self.modAx = self.clear_model_image()
        self.uvAx = self.clear_uvcov()
        self.beamAx = self.clear_beam()
        self.fftAx = self.clear_model_fft()
        self.gridAx = self.clear_obs_fft()
        self.obsAx = self.clear_obs_image()
        self.show()
        
        # Add the information panel
        self.infoPanel = InformationPanel(self)
        self.infoPanel.grid(column=0, row=1, columnspan=3, padx=5, pady=5,
                            sticky="EW")

    def _plot_image(self, ax, imgArr, title=""):
        ax.clear()
        ax.imshow(imgArr, cmap=plt.cm.cubehelix,
                  interpolation="nearest", origin="lower")
        plt.setp(ax.get_yticklabels(), visible=False)
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.set_title(title)
        ax.set_aspect('equal')

    def _plot_fft(self, ax, fftArr, limit=None, title=""):
        if limit:
            extent=[-limit, limit, -limit, limit]
        else:
            extent=None
        ax.clear()
        ax.imshow(np.abs(fftArr), norm=LogNorm(), cmap=plt.cm.cubehelix,
                  interpolation="nearest", origin="lower", extent=extent)
        ax.set_title(title)
        ax.set_xlabel(u"u (k$\lambda$)")
        ax.set_ylabel(u"v (k$\lambda$)")
        ax.set_aspect('equal')#, 'datalim')
        
    def _plot_uvcov(self, ax, arrsSelected):
        colLst=["r", "b", "g", "m", "c", "y", "k"]
        
        ax.clear()
        for i, e in enumerate(arrsSelected):
            u = e["uArr_lam"]
            v = e["vArr_lam"]
            ax.scatter(x=u/1000, y=v/1000, marker=".", edgecolor='none', s=2,
                       color=colLst[i%len(colLst)])
            ax.scatter(x=-u/1000, y=-v/1000, marker=".", edgecolor='none', s=2,
                       color=colLst[i%len(colLst)])
        ax.set_xlabel(u"u (k$\lambda$)")
        ax.set_ylabel(u"v (k$\lambda$)")
        ax.set_aspect('equal', 'datalim')
        ax.margins(0.02)
        
    def plot_model_image(self, imgArr):
        self._plot_image(self.modAx, imgArr)
        self.pltDict["modelImg"] = 1
        
    def clear_model_image(self):
        if self.modAx:
            self.modAx.clear()
        else:
            self.modAx = self.fig.add_subplot(231)
        plt.setp(self.modAx.get_yticklabels(), visible=False)
        plt.setp(self.modAx.get_xticklabels(), visible=False)
        self.pltDict["modelImg"] = 0
        return self.modAx

    def plot_model_fft(self, fftArr, limit=None):
        self._plot_fft(self.fftAx, fftArr, limit=limit)
        self.pltDict["modelFFT"] = 1
    
    def clear_model_fft(self):
        if self.fftAx:
            self.fftAx.clear()
        else:
            self.fftAx = self.fig.add_subplot(234)
        plt.setp(self.fftAx.get_yticklabels(), visible=False)
        plt.setp(self.fftAx.get_xticklabels(), visible=False)
        self.pltDict["modelFFT"] = 0
        return self.fftAx
    
    def plot_uvcov(self, arrsSelected):
        self._plot_uvcov(self.uvAx, arrsSelected)
        self.pltDict["uvCov"] = 1
    
    def clear_uvcov(self):
        if self.uvAx:
            self.uvAx.clear()
        else:
            self.uvAx = self.fig.add_subplot(232)
        plt.setp(self.uvAx.get_yticklabels(), visible=False)
        plt.setp(self.uvAx.get_xticklabels(), visible=False)
        self.pltDict["uvCov"] = 0
        return self.uvAx
        
    def plot_beam(self, beamArr):
        self._plot_image(self.beamAx, beamArr)
        self.pltDict["beam"] = 1
    
    def clear_beam(self):
        if self.beamAx:
            self.beamAx.clear()
        else:
            self.beamAx = self.fig.add_subplot(233)
        plt.setp(self.beamAx.get_yticklabels(), visible=False)
        plt.setp(self.beamAx.get_xticklabels(), visible=False)
        self.pltDict["beam"] = 0
        return self.beamAx

    def plot_obs_fft(self, fftArr, limit=None):
        self._plot_fft(self.gridAx, fftArr, limit=limit)
        self.pltDict["obsFFT"] = 1

    def clear_obs_fft(self):
        if self.gridAx:
            self.gridAx.clear()
        else:
            self.gridAx = self.fig.add_subplot(235)
        plt.setp(self.gridAx.get_yticklabels(), visible=False)
        plt.setp(self.gridAx.get_xticklabels(), visible=False)
        self.pltDict["obsFFT"] = 0
        return self.gridAx
    
    def plot_obs_image(self, obsImgArr):
        self._plot_image(self.obsAx, obsImgArr)
        self.pltDict["obsImg"] = 1
    
    def clear_obs_image(self):
        if self.obsAx:
            self.obsAx.clear()
        else:
            self.obsAx = self.fig.add_subplot(236)
        plt.setp(self.obsAx.get_yticklabels(), visible=False)
        plt.setp(self.obsAx.get_xticklabels(), visible=False)
        self.pltDict["obsImg"] = 0
        return self.obsAx
        
    def show(self):
        self.fig.subplots_adjust(left=0.05, right=0.97, top=0.97, bottom=0.07,
                                  wspace=0.16, hspace=0.16)
        self.figCanvas.show()

    def clear_by_dict(self, stateDict):
        """Use a dictionary toclear downstream plots based on state."""

        # Track which plots need to be cleared so we don't make multiple calls
        self.clearDict = {"modelImg": 0, "modelFFT": 0, "uvCov": 0,
                          "beam": 0, "obsFFT": 0, "obsImg": 0}

        # Clear the downstream plots base on process status and plot status
        if not stateDict["statusObsDone"]:
            self.clearDict["obsImg"] =  self.pltDict["obsImg"]
        if not stateDict["statusuvGrid"]:
            self.clearDict["obsFFT"] =  self.pltDict["obsFFT"]
            self.clearDict["obsImg"] =  self.pltDict["obsImg"]
            self.clearDict["beam"] =  self.pltDict["beam"]
        if not stateDict.get("statusModelFFT", 1):
            self.clearDict["modelFFT"] =  self.pltDict["modelFFT"]
            self.clearDict["obsFFT"] =  self.pltDict["obsFFT"]
            self.clearDict["obsImg"] =  self.pltDict["obsImg"]
            self.clearDict["beam"] =  self.pltDict["beam"]
        if not stateDict.get("statusSelection", 1):
            self.clearDict["uvCov"] =  self.pltDict["uvCov"]
            self.clearDict["obsFFT"] =  self.pltDict["obsFFT"]
            self.clearDict["obsImg"] =  self.pltDict["obsImg"]
            self.clearDict["beam"] =  self.pltDict["beam"]
        if not stateDict.get("statusuvCalc", 1):
            self.clearDict["uvCov"] =  self.pltDict["uvCov"]
            self.clearDict["obsFFT"] =  self.pltDict["obsFFT"]
            self.clearDict["beam"] =  self.pltDict["beam"]
            self.clearDict["obsImg"] =  self.pltDict["obsImg"]
        if not stateDict.get("statusBeam", 1):
            self.clearDict["beam"] =  self.pltDict["beam"]
            self.clearDict["obsImg"] =  self.pltDict["obsImg"]
        if not stateDict.get("statusModel", 1):
            self.clearDict["modelImg"] =  self.pltDict["modelImg"]
            self.clearDict["modelFFT"] =  self.pltDict["modelFFT"]
            self.clearDict["obsImg"] =  self.pltDict["obsImg"]
            self.clearDict["beam"] =  self.pltDict["beam"]

        # Clear the relevant plots and show 
        if self.clearDict["modelImg"]:
            self.clear_model_image()
        if self.clearDict["modelFFT"]:
            self.clear_model_fft()
        if self.clearDict["uvCov"]:
            self.clear_uvcov()
        if self.clearDict["beam"]:
            self.clear_beam()
        if self.clearDict["obsFFT"]:
            self.clear_obs_fft()
        if self.clearDict["obsImg"]:
            self.clear_obs_image()
        self.show()


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
    
