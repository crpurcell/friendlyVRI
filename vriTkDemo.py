#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     vriTKdemo.py                                                      #
#                                                                             #
# PURPOSE:  A virtual interferometer application written in Tkinter.          #
#                                                                             #
# REQUIRED: Requires numpy, tkinter, matplotlib and PIL/PILLOW                #
#                                                                             #
# CREDITS:  Cormac R. Purcell (cormac.purcell at mq.edu.au)                   #
#           Roy Truelove  (Macquarie University)                              #
#                                                                             #
# MODIFIED: 07-Jul-2018 by C.Purcell                                          #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
# App (class)        ... class containing the main application logic          #
#     _applicationExit                                                        #
#     _show_textfile                                                          #
#     _show_lone_figure                                                       #
#     _update_status                                                          #
#     _on_select_config                                                       #
#     _on_sel_change                                                          #
#     _on_obsparm_change                                                      #
#     _on_pixscale_change                                                     #
#     _on_plot_modFFT                                                         #
#     _on_plot_uvcov                                                          #
#     _on_plot_elevation                                                      #
#     _on_plot_scalehist                                                      #
#     _on_load_model                                                          #
#     _on_do_observation                                                      #
#     _on_show_results                                                        #
#                                                                             #
# ModelSelector      ... class defining the model selection interface         #
#     _handler_image_selected                                                 #
#     _handler_capture_photo                                                  #
#     _round_pixscale                                                         #
#     _round_decscale                                                         #
#                                                                             #
# ArraySelector      ... class defining the array selection interface         #
#     _handler_add_button                                                     #
#     _handler_clear_button                                                   #
#     _handler_clear_all_button                                               #
#                                                                             #
# ArrayScanner       ... class defining the array scanning interface          #
#     _handler_scan_button                                                    #
#     _handler_flat_button                                                    #
#     _handler_save_button                                                    #
#     _show_control_window                                                    #
#                                                                             #
# ObsPanel          ... class exposing final observation inputs and plots     #
#     _handler_browse_button                                                  #
#     _handler_capture_photo                                                  #
#     _round_scale                                                            #
#                                                                             #
# StatusFrame        ... class defining status indicators and action buttons  #
#     _draw_checkbox                                                          #
#     set_state_by_dict                                                       #
#     set_state                                                               #
#                                                                             #
# InformationPanel   ... class showing derived properties of observation      #
#     update                                                                  #
#                                                                             #
# InformationPanelGraphic ...class showing derived properties as graphic      #
#     _draw_block                                                             #
#     _world2canvas                                                           #
#     _set_scale                                                              #
#     update                                                                  #
#                                                                             #
# PlotFrame          ... class defining the plotting window                   #
#     _show_control_window                                                    #
#     plot_image                                                              #
#     plot_fft                                                                #
#     plot_uvcov                                                              #
#     show                                                                    #
#     clear_by_state                                                          #
#                                                                             #
# MPLnavToolbar      ... subclass the MPL nav toolbar to disable readout      #
#     set_message                                                             #
#                                                                             #
#=============================================================================#
#                                                                             #
# The MIT License (MIT)                                                       #
#                                                                             #
# Copyright (c) 2017 Cormac R. Purcell and Roy Truelove                       #
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
    from tkinter.scrolledtext import ScrolledText as tkScrolledText
import numpy as np
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from distutils.version import LooseVersion
if LooseVersion(mpl.__version__) < LooseVersion("2.2"):
    from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg \
        as NavigationToolbar2Tk
else:
    from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk

# Webcam library
try:
    import cv2
    hasCV2 = True
except ImportError:
    hasCV2 = False

# Disable cv2 use on Mac OS because of buggy implementation
if sys.platform=="darwin":
    hasCV2 = False
    
from Imports.util_tk import *
from Imports.vriCalc import *
from Imports.util_scan import *


#-----------------------------------------------------------------------------#
class App(ttk.Frame):
    """Class defining the Virtual Radio Interferometer application.

    This class creates a root window used to set observation parameters, choose
    and accumulate array configurations and load a model image. A secondary 
    Toplevel window is used to display the input, output and intermediate
    Fourier transforms and images."""
    
    def __init__(self, parent, bgColour=None, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        self.bgColour=bgColour
        self.parent.title("Friendly VRI: Control Window")
        self.obsManager = None
        
        # Create the notebook and grid
        self.nb = ttk.Notebook(self, name="notebook")
        self.nb.enable_traversal()
        self.nb.grid(row=0, column=0, padx=5, pady=5, sticky="NSEW")
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

        # Add the panels in order
        self._add_model_tab()
        self._add_scan_tab()
        self._add_array_tab()
        self._add_observe_tab()
        self._add_results_tab()
        
        # Create the display window
#        self.dispWin = tk.Toplevel(self)
#        self.dispWin.title("Friendly VRI: Plot Window")
#        self.dispWin.protocol("WM_DELETE_WINDOW", self._applicationExit)
#        self.dispWin.columnconfigure(0, weight=1)
#        self.dispWin.rowconfigure(0, weight=1)

        # Draw the display interface
#        self.pltFrm = PlotFrame(self.dispWin, bgColour=self.bgColour)
#        self.pltFrm.grid(row=0, column=0, padx=0, pady=0, sticky="NSEW")

        # Set focus to the main window and bring to the fore
        self.focus_force()
        self.lift()
        
        # DEBUG
        if False:
            self.testWin = tk.Toplevel(self)
            self.testWin.title(" TEST WINDOW ")
            self.testWin.protocol("WM_DELETE_WINDOW", self._applicationExit)
            self.testWin.resizable(True, True)
            self.testWin.columnconfigure(0, weight=1)
            self.testWin.rowconfigure(0, weight=1)
            self.a = ObsPanel(self.testWin)
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
        #self.parent.bind("<<scan_array>>",
        #          lambda event : self._on_scan_array(event))
        self.parent.bind("<<array_scanned>>",
                  lambda event : self._on_scan_array(event))
        self.parent.bind("<<config_in_selected>>",
                  lambda event : self._on_select_config(event))
        self.parent.bind("<<selection_changed>>",
                  lambda event : self._on_sel_change(event))
#        self.parent.bind("<<obsparm_changed>>",
#                  lambda event : self._on_obsparm_change(event))    
        self.parent.bind("<<pixscale_changed>>",
                  lambda event : self._on_pixscale_change(event))
        self.parent.bind("<<plot_modFFT>>",
                  lambda event : self._on_plot_modFFT(event))
        self.parent.bind("<<plot_uvcoverage>>",
                  lambda event : self._on_plot_uvcov(event))
        self.parent.bind("<<plot_elevation>>",
                  lambda event : self._on_plot_elevation(event))
        self.parent.bind("<<plot_scalehist>>",
                  lambda event : self._on_plot_scalehist(event))
        self.parent.bind("<<do_observation>>",
                  lambda event : self._on_do_observation(event))
        self.parent.bind("<<show_results>>",
                  lambda event : self._on_show_results(event))
        self.parent.bind("<<load_model_image>>",
                       lambda event : self._on_load_model(event))
        
        # Force a minimum size on the windows & set resize properties
        self.parent.update()
        self.parent.minsize(self.parent.winfo_width(),
                            self.parent.winfo_height())
        #self.parent.resizable(False, False)
#        self.dispWin.update()
#        self.dispWin.minsize(self.dispWin.winfo_width(),
#                             self.dispWin.winfo_height())
#        self.dispWin.resizable(True, True)        

    def _add_model_tab(self):
        tabFrm = tk.Frame(self.nb)
        tabFrm.rowconfigure(0, weight=1)
        tabFrm.columnconfigure(0, weight=1)
        self.modelSelector = ModelSelector(tabFrm, bgColour=self.bgColour)
        self.modelSelector.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrm, text="    1) Choose Model Image    ", padding=5)

    def _add_scan_tab(self):
        tabFrm = tk.Frame(self.nb)
        tabFrm.rowconfigure(0, weight=1)
        tabFrm.columnconfigure(0, weight=1)
        self.arrayScanner = ArrayScanner(tabFrm, bgColour=self.bgColour)
        self.arrayScanner.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrm, text="    2) Create New Telescope Array    ",
                    padding=5)
        
    def _add_array_tab(self):
        tabFrm = tk.Frame(self.nb)
        tabFrm.rowconfigure(0, weight=1)
        tabFrm.columnconfigure(0, weight=1)
        self.selector = ArraySelector(tabFrm, bgColour=self.bgColour)
        self.selector.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrm, text="    3) Select Array Configurations    ",
                    padding=5)

    def _add_observe_tab(self):
        tabFrm = tk.Frame(self.nb)
        tabFrm.rowconfigure(0, weight=1)
        tabFrm.columnconfigure(0, weight=1)
        self.inputs = ObsPanel(tabFrm)
        self.inputs.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrm, text="    3) Do Observations    ", padding=5)

    def _add_results_tab(self):
        tabFrm = tk.Frame(self.nb)
        tabFrm.rowconfigure(0, weight=1)
        tabFrm.columnconfigure(0, weight=1)
        self.pltFrm = PlotFrame(tabFrm, bgColour=self.bgColour)
        self.pltFrm.grid(column=0, row=0, sticky="NSEW")
        self.nb.add(tabFrm,  text="    4) View Results    ", padding=5)
        
    def _applicationExit(self):
        """Exit the application cleanly if the window is closed."""
        
        self.parent.destroy()
        
    def _show_textfile(self, fileName, title=""):
        """Show a text file in a new window."""
        
        self.helpWin = tk.Toplevel(self, background=self.bgColour)
        self.helpWin.title(title)
        self.helpTxt = tkScrolledText(self.helpWin, width=80,
                                      font=fontFixed)
        self.helpTxt.config(state="normal")
        with open(fileName,'r') as f:
            text = f.read()
        self.helpTxt.insert('1.0', text)
        self.helpTxt.config(state="disabled")
        self.helpTxt.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW")
        self.closeBtn = ttk.Button(self.helpWin, text='Close',
                                   command=self.helpWin.destroy)
        self.closeBtn.grid(column=0, row=1, padx=5, pady=5, sticky="E")
        self.helpWin.rowconfigure(0, weight=1)
        self.helpWin.columnconfigure(0, weight=1)
        
    def _show_lone_figure(self, fig, title="Plot Window"):
        """Show a matplotlib figure in a new window."""
        
        self.figWin = tk.Toplevel(self, background=self.bgColour)
        self.figWin.title(title)
        figCanvas = FigureCanvasTkAgg(fig, master=self.figWin)
        loneCan = figCanvas.get_tk_widget()
        loneCan.configure(highlightthickness=0)
        loneCan.configure(background=self.bgColour)
        loneCan.grid(column=0, row=0, columnspan=2, padx=0, pady=0,
                     sticky="NSEW")
        tbarFrm = ttk.Frame(self.figWin)
        toolbar = MPLnavToolbar(figCanvas, tbarFrm)
        tbarFrm.grid(column=0, row=1, padx=5, pady=5, sticky="W")
        closeBtn = ttk.Button(self.figWin, text='Close',
                              command=self.figWin.destroy)
        closeBtn.grid(column=1, row=1, padx=5, pady=5, sticky="E")
        self.figWin.rowconfigure(0, weight=1)
        self.figWin.columnconfigure(0, weight=1)

    def _update_status(self):
        """Update the status of the user interface, including the checkbox
        indicators, plot axes and information panels. The status of each step
        is queried from the instance of the observationManager class."""

        # Query the status and set the indicators
        stateDict = self.obsManager.get_status()
        self.inputs.statFrm.set_state_by_dict(stateDict)
        
        # Clear the plots based on the status
        self.pltFrm.clear_by_state(stateDict)
        
        # Query the scales and update the information panel
#        parmDict = self.obsManager.get_scales()
#        self.pltFrm.infoPanel.update(parmDict)

        # Force an update of the GUI
        root.update()
        
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
        text = u"{:.4f}\u00B0".format(d["latitude_deg"])
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
                                      self.inputs.dec_deg.get())

        # Reset the array and hour-angle selections
        self.obsManager.clear_all_selections()
        for selection in self.selector.configOutTab.get_all_text():
            key = "_".join(selection[:2])
            haStart = float(selection[2])
            haEnd = float(selection[3])
            sampRate = float(selection[4])
            self.obsManager.select_array(key, haStart, haEnd, sampRate)

        # Update the elevation and scale plots
        self.event_generate("<<plot_elevation>>")
        self.event_generate("<<plot_scalehist>>")
            
        # Update the status
        self._update_status()

    def _on_obsparm_change(self, event=None):
        """When the observation parameters change in the GUI, reset the 
        calculations."""

        # Reset the common parameters
        try:
            self.obsManager.set_obs_parms(float(self.inputs.freq_MHz.get()),
                                    float(self.inputs.dec_deg.get()))
        except Exception:
            pass

        # Update the status
        self._update_status()
        
    def _on_pixscale_change(self, event=None):
        """When the pixel scale is changed in the GUI, clear the FFT plot."""

        # Only operate if a model is already loaded
        stateDict = self.obsManager.get_status()
        if not stateDict["statusModel"]:
            return

        # Reset the pixel scale
        try:
            pixScale_asec = float(self.modelSelector.pixScale_asec.get())
            self.obsManager.set_pixscale(pixScale_asec)
        except Exception:
            self.modelSelector.extent.set("")
            return
        
        # Set the extent label in the load frame
        pixScaleImg_deg = self.obsManager.pixScaleImg_asec / 3600.0
        nX = self.obsManager.nX
        nY = self.obsManager.nY
        text = " %d x %d pix  /  %s x %s" % (nX, nY,
                                             ang2str(nX * pixScaleImg_deg),
                                             ang2str(nY * pixScaleImg_deg))
        self.modelSelector.extent.set(text)
        
        # Only update the status if this is the first change after processing
        if stateDict["statusModelFFT"]:
            self._update_status()
        
    def _on_plot_modFFT(self, event=None):
        """Show the FFT of the model image."""
            
        # Invert the model image
        self.obsManager.invert_model()
        
        # Plot the model FFT
        parmDict = self.obsManager.get_scales()
        lim_kl = parmDict["fftScale_lam"]/1e3
        self.pltFrm.plot_fft("modelFFT", self.obsManager.modelFFTarr,
                             limit=lim_kl, title="Model FFT")
        
        # Update the status
        self._update_status()
        
    def _on_plot_uvcov(self, event=None):
        """Plot the uv-coverage for all selected array configurations"""

        # Calculate the uv-coverage for the selected observation parameters
        self.obsManager.calc_uvcoverage()

        # Plot the uv-coverage in the display window
        self.pltFrm.plot_uvcov("uvCov", self.obsManager.arrsSelected,
                               title="uv-Coverage")
        
        # Update the status
        self._update_status()
        
    def _on_plot_elevation(self, event=None):
        """Create a plot showing the elevation of the source as seen from 
        the currently selected telescopes."""
        
        # Query the selected array configurations
        selTab = self.obsManager.get_selected_arrays()
        telescopeLst = []
        if not selTab is None:
            telescopeLst = set(selTab["telescope"])

        # Calculate each of the elevation curves
        elCurveDict = {}
        for i, e in enumerate(telescopeLst):
            label = e.decode("utf-8")
            haArr_hr, elArr_deg = self.obsManager.calc_elevation_curve(e)
            elCurveDict[label] = [haArr_hr, elArr_deg]

        # Call the plotting function. Will initialise the axes if blank.
        self.inputs.plot_elevation(elCurveDict)

    def _on_plot_scalehist(self, event=None):
        """Create a plot showing a histogram of the baseline lengths and
        compare to the power in the source model."""
        
        # Query the selected array configurations
        selTab = self.obsManager.get_selected_arrays()
        configLst = []
        if not  selTab is None:
            configLst = set(selTab["key"])

        # Get the baseline lengths
        lBaseArr_lam = []
        for k in configLst:
            lBase_m = self.obsManager.get_baseline_lengths(key=k)
            lBase_lam = lBase_m / self.obsManager.lambda_m
            lBaseArr_lam = np.append(lBaseArr_lam, lBase_lam)
        lBaseArr_lam = np.unique(lBaseArr_lam)
        
        # Call the plotting function. Will initialise the axes if blank.
        self.inputs.plot_scales(lBaseArr_lam)
        
    def _on_load_model(self, event=None):
        """Load a model image into the observation manager"""
        
        # Load the model image
        modelFile = self.modelSelector.modelPath
        pixScale_asec = self.modelSelector.pixScale_asec.get()
        self.obsManager.load_model_image(modelFile, pixScale_asec)

        # Set the extent label in the load frame
        pixScaleImg_deg = self.obsManager.pixScaleImg_asec / 3600.0
        nX = self.obsManager.nX
        nY = self.obsManager.nY
        text = " %d x %d pix  /  %s x %s" % (nX, nY,
                                             ang2str(nX * pixScaleImg_deg),
                                             ang2str(nY * pixScaleImg_deg))
        self.modelSelector.extent.set(text)
        
        # Plot the model image
        self.pltFrm.plot_image("modelImg", self.obsManager.modelImgArr,
                               title="Model Image")
        
        # Update the status
        self._update_status()
        
    def _on_do_observation(self, event=None):
        """Perform the bulk of the observing steps"""
        
        # Calculate the uv-coverage if not cached
        stateDict = self.obsManager.get_status()
        if not stateDict['statusuvCalc']:
            self.obsManager.calc_uvcoverage()

            # Plot the uv-coverage
            self.pltFrm.plot_uvcov("uvCov", self.obsManager.arrsSelected,
                                   title="uv-Coverage")
            
            # Update the status
            self._update_status()
        
        # Calculate the Fourier transform of the model if not cached
        stateDict = self.obsManager.get_status()
        if not stateDict['statusModelFFT']:
            self.obsManager.invert_model()
        
            # Plot the model FFT
            parmDict = self.obsManager.get_scales()
            lim_kl = parmDict["fftScale_lam"]/1e3
            self.pltFrm.plot_fft("modelFFT", self.obsManager.modelFFTarr,
                                 limit=lim_kl, title="Model FFT")
            
            # Update the status
            self._update_status()
            
        # Grid the uv-coverage to make a mask
        self.obsManager.grid_uvcoverage()

        # Show the observed FFT
        parmDict = self.obsManager.get_scales()
        lim_kl = parmDict["fftScale_lam"]/1e3
        self.pltFrm.plot_fft("obsFFT", self.obsManager.obsFFTarr, limit=lim_kl,
                             title="Observed FFT")
        #self.pltFrm.plot_fft("obsFFT", self.obsManager.uvCntArr, limit=lim_kl,
        #                     title="Observed FFT")
        
        # Update the status
        self._update_status()
        
        # Calculate the PSF
        self.obsManager.calc_beam()
        
        # Show the synthesised beam      
        self.pltFrm.plot_image("beam", np.abs(self.obsManager.beamArr),
                               title="Synthesised Beam", pRng=(-0.1, 0.5))
        
        # Update the status
        self._update_status()
        
        # Invert the observed FFT
        self.obsManager.invert_observation()
        
        # Show the observed image
        self.pltFrm.plot_image("obsImg", np.abs(self.obsManager.obsImgArr),
                               title="Observed Image")

        # Update the status
        self._update_status()

    def _on_show_results(self, event=None):
        """Raise the focus of the plotting window."""
        pass
#        self.dispWin.focus_force()
#        self.dispWin.lift()
        
    def _on_scan_array(self, event=None):
        """Load the latest scan of the antenna plate."""
        
        # Load the custom array file into the observation manager
        tmpFile = "arrays/custom.config"
        self.obsManager._load_one_array(tmpFile)
        
        # Append to the list of available array configurations
        vals = self.obsManager.arrsAvailable.values()[-1]
        configLst = zip([vals["telescope"]],
                        [vals["config"]])
        self.selector.configInTab.insert_rows(configLst,
                                              ("Telescope", "Array"))
        
        # Clean up
        os.remove(tmpFile)
        
        
#-----------------------------------------------------------------------------#
class ModelSelector(ttk.Frame):
    """Interface to preview and select a model image or capture an image 
    using the webcam."""
    
    def __init__(self, parent, bgColour=None, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        if bgColour is None:
            bgColour = ttk.Style().lookup("TFrame", "background")

        # Set the expansion properties
        # padding cols = [3, 7]
        self.rowconfigure(0, weight=1)
        self.columnconfigure(3, weight=1)
        self.columnconfigure(5, weight=1)

        # Preview panel showing the model images
        self.imgPreview = ImageBrowser(self, virtEvent="<<model_selected>>",
                                       path="models")
        self.imgPreview.grid(column=0, row=0, columnspan=6, sticky="NSEW")
        self.imgPreview.bind("<<model_selected>>",
                         lambda event : self._handler_image_selected(event))

        # Webcam button
        self.cameraPhoto = tk.PhotoImage(file='Imports/camera.gif')
        self.cameraBtn = ttk.Button(self, image=self.cameraPhoto, width=16,
                                    compound="left", text=" Capture Photo",
                                    command=self._handler_capture_photo)
        self.cameraBtn.grid(column=0, row=1, rowspan=2, padx=5, pady=5,
                            sticky="NSW")
        
        # Image details
        self.fileLab = ttk.Label(self, text="Chosen Image:")
        self.fileLab.grid(column=1, row=1, padx=(15,5), pady=5, sticky="E")
        self.modelFile = tk.StringVar()
        self.fileEnt = ttk.Entry(self, width=33, textvariable=self.modelFile)
        self.fileEnt.configure(state="readonly")
        self.fileEnt.grid(column=2, row=1, columnspan=1, padx=5, pady=5,
                          sticky="EW")
        self.extLab = ttk.Label(self, text="Image Scale:")
        self.extLab.grid(column=1, row=2, padx=(15,5), pady=5, sticky="E")
        self.extent = tk.StringVar()
        self.extent.set("")
        self.extentLab = ttk.Label(self, textvariable=self.extent)
        self.extentLab.grid(column=2, row=2, columnspan=1, padx=5, pady=5,
                            sticky="E")
        
        # Pixel scale slider
        self.pixScaLab = ttk.Label(self,
                                   text="Pixel Scale (arcseconds):")
        self.pixScaLab.grid(column=4, row=1, padx=5, pady=5, sticky="W")
        self.pixScale_asec = tk.DoubleVar()
        self.pixScale_asec.set(1.0)
        self.pixScaValLab = ttk.Label(self, textvariable=self.pixScale_asec,
                                      width=5, anchor="e")
        self.pixScaValLab.grid(column=5, row=1, padx=5, pady=5, sticky="W")
        self.pixScale = ttk.Scale(self, from_=0.1, to=10.0, length=450,
                                  variable=self.pixScale_asec,
                                  command=self._round_pixscale)
        self.pixScale.grid(column=4, row=2, columnspan=2, padx=5, pady=5,
                           sticky="EW")
        self.pixScale.bind('<Any-ButtonRelease-1>',
                        lambda e: self.event_generate("<<pixscale_changed>>"))
                
    def _handler_image_selected(self, event=None):
        
        modelPath = self.imgPreview.imagePath
        if not len(modelPath)==0:
            if os.path.exists(modelPath):
                s = os.path.split(modelPath)
                if len(s)==2:
                    self.modelFile.set(s[-1])
                    self.modelPath = modelPath
                    self.event_generate("<<load_model_image>>")
        
    def _handler_capture_photo(self):
        """Capture a photo using the webcam."""
        
        try:
            cam = cv2.VideoCapture()
            cam.open(0)
            #cam.set(3, int(self.resX.get()))
            #cam.set(4, int(self.resX.get()))
            for i in range(10):
                success, img = cam.read()
            cam.release()
            if success:
                cv2.imwrite("models/webcam.png", img)
                self.modelFile.set("webcam.png")
                self.modelPath = "models/webcam.png"
                self.imgPreview.display_image(self.modelPath)

                
                #self.event_generate("<<load_model_image>>")
        except Exception:
            pass
    
    def _round_pixscale(self, e=None):
        value = self.pixScale.get()
        value *= 10
        if int(value) != value:
            self.pixScale.set(round(value)/10)
            
    def _round_decscale(self, e=None):
        value = self.decScale.get()
        if int(value) != value:
            self.decScale.set(round(value))


#-----------------------------------------------------------------------------#
class ArraySelector(ttk.Frame):
    """Two multi-column listboxes an hour-angle slider and a 'Select' button
    that make up the array selection interface."""
    
    def __init__(self, parent, bgColour=None, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        if bgColour is None:
            bgColour = ttk.Style().lookup("TFrame", "background")
        
        # Set the expansion properties
        self.columnconfigure(0, weight=10)
        self.columnconfigure(7, weight=20)
        self.columnconfigure(8, weight=20)
        self.rowconfigure(2, weight=1)
        self.rowconfigure(3, weight=1)

        # Listbox showing the available array configurations
        self.configInTab = ScrolledTreeTab(self,
                                           virtEvent="<<config_in_selected>>")
        self.configInTab.name_columns(("Telescope", "Array"))
        self.configInTab.grid(column=0, row=0, columnspan=1, rowspan=7,
                              padx=5, pady=5, sticky="NSEW")

        # Hour angle slider
        self.haLab = ttk.Label(self, text="Hour Angle Range (hours):")
        self.haLab.grid(column=1, row=0, columnspan=5, padx=0, pady=5,
                        sticky="NW")
        self.haScale = DoubleScale(self, from_=-12.0, to=12.0,
                                   initLeft=-1.0, initRight=+1.0,
                                   tickIntMajor=6, tickIntMinor=1, width=400)
        self.haScale.grid(column=1, row=1, columnspan=5, padx=0, pady=5,
                          sticky="N")

        # Sampling cadence
        self.sampRtLab = ttk.Label(self, text="Sampling Cadence (s):")
        self.sampRtLab.grid(column=1, row=2, columnspan=2, padx=5, pady=5,
                            sticky="E")
        sampRtLst_s = ["10", "30", "60", "100", "300", "600", "1200", "1800",
                       "3600"]
        self.sampRt_s = tk.StringVar()
        self.sampRtComb = ttk.Combobox(self, state="readonly",
                                       textvariable=self.sampRt_s,
                                       values=sampRtLst_s, width=5)
        self.sampRtComb.current(4)
        self.sampRtComb.grid(column=3, row=2, columnspan=2, padx=5, pady=5,
                             sticky="W")

        # Fancy add button with strike-through arrow
        self.canvas = tk.Canvas(self, width=400, height=30,
                               background=bgColour, highlightthickness=0)
        self.canvas.create_line(10,15,390,15, width=2, arrow=tk.LAST,
                                arrowshape=(10,15,5), fill="black")
        self.addBtn = ttk.Button(self, text="Add", width=10,
                                 command=self._handler_add_button)
        self.canvas.create_window(200, 15, window=self.addBtn)
        self.canvas.grid(column=1, row=3, columnspan=5, padx=0, pady=5)
        
        # Scatter plot of antenna locations
        self.antPosPlot = ScatterPlot(self, width=440, height=400,
                                      axPad=(100,25,70,25), aspect="equal",
                                      pntSize=2.5)
        self.antPosPlot.grid(column=1, row=4, columnspan=5, rowspan=1,
                             padx=5, pady=5)
        txt = "Select an array configuration from the table\n"
        txt += "and click the 'Add' button to build up a list.\n"
        self.antPosPlot.display_message(txt)

        # Antenna Diameter and array latitude labels
        self.antDlab = ttk.Label(self, text=u"Antenna \u0398:   ")
        self.antDlab.grid(column=1, row=5, padx=0, pady=0, sticky="E")
        self.antD_m = tk.StringVar()
        self.antDval = ttk.Label(self, textvariable=self.antD_m)
        self.antDval.grid(column=2, row=5, padx=0, pady=0, sticky="W")
        self.antLlab = ttk.Label(self, text=u"Telescope \u03c6:   ")
        self.antLlab.grid(column=1, row=6, padx=0, pady=0, sticky="E")
        self.antL_deg = tk.StringVar()
        self.antLval = ttk.Label(self, textvariable=self.antL_deg)
        self.antLval.grid(column=2, row=6, padx=0, pady=0, sticky="W")
        
        # Minimum and maximum baseline labels
        self.minBaseLab = ttk.Label(self, text=u"Min Baseline:   ")
        self.minBaseLab.grid(column=4, row=5, padx=0, pady=0, sticky="E")
        self.minBase_km = tk.StringVar()
        self.minBaseVal = ttk.Label(self, textvariable=self.minBase_km)
        self.minBaseVal.grid(column=5, row=5, padx=0, pady=0, sticky="W")

        self.maxBaseLab = ttk.Label(self, text=u"Max Baseline:   ")
        self.maxBaseLab.grid(column=4, row=6, padx=0, pady=0, sticky="E")
        self.maxBase_km = tk.StringVar()
        self.maxBaseVal = ttk.Label(self, textvariable=self.maxBase_km )
        self.maxBaseVal.grid(column=5, row=6, padx=0, pady=0, sticky="W")
        
        # Listbox showing the selected array configurations
        self.configOutTab = ScrolledTreeTab(self,
                                        virtEvent="<<config_out_selected>>")
        self.configOutTab.name_columns(("Telescope", "  Array  ",
                                        "HA-Start", " HA-End ", "Cadence"))
        self.configOutTab.grid(column=7, row=0, columnspan=2, rowspan=6,
                               padx=5, pady=5, sticky="NSEW")
        
        # Clear selected and all buttons
        self.delBtn = ttk.Button(self, text="Clear Selected", width=20,
                                 command=self._handler_clear_button)
        self.delBtn.grid(column=7, row=6, padx=5, pady=5, sticky="EW" )
        self.delAllBtn = ttk.Button(self, text="Clear All", width=20,
                                    command=self._handler_clear_all_button)
        self.delAllBtn.grid(column=8, row=6, padx=5, pady=5, sticky="EW" )
        
    def _handler_add_button(self):
        """Add the selected configuration to the list box"""
        
        selConf = self.configInTab.get_text_selected()
        if selConf is not None:
            selConf = list(selConf)
            selConf.append(self.haScale.valueLeft.get())
            selConf.append(self.haScale.valueRight.get())
            selConf.append(self.sampRt_s.get())
            self.configOutTab.insert_rows(
                [selConf], ("Telescope", "  Array  ", "HA-Start", " HA-End ",
                            "Cadence"))
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
class ArrayScanner(ttk.Frame):
    """Interface to scan the physical model of the array using webcam 2."""

    def __init__(self, parent, bgColour=None, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        if bgColour is None:
            bgColour = ttk.Style().lookup("TFrame", "background")
        self.X_m = None
        self.Y_m = None
        self.shape = None
        self.doFlat = None
        
        # Set the expansion properties
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        # Scanner settings
        self.scanFrm = ttk.Labelframe(self, text=" Scanner Settings ")
        self.scanFrm.grid(column=0, row=0, padx=5, pady=5, sticky="NSEW")
        # Resolution
        self.resLab = ttk.Label(self.scanFrm, text="Resolution:")
        self.resLab.grid(column=0, row=1, padx=5, pady=5, sticky="E")
        self.resX = tk.StringVar()
        self.resY = tk.StringVar()
        self.resX.set("1280")
        self.resY.set("800")
        self.resXEnt = ttk.Entry(self.scanFrm, width=5,
                                  textvariable=self.resX)
        self.resXEnt.grid(column=1, row=1, padx=5, pady=5, sticky="EW")
        self.resYEnt = ttk.Entry(self.scanFrm, width=5,
                                  textvariable=self.resY)
        self.resYEnt.grid(column=2, row=1, padx=5, pady=5, sticky="EW")
        # Crop size
        self.cropLab = ttk.Label(self.scanFrm, text="Crop Size:")
        self.cropLab.grid(column=0, row=2, padx=5, pady=5, sticky="E")
        self.cropX = tk.StringVar()
        self.cropY = tk.StringVar()
        self.cropX.set("30")
        self.cropY.set("30")
        self.cropXEnt = ttk.Entry(self.scanFrm, width=5,
                                  textvariable=self.cropX)
        self.cropXEnt.grid(column=1, row=2, padx=5, pady=5, sticky="EW")
        self.cropYEnt = ttk.Entry(self.scanFrm, width=5,
                                  textvariable=self.cropY)
        self.cropYEnt.grid(column=2, row=2, padx=5, pady=5, sticky="EW")
        # Smoothing kernel
        self.kernLab = ttk.Label(self.scanFrm, text="Kernel Size:")
        self.kernLab.grid(column=0, row=3, padx=5, pady=5, sticky="E")
        self.kernSize = tk.StringVar()
        self.kernSize.set("0")
        self.kernEnt = ttk.Entry(self.scanFrm, width=5,
                                 textvariable=self.kernSize)
        self.kernEnt.grid(column=1, row=3, columnspan=2, padx=5, pady=5,
                          sticky="EW")
        # Detection threshold
        self.sigmaLab = ttk.Label(self.scanFrm, text="Threshold (sigma):")
        self.sigmaLab.grid(column=0, row=4, padx=5, pady=5, sticky="E")
        self.sigma = tk.StringVar()
        self.sigma.set("3")
        self.sigmaEnt = ttk.Entry(self.scanFrm, width=5,
                                  textvariable=self.sigma)
        self.sigmaEnt.grid(column=1, row=4, columnspan=2, padx=5, pady=5,
                          sticky="EW")
        # Min number pixels in island
        self.minPixLab = ttk.Label(self.scanFrm, text="Min # Pixels:")
        self.minPixLab.grid(column=0, row=5, padx=5, pady=5, sticky="E")
        self.minPix = tk.StringVar()
        self.minPix.set("200")
        self.minPixEnt = ttk.Entry(self.scanFrm, width=5,
                                   textvariable=self.minPix)
        self.minPixEnt.grid(column=1, row=5, columnspan=2, padx=5, pady=5,
                          sticky="EW")
        # Flat Field button
        self.flatBtn = ttk.Button(self.scanFrm, text="Flat Field", width=20,
                                  command=self._handler_flat_button)
#        self.flatBtn.grid(column=0, row=6, columnspan=3, padx=5, pady=5,
#                          sticky="SEW" )
        
        # Scan button
        self.scanBtn = ttk.Button(self.scanFrm, text="Scan Array", width=20,
                                  command=self._handler_scan_button)
        self.scanBtn.grid(column=0, row=7, columnspan=3, padx=5, pady=5,
                          sticky="SEW" )
        self.scanFrm.columnconfigure(0, weight=1)
        self.scanFrm.rowconfigure(6, weight=1)

        # Array settings
        self.arrFrm = ttk.Labelframe(self, text=" Custom Array Settings ")
        self.arrFrm.grid(column=0, row=1, padx=5, pady=5, sticky="NSEW")
        # Telescope name
        self.myScopeLab = ttk.Label(self.arrFrm, text="Telescope Name:")
        self.myScopeLab.grid(column=0, row=0, padx=5, pady=5, sticky="E")
        self.myScope = tk.StringVar()
        self.myScope.set("Custom_1")
        self.myScopeEnt = ttk.Entry(self.arrFrm, width=20, state="disabled",
                                    textvariable=self.myScope)
        self.myScopeEnt.grid(column=1, row=0, padx=5, pady=5, sticky="EW")
        # Array config
        self.myArrayLab = ttk.Label(self.arrFrm, text="Array Name:")
        self.myArrayLab.grid(column=0, row=1, padx=5, pady=5, sticky="E")
        self.myArray = tk.StringVar()
        self.myArray.set("Array")
        self.myArrayEnt = ttk.Entry(self.arrFrm, width=10,
                                    textvariable=self.myArray)
        self.myArrayEnt.grid(column=1, row=1, padx=5, pady=5, sticky="EW")
        # Latitude
        self.myLatLab = ttk.Label(self.arrFrm, text="Latitude (deg):")
        self.myLatLab.grid(column=0, row=2, padx=5, pady=5, sticky="E")
        latLst = [str(x) for x in range(90, -91, -10)]
        self.myLat = tk.StringVar()
        self.myLatComb = ttk.Combobox(self.arrFrm, textvariable=self.myLat,
                                      values=latLst, width=10)        
        self.myLatComb.current(11)
        self.myLatComb.grid(column=1, row=2, padx=5, pady=5, sticky="EW")
        # Antenna diameter
        self.myAntDiamLab = ttk.Label(self.arrFrm, text=u"Antenna \u0398 (m):")
        self.myAntDiamLab.grid(column=0, row=3, padx=5, pady=5, sticky="E")
        dLst = ["22"]
        self.myAntDiam = tk.StringVar()
        self.myAntDiamComb = ttk.Combobox(self.arrFrm, width=12, values=dLst,
                                         textvariable=self.myAntDiam)
        self.myAntDiamComb.current(0)
        self.myAntDiamComb.grid(column=1, row=3, padx=5, pady=5, sticky="EW")
        # Plate scale
        self.plateScaleLab = ttk.Label(self.arrFrm, text=u"Array Scale (m):")
        self.plateScaleLab.grid(column=0, row=4, padx=5, pady=5, sticky="E")
        self.plateScale = tk.StringVar()
        self.plateScale.set("2000.0")
        self.plateScaleEnt = ttk.Entry(self.arrFrm, width=12,
                                       textvariable=self.plateScale)
        self.plateScaleEnt.grid(column=1, row=4, padx=5, pady=5, sticky="EW")
        # Save button
        self.saveBtn = ttk.Button(self.arrFrm, text="Save Array", width=20,
                                  command=self._handler_save_button)
        self.saveBtn.configure(state="disabled")
        self.saveBtn.grid(column=0, row=5, columnspan=2, padx=5, pady=5,
                          sticky="SEW" )
        # Show control window
        self.showBtn = ttk.Button(self.arrFrm, text = "Show Control Window",
                                  width=20, command=self._show_control_window)
        self.showBtn.grid(column=0, row=6,  columnspan=3, padx=5, pady=5,
                          sticky="SEW")
        self.arrFrm.columnconfigure(0, weight=1)
        self.arrFrm.rowconfigure(5, weight=1)

        # Figure showing the scan of the antenna models
        self.fig = Figure(figsize=(8.0, 7.0), facecolor=bgColour)
        self.figCanvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas = self.figCanvas.get_tk_widget()
        self.canvas.configure(highlightthickness=0)
        self.canvas.configure(background=bgColour)
        self.canvas.grid(column=1, row=0, columnspan=1, rowspan=2,
                         padx=0, pady=0, sticky="NSEW")
        self.ax = self.fig.add_subplot(111)
        plt.setp(self.ax.get_yticklabels(), visible=False)
        plt.setp(self.ax.get_xticklabels(), visible=False)

    def _handler_scan_button(self):

        # Capture an image of the array via webcam 2
        cam = cv2.VideoCapture()
        cam.open(1)
        cam.set(3, int(self.resX.get()))
        cam.set(4, int(self.resX.get()))
        for i in range(10):
            success, img = cam.read()
        cam.release()
        if success:

            # Resize the image to make it easier to manage
            img = cv2.resize(img, (0,0), fx=0.5, fy=0.5)

            # Save scan to a png file
            cv2.imwrite("arrays/scan.png", img)
            
            # Detect the islands in the image
            imgName = "arrays/scan.png"
            if self.doFlat:
                flatImgName = "arrays/scan.png"
            else:
                flatImgName = None
            self.X_m, self.Y_m, self.shape = \
                scan_to_pixcoords("arrays/scan.png",
                                  eSize=float(self.kernSize.get()),
                                  threshold_sigma=float(self.sigma.get()),
                                  minPix=int(self.minPix.get()),
                                  cropX=img.shape[1]-int(self.cropX.get()),
                                  cropY=img.shape[0]-int(self.cropY.get()),
                                  flatImgName = flatImgName, 
                                  ax=self.ax)

            # Show the scanned figure and allow saving
            self.figCanvas.draw()

            # Enable/disable saving
            if len(self.X_m)>0:
                self.saveBtn.configure(state="enabled")
            else:
                self.saveBtn.configure(state="disabled")
    
    def _handler_flat_button(self):

        # Capture an image of the array via webcam 2
        cam = cv2.VideoCapture()
        cam.open(0)
        cam.set(3, int(self.resX.get()))
        cam.set(4, int(self.resX.get()))
        for i in range(10):
            success, img = cam.read()
        cam.release()
        if success:
            # Resize the image to make it easier to manage
            img = cv2.resize(img, (0,0), fx=0.5, fy=0.5)

            # Save scan to a png file
            cv2.imwrite("arrays/flat.png", img)
            self.doFlat = True
            
    def _handler_save_button(self):
        
        # Write a temporary array definition file
        tmpFile = "arrays/custom.config"
        write_arrayfile(tmpFile,
                        X_m=self.X_m,
                        Y_m=self.Y_m,
                        Nx=self.shape[-1],
                        Ny=self.shape[-2],
                        scale_m=float(self.plateScale.get()),
                        telescope=self.myScope.get(),
                        config=self.myArray.get(),
                        latitude_deg=float(self.myLat.get()),
                        diameter_m=float(self.myAntDiam.get()))

        # Increment the number of the telescope
        pre, i = self.myScope.get().rsplit("_")
        self.myScope.set(pre + "_" + str(int(i)+1))
        self.event_generate("<<array_scanned>>")

    def _show_control_window(self):
        """Set focus back to the main control window."""
        
        root.focus_force()
        root.lift()
        
#-----------------------------------------------------------------------------#
class ObsPanel(ttk.Frame):
    """Input settings for the observation (Model image, Declination of source,
    observing frequency, Robust weighting factor."""

    def __init__(self, parent, bgColour=None, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        if bgColour is None:
            bgColour = ttk.Style().lookup("TFrame", "background")

        # Set the expansion properties
        self.columnconfigure(2, weight=1)
        self.columnconfigure(5, weight=1)
        self.rowconfigure(0, weight=1)

        # Figure showing the elevation of the source
        self.elFig = Figure(figsize=(8.0, 7.0), facecolor=bgColour)
        self.elFigCanvas = FigureCanvasTkAgg(self.elFig, master=self)
        self.elCanvas = self.elFigCanvas.get_tk_widget()
        self.elCanvas.configure(highlightthickness=0)
        self.elCanvas.configure(background=bgColour)
        self.elCanvas.grid(column=0, row=0, columnspan=3, rowspan=1,
                           padx=0, pady=0, sticky="NSEW")
        self.elAx = self.elFig.add_subplot(111)
        plt.setp(self.elAx.get_yticklabels(), visible=False)
        plt.setp(self.elAx.get_xticklabels(), visible=False)
        
        # Source Declination slider
        self.decSrcLab = ttk.Label(self,
                                   text="Source Declination (degrees):")
        self.decSrcLab.grid(column=0, row=1, padx=5, pady=5, sticky="E")
        self.dec_deg = tk.DoubleVar()
        self.dec_deg.set(-20.0)
        self.decValLab = ttk.Label(self, textvariable=self.dec_deg, width=5,
                                   anchor="e")
        self.decValLab.grid(column=1, row=1, padx=5, pady=5, sticky="E")
        self.decScale = ttk.Scale(self, from_=-90, to=90, variable=self.dec_deg,
                                  command=self._round_decscale)
        self.decScale.grid(column=0, row=2, columnspan=3, padx=(5, 25),
                           pady=5, sticky="EW")
#        self.decScale.bind('<Any-ButtonRelease-1>',
#                        lambda e: self.event_generate("<<obsparm_changed>>"))
        self.decScale.bind('<Any-ButtonRelease-1>',
                        lambda e: self.event_generate("<<selection_changed>>"))
        
        # Figure showing a histogram of the baseline length
        self.blFig = Figure(figsize=(8.0, 7.0), facecolor=bgColour)
        self.blFigCanvas = FigureCanvasTkAgg(self.blFig, master=self)
        self.blCanvas = self.blFigCanvas.get_tk_widget()
        self.blCanvas.configure(highlightthickness=0)
        self.blCanvas.configure(background=bgColour)
        self.blCanvas.grid(column=3, row=0, columnspan=3, rowspan=1,
                           padx=0, pady=0, sticky="NSEW")
        self.blAx = self.blFig.add_subplot(111)
        self.blAxTop = self.blAx.twiny()
        #self.blAx.xaxis.set_major_locator(MaxNLocator(4))
        plt.setp(self.blAx.get_yticklabels(), visible=False)
        plt.setp(self.blAx.get_xticklabels(), visible=False)
        
        # Frequency slider
        self.freqLab = ttk.Label(self, text="Observing Frequency (MHz):")
        self.freqLab.grid(column=3, row=1, padx=(25, 5), pady=5, sticky="E")
        self.freq_MHz = tk.StringVar()
        self.freq_MHz.set(1420.0)
        self.freqValLab = ttk.Label(self, textvariable=self.freq_MHz, width=12,
                                   anchor="e")
        self.freqValLab.grid(column=4, row=1, padx=5, pady=5, sticky="E")
        self.freqScale = ttk.Scale(self, from_=500, to=110000,
                                   variable=self.freq_MHz,
                                   command=self._round_freqscale)
        self.freqScale.grid(column=3, row=2, columnspan=3, padx=(25,5),
                            pady=5, sticky="EW")
#        self.freqScale.bind('<Any-ButtonRelease-1>',
#                        lambda e: self.event_generate("<<obsparm_changed>>"))
        self.freqScale.bind('<Any-ButtonRelease-1>',
                        lambda e: self.event_generate("<<selection_changed>>"))

        # Add the controls and status icons
        self.statFrm = StatusFrame(self, boxWidth=22, gapWidth=165)
        self.statFrm .grid(column=0, row=4, columnspan=6, padx=5, pady=5,
                           sticky="EW")
        
    def _round_freqscale(self, e=None):
        value = self.freqScale.get()
        if int(value) != value:
            self.freqScale.set(round(value))
        
    def _round_decscale(self, e=None):
        value = self.decScale.get()
        if int(value) != value:
            self.decScale.set(round(value))

    def plot_elevation(self, elCurveDict=None):
        """Create a plot showing the elevation of the source as seen from 
        the currently selected telescopes."""
        
        # Clear the axis. Return if nothing to plot
        self.elAx.clear()
        self.elAx.set_title("Elevation of Source")
        if elCurveDict is None or len(elCurveDict)==0:
            self.elFigCanvas.draw()
            return
        
        # Re-plot each of the elevation curves
        colLst=["r", "b", "g", "m", "c", "y", "k"]
        for i, label in enumerate(elCurveDict.keys()):
            haArr_hr, elArr_deg = elCurveDict[label]
            self.elAx.plot(haArr_hr, elArr_deg, color=colLst[i%len(colLst)],
                           label=label)
            
        # Format labels and legend
        self.elAx.set_xlim(-12.0, 12.0)
        self.elAx.set_ylim(0.0, 90.0)
        self.elAx.set_xlabel("Hour Angle (hours)")
        self.elAx.set_ylabel("Elevation (degrees)")
        self.elAx.margins(0.02)
        leg = self.elAx.legend(shadow=False)
        for t in leg.get_texts():
            t.set_fontsize('small')

        # Show the plot
        self.elFigCanvas.draw()
        
    def plot_scales(self, lBaseArr_lam=None):
        """Create a plot showing the distribution of baseline lengths and
        equivalent scales."""
        
        # Clear the axis. Return if nothing to plot
        self.blAx.clear()
        self.blAxTop.clear()
        if lBaseArr_lam is None or len(lBaseArr_lam)==0:
            self.blFigCanvas.draw()
            return
        
        # Re-plot a histogram of baseline lengths
        xLst = []
        yLst = []
        for l in lBaseArr_lam:
            xLst += [l/1e3, l/1e3, l/1e3]
            yLst += [-0.1, 1, -0.1]
        self.blAx.plot(xLst, yLst, lw=1.0, alpha=0.4)
        self.blAx.set_ylim(0, 1.4)
        self.blAx.set_xlim(0.0, np.max(lBaseArr_lam)*1.05e-3)
        self.blAxTop.set_xlim(self.blAx.get_xlim()) 

        # Angle 2 str
        def ang2str(ang_deg):
            try:
                ang_deg = float(ang_deg)
                ang_arcsec = ang_deg*3600.0
                if ang_arcsec<60.0:
                    text = u'{:.1f}"'.format(ang_arcsec)
                elif ang_arcsec>=60.0 and ang_arcsec<3600.0:
                    text = u"{:.1f}'".format(ang_deg*60.0)
                else:
                    text = u"{:.1f}\u00B0".format(ang_deg)
                return text
            except Exception:
                return "*"
        
        # Conversion between uvdist & scale
        def tick_function(X):
            np.seterr(divide='ignore', invalid='ignore')
            altUnit = np.degrees(1.22  / (X * 1e3))*3600.0
            return [ang2str(z) for z in altUnit]
        
        # Format labels and legend
        tickLocs = self.blAx.xaxis.get_majorticklocs()
        self.blAxTop.set_xticks(tickLocs)
        self.blAxTop.set_xlim(self.blAx.get_xlim()) 
        self.blAxTop.set_xticklabels(tick_function(tickLocs))        
        self.blAx.set_xlabel("Baseline Length (kilo-lambda)")
        self.blAxTop.set_xlabel("Angular Scale")
        self.blAx.set_ylabel("Arbritrary Units")
        #self.blAx.margins(0.02)
        #self.blAxTop.margins(0.02)

        # Show the plot
        self.blFigCanvas.draw()

        
#-----------------------------------------------------------------------------#
class StatusFrame(ttk.Frame):
    """Frame presenting status indicators lights and action buttons for the
    steps in the observation process."""
    
    def __init__(self, parent, bgColour=None, boxWidth=20, gapWidth=30,
                 yPad=25):
        ttk.Frame.__init__(self, parent)
        self.parent = parent
        if bgColour is None:
            bgColour = ttk.Style().lookup("TFrame", "background")

        # Properties of the status boxes
        self.tagLst = ['statusSelection', 'statusModel', 'statusuvCalc',
                       'statusModelFFT', 'statusuvGrid', 'statusBeam',
                       'statusObsDone']
        self.labLst = ['Array', 'Model', 'uv-Coverage', 'Model FFT', 'uv-Grid',
                       'Beam', 'Observation']
        self.state = [0] * len(self.tagLst)
        
        # Calculate canvas dimensions and internal grid coordinates
        self.boxWidth = boxWidth
        self.gapWidth = gapWidth
        self.yPad = yPad
        self.nBox = len(self.tagLst)
        self.width = self.nBox * boxWidth + self.nBox * gapWidth
        self.height = boxWidth + yPad * 4.8
        self.Y1 = yPad * 0.75
        self.Y2 = boxWidth * 2.0
        self.Y3 = self.Y2 + boxWidth / 2.0 + yPad
        self.Y4 = self.Y3 + yPad * 1.5
        self.X = []
        self.dX = gapWidth + boxWidth
        for i in range(self.nBox):
            self.X.append(gapWidth/2.0 + boxWidth/2.0+ i*(boxWidth + gapWidth))
            
        # Insert the canvas
        self.canvas = tk.Canvas(self, width=self.width, height=self.height,
                                background=bgColour, highlightthickness=0)
        self.canvas.grid(row=0, column=0, padx=0, pady=0)

        # Draw the boxes and labels in order
        for x, tag, text, in zip(self.X, self.tagLst, self.labLst):
            self.canvas.create_text(x, self.Y1, text=text)
            self._draw_checkbox(x, self.Y2, self.boxWidth, tag)

        # Draw the inputs bracket and text
        self.canvas.create_line(self.X[0] - gapWidth / 2.5, self.Y2,
                                self.X[0] - gapWidth / 2.5, self.Y3,
                                self.X[1] + gapWidth / 2.5, self.Y3,
                                self.X[1] + gapWidth / 2.5, self.Y2,
                                width=2, fill="black", joinstyle=tk.MITER)
        self.canvas.create_line(self.X[0] + self.dX/2.0, self.Y3,
                                self.X[0] + self.dX/2.0, self.Y4-16,
                                width=2, fill="black",
                                joinstyle=tk.MITER)
        self.canvas.create_text(self.X[0] + self.dX/2.0, self.Y4,
                                text="INPUTS")
        
        # Draw plot uv-coverage button & line
        self.canvas.create_line(self.X[2], self.Y4,
                                self.X[2], self.Y3 + 3.0,
                                width=2, fill="grey", joinstyle=tk.MITER)
        self.canvas.create_line(self.X[2], self.Y3 - 3.0,
                                self.X[2], self.Y2 + boxWidth / 2.0,
                                width=2, fill="grey", joinstyle=tk.MITER,
                                arrow=tk.LAST)
        self.pltuvCovBtn = ttk.Button(self, text="Calculate", width=16,
                    command=lambda: self.event_generate("<<plot_uvcoverage>>"))
        self.pltuvCovBtn.configure(state="disabled")
        self.canvas.create_window(self.X[2], self.Y4, window=self.pltuvCovBtn)

        # Draw plot model FFT button & line
        self.canvas.create_line(self.X[3], self.Y4,
                                self.X[3], self.Y3 + 3.0,
                                width=2, fill="grey", joinstyle=tk.MITER)
        self.canvas.create_line(self.X[3], self.Y3 - 3.0,
                                self.X[3], self.Y2 + boxWidth / 2.0,
                                width=2, fill="grey", joinstyle=tk.MITER,
                                arrow=tk.LAST)
        self.pltModFFTbtn = ttk.Button(self, text="Calculate", width=16,
                        command=lambda: self.event_generate("<<plot_modFFT>>"))
        self.pltModFFTbtn.configure(state="disabled")
        self.canvas.create_window(self.X[3], self.Y4, window=self.pltModFFTbtn)

        # Draw the observe button & bracket
        self.canvas.create_line(self.X[2] - gapWidth / 2.5, self.Y2,
                                self.X[2] - gapWidth / 2.5, self.Y3,
                                self.X[-1] + gapWidth / 2.5, self.Y3,
                                self.X[-1] + gapWidth / 2.5, self.Y2,
                                width=2, fill="black", joinstyle=tk.MITER)
        self.canvas.create_line(self.X[6], self.Y3,
                                self.X[6], self.Y4,
                                width=2, fill="black",
                                joinstyle=tk.MITER)
        self.obsBtn = ttk.Button(self, text = "Do Observation", width=16,
                    command=lambda: self.event_generate("<<do_observation>>"))
        self.obsBtn.configure(state="disabled")
        obsBtnW = self.canvas.create_window(self.X[6], self.Y4,
                                            window=self.obsBtn)

        # Draw the show results button
#        self.showBtn = ttk.Button(self, text = "Show Plot Window", width=16,
#                    command=lambda: self.event_generate("<<show_results>>"))
#        showBtnW = self.canvas.create_window(self.X[6], self.Y4,
#                                             window=self.showBtn)
        
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
        self._draw_checkbox(self.X[indx], self.Y2, self.boxWidth,
                            tag, state=self.state[indx])

        # Enable or dissable buttons based on state
        if self.state[0]:
            self.pltuvCovBtn.configure(state="enabled")
        else:
            self.pltuvCovBtn.configure(state="disabled")
        if self.state[1]:
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
class InformationPanelGraphic(ttk.Frame):
    """Frame presenting the information on the selected arrays, input model
    and results"""
    
    def __init__(self, parent, width=500, rowHeight=20, axPad=(70,25),
                 tickLen=10, nYticks=3, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        bgColour = ttk.Style().lookup("TFrame", "background")
        bgColour = "lightgreen"
        
        # Canvas & graphic layout parameters
        self.axPad = axPad
        self.width = width
        self.height = rowHeight * 5
        self.canMin = self.axPad[0]
        self.canMax = self.width - self.axPad[1]

        # Plot uv-limits
        self.plotMin = None
        self.plotMax = None
        self.uvCovMin = None
        self.uvCovMax = None        
        self.imgCovMin = None
        self.imgCovMax = None
        
        # Insert the canvas
        self.canvas = tk.Canvas(self, background=bgColour, width=self.width,
                                height=self.height, highlightthickness=0)
        self.canvas.grid(row=0, column=0, padx=0, pady=0)

        # Layout grid coordinates
        self.rowHeight = rowHeight
        self.Y = np.arange(1, 5) * self.rowHeight

        # Main labels
        self.canvas.create_text(self.canMin-10, self.Y[1], anchor="e",
                                text="Array:")
        self.canvas.create_text(self.canMin-10, self.Y[2], anchor="e",
                                text="Model:")
        # Axis line 
        self.canvas.create_line(self.canMin, self.Y[1] + self.rowHeight/2.0,
                                self.canMax, self.Y[1] + self.rowHeight/2.0,
                                width=1, fill="black",
                                joinstyle=tk.MITER)
        
    def _draw_block(self, uvMin, uvMax, y1, y2, fill="blue", tag=""):
        x1 =  self._world2canvas(uvMin)
        x2 =  self._world2canvas(uvMax)
        item = self.canvas.create_rectangle(x1, y1, x2, y2, outline="black",
                                            fill=fill, width=1.0,
                                            tag=("bar", tag))
        
    def _world2canvas(self, x):
        """Convert an array of world coordinates to canvas coordinates"""

        canRng = float(self.canMax - self.canMin)
        plotRng = float(self.plotMax - self.plotMin)
        l = (x-self.plotMin)*canRng/plotRng + self.canMin
        
        return l
    
    def _set_scale(self):
        if self.uvCovMin is None and self.imgCovMin is None:
            self.plotMin = None
        else:
            self.plotMin = min(i for i in [self.uvCovMin, self.imgCovMin]
                               if i is not None)
        if self.uvCovMax is None and self.imgCovMax is None:
            self.plotMax = None
        else:
            self.plotMax = max(i for i in [self.uvCovMax, self.imgCovMax]
                               if i is not None)
        #print 
        #print "ARRAY RNG", self.uvCovMin, self.uvCovMax
        #print "IMAGE RNG",  self.imgCovMin, self.imgCovMax
        #print "PLOT RNG", self.plotMin, self.plotMax 
        #print 
        
    def update(self, parmDict):
        """Update the graphic using a dictionary."""
        
        p = parmDict

        # uv-coverage limits
        if not (p["uvRngMin_lam"] is None or p["uvRngMax_lam"] is None):
            self.uvCovMin = p["uvRngMin_lam"]
            self.uvCovMax = p["uvRngMax_lam"]
        else:
            self.uvCovMin = None
            self.uvCovMax = None
        self._set_scale()            

        # Model coverage limits
        if not (p["pixScaleFFTX_lam"] is None or p["pixScaleFFTY_lam"] is None
                or p["fftScale_lam"] is None):     
            self.imgCovMin = min(i for i in [p["pixScaleFFTX_lam"],
                                             p["pixScaleFFTY_lam"]]
                                 if i is not None)
            self.imgCovMax = p["fftScale_lam"]
        else:
            self.imgCovMin = None
            self.imgCovMax = None            
        self._set_scale()

        # Re-plot the bars
        itemLst = self.canvas.find_withtag("bar")
        for item in itemLst:
            self.canvas.delete(item)
            
        if self.uvCovMin and self.uvCovMax:
            self._draw_block(self.uvCovMin, self.uvCovMax,
                             self.Y[1]-self.rowHeight/2,
                             self.Y[1]+self.rowHeight/2)
            
        if self.imgCovMin and self.imgCovMax:
            self._draw_block(self.imgCovMin, self.imgCovMax,
                             self.Y[2]-self.rowHeight/2,
                             self.Y[2]+self.rowHeight/2, fill="red")

            
#-----------------------------------------------------------------------------#
class PlotFrame(ttk.Frame):
    """Frame showing the plots produced by the virtual interferometer."""

    def __init__(self, parent, bgColour=None, *args, **kwargs):
        ttk.Frame.__init__(self, parent, *args, **kwargs)
        self.parent = parent
        if bgColour is None:
            bgColour = ttk.Style().lookup("TFrame", "background")

        # Dictionary tracking the 6 plot axes and their state (active/clear)
        #              {"AxisName": [axis, location, state],}
        self.axDict  = {"modelImg": [None, "231", 0],
                        "modelFFT": [None, "234", 0],
                        "uvCov":    [None, "232", 0],
                        "beam":     [None, "233", 0],
                        "obsFFT":   [None, "235", 0],
                        "obsImg":   [None, "236", 0]}
        
        # Create the blank figure canvas and grid its tk canvas
        self.fig = Figure(figsize=(13.0, 8.0), facecolor=bgColour)
        self.figCanvas = FigureCanvasTkAgg(self.fig, master=self)
        self.canvas = self.figCanvas.get_tk_widget()
        self.canvas.configure(highlightthickness=0)
        self.canvas.configure(background=bgColour)
        self.canvas.grid(column=0, row=0, columnspan=2, padx=0, pady=0,
                         sticky="NSEW")
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        
        # Add the axes for the plots - default plot is empty
        self.axDict["modelImg"][0] = \
                            self.plot_image("modelImg", title="Model Image")
        self.axDict["modelFFT"][0] = \
                            self.plot_fft("modelFFT", title="Model FFT")
        self.axDict["uvCov"][0] = \
                            self.plot_uvcov("uvCov", title="uv-Coverage")
        self.axDict["beam"][0] = \
                            self.plot_image("beam", title="Synthesised Beam")
        self.axDict["obsFFT"][0] = \
                            self.plot_fft("obsFFT", title="Observed FFT")
        self.axDict["obsImg"][0] = \
                            self.plot_image("obsImg", title="Observed Image")
        
        # Add the information panel
#        self.infoPanel = InformationPanel(self)
#        self.infoPanel.grid(column=0, row=1, columnspan=1, rowspan=2,
#                            padx=5, pady=5, sticky="W")
        
        # Add the matplotlib toolbar
        self.tbarFrm = ttk.Frame(self)
        self.toolbar = MPLnavToolbar(self.figCanvas, self.tbarFrm)
        self.tbarFrm.grid(column=1, row=1, padx=5, pady=5, sticky="NE")

        # Add the show control window Button
#        self.showBtn = ttk.Button(self, text = "Show Control Window", width=20,
#                                  command=self._show_control_window)
#         self.showBtn.grid(column=1, row=2, padx=5, pady=5, sticky="SE")

        # Show the plot
        self.show()
        
    def _show_control_window(self):
        """Set focus back to the main control window."""
        
        root.focus_force()
        root.lift()

    def plot_image(self, axName, imgArr=None, title="", pRng=None):
        """Plot an image with a scalebar (TBD)."""
        
        ax = self.axDict[axName][0]
        loc =  self.axDict[axName][1]
        if ax:
            ax.clear()
            self.axDict[axName][2] = 0
        else:
            ax = self.fig.add_subplot(loc)
        plt.setp(ax.get_yticklabels(), visible=False)
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.set_title(title)
        
        # Catch blank arrays
        if imgArr is None:
            return ax
        if imgArr.max()==imgArr.min():
            return ax

        # Set the colour clip to fractions of range, if requested
        zMin = None
        zMax = None
        if pRng is not None:
            zMin = np.nanmin(imgArr)
            zMax = np.nanmax(imgArr)
            zRng = zMin - zMax
            zMin -= zRng * pRng[0]
            zMax += zRng * pRng[1]
        if zMax==zMin:
            zMin = None
            zMax = None
        
        # Show the image array
        ax.imshow(imgArr, cmap=plt.cm.cubehelix, interpolation="nearest",
                  origin="lower", vmin=zMin, vmax=zMax)
        ax.set_aspect('equal')
        self.axDict[axName][2] = 1
            
        return ax
    
    def plot_fft(self, axName, fftArr=None, limit=None , title=""):
        """Plot a 2D Fourier transform image."""
        
        ax = self.axDict[axName][0]
        loc =  self.axDict[axName][1]
        if ax:
            ax.clear()
            self.axDict[axName][2] = 0
        else:
            ax = self.fig.add_subplot(loc)
        plt.setp(ax.get_yticklabels(), visible=False)
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.set_title(title)
        
        # Catch blank arrays
        if fftArr is None:
            return ax
        if fftArr.max()==fftArr.min():
            return ax
        
        # Show the image array
        if limit:
            extent=[-limit, limit, -limit, limit]
        else:
            extent=None
        ax.imshow(np.abs(fftArr), norm=LogNorm(), cmap=plt.cm.cubehelix,
                  interpolation="nearest", origin="lower", extent=extent)
        ax.set_xlabel(u"u (k$\lambda$)")
        ax.set_ylabel(u"v (k$\lambda$)")
        ax.set_aspect('equal')
        plt.setp(ax.get_yticklabels(), visible=True)
        plt.setp(ax.get_xticklabels(), visible=True)
        self.axDict[axName][2] = 1
            
        return ax

    def plot_uvcov(self, axName, arrsSelected=None, title=""):
        """Plot the uv-coverage of the selected arrays."""
        
        colLst=["r", "b", "g", "m", "c", "y", "k"]
        ax = self.axDict[axName][0]
        loc =  self.axDict[axName][1]
        if ax:
            ax.clear()
            self.axDict[axName][2] = 0
        else:
            ax = self.fig.add_subplot(loc)
        plt.setp(ax.get_yticklabels(), visible=False)
        plt.setp(ax.get_xticklabels(), visible=False)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.set_title(title)

        # Set the z-order of the plots based on the uv-coverage scale
        if arrsSelected is not None:
            oLst = []
            sLst = []
            for i, e in enumerate(arrsSelected):
                oLst.append(i)
                sLst.append(e["scaleMin_deg"])
            multiLst = list(zip(sLst, oLst))
            multiLst.sort()
            sLst, oLst = zip(*multiLst)
            zLst = range(len(oLst))
            multiLst = list(zip(oLst, zLst))
            multiLst.sort()
            oLst, zLst = zip(*multiLst)
            
        if arrsSelected is not None:
            for i, e in enumerate(arrsSelected):
                u = e["uArr_lam"]
                v = e["vArr_lam"]
                ax.scatter(x=u/1000, y=v/1000, marker=".", edgecolor='none',
                           s=2, color=colLst[i%len(colLst)], zorder=zLst[i])
                ax.scatter(x=-u/1000, y=-v/1000, marker=".", edgecolor='none',
                           s=2, color=colLst[i%len(colLst)], zorder=zLst[i])
            ax.set_xlabel(u"u (k$\lambda$)")
            ax.set_ylabel(u"v (k$\lambda$)")
            ax.set_aspect('equal', 'datalim')
            ax.margins(0.02)
            plt.setp(ax.get_yticklabels(), visible=True)
            plt.setp(ax.get_xticklabels(), visible=True)
            self.axDict[axName][2] = 1
            
        return ax

    def show(self):
        """Show the plots in the plotting window."""
        
        self.fig.subplots_adjust(left=0.07, right=0.97, top=0.95, bottom=0.07,
                                 wspace=0.27, hspace=0.24)
        self.toolbar.update()
        self.figCanvas.draw()
        

    def clear_by_state(self, stateDict):
        """Use a dictionary to clear downstream plots based on state."""

        # Track which plots need to be cleared so we don't make multiple calls
        self.clearDict = {"modelImg": 0, "modelFFT": 0, "uvCov": 0,
                          "beam": 0, "obsFFT": 0, "obsImg": 0}

        # Clear the downstream plots base on process status and plot status
        if not stateDict["statusObsDone"]:
            self.clearDict["obsImg"] =  self.axDict["obsImg"][2]
        if not stateDict["statusuvGrid"]:
            self.clearDict["obsFFT"] =  self.axDict["obsFFT"][2]
            self.clearDict["obsImg"] =  self.axDict["obsImg"][2]
            self.clearDict["beam"] =  self.axDict["beam"][2]
        if not stateDict.get("statusModelFFT", 1):
            self.clearDict["modelFFT"] =  self.axDict["modelFFT"][2]
            self.clearDict["obsFFT"] =  self.axDict["obsFFT"][2]
            self.clearDict["obsImg"] =  self.axDict["obsImg"][2]
            self.clearDict["beam"] =  self.axDict["beam"][2]
        if not stateDict.get("statusSelection", 1):
            self.clearDict["uvCov"] =  self.axDict["uvCov"][2]
            self.clearDict["obsFFT"] =  self.axDict["obsFFT"][2]
            self.clearDict["obsImg"] =  self.axDict["obsImg"][2]
            self.clearDict["beam"] =  self.axDict["beam"][2]
        if not stateDict.get("statusuvCalc", 1):
            self.clearDict["uvCov"] =  self.axDict["uvCov"][2]
            self.clearDict["obsFFT"] =  self.axDict["obsFFT"][2]
            self.clearDict["beam"] =  self.axDict["beam"][2]
            self.clearDict["obsImg"] =  self.axDict["obsImg"][2]
        if not stateDict.get("statusBeam", 1):
            self.clearDict["beam"] =  self.axDict["beam"][2]
            self.clearDict["obsImg"] =  self.axDict["obsImg"][2]
        if not stateDict.get("statusModel", 1):
            self.clearDict["modelImg"] =  self.axDict["modelImg"][2]
            self.clearDict["modelFFT"] =  self.axDict["modelFFT"][2]
            self.clearDict["obsImg"] =  self.axDict["obsImg"][2]
            self.clearDict["beam"] =  self.axDict["beam"][2]

        # Clear the relevant plots and show 
        if self.clearDict["modelImg"]:
            self.plot_image("modelImg")
        if self.clearDict["modelFFT"]:
            self.plot_fft("modelFFT")
        if self.clearDict["uvCov"]:
            self.plot_uvcov("uvCov")
        if self.clearDict["beam"]:
            self.plot_image("beam")
        if self.clearDict["obsFFT"]:
            self.plot_fft("obsFFT")
        if self.clearDict["obsImg"]:
            self.plot_image("obsImg")
        self.show()

#-----------------------------------------------------------------------------#
#class MPLnavToolbar(NavigationToolbar2TkAgg):
class MPLnavToolbar(NavigationToolbar2Tk):
    """Subclass the MPL navigation toolbar to disable the coord readout."""
    
    def set_message(self, msg):
        pass
    
    
#-----------------------------------------------------------------------------#
if __name__ == "__main__":
    
    root = tk.Tk()
    
    # Force platform specific colours and fonts
    if sys.platform=="darwin":
        bgColour = "#ececec"
        fontSize = 12
    else:
        bgColour = ttk.Style().lookup("TFrame", "background")
        fontSize = 12       
    ttk.Style().configure("TFrame", background=bgColour)
    ttk.Style().configure("TLabelframe", background=bgColour)
    ttk.Style().configure("TLabel", background=bgColour)
    fontDefault = tkFont.nametofont("TkDefaultFont")
    fontDefault.configure(size=fontSize)
    fontFixed = tkFont.nametofont("TkFixedFont")
    fontFixed.configure(size=fontSize)
    root.option_add("*Font", fontDefault)
    
    # Hack to hide dot files in the Linux tk file dialog
    # https://mail.python.org/pipermail/tkinter-discuss/2015-August/003762.html
    try:
        # call a dummy dialog with an impossible option to initialize the file
        # dialog without really getting a dialog window; this will throw a
        # TclError, so we need a try...except :
        try:
            root.tk.call('tk_getOpenFile', '-foobarbaz')
        except:
            pass
        # now set the magic variables accordingly
        root.tk.call('set', '::tk::dialog::file::showHiddenBtn', '1')
        root.tk.call('set', '::tk::dialog::file::showHiddenVar', '0')
    except:
        pass
    
    # Attempt to compensate for high-DPI displays (not working)
    #root.tk.call('tk', 'scaling', 4.0)
    #root.tk.call('tk', 'scaling', '-displayof', '.', 50)
    
    # Grid the main window and start mainloop
    app = App(root, bgColour).pack(side="top", fill="both", expand=True)
    root.geometry("1428x768")    
    root.update()
    root.minsize(root.winfo_width(), root.winfo_height())

    root.mainloop()
