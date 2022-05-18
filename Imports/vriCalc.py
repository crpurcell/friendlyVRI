#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     vriCalc.py                                                        #
#                                                                             #
# PURPOSE:  Back-end for a virtual interferometer application.                #
#                                                                             #
# REQUIRED: Requires numpy and pillow                                         #
#                                                                             #
# CREDITS:  Cormac R. Purcell (cormac.purcell at mq.edu.au)                   #
#           Roy Truelove      (Macquarie University)                          #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  observationManager (class)                                                 #
#      _reset_uv_vars       ... reset variables associated with uv-coverage   #
#      _reset_model_vars    ... reset variables associated with model image   #
#      _load_one_array      ... load a single array configuration             #
#      _load_all_arrays     ... create the table of available array configs   #
#      get_available_arrays ... list the available arrays                     #
#      select_array         ... select an array config & HA-range             #
#      get_selected_arrays  ... list the selected cofigurations & HA-ranges   #
#      clear_all_selections ... clear the current selections                  #
#      set_obs_parms        ... set the common parameters (freq, dec)         #
#      get_obs_parms        ... get a dict of common parameters               #
#      calc_uvcoverage      ... calculate the uv-coverage for selected arrays #
#      load_model_image     ... load a model image                            #
#      set_pixscale         ... set a new pixel scale for the model           #
#      invert_model         ... calculate the FFT of the model image          #
#      grid_uvcoverage      ... grid the uv-coverage onto the image grid      #
#      calc_beam            ... calculate the beam image                      #
#      invert_observation   ... apply the uv-coverage to the model image      #
#      get_status           ... set the status flags of each step             #
#      get_array_params     ... get the parameters of an available array      #
#      get_baseline_lengths ... get the baseline lengths for an array         #
#      calc_elevation_curve ... calculate elevation curves                    #
#                                                                             #
#  antArray (class)                                                           #
#      _init_variables      ... define variables describing telescope & array #
#      _load_arrayfile      ... load an array definition file into the class  #
#      _calc_baselines      ... calculate the baseline vectors                #
#                                                                             #
#  od2list                  ... convert an ordered dictionary to a list       #
#  sort_nicely              ... sort words as a human would (1 before 10)     #
#  ang2str                  ... convert a degree angle to a string+units      #
#                                                                             #
#=============================================================================#
#                                                                             #
# The MIT License (MIT)                                                       #
#                                                                             #
# Copyright (c) 2017 - 2022 Cormac R. Purcell and Roy Truelove                #
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

import re
import glob
import copy
import traceback
from collections import OrderedDict as od
import numpy as np
from PIL import Image

# Speed of light
C = 2.99792458e8

#-----------------------------------------------------------------------------#
class observationManager:
    """Class used to simulate an observation of an astronomical source with
    a combination of radio-interferometer array configurations.

    The user selects one or more array configurations from the available list
    and chooses a range of observed hour-angles for each. The sampling cadence
    and observing frequency are common to all array configurations.
    
    The observationManager class keeps track of the observing parameters and
    calculates the uv-coverage, which is then applied to a model image to
    simulate the observation.

    The current version assumes a noiseless image, does not include sensitivity
    calculations and uses natural weighting only."""
    
    def __init__(self, arrayDir="arrays", verbose=True, debug=True):
        """Constructor for the observationManager class.
 
        Initialise the array configurations and default observing parameters.
        The array configurations are defined in ASCII files in the arrayDir
        subdirectory [default = 'arrays/']."""
        
        # Ordered dictionary of available arrays as antArray objects.
        # Unique key = "telescope_config" with entries:
        #  "row"       ... <row_number>
        #  "telescope" ... <telescope_name>
        #  "config"    ... <config_name>
        #  "antArray"  ... <antennaArray_object>
        self.arrsAvailable = od()
        
        # Telescope parameter lookup tables
        self.telescopeLatDict = {}
        self.telescopeDiamDict = {}
        
        # Table of selected arrays. Each entry is a dictionary that will
        # contain the following keys:
        #  "key"          ... unique key composed of <telescope>_<config>
        #  "telescope"    ... name of telescope
        #  "config"       ... name of array configuration
        #  "row"          ... row number in the arrsAvailable table
        #  "haStart"      ... starting Hour angle in hours
        #  "haEnd"        ... ending Hour angle in hours
        #  "sampRate"     ... sampling cadence in seconds
        #  "haArr_rad"    ... vector of hour angles
        #  "uArr_lam"     ... array of u values
        #  "vArr_lam"     ... array of v values
        #  "scaleMin_deg" ... minimum uv-coverage scale
        #  "scaleMax_deg" ... maximum uv-coverage scale
        #  "priBeam_deg"  ... primary beam FWHM at current frequency
        self.arrsSelected = []

        # Observing parameters common to all array configurations.
        self.freq_Hz = 1420e6
        self.lambda_m = C/self.freq_Hz
        self.dec_deg = 20.0

        # Ensemble uv-coverage parameters
        self.scaleMin_deg = None
        self.scaleMax_deg = None
        self.uvRngMin_lam = None
        self.uvRngMax_lam = None
        self.priBeamMax_deg = None
        
        # Model image and parameters
        self.modelImgArr = None
        self.nX = None
        self.nY = None
        self.pixScaleImg_asec = None

        # Model FFT image and parameters
        self.modelFFTarr = None
        self.fftScale_lam = None
        self.pixScaleFFTX_lam = None
        self.pixScaleFFTY_lam = None
        
        # uv-mask and synthesised beam
        self.uvMaskArr = None
        self.uvCntArr = None
        self.beamArr = None
        
        # Observed FFT and final image
        self.obsFFTarr = None
        self.obsImgArr = None

        # Flags
        self.statusSelection = False
        self.statusuvCalc = False
        self.statusModel = False
        self.statusModelFFT = False
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False
        self.verbose = verbose
        self.debug = debug
        
        # Read the array definition files from the array directory
        self._load_all_arrays(arrayDir)
    
    def _reset_uv_vars(self):
        """Reset the variables associated with uv-coverage."""
        
        self.uvMaskArr = None
        self.beamArr = None
        self.scaleMin_deg = None
        self.scaleMax_deg = None
        self.priBeamMax_deg = None
        self.obsFFTarr = None
        self.obsImgArr = None

    def _reset_model_vars(self):
        """Reset the variables associated with the model."""
        
        self.pixScaleImg_asec = None
        self.modelFFTarr = None
        self.fftScale_lam = None
        self.pixScaleFFTX_lam = None
        self.pixScaleFFTY_lam = None
        self.uvMaskArr = None
        self.beamArr = None
        self.obsFFTarr = None
        self.obsImgArr = None
        
    def _load_one_array(self, arrayFile):
        """Load a single ASCII array file into memory."""
        
        # Create an antArray object for each and store in a dictionary
        ar = antArray(arrayFile)
        if ar.telescope is not None:
            key = ar.telescope + "_" + ar.config
            value = {"row": len(self.arrsAvailable),
                     "telescope": ar.telescope,
                     "config": ar.config,
                     "antArray": ar}
            self.arrsAvailable[key] = value

            # Populate the telescope lookup tables
            self.telescopeLatDict[ar.telescope] = ar.latitude_deg
            self.telescopeDiamDict[ar.telescope] = ar.diameter_m
        
    def _load_all_arrays(self, arrayDir, pattern="*.config"):
        """Read and parse each of the ASCII files defining the antenna
        coordinates and telescope parameters (telescope, configuration, 
        latitude, antenna diameter). By default all files with extension
        '.config' in the arrayDir directory are read. The information on 
        each array configuration is stored in an ordered dictionary of 
        antArray objects."""
        
        # Read all files ending with '.config' in the array directory
        arrayFileLst = glob.glob(arrayDir + "/" + pattern)
        sort_nicely(arrayFileLst)
        for i in range(len(arrayFileLst)):
            arrayFile = arrayFileLst[i]
            
            # Create an antArray object for each and store in a dictionary
            self._load_one_array(arrayFileLst[i])
            
        if self.verbose:
            print("Successfully loaded %d array configurations." 
                  % len(self.arrsAvailable))

    def get_available_arrays(self):
        """Return a list of available array configurations from the
        'arrsAvailable' dict. Returns a numpy structured array with column
        names ['key', 'telescope', 'config', 'minBaseline', 'maxBaseline'].
        """
        
        if len(self.arrsAvailable)>0:
            
            # Create a numpy structured array
            dtype=[('key', 'a100'), ('telescope', 'a50'), ('config', 'a50'),
                   ('minBaseline', 'f2'), ('maxBaseline', 'f2')]
            a = np.zeros(len(self.arrsAvailable), dtype=dtype)

            # Fill with a summary of available array configurations and return
            for i, k in enumerate(self.arrsAvailable.keys()):
                v = self.arrsAvailable[k]
                a[i] = (k, v['telescope'], v['config'],
                        np.min(v['antArray'].lBase_m),
                        np.max(v['antArray'].lBase_m))
            return a
        
    def select_array(self, key, haStart=-6.0, haEnd=6.0, sampRate_s=300.0,
                     byRow=False):
        """Select an array configuration, hour-angle range and sampling rate
        to include in the observation. The key must be the name of one of the
        available array configurations (format = <telescope>_<array>). 
        Available keys can be queried by invoking the get_available_arrays()
        method."""

        # Clear downstream uv-coverage variables
        self._reset_uv_vars()
        
        # Set the downstream status flags
        self.statusSelection = False
        self.statusuvCalc = False
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False

        # Limit the values of haStart and haEnd
        try:
            haStart = max(-12.0, float(haStart))
            haEnd = min(12.0, float(haEnd))
        except Exception:
            if self.verbose:
                print("Setting the hour-angle range failed!")
            if self.debug:
                print(traceback.format_exc())
                return
            
        # Append selection to the observation table
        if key in self.arrsAvailable:
            try: 
                self.arrsSelected.append(
                    {"key": key,
                     "telescope": self.arrsAvailable[key]["telescope"],
                     "config": self.arrsAvailable[key]["config"],
                     "row": int(self.arrsAvailable[key]["row"]),
                     "haStart": haStart,
                     "haEnd": haEnd,
                     "sampRate": sampRate_s,
                     "haArr_rad": None,
                     "uArr_lam": None,
                     "vArr_lam": None,
                     "scaleMin_deg": None,
                     "scaleMax_deg": None,
                     "priBeam_deg": None}
                )
                if self.verbose:
                    print("Array '{:s}' HA = {:.2f} to {:.2f} by {:.0f}s".\
                          format(key, float(haStart), float(haEnd),
                                 sampRate_s))
            except Exception:
                if self.verbose:
                    print("Selection failed!")
                if self.debug:
                    print(traceback.format_exc())
                return
        else:
            if self.verbose:
                print("Warning: array key '%s' not found!" % key)
                print("\nAvailable array configurations are:")
                print(self.arrsAvailable.keys())

        # Update the selection flag
        if len(self.arrsSelected)>0:
            self.statusSelection = True
                
    def get_selected_arrays(self):
        """Return a list of array configurations selected for the current
        observation (from the 'arrsSelected' list). Returns a numpy structured
        array with column names ['#', 'key', 'telescope', 'config', 'haStart',
        'haEnd']."""

        if len(self.arrsSelected)>0:
            
            # Create a numpy structured array
            dtype=[('#', 'i4'), ('key', 'a100'), ('telescope', 'a50'),
                   ('config', 'a50'), ('haStart', 'f4'), ('haEnd', 'f4'),
                   ('sampRate', 'f4')]
            a = np.zeros(len(self.arrsSelected), dtype=dtype)
            
            # Fill with a summary of the arrsSelected table
            for i, e in enumerate(self.arrsSelected):
                a[i] = (i, e["key"], e["telescope"], e["config"],
                        e["haStart"], e["haEnd"], e["sampRate"])
            return a
        
    def clear_all_selections(self):
        """Clear all array configuration selections from the observation."""

        # Clear all uv-coverage variables
        self.arrsSelected = []
        self._reset_uv_vars()

        # Set the downstream status flags
        self.statusSelection = False
        self.statusuvCalc = False
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False
            
    def set_obs_parms(self, freq_MHz=1420.0, dec_deg=20.0):
        """Set the common observation parameters:
            Observing frequency (MHz) = 1420.0
            Source declination (deg)  =   20.0"""
        
        self.freq_Hz = float(freq_MHz)*1e6
        self.lambda_m = C/self.freq_Hz
        self.dec_deg = float(dec_deg)

        # Set the downstream status flags
        self.statusuvCalc = False
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False

    def get_obs_parms(self):
        """Return the common observation parameters as a dictionary."""

        d = {"freq_Hz": self.freq_Hz,
             "dec_deg": self.dec_deg}
        
        return d
    
    def calc_uvcoverage(self, redo=False):
        """Calculate the uv-coverage for each of the selected array configs."""

        # Check for selected array configurations
        if len(self.arrsSelected)==0:
            if self.verbose:
                print("No array configurations selected.")
                return

        # Reset the downstream variables
        self._reset_uv_vars()
        
        # Set the downstream status flags
        self.statusuvCalc = False
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False

        # Unit conversions of common parameters
        self.lambda_m = C/self.freq_Hz
        dec_rad = np.radians(self.dec_deg)
        
        # Loop through each selected array configuration
        scaleMinLst_deg = []
        scaleMaxLst_deg = []
        priBeamLst_deg = []
        for e in self.arrsSelected:
            if self.verbose:
                print("\nCalculating uv-coverage for %s" % \
                      [e["telescope"], e["config"], e["haStart"], e["haEnd"],
                       e["sampRate"]])
            try:
                
                # Calculate the hour-angle sampling vector
                sampRate_deg = e["sampRate"] * 15.0 / 3600.0
                sampRate_hr = e["sampRate"] / 3600.0
                nSamps = int((e["haEnd"] - e["haStart"])/sampRate_hr +1)
                haArr_hr = np.linspace(e["haStart"], e["haEnd"], nSamps)
                haArr_rad = np.radians(haArr_hr * 15.0)
                e["haArr_rad"] = haArr_rad

                # Calculate the elevation curve and mask HA-range
                dummy, elArr_deg = \
                        self.calc_elevation_curve(e["telescope"], haArr_hr)
                haArr_rad[elArr_deg<=0] = np.nan
                    
                # Fill the uv-plane with samples over the hour-angle range
                ar = od2list(self.arrsAvailable)[e["row"]]["antArray"]
                latitude_rad = np.radians(ar.latitude_deg)
                u_m = np.zeros((ar.nBase, nSamps))
                v_m = np.zeros((ar.nBase, nSamps))
                for i in range(ar.nBase):
                    u_m[i, :] = (ar.Bx_m[i] * np.sin(haArr_rad) + 
                                 ar.By_m[i] * np.cos(haArr_rad))
                    v_m[i, :] = (-ar.Bx_m[i] * np.sin(dec_rad) *
                                 np.cos(haArr_rad) +
                                 ar.By_m[i] * np.sin(dec_rad) *
                                 np.sin(haArr_rad) +
                                 ar.Bz_m[i] * np.cos(dec_rad))
                e["uArr_lam"] = u_m/self.lambda_m
                e["vArr_lam"] = v_m/self.lambda_m
                
                # Calculate the max & min scales from the uv-coverage
                if np.all(haArr_rad!=haArr_rad):
                    e["scaleMin_deg"] = np.nan
                    e["scaleMax_deg"] = np.nan
                    e["priBeam_deg"] = np.nan
                else:
                    lArr_lam = np.sqrt(e["uArr_lam"]**2.0 + e["vArr_lam"]**2.0)
                    e["scaleMin_deg"] = np.degrees(1.0/np.nanmax(lArr_lam))
                    e["scaleMax_deg"] = np.degrees(1.0/np.nanmin(lArr_lam))
                    e["priBeam_deg"] = \
                                np.degrees(1.22*self.lambda_m/ar.diameter_m)
                scaleMinLst_deg.append(e["scaleMin_deg"])
                scaleMaxLst_deg.append(e["scaleMax_deg"])
                priBeamLst_deg.append(e["priBeam_deg"])
                if self.verbose:
                    print("Scale Range = %s to %s" %
                          (ang2str(e["scaleMin_deg"]),
                                   ang2str(e["scaleMax_deg"])))
                    print("Primary Beam FWHM = %s" %
                          (ang2str(e["priBeam_deg"])))
                    
            except Exception:        
                if self.verbose:
                    print("uv-coverage calculation failed for:")
                    print(e["telescope"], e["config"])
                if self.debug:
                    print(traceback.format_exc())
                return

        # Remember the range of scales and largest primary beam
        self.scaleMin_deg = np.min(scaleMinLst_deg)
        self.scaleMax_deg = np.max(scaleMaxLst_deg)
        self.uvRngMin_lam = 1.0/np.radians(np.max(scaleMaxLst_deg))
        self.uvRngMax_lam = 1.0/np.radians(np.min(scaleMinLst_deg))
        self.priBeamMax_deg = np.max(priBeamLst_deg)
        if self.verbose:
            print ("\nuv-Coverage Parameters:")
            print(u"Minimum uv-Spacing = %.3f k\u03bb" %
                  (self.uvRngMin_lam/1e3))
            print(u"Maximum uv-Spacing = %.3f k\u03bb" %
                  (self.uvRngMax_lam/1e3))
            print("Minimum scale = %s" % (ang2str(self.scaleMin_deg)))
            print("Maximum scale = %s" % (ang2str(self.scaleMax_deg)))
            print("Field of View = %s" % (ang2str(self.priBeamMax_deg)))

        # Update the status flag for the uvcoverage calculation
        self.statusuvCalc = True
            
    def load_model_image(self, modelFile, pixScaleImg_asec=0.5):
        """Load a model from a standard image file.  Valid file formats are
        raster images supported by the Python Imaging Library (PIL),
        e.g.: PNG, JPEG, TIFF, GIF and BMP. Images can be rectangular, but are
        assumed to have square pixels."""

        # Reset the downstream model and the current model variables
        self._reset_model_vars()
        self.modelImgArr = None
        self.nX = None
        self.nY = None
        
        # Set the downstream status flags
        self.statusModel = False
        self.statusModelFFT = False
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False
        
        # Open the image, convert to luminance and then a numpy array
        try:
            imgPIL = Image.open(modelFile).convert("L")
            self.modelImgArr = np.flipud(np.asarray(imgPIL))
            self.pixScaleImg_asec = pixScaleImg_asec
            self.nY, self.nX = self.modelImgArr.shape
        except Exception:
            if self.verbose:
                print("Failed to load model image.")
            if self.debug:
                print(traceback.format_exc())
            return
        
        # Print the model image parameters
        if self.verbose:
            print ("\nModel Image Parameters:")
            print("Pixel scale = %s" % (ang2str(self.pixScaleImg_asec/3600.0)))
            print("Image size = %d x %d pixels [%s x %s]" % (self.nX, self.nY,
                                ang2str(self.nX*self.pixScaleImg_asec/3600.0),
                                ang2str(self.nY*self.pixScaleImg_asec/3600.0)))
        
        # Set the status of the model image & FFT flags to True
        self.statusModel = True

    def set_pixscale(self, pixScaleImg_asec=0.5):
        """Set a new value for the size of the pixels in the model image."""

        # Reset the downstream model variables
        self._reset_model_vars()
        
        # Set the downstream status flags
        self.statusModelFFT = False
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False

        # Set the new pixel scale
        self.pixScaleImg_asec = pixScaleImg_asec
        
        # Print the model image parameters
        if self.verbose:
            print ("\nModel Image Parameters:")
            print("Pixel scale = %s" % (ang2str(self.pixScaleImg_asec/3600.0)))
            print("Image size = %d x %d pixels [%s x %s]" % (self.nX, self.nY,
                                ang2str(self.nX*self.pixScaleImg_asec/3600.0),
                                ang2str(self.nY*self.pixScaleImg_asec/3600.0)))

    def invert_model(self):
        """Calculate the 2D Fast Fourier Transform of the model image."""
        
        # First check that a model has been loaded
        if not self.statusModel:
            if self.verbose:
                print("A model image has not been loaded.")
            return
        
        # Set the downstream status flags
        self.statusModelFFT = False
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False
        
        # Calculate the 2D FFT and scaling factors.
        # The shape of FFT array is same as the model image.
        try:
            self.modelFFTarr = np.fft.fft2(self.modelImgArr)
            self.modelFFTarr = np.fft.fftshift(self.modelFFTarr)
            pixScaleImg_lam = np.radians(self.pixScaleImg_asec/3600.0)
            self.fftScale_lam = 1.0/pixScaleImg_lam
            self.pixScaleFFTX_lam = 2.0*self.fftScale_lam/self.nX
            self.pixScaleFFTY_lam = 2.0*self.fftScale_lam/self.nY            
        except Exception:
            if self.verbose:
                print("Failed to calculate the FFT of the model image.")
            if self.debug:
                print(traceback.format_exc())
            return        
        
        # Print the model FFT parameters
        if self.verbose:
            print ("\nModel FFT Parameters:")
            print(u"Pixel scale = %.3f x %.3f k\u03bb" % \
                  (self.pixScaleFFTX_lam/1000.0, self.pixScaleFFTY_lam/1000.0))
            print(u"Image limits = -%.3f to +%.3f k\u03bb" % \
                  (self.fftScale_lam/1000.0, self.fftScale_lam/1000.0))
        
        # Set the status of the model image & FFT flags to True
        self.statusModelFFT = True
        
    def grid_uvcoverage(self):
        """Grid the uv-coverage to use as a mask for the model FFT image."""

        # First check uv-coverage and model are available
        if not self.statusModelFFT or not self.statusuvCalc:
            if self.verbose:
                print("Model FFT or uv-Coverage unavailable!")
            return
        
        # Set the status flags for next calculations
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False
        
        # Grid the uv-coverage
        try:
            self.uvMaskArr = np.zeros(self.modelFFTarr.shape, dtype=np.int32)
            self.uvCntArr = np.zeros(self.modelFFTarr.shape, dtype=np.int32)
            for i, e in enumerate(self.arrsSelected):
                u_lam = e["uArr_lam"].flatten()
                v_lam = e["vArr_lam"].flatten()
                u_pix = (u_lam+self.fftScale_lam)/self.pixScaleFFTX_lam
                v_pix = (v_lam+self.fftScale_lam)/self.pixScaleFFTY_lam
                u2_pix = (-u_lam+self.fftScale_lam)/self.pixScaleFFTX_lam
                v2_pix = (-v_lam+self.fftScale_lam)/self.pixScaleFFTY_lam
                for j in range(len(u_pix)):
                    try:
                        self.uvMaskArr[int(v_pix[j]), int(u_pix[j])] = 1
                        self.uvMaskArr[int(v2_pix[j]), int(u2_pix[j])] = 1
                        self.uvCntArr[int(v_pix[j]), int(u_pix[j])] += 1
                        self.uvCntArr[int(v2_pix[j]), int(u2_pix[j])] += 1
                    except Exception:
                        # Ignore if visibility falls outside of the FFT image
                        pass
        except Exception:
                if self.verbose:
                    print("Gridding failed!")
                if self.debug:
                    print(traceback.format_exc())
                return

        # Apply the gridded uv-coverage to the model FFT
        try:
            self.obsFFTarr = self.modelFFTarr.copy()*self.uvMaskArr
        except Exception:
                if self.verbose:
                    print("Masking failed!")
                if self.debug:
                    print(traceback.format_exc())
                return
                    
        # Print the percentage coverage
        if self.verbose:
            nPix = self.uvMaskArr.shape[0] * self.uvMaskArr.shape[1]
            nUsedPix = np.sum(self.uvMaskArr)
            pc = nUsedPix*100.0/nPix
            print("{:.2f} % of pixels used in observed FFT image.".\
                  format(pc))                

        # Set the status of the uv-grid flag
        self.statusuvGrid = True
        
    def calc_beam(self):
        """Calculate the beam image for the gridded uv-coverage."""

        # Check the gridding status
        if not self.statusuvGrid:
            if self.verbose:
                print("The uv-coverage has not been gridded!")
            return

        # Set the downstream flag
        self.statusBeam = False
        
        # Calculate the beam image
        try:
            self.beamArr = np.fft.ifft2(self.uvMaskArr)
            self.beamArr = np.fft.ifftshift(self.beamArr)
        except Exception:
            if self.verbose:
                print("Failed to calculate the beam image!")
            if self.debug:
                print(traceback.format_exc())
            return
        
        # Update the status
        self.statusBeam = True
        
    def invert_observation(self):
        """Apply the gridded uv-coverage to the model FFT image, and
        invert to produce the final observed image."""
        
        # First check that tge gridding has been done
        if not self.statusuvGrid:
            if self.verbose:
                print("The uv-coverage has not been gridded!")
            return
        
        # Set the downstream flag
        self.statusObsDone = False
        
        try:
            
            # Invert to produce the final image
            self.obsImgArr = np.fft.ifft2(np.fft.ifftshift(self.obsFFTarr))
        except Exception:
            if self.verbose:
                print("Failed produce the observed image!")
            if self.debug:
                print(traceback.format_exc())
            return

        # Update the status 
        self.statusObsDone = True
        if self.verbose:
            print("Observation complete!")

    def get_status(self):
        """Return the values of the status flags as a dictionary."""

        d = {"statusSelection": self.statusSelection,
             "statusuvCalc": self.statusuvCalc,
             "statusModel": self.statusModel,
             "statusModelFFT": self.statusModelFFT,
             "statusuvGrid": self.statusuvGrid,
             "statusBeam": self.statusBeam,
             "statusObsDone": self.statusObsDone}
        return d
            
    def get_scales(self):
        """Return the values of the status flags as a dictionary."""
        
        d = {"scaleMin_deg": self.scaleMin_deg,
             "scaleMax_deg": self.scaleMax_deg,
             "uvRngMin_lam": self.uvRngMin_lam,
             "uvRngMax_lam": self.uvRngMax_lam,
             "priBeamMax_deg": self.priBeamMax_deg,
             "pixScaleFFTX_lam": self.pixScaleFFTX_lam,
             "pixScaleFFTY_lam": self.pixScaleFFTY_lam,
             "fftScale_lam": self.fftScale_lam}
        return d
             
    def get_array_params(self, row=None, key=None):
        """Return the basic parameters of an antenna array from one entry
        in the 'arrsAvailable' dictionary. The query can take either a row
        number or the key=['<telescope>_<configuration>']."""

        d = {"diameter_m": None,
             "latitude_deg": None,
             "baseMax_m": None,
             "baseMin_m": None,
             "x": None,
             "y": None}
        
        if row is not None:
            arrsAvailLst = od2list(self.arrsAvailable)
            if row>=len(arrsAvailLst):
                return d
            d["diameter_m"] = arrsAvailLst[row]["antArray"].diameter_m
            d["latitude_deg"] = arrsAvailLst[row]["antArray"].latitude_deg
            d["baseMin_m"] = np.nanmin(arrsAvailLst[row]["antArray"].lBase_m)
            d["baseMax_m"] = np.nanmax(arrsAvailLst[row]["antArray"].lBase_m)
            d["x"] = arrsAvailLst[row]["antArray"].eastArr_m.copy()
            d["y"] = arrsAvailLst[row]["antArray"].northArr_m.copy()
        if key is not None:
            if not key in self.arrsAvailable:
                return d
            antArrayRow = self.arrsAvailable[key]["antArray"]
            d["diameter_m"] = antArrayRow.diameter_m
            d["latitude_deg"] = antArrayRow.latitude_deg
            d["baseMin_m"] = np.nanmin(antArrayRow.lBase_m)
            d["baseMax_m"] = np.nanmax(antArrayRow.lBase_m)
            d["x"] = antArrayRow.eastArr_m.copy()
            d["y"] = antArrayRow.northArr_m.copy()
        return d
    
    def get_baseline_lengths(self, row=None, key=None, doKiloLam=True):
        """Return the baseline lengths for an antenna array from one entry
        in the 'arrsAvailable' dictionary. The query can take either a row
        number or the key=['<telescope>_<configuration>']."""

        if row is not None:
            arrsAvailLst = od2list(self.arrsAvailable)
            if row>=len(arrsAvailLst):
                return None
            lBase_m =  arrsAvailLst[row]["antArray"].lBase_m
        if key is not None:
            key = key.decode("utf-8")
            if not key in self.arrsAvailable:
                return None
            antArrayRow = self.arrsAvailable[key]["antArray"]
            lBase_m = antArrayRow.lBase_m
        return lBase_m

    def calc_elevation_curve(self, telescope, haArr_hr=None):
        """Calculate the elevation curves for a telescope and the current
        source declination over a vector of hour-angles."""
        
        if haArr_hr is None:
            haArr_hr = np.linspace(-12.0, +12.0, 600)
        haArr_rad = np.radians(haArr_hr * 15.0)
        dec_rad = np.radians(self.dec_deg)

        telescope = np.bytes_(telescope).decode("utf-8")  # Python 3 fix
        if telescope in self.telescopeLatDict.keys():
            latitude_rad = np.radians(self.telescopeLatDict[telescope])
            elArr_rad  = (np.sin(latitude_rad) * np.sin(dec_rad) +
                np.cos(latitude_rad) * np.cos(dec_rad) * np.cos(haArr_rad))
            elArr_deg = np.degrees(elArr_rad)
            elArr_deg[elArr_deg<0] = 0.0
        else:
            return None, None
            
        return haArr_hr, elArr_deg


#-----------------------------------------------------------------------------#
class antArray:
    """Class defining the parameters of an interferometer array configuration.
    The antenna coordinates and other parameters are read from a formatted
    ASCII file on initialisation. In the file, general parameters should be
    encoded as 'key = value' pairs. Antenna positions should be specified as
    two comma-separated columns in 'EAST, NORTH' coordinate format. 
    The '#' character comments out lines."""
    
    def __init__(self, arrayFile):
        """Constructor for the antArray class. Object is always constructed
        but has 'None' values if the array file does not load successfully."""
        
        # Initialise the class variables
        self._init_variables()
        
        # Read and parse the array file
        success = self._load_arrayfile(arrayFile)
        if not success:
            self._init_variables()
            return
        
        # Calculate the baseline vector components
        success = self._calc_baselines()
        if not success:
            self._init_variables()
            return

    def _init_variables(self):
        """Initialise the class variables"""
        
        self.telescope = None
        self.config = None
        self.latitude_deg = None
        self.diameter_m = None
        self.nAnt = None
        self.nBase = None
        self.eastArr_m = None
        self.northArr_m = None
        self.xArr_m = None
        self.yArr_m = None
        self.zArr_m = None
        self.Bx_m = None
        self.By_m = None
        self.Bz_m = None
        self.lBase_m = None
        
    def _load_arrayfile(self, arrayFile):
        """Read and parse the ASCII file defining the array parameters and 
        antenna coordinates."""
        
        # Temporary variables
        keyValDict = {}
        eastLst = []
        northLst = []
        
        # Compile necessary regular expressions
        spaces = re.compile('\s+')
        commaAndSpaces = re.compile(',\s+')
        comment = re.compile('#.*')
        quotes = re.compile('\'[^\']*\'')
        keyVal = re.compile('^.+=.+')
        twoCSVcols = re.compile('^.+,.+')

        # Open the file and scan, line by line
        try:
            FH = open(arrayFile, "r")
        except Exception:
            if self.verbose:
                print("Failed to open file '%s'." % arrayFile)
            if self.debug:
                print(traceback.format_exc())
            return 0
        for line in FH:
            line = line.rstrip("\n\r")
            if not comment.match(line):
                line = comment.sub('', line)           # internal comments 
                line = line.replace("'", '')           # remove quotes
                line = commaAndSpaces.sub(',', line)   # kill ambiguous spaces
                
                # Capture key=value pairs
                if keyVal.match(line):
                    keyword, value = line.split('=',1)
                    value = value.strip()              # kill ext whitespace 
                    keyword = keyword.strip()       
                    value = spaces.sub('', value)      # shrink int whitespace
                    keyword = spaces.sub('', keyword)    
                    if value:
                        keyValDict[keyword] = value

                # Capture antenna coordinate entries
                if twoCSVcols.match(line):
                    east, north = line.split(',')
                    eastLst.append(east)
                    northLst.append(north)
        FH.close()

        # Set the class variables, using defaults if necessary
        try:
            self.telescope = keyValDict.get("telescope", "UNKNOWN")
            self.config = keyValDict.get("config", "UNKNOWN")
            self.latitude_deg = float(keyValDict.get("latitude_deg", 20.0))
            self.latitude_rad = np.radians(self.latitude_deg)
            self.diameter_m = float(keyValDict.get("diameter_m", 22.0))
            self.eastArr_m = np.array(eastLst, dtype="f4")
            self.northArr_m = np.array(northLst, dtype="f4")
            self.nAnt = int(len(self.eastArr_m))
            self.nBase = int(self.nAnt*(self.nAnt-1)/2)
        except Exception:
            if self.verbose:
                print("Failed to set telescope parameter!")
            if self.debug:
                print(traceback.format_exc())
            return 0

        # Catch array definition files with < 2 antennas
        if len(self.eastArr_m)<2:
            return 0
        
        # Calculate the antenna coordinates in Earth-centred coordinate frame
        # Technically, we should have terms for the distance from the
        # centre of the Earth, but if the elevation is the same for all
        # antennas, these cancel out when calculating the baseline vectors.
        try:
            self.xArr_m = -self.northArr_m*np.sin(self.latitude_rad)
            self.yArr_m = self.eastArr_m
            self.zArr_m = self.northArr_m*np.cos(self.latitude_rad)
        except Exception:
            if self.debug:
                print(traceback.format_exc())
            return 0

        return 1
        
    def _calc_baselines(self):
        """Calculate the baselines vectors in metres from the antenna
        coordinates in the Earth-centred system."""

        try:
            self.Bx_m = np.zeros((self.nBase))
            self.By_m = np.zeros((self.nBase))
            self.Bz_m = np.zeros((self.nBase))

            # Loop through the unique antenna pairs
            n = 0
            for i in range(self.nAnt-1):
                for j in range(i+1, self.nAnt):
                    self.Bx_m[n] = self.xArr_m[j] - self.xArr_m[i]
                    self.By_m[n] = self.yArr_m[j] - self.yArr_m[i]
                    self.Bz_m[n] = self.zArr_m[j] - self.zArr_m[i]
                    n += 1

            # Calculate vector of baseline lengths
            self.lBase_m = np.sqrt(self.Bx_m**2.0 + self.By_m**2.0
                                   + self.Bz_m**2.0)
        except Exception:
            if self.debug:
                print(traceback.format_exc())
            return 0
        
        return 1

    
#-----------------------------------------------------------------------------#
def od2list(od):
    """Convert an ordered dictionary to a list."""
        
    return list(od.values())


#-----------------------------------------------------------------------------#
def sort_nicely( l ):
    """Sort in human order. Code from StackOverflow:
    http://stackoverflow.com/questions/4836710/
    does-python-have-a-built-in-function-for-string-natural-sort"""
    
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 
    l.sort( key=alphanum_key ) 

    
#-----------------------------------------------------------------------------#
def ang2str(angle_deg):
    """Convert an angle in degrees to a unicode string with appropriate units.
    """
    try:
        angle_deg = float(angle_deg)
        angle_arcsec = angle_deg*3600.0
        if angle_arcsec<60.0:
            text = u'{:.2f}"'.format(angle_arcsec)
        elif angle_arcsec>=60.0 and angle_arcsec<3600.0:
            text = u"{:.2f}'".format(angle_deg*60.0)
        else:
            text = u"{:.2f}\u00B0".format(angle_deg)
        return text
    except Exception:
        return ""
