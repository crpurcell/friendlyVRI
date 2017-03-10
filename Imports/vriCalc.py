#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     vriCalc.py                                                        #
#                                                                             #
# PURPOSE:  Back-end for a virtual interferometer application.                #
#                                                                             #
# MODIFIED: 10-Mar-2017 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  observationManager (class)                                                 #
#      _load_all_arrays     ... create the table of available array configs   #
#      get_available_arrays ... list the available arrays                     #
#      select_array         ... select an array config & HA-range             #
#      get_selected_arrays  ... list the selected cofigurations & HA-ranges   #
#      clear_all_selections ... clear the current selections                  #
#      set_obs_parms        ... set the common parameters (freq, samp, dec)   #
#      get_obs_parms        ... get a dict of common parameters               #
#      od2list              ... convert an orderd dictionary to a list        #
#      calc_uvcoverage      ... calculate the uv-coverage for selected arrays #
#      load_model_image     ... load a model image                            #
#      invert_model         ... calculate the FFT of the model image          #
#      grid_uvcoverage      ... grid the uv-coverage onto the image grid      #
#      calc_beam            ... calculate the beam image                      #
#      invert_observation   ... apply the uv-coverage to the model image      #
#      get_status           ... set the status flags of each step             #
#      get_ant_coordinates  ... get the antenna coords of an array            #
#      calc_elevation_curve ... calculate elevation curves                    #
#                                                                             #
#  antArray (class)                                                           #
#      _init_variables      ... define variables describing telescope & array #
#      _load_arrayfile      ... load an array definition file into the class  #
#      _calc_baselines      ... calculate the baseline vectors                #
#                                                                             #
#  sort_nicely              ... sort words as a human would (1 before 10)     #
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

import re
import glob
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
    
    def __init__(self, arrayDir="arrays", verbose=False, debug=False):
        """Constructor for the observationManager class.
 
        Initialise the array configurations and default observing parameters. 
        The array configurations are defined in ASCII files in the arrayDir
        subdirectory [default = 'arrays/']."""
        
        # Ordered dictionary of available arrays as antArray objects.
        # Unique key = "telescope_config" with entries:
        #  { "row": <rownum>,
        #    "telescope": <telescope_name>,
        #    "config": <config_name>,
        #    "antArray": <antennaArray_object> }
        self.arrsAvailable = od()
        
        # Telescope parameter lookup tables
        self.telescopeLatDict = {}
        self.telescopeDiamDict = {}
        
        # Table of selected arrays. Each entry is a dictionary that will
        # contain the HA array, uv-coverage arrays and scale information.
        self.arrsSelected = []

        # Observing parameters common to all array configurations.
        self.freq_Hz = 1420e6
        self.lambda_m = C/self.freq_Hz
        self.sampRate_s = 60.0
        self.dec_deg = 20.0

        # uv-coverage parameters
        self.scaleMin_deg = None
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
        
        # uv-mask and beam
        self.uvMaskArr = None
        self.beamArr = None
        
        # Observed FFT and final image
        self.obsFFTarr = None
        self.obsImgArr = None

        # Status flags
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
    
    def _load_all_arrays(self, arrayDir, pattern="*.config"):
        """Read and parse each of the ASCII files defining the antenna
        coordinates and telescope parameters (telescope, configuration, 
        latitude, antenna diameter). By default all files with extension
        '.config' in the arrayDir directory are read. The information on 
        each array configuration is stored in an orderd dictionary of 
        antArray objects."""
        
        # Read all files ending with '.config' in the array directory
        arrayFileLst = glob.glob(arrayDir + "/" + pattern)
        sort_nicely(arrayFileLst)
        for i in range(len(arrayFileLst)):
            arrayFile = arrayFileLst[i]
            
            # Create an antArray object for each and store in a dictionary
            ar = antArray(arrayFileLst[i])
            if ar.telescope is not None:
                key = ar.telescope + "_" + ar.config
                value = {"row": i,
                         "telescope": ar.telescope,
                         "config": ar.config,
                         "antArray": ar}
                self.arrsAvailable[key] = value

                # Populate the telescope lookup tables
                self.telescopeLatDict[ar.telescope] = ar.latitude_deg
                self.telescopeDiamDict[ar.telescope] = ar.diameter_m
            
        if self.verbose:
            print("Successfully loaded %d array configurations." 
                  % len(self.arrsAvailable))

    def get_available_arrays(self):
        """Get a list of available array configurations from the 
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
        
    def select_array(self, key, haStart=-6.0, haEnd=6.0, byRow=False):
        """Select an array configuration and hour-angle range to include in
        the observation. The key must be the name of one of the available
        array configurations (format = <telescope>_<array>). Available keys
        can be queried by invoking the get_available_arrays() method."""

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
            
        # Append to the observation table
        if key in self.arrsAvailable:
            try: 
                self.arrsSelected.append(
                    {"key": key,
                     "telescope": self.arrsAvailable[key]["telescope"],
                     "config": self.arrsAvailable[key]["config"],
                     "row": int(self.arrsAvailable[key]["row"]),
                     "haStart": haStart,
                     "haEnd": haEnd,
                     "haArr_rad": None,
                     "uArr_lam": None,
                     "vArr_lam": None,
                     "scaleMin_deg": None,
                     "scaleMax_deg": None,
                     "priBeam_deg": None}
                )
                if self.verbose:
                    print("Selected array '{:s}' HA = {:.2f} to {:.2f}".\
                          format(key, float(haStart), float(haEnd)))
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
        """Get a list of array configurations selected for the current
        observation (from the 'arrsSelected' list). Returns a numpy structured
        array with column names ['#', 'key', 'telescope', 'config', 'haStart',
        'haEnd']."""

        if len(self.arrsSelected)>0:
            
            # Create a numpy structured array
            dtype=[('#', 'i4'), ('key', 'a100'), ('telescope', 'a50'),
                   ('config', 'a50'), ('haStart', 'f4'), ('haEnd', 'f4')]
            a = np.zeros(len(self.arrsSelected), dtype=dtype)
            
            # Fill with a summary of the arrsSelected table
            for i, e in enumerate(self.arrsSelected):
                a[i] = (i, e["key"], e["telescope"], e["config"],
                        e["haStart"], e["haEnd"])
            return a
        
    def clear_all_selections(self):
        """Clear all array configuration selections from the observation."""
        
        self.arrsSelected = []
        
        # Set the downstream status flags
        self.statusSelection = False
        self.statusuvCalc = False
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False
            
    def set_obs_parms(self, freq_MHz=1420.0, sampRate_s=60.0, dec_deg=20.0):
        """Set the common observation parameters: 
            Observing frequency (MHz) = 1420.0  
            Sampling rate (s)         =   60.0
            Source declination (deg)  =   20.0"""
        
        self.freq_Hz = float(freq_MHz)*1e6
        self.sampRate_s = float(sampRate_s)
        self.dec_deg = float(dec_deg)

        # Set the downstream status flags
        self.statusuvCalc = False
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False

    def get_obs_parms(self):
        """Return the common observation parameters as a dictionary."""

        d = {"freq_Hz": self.freq_Hz,
             "sampRate_s": self.sampRate_s,
             "dec_deg": self.dec_deg}
        
        return d

    def od2list(self, od):
        """Convert an ordered dictionary to a list."""
        
        return list(od.values())
    
    def calc_uvcoverage(self, redo=False):
        """Calculate the uv-coverage for each of the selected array configs."""
        
        # Set the downstream status flags
        self.statusuvCalc = False
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False

        # Check for selected array configurations
        if len(self.arrsSelected)==0:
            if self.verbose:
                print("No array configurations selected.")
                return

        # Unit conversions of common parameters
        self.lambda_m = C/self.freq_Hz
        dec_rad = np.radians(self.dec_deg)
        sampRate_deg = self.sampRate_s * 15.0 / 3600.0
        sampRate_hr = self.sampRate_s / 3600.0
        
        # Loop through each selected array configuration
        scaleMinLst_deg = []
        priBeamLst_deg = []
        for e in self.arrsSelected:
            if self.verbose:
                print("\nCalculating uv-coverage for %s" % \
                      [e["telescope"], e["config"], e["haStart"], e["haEnd"]])
            try:
                
                # Calculate the hour-angle sampling vector
                nSamps = int((e["haEnd"] - e["haStart"])/sampRate_hr +1)
                haArr_hr = np.linspace(e["haStart"], e["haEnd"], nSamps)
                haArr_rad = np.radians(haArr_hr * 15.0)
                e["haArr_rad"] = haArr_rad
            
                # Fill the uv-plane with samples over the hour-angle range
                ar = self.od2list(self.arrsAvailable)[e["row"]]["antArray"]
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
                lArr_lam = np.sqrt(e["uArr_lam"]**2.0 + e["vArr_lam"]**2.0)
                e["scaleMin_deg"] = np.degrees(1.0/np.nanmax(lArr_lam))
                e["scaleMax_deg"] = np.degrees(1.0/np.nanmin(lArr_lam))
                e["priBeam_deg"] = np.degrees(1.22*self.lambda_m/ar.diameter_m)
                scaleMinLst_deg.append(e["scaleMin_deg"])
                priBeamLst_deg.append(e["priBeam_deg"])
                if self.verbose:
                    print("Scale Range = {:.2f} - {:.2f} arcseconds.".\
                          format(e["scaleMin_deg"]*3600.0,
                                 e["scaleMax_deg"]*3600.0))
                    print("Primary Beam = {:.2f} arcmin.".\
                          format(e["priBeam_deg"]*60))
            except Exception:        
                if self.verbose:
                    print("uv-coverage calculation failed for:")
                    print(e["telescope"], e["config"])
                if self.debug:
                    print(traceback.format_exc())
                return

        # Remember the minumum resolution and largest primary beam
        self.scaleMin_deg = np.min(scaleMinLst_deg)
        self.priBeamMax_deg = np.max(priBeamLst_deg)
        if self.verbose:
            print("\nResolution = {:.2f} arcseconds.".\
                  format(self.scaleMin_deg*3600.0))
            print("Field of View =  {:.2f} arcminutes.".\
                  format(self.priBeamMax_deg*60.0))

        # Update the status flag for the uvcoverage calculation
        self.statusuvCalc = True
            
    def load_model_image(self, modelFile, pixScaleImg_asec=0.5):
        """Load a model image from a standard image file.  Valid file formats
        are raster images supported by the Python Imaging Library (PIL), 
        e.g.: PNG, JPEG, TIFF, GIF and BMP. Images can be rectangular, but are
        assumed to have square pixels."""

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
            print("Pixel scale = %.3f arcsec " % self.pixScaleImg_asec)
            print("Image size = %d x %d pixels [%.1f x %.1f arcsec]" % \
                  (self.nX, self.nY, self.nX*self.pixScaleImg_asec,
                   self.nY*self.pixScaleImg_asec))
        
        # Set the status of the model image & FFT flags to True
        self.statusModel = True
        
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
                print("Failed to take the FFT of the model image.")
            if self.debug:
                print(traceback.format_exc())
            return        
        
        # Print the model FFT parameters
        if self.verbose:
            print ("\nModel FFT Parameters:")
            print("Pixel scale = %.3f x %.3f kilo-lambda" % \
                (self.pixScaleFFTX_lam/1000.0, self.pixScaleFFTY_lam/1000.0))
            print("Image limits: -%.3f to +%.3f kilo-lambda" % \
                (self.fftScale_lam/1000.0, self.fftScale_lam/1000.0))
        
        # Set the status of the model image & FFT flags to True
        self.statusModelFFT = True
        
    def grid_uvcoverage(self):
        """Grid the uv-coverage to use as a mask for the model FFT image."""

        # First check that a model has been loaded
        if not self.statusModelFFT:
            if self.verbose:
                print("The FFT of the model image is not available!")
            return
        
        # Set the status flags for next calculations
        self.statusuvGrid = False
        self.statusBeam = False
        self.statusObsDone = False
        
        # Grid the uv-coverage
        try:
            self.uvMaskArr = np.zeros(self.modelFFTarr.shape, dtype=np.int32)
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
                    except Exception:
                        # Ignore if visibility falls outside of the FFT image
                        pass
        except Exception:
                if self.verbose:
                    print("Gridding failed!")
                if self.debug:
                    print(traceback.format_exc())
                return
                    
        # Print the percentage coverage
        if self.verbose:
            nPix = self.uvMaskArr.shape[0] * self.uvMaskArr.shape[1]
            nUsedPix = np.sum(self.uvMaskArr)
            pc = nUsedPix*100.0/nPix
            print("{:.1f} % of pixels used in observed FFT image.".\
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
        
        # First check for that a model has been loaded
        if not self.statusuvGrid:
            if self.verbose:
                print("The uv-coverage has not been gridded!")
            return
        
        # Set the downstream flag
        self.statusObsDone = False
        
        try:
            # Apply the gridded uv-coverage to the model FFT
            self.obsFFTarr = self.modelFFTarr.copy()*self.uvMaskArr

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
            
    def get_ant_coordinates(self, row=None, key=None):
        """Return the antenna coordinates from one entry in the 
        'arrsAvailable' dictionary. The query can either take row number 
        or the key=['<telescope>_<configuration>']."""

        if row is None and key is None:
            return None, None
        
        if row is not None:
            arrsAvailLst = self.od2list(self.arrsAvailable)
            if row>=len(arrsAvailLst):
                return None, None           
            x = arrsAvailLst[row]["antArray"].eastArr_m.copy()
            y = arrsAvailLst[row]["antArray"].northArr_m.copy()
            return x, y
        if key is not None:
            if not key in self.arrsAvailable:
                return None, None          
            x = self.arrsAvailable[key]["antArray"].eastArr_m.copy()
            y = self.arrsAvailable[key]["antArray"].northArr_m.copy()
            return x, y
        
    def calc_elevation_curve(self, telescope, haArray=None):
        """Calculate the elevation curves for a telescope and the current
        source declination over a vector of hour-angles."""
        
        pass
#        haArr_hr = np.linspace(-6.0, +6.0, 100)
#        haArr_rad = np.radians(haArr_hr * 15.0)
#        dec_rad = np.radians(self.dec_deg)
#        elLst = []        
#        telescopeLst = self.telescopeLatDict.keys()
#        for telescope in telescopeLst:                    
#            latitude_rad = np.radians(self.telescopeLatDict[telescope])
#            el  = (np.sin(latitude_rad) * np.sin(dec_rad) +
#                   np.cos(latitude_rad) * np.cos(dec_rad) * np.cos(haArr_rad))
#            elLst.append(el)
#
#        return haArr_rad, elLst



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

            # Calculate vector of baseline length
            self.lBase_m = np.sqrt(self.Bx_m**2.0 + self.By_m**2.0
                                   + self.Bz_m**2.0)
        except Exception:
            if self.debug:
                print(traceback.format_exc())
            return 0
        
        return 1
        
#-----------------------------------------------------------------------------#
def sort_nicely( l ):
    """Sort in human order"""
    
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 
    l.sort( key=alphanum_key ) 
