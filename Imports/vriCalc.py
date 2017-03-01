#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     vri_calc.py                                                       #
#                                                                             #
# PURPOSE:  Back-end for the virtual interferometer application.              #
#                                                                             #
# MODIFIED: 01-Mar-2017 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#  observationManager                                                         #
#  antArray                                                                   #
#  sort_nicely                                                                #
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
from collections import OrderedDict as od
import numpy as np

# Speed of light
C = 2.99792458e8

#-----------------------------------------------------------------------------#
class observationManager:
    """Base class used to perform calculations for the simple virtual radio-
    interferometer application. The user selects one or more array
    configurations from the available list, chooses an hour-angle range for
    each and a common observing frequency. The class calculates the
    uv-coverage and applies to a model image to derive the dirty map."""
    
    def __init__(self, arrayDir="arrays"):

        # ordered dictionary of available arrays as antArray objects.
        # Unique key = "telescope_config" with entries:
        #  { "row": <rownum>,
        #    "telescope": <telescope_name>,
        #    "config": <config_name>,
        #    "antArray": <antennaArray_object> }
        self.arrsAvailable = od()
        
        # Telescope parameter lookup tables
        self.telescopeLatDict = {}
        self.telescopeDiamDict = {}
        
        # Table of selected arrays, including uv-coverage
        self.arrsSelected = []

        # Observing parameters        
        self.freq_Hz = 1420e6
        self.sampRate_s = 60.0
        self.dec_deg = 20.0
        
        # Model image and parameters
        self.modelImg = None
        self.modelFFT = None
        self.pixScale_asec = None

        # uv-mask and beam
        
        # Read the array definition files from the array directory
        self._load_all_arrays(arrayDir)
    
    def _load_all_arrays(self, arrayDir, pattern="*.config"):
        """Read and parse each of the ASCII files defining the antenna
        coordinates and telescope parameters (latitude, antenna diameter).
        Store the information in an orderd dict of 'antArray' objects."""
        
        # Read all files ending with '.config' in the array directory
        arrayFileLst = glob.glob(arrayDir + "/" + pattern)
        sort_nicely(arrayFileLst)        
        for i in range(len(arrayFileLst)):
            arrayFile = arrayFileLst[i]
            
            # Create an antArray object for each and store in dictionary
            ar = antArray(arrayFileLst[i])
            key = ar.telescope + "_" + ar.config
            value = {"row": i,
                     "telescope": ar.telescope,
                     "config": ar.config,
                     "antArray": ar}
            self.arrsAvailable[key] = value

            # Populate the telescope lookup tables
            self.telescopeLatDict[ar.telescope] = ar.latitude_deg
            self.telescopeDiamDict[ar.telescope] = ar.diameter_m
            
    def select_array(self, key, haStart=-6.0, haEnd=6.0, byRow=False):
        """Select an array configuration to include in the observation. The
        hour-angle range is unique to each selected configuration"""
        
        self.arrsSelected.append(
            {"row": int(self.arrsAvailable[key]["row"]),
             "haStart": float(haStart),
             "haEnd": float(haEnd),
             "haArr_rad": None,
             "uArr_lam": None,
             "vArr_lam": None}
        )
        
    def clear_selection(self):
        """Clear the selections from the arrsSelected table."""
        
        self.arrsSelected = []
            
    def set_obs_parms(self, freq_MHz=1420.0, sampRate_s=60.0, dec_deg=20.0):
        """Set the common observation parameters."""
        
        self.freq_Hz = float(freq_MHz)*1e6
        self.sampRate_s = float(sampRate_s)
        self.dec_deg = float(dec_deg)

    def get_array_by_row(self, row):
        """Access the ordered dict of available configurations by row index."""
        
        return list(self.arrsAvailable.values())[row]
        
    def calc_selected_uvcoverage(self, redo=False):
        """Calculate the uv-coverage for each of the selected array configs."""

        # Unit conversions of common parameters
        lambda_m = C/self.freq_Hz
        dec_rad = np.radians(self.dec_deg)
        sampRate_deg = self.sampRate_s * 15.0 / 3600.0
        sampRate_hr = self.sampRate_s / 3600.0
        
        # Loop through each selected array configuration
        for e in self.arrsSelected:
            
            # Calculate the hour-angle sampling vector
            nSamps = int((e["haEnd"] - e["haStart"])/sampRate_hr +1)
            haArr_hr = np.linspace(e["haStart"], e["haEnd"], nSamps)
            haArr_rad = np.radians(haArr_hr * 15.0)
            e["haArr_rad"] = haArr_rad
            
            # Fill the uv-plane with samples over the hour-angle range
            ar = self.get_array_by_row(e["row"])["antArray"]
            latitude_rad = np.radians(ar.latitude_deg)
            u_m = np.zeros((ar.nBase, nSamps))
            v_m = np.zeros((ar.nBase, nSamps))
            for i in range(ar.nBase):
                u_m[i, :] = (ar.Bx_m[i] * np.sin(haArr_rad) + 
                             ar.By_m[i] * np.cos(haArr_rad))
                v_m[i, :] = (-ar.Bx_m[i] * np.sin(dec_rad) * np.cos(haArr_rad) +
                             ar.By_m[i] * np.sin(dec_rad) * np.sin(haArr_rad) +
                             ar.Bz_m[i] * np.cos(dec_rad))
            e["uArr_lam"] = u_m/lambda_m
            e["vArr_lam"] = v_m/lambda_m

    def grid_uvcoverage(self):
        """Grid the uv-coverage to use as a mask for the Fourier domain image
        of the input model."""
        
        pass

    def get_ant_coordinates(self, row):
        """Return the antenna coordinates of the available array in the
        requested row number."""

        arrsAvailLst = list(self.arrsAvailable.values())
        x = arrsAvailLst[row]["antArray"].eastArr_m.copy()
        y = arrsAvailLst[row]["antArray"].northArr_m.copy()

        return x, y

    def get_elevation(self):
        """Return the elevation curves of each selected telescope for the 
        current source declination."""
        
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
    """Class containing the layout of an interferometer array configuration.
    The antenna coordinates and other parameters are read from a formatted
    ASCII file on initialisation. In the file, general parameters should be
    encoded as 'key = value' pairs. Antenna positions should be specified as
    two comma-separated columns in 'EAST, NORTH' coordinate format. 
    The '#' character comments out lines."""
    
    def __init__(self, arrayFile):

        # Class variables
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
        
        # Read and parse the array file
        self._load_arrayfile(arrayFile)

        # Calculate the baseline vector components
        self._calc_baselines()
        
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
        FH = open(arrayFile, "r")
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
        
        # Set the class variables
        self.telescope = keyValDict.get("telescope", "UNKNOWN")
        self.config = keyValDict.get("config", "UNKNOWN")
        self.latitude_deg = float(keyValDict.get("latitude_deg", 20.0))
        self.latitude_rad = np.radians(self.latitude_deg)
        self.diameter_m = float(keyValDict.get("diameter_m", 22.0))
        self.eastArr_m = np.array(eastLst, dtype="f4")
        self.northArr_m = np.array(northLst, dtype="f4")
        self.nAnt = len(self.eastArr_m)
        self.nBase = self.nAnt*(self.nAnt-1)/2

        # Calculate the antenna coordinates in Earth-centred coordinate frame
        # Technically, we should have terms for the distance from the
        # centre of the Earth, but if the elevation is the same for all
        # antennas, these cancel out when calculating the baseline vectors.
        self.xArr_m = -self.northArr_m*np.sin(self.latitude_rad)
        self.yArr_m = self.eastArr_m
        self.zArr_m = self.northArr_m*np.cos(self.latitude_rad)

    def _calc_baselines(self):
        """Calculate the baselines vectors in metres from the antenna
        coordinates in the Earth-centred system."""

        # Calculate the baselines vectors
        self.Bx_m = np.zeros((self.nBase))
        self.By_m = np.zeros((self.nBase))
        self.Bz_m = np.zeros((self.nBase))
        n = 0
        for i in range(self.nAnt-1):
            for j in range(i+1, self.nAnt):
                self.Bx_m[n] = self.xArr_m[j] - self.xArr_m[i]
                self.By_m[n] = self.yArr_m[j] - self.yArr_m[i]
                self.Bz_m[n] = self.zArr_m[j] - self.zArr_m[i]
                n += 1
        self.lBase_m = np.sqrt(self.Bx_m**2. + self.By_m**2. + self.Bz_m**2.)

        
#-----------------------------------------------------------------------------#
def sort_nicely( l ):
    """Sort in human order"""
    
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)] 
    l.sort( key=alphanum_key ) 
