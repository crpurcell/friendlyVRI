#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_scan.py                                                      #
#                                                                             #
# PURPOSE:  Routines for processing a model  array captured via webcam.       #
#                                                                             #
# REQUIRED: Requires numpy, scipy, matplotlib and PIL/PILLOW                  #
#                                                                             #
# CREDITS:  Cormac R. Purcell (cormac.purcell at mq.edu.au)                   #
#           Roy Truelove (Macquarie University)                               #
#                                                                             #
# MODIFIED: 29-Jun-2017 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
#                                                                             #
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

import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator
from PIL import Image


def scan_to_coords(imgName, scale_m=1000.0, eSize=13, threshold_sigma=5.0,
                   crop=1024, ax=None):
    
    # Open the image, convert to luminance greyscale and then a numpy array
    imgPIL = Image.open(imgName).convert("L")
    imgArr = 256 - np.flipud(np.asarray(imgPIL))

    # Crop the image
    Ny, Nx = imgArr.shape
    crop1 = min([crop, Nx, Ny])
    dx = max(0, Nx - crop1)
    dy = max(0, Ny - crop1)
    imgArr = imgArr[dy/2:Ny-dy/2, dx/2:Nx-dx/2]
    
    # Remove large-scale background using morphological opening
    foot = generate_footprint(eSize)
    imgErode = ndimage.grey_erosion(imgArr, size=(eSize, eSize), footprint=foot)
    imgOpen = ndimage.grey_dilation(imgErode, size=(eSize, eSize),
                                    footprint=foot)
    imgBgArr = imgArr - imgOpen
    
    # Determine the finding threshold
    zMax = np.nanmax(imgBgArr)
    zMin = np.nanmin(imgBgArr)
    rms = np.std(imgBgArr)
    zMed = np.median(imgBgArr)
    threshold = zMed + rms * threshold_sigma
    
    # Set everything below the threshold to zero:
    imgMskArr = np.copy(imgBgArr)
    imgMskArr[imgBgArr<threshold] = 0

    # Find the objects and extract subimages
    imgLabeled, Nobjects = ndimage.label(imgMskArr)
    islands = ndimage.find_objects(imgLabeled)

    # Find the centroids of the islands
    pixScale_m = scale_m/float(crop1)
    E_m = []
    N_m = []
    for island in islands:
        dy, dx  = island
        x, y = dx.start, dy.start
        cx, cy = centroid(imgMskArr[island])
        xCent_pix = x + cx
        yCent_pix = y + cy
        xCent_m = (xCent_pix - float(crop1)/2.0) * pixScale_m
        yCent_m = (yCent_pix - float(crop1)/2.0) * pixScale_m
        E_m.append(xCent_m)
        N_m.append(yCent_m)

    # Plot the positions of the antennae
    if not ax==None:
        
        ax.imshow(imgBgArr, interpolation="nearest", cmap="gray_r",
                  origin='lower', extent=[-scale_m/2000.0, scale_m/2000.0,
                                          -scale_m/2000.0, scale_m/2000.0])

        # Annotate the detected antennae
        for x, y in zip(E_m, N_m):
            ax.plot(x/1000.0, y/1000.0, 'kx', ms=10, color="green",
                    linewidth=5)
            print "%.1f, %.1f" % (x, y)
        
        # Format labels and legend
        ax.set_aspect('equal')
        ax.set_xlim(-scale_m/2000.0, scale_m/2000.0)
        ax.set_ylim(-scale_m/2000.0, scale_m/2000.0)
        ax.xaxis.set_major_locator(MaxNLocator(4))
        ax.yaxis.set_major_locator(MaxNLocator(4))
        ax.margins(0.02)
        ax.set_xlabel("East Offset (km)")
        ax.set_ylabel("North Offset (km)")

        
    return np.array(E_m), np.array(N_m), None

        
#-----------------------------------------------------------------------------#
def centroid(data):
    h,w = np.shape(data)   
    x = np.arange(0,w)
    y = np.arange(0,h)

    X,Y = np.meshgrid(x,y)

    cx = np.sum(X*data)/np.sum(data)
    cy = np.sum(Y*data)/np.sum(data)

    return cx, cy


#-----------------------------------------------------------------------------#
def generate_footprint(size):
    f = np.ones((size, size))
    yi, xi =np.indices(f.shape, dtype='float')
    d = np.sqrt((yi-size/2.0+0.5)**2.0 + (xi-size/2.0+0.5)**2.0)
    b = np.where(d<(size/2.0), 1, 0)
    return b


#-----------------------------------------------------------------------------#
def write_arrayfile(fileName, E_m, N_m, telescope="MyTelescope_1",
                    config="Config_1", latitude_deg=-20.0, diameter_m=22.0):

    FH = open(fileName, "w")
    FH.write("telescope = %s\n" % telescope)
    FH.write("config =  %s\n" % config)
    FH.write("latitude_deg = %f\n" % latitude_deg)
    FH.write("diameter_m = %f\n" % diameter_m)

    
    for i in range(len(E_m)):
        FH.write("%f, %f\n" % (E_m[i], N_m[i]))
    FH.close()

    
