#!/usr/bin/env python
#=============================================================================#
#                                                                             #
# NAME:     util_plop.py                                                      #
#                                                                             #
# PURPOSE:  Matplotlib plotting functions for the VRI application.            #
#                                                                             #
# MODIFIED: 03-Mar-2017 by C. Purcell                                         #
#                                                                             #
# CONTENTS:                                                                   #
#                                                                             #
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

import numpy as np
#import StringIO
from io import StringIO
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Polygon
from matplotlib.ticker import FuncFormatter
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg


#-----------------------------------------------------------------------------#
def plot_uvcov_ax(ax, arrsSelected):
    
    # Plot the uv-coverage
    ax.cla()
    for e in arrsSelected:
        u = e["uArr_lam"]
        v = e["vArr_lam"]
        ax.scatter(x=u/1000, y=v/1000, marker=".", color="r",
                   edgecolor='none', s=2)
        ax.scatter(x=-u/1000, y=-v/1000, marker=".", color="b",
                   edgecolor='none', s=2)
    ax.set_xlabel(u"u (k$\lambda$)")
    ax.set_ylabel(u"v (k$\lambda$)")
    ax.set_aspect('equal', 'datalim')
    ax.margins(0.02)
    

#-----------------------------------------------------------------------------#
def plot_image_ax(ax, imgArr):
    
    ax.cla()
    ax.imshow(np.abs(imgArr), cmap=plt.cm.cubehelix,
              interpolation="nearest", origin="lower")
    #ax.set_title("Model Image")
    ax.set_xlabel(u"X (pixels)")
    ax.set_ylabel(u"Y (pixels)")
    ax.set_aspect('equal', 'datalim')
    
#-----------------------------------------------------------------------------#
def plot_fft_ax(ax, imgArr, extent=[]):
    
    ax.cla()
    ax.imshow(np.abs(imgArr), norm=LogNorm(), cmap=plt.cm.cubehelix,
              interpolation="nearest", origin="lower")
    #ax.set_title("FFT Image")
    ax.set_xlabel(u"u (kilo-lambda)")
    ax.set_ylabel(u"v (kilo-lambda)")
    ax.set_aspect('equal', 'datalim')
    
