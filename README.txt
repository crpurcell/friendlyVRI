#-----------------------------------------------------------------------------#
#
# A friendly Virtual Radio Interferometer tool.
# by Cormac R. Purcell and Roy Truelove (Macquarie University, Sydney)
#
# Copyright (c) 2017 Cormac R. Purcell and Roy Truelove
#
#-----------------------------------------------------------------------------#

This application is built to help astronomers investigate the effect of
combining different array configurations when observing an astronomical object
using an interferometer.

The graphical interface is written using the tkinter library for maximum
portability.

TODO:

* ScatterPlot widget:
DONE  - equalise X & y scales
DONE  - Format of text labels and margins
  - Alternative vertical text for older Tkinter.

* TreeTables:
DONE  - Format of column spacing & text.

* DoubleSlider:
NA  - Bind entry boxes
DONE  - Limit parameter ranges.
DONE  - Change to asymmetric handle shape.
NA  - Snap to grid.
DONE  - Limit significant figures.

* Calculations and logic:
DONE - flags showing & controlling steps done
DONE - catch obvious errors
 - mask uv-coverage for elevation limits.
DONE - status lights
 - Create ATCA array configuration files.
 - Limit max & min values in vriCalc
 - Limt max & min value in vriTk
 - In vriTk, check for numbers & disallow non-numeric characters.

* Model browser & observation control
DONE  - relative path in file entry box ( + make longer )
DONE  - move all controlls to controller window.
IN PROG  - Add an information panel with:
DONE    -- max & min size scale
    -- synthesised beam size
DONE    -- image sampling parameters
NA    -- warning about miss-match between model & observing parameters
    
* Nice to have:
  - graphical comparison of model size scale (structure fn?, dispersion)
    versus uv-coverage (histogram of baselines?).
  - shadowing calculation.
  - robust weighting calculation
  - fitted synthesized beam size
  - Elevation plot for selected telescopes & uv-coverage
  - Initial instructions on the array plot window.
  - Array selection table: selection & plot should change by arrow. <Ret>
       should choose>.

* Plotting
 - scalebar on images rather than axes
 - primary beam circles on images.
 - fix colours on uv-coverage plot
 - legend for uv-coverage plot
 - zoom & pan bindings
 - images resize when window is maximized
 - options to save publication quality plots
 - zoom in on synthesized beam (1/3rd of plot)

#-----------------------------------------------------------------------------#
# Using the observationManager from the ipython shell

# Load the observation manager
from Imports.vriCalc import observationManager
obsMan = observationManager(verbose=True, debug=True)
obsMan.get_available_arrays()

# Select array configurations and hour-angle ranges.
obsMan.select_array('VLA_A')
obsMan.select_array('VLA_B')
obsMan.select_array('VLA_C')
obsMan.select_array('VLA_D')
obsMan.get_selected_arrays()

# Calculate the uv-coverage
obsMan.calc_uvcoverage()

# Load the model
obsMan.load_model_image("models/Lenna.png")

# Grid the uv-coverage onto the model grid
obsMan.grid_uvcoverage()

# Create the beam image
obsMan.calc_beam()

# Apply the uv-coverage and create observed image
obsMan.invert_observation()

