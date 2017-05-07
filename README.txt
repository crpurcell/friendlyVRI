#-----------------------------------------------------------------------------#
#                                                                             #
# The friendly Virtual Radio Interferometer tool,                             #
# by Cormac R. Purcell and Roy Truelove (Macquarie University, Sydney).       #
#                                                                             #
# Copyright (c) 2017 Cormac R. Purcell and Roy Truelove.                      #
#                                                                             #
# Released under the MIT licence (see LICENCE.txt).                           #
#                                                                             #
#-----------------------------------------------------------------------------#

ABOUT:
This application is built to help astronomers investigate the effect of
combining different array configurations when observing an astronomical object
using a radio interferometer. The graphical interface is written using the
Python tkinter library for maximum portability.

CONTACT:
Questions or comments should be directed to 'cormac.purcell (at) mq.edu.au'.


#-----------------------------------------------------------------------------#
DEVELOPMENT TO-DO LIST:

* DoubleSlider:
  - Bind jump to nearest click
  - (Bind simultaneous drag)

* Calculations and logic:
 - Deal cleanly with telescopes that can can never see source (below horizon)
 - Create ATCA array configuration files.
 - In vriTk, check for numbers & disallow non-numeric characters.

* Model browser & observation control
 - Convert Information panel to graphic
    
* Nice to have:
  - Array selection table: selection & plot should change by arrow. <Ret>
       should choose>.
  - (shadowing calculation.)
  - (robust weighting calculation)
  - (fitted synthesized beam size)

* Plotting
 - Modal window for other plots.
 - scalebar on images
 - primary beam circles on images.
 - legend for uv-coverage plot
 - zoom & pan bindings
 - options to save publication quality plots
 - options to save FITS files
 - zoom & auto colour-scale for synthesised beam
 - gamma slider for each MPL image figure


#-----------------------------------------------------------------------------#
THE CALCULATION MODULE
The calculations underlying the graphical applicatiion are seperated into the
file 'vriCalc.py' to facilitate use with altrernative interfaces. For example,
observations may be simulated from the basic ipython shell as follows:

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
