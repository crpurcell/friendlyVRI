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
Python tkinter library for maximum portability. For a quick-start guide see
the file 'HELP.txt'.

INSTALLATION REQUIREMENTS:
The friendlyVRI tool is written in Python and requires the numpy, matplotlib
tkinter and PIL (or PILLOW) modules. If the python OpenCV module (cv2) is
available a webcam can be used to capture a model image.

CONTACT:
Questions or comments should be directed to 'cormac.purcell (at) mq.edu.au'.


#-----------------------------------------------------------------------------#
THE CALCULATION MODULE
The calculations underlying the graphical application are seperated into the
file 'vriCalc.py' to facilitate use with alternative interfaces. For example,
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

# Set the observing frequency (MHz) and source declination (deg).
obsMan.set_obs_parms(1420.0, 20.0)

# Calculate the uv-coverage
obsMan.calc_uvcoverage()

# Load the model and set the pixel scale in arcsec
obsMan.load_model_image("models/galaxy_lobes.png")
obsMan.set_pixscale(1.0)

# Calculate the FFT of the model image
obsMan.invert_model()

# Grid the uv-coverage onto the same pixels as the FFT as the model image
obsMan.grid_uvcoverage()

# Create the beam image
obsMan.calc_beam()

# Apply the uv-coverage and create observed image
obsMan.invert_observation()


#-----------------------------------------------------------------------------#
DEVELOPMENT TO-DO LIST:

* Calculations and logic:
  - robust weighting calculation
    
* Nice to have:
  - Convert Information panel in plotting window to table.
  - Array selection table: <Ret> should choose (currently mouse only).
  - Fit Gaussian to synthesized beam to evaluate size & PA

* Plotting
 - legend for uv-coverage plot.
 - scalebar on images.
 - primary beam circles on images.
 - options to save publication quality plots (individualy & all 6 panels).
 - gamma slider for each MPL image figure.
 - option to save FITS files.
