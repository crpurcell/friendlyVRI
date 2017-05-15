
::

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

ABOUT
=======

This application shows the effect of combining different array
configurations when observing an astronomical object using a radio
interferometer.


INSTALLATION
==============

The Friendly VRI tool is written in Python and requires the numpy, matplotlib
tkinter and PIL (or PILLOW) modules. If the python OpenCV module (cv2) is
available a webcam can be used to capture a model image. This feature is
currently disabled on Mac OS as the default OpenCV module is faulty.

See the file `INSTALL <INSTALL.txt>`_ for full instructions.


USAGE INSTRUCTIONS
===================

The interface is designed to be relatively intuitive for non-experts. The
control window allows you to:
* Plot the layout of the antennas in an ``array configuration`` for a particular
  telescope.
* Create a list of observations to be done using different array
  configurations, over different time ranges (hour angle ranges) and with
  different sampling cadence.
* Load in a model image and set its angular scale on the sky.
* Apply the uv-coverage of the observations to the model to simulate an
  observation.
The plotting window shows inputs, outputs an intermediate steps in the process:
model image, fast Fourier transform (FFT) of the model, plot of uv-coverage,
filtered model FFT, synthesised beam and final observed image.

See the file `HELP <HELP.txt>`_ or the help menu for step-by-step instructions.


CAVEATS
=========

This software is designed to be a simple 'quick-look' tool and
ignores lots of effects such as non-co-planar arrays, time averaging,
multi-frequency synthesis etc. It also sets the sampling grid in the
Fourier domain from the extent and pixel spacing of the input image,
which can result under-sampled synthesised beams and artefacts. Tasks
to perform robust simulations of an observation exist in CASA, MIRIAD
and AIPS, but are more difficult to use.


CONTACT:
==========
Questions or comments should be directed to 'cormac.purcell (at) mq.edu.au'.


THE CALCULATION MODULE
-----------------------------------------------------------------------------

The calculations underlying the graphical application are separated into the
file ``vriCalc.py`` to facilitate use with alternative interfaces. For example,
observations may be simulated from the basic ipython shell as follows:

.. code:: python 

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
