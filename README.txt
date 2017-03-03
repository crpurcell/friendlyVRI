#-----------------------------------------------------------------------------#
#
# A friendly Virtual Radio Interferometer tool.
# by Cormac R. Purcell and Roy Truelove (Macquarie University, Sydney)
#
# Copyright (c) 2017 Cormac R. Purcell and Roy Truelove
#
#-----------------------------------------------------------------------------#

This application is built to help astronomers investigate the effect of
combining different array configurations when observing an astronomical object.

The graphical interface is written using the tkinter library for maximum
portability.

TODO:

* ScatterPlot widget:
  - Format of text labels and margins
  - Alternative vertical text for buggy python3
  - equalise X & y scales
* TreeTables:
  - Format of column spacing & text.
* DoubleSlider:
  - Change to asymmetric handle shape.
  - Snap to grid.
  - Bind entry boxes & limit significant figures.
* Observing parameters:
  - Limit max & min values for all.
  - Check for numbers & disallow non-numeric characters,
* Model browser & observation control
  - relative path in file entry box ( + make longer )
  - move all controlls to controller window.
  - Add an information panel with:
    -- max & min size scale
    -- synthesised beam size
    -- ideal sampling parameters
    -- image sampling parameters
    -- warning about miss-match between model & observing parameters
  - Nice to have:
    -- graphical comparison of model size scale (structure fn?, dispersion)
       versus uv-coverage (histogram of baselines?).
    -- shadowing calculation.
    -- robust weighting calculation
    -- fitted synthesized beam size
    -- Elevation plot for selected telescopes & uv-coverage 
* Plotting
 - scalebar on images rather than axes
 - primary beam circles on images.
 - fix colours on uv-coverage plot
 - legend for uv-coverage plot
 - zoom & pan bindings
 - images resize when window is maximized
 - options to save publication quality plots
 - zoom in on synthesized beam (1/3rd of plot)

* Calculations and logic:
 - mask uv-coverage for elevation limits.
 - status lights
 - robust weighting calculation.
