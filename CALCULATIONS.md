# The ObservationManager Class

The calculations underlying the graphical application are separated into the
file ``vriCalc.py`` to facilitate use with alternative interfaces. For example,
observations may be simulated from the basic ipython shell as follows:

```python 
	  
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
```