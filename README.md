# GLACiAR2

Overview 
----------

``GLACiAR2`` is an open-source python3 tool for injection-recovery simulations that determine completeness in galaxy surveys.
It's an updated version of GLACiAR developped by Carrasco et al. 2018 (https://arxiv.org/abs/1805.08985, https://github.com/danielacarrasco/GLACiAR).

Requirements
----------

Python and SExtractor are required to run the program. We also recommend DS9 to visualise the images.

We suggest users to download Anaconda https://www.anaconda.com/download/ which includes all the softwares needed to run ``GLACiAR``. 

``GLACiAR2`` uses a module that may not be included in most python installations, which is pysysp. User may need to edit the source code of pysysp so that it is compatible with python3. 

Running GLACiAR 
----------

1. Download the source code from GitHub
2. Modify or create a parameter file (e.g. ``parameters.yaml``).
3. Create a directory that contains the images for all different fields and bands required. In this directory, a subfolder ``Results`` will be created by ``GLACiAR2``. All the files created will be saved here. 
4. Make sure that the throughput tables for all bands and the PSF images are in the ``Files`` folder.
5. Modify the SExtractor file ``parameters.sex`` in the folder ``SExtractor_files`` if needed.
5. Modify '`dropouts.py'` if needed.
6. Run ``python completeness.py -s source-extractor parameters.yaml``. The -s option is for inputting the command that runs ``SExtractor`` (default = 'source-extractor').

Parameters
----------

Modify or create a yaml parameter file (e.g. 'parameters.yaml').

**Parameters files**
- *LF_shape:* List of the underlying distributions of the injected galaxies. The choices are the following (default = *schechter_flat*).
    - *flat* - All input magnitude bins have the same number of injected galaxies, equal to *n_galaxies* x *n_iterations*. Each M<sub>UV</sub> is sampled from a uniform distribution.
    - *schechter_flat* - All input magnitude bins have the same number of injected galaxies but each M<sub>UV</sub> is sampled from the Schechter function specified in ``Files/LF_Schechter_params.txt``
    - *schechter* The number of galaxies in each magnitude bin follows the specified Schechter function. Each M<sub>UV</sub> is also sampled from the Schechter function.
    - *linear* The number of galaxies in each magnitude bin follows a linear function with a slope *lin_slope*. ach M<sub>UV</sub> is also sampled from the linear function.
    - *exp* The number of galaxies in each magnitude bin follows an exponential function with an exponential base *exp_base*. Each M<sub>UV</sub> is also sampled from the exponential function.  
For the latter three options, the number of the injected galaxies in the brightest input magnitude bin is *n_galaxies* x *n_iterations*. Use these features with caution because the number of injected galaxies can get extremely large.
- *n_galaxies:* When *LF_shape* is *flat* or *schechter_flat*, it is the number of galaxies to be injected per image per iteration. Otherwise, it is the number of galaxies to be injected per iteration at the brightest magnitude bin -- see above (default=100).
- *n_iterations:* Number of iterations, i.e., the number of times the simulation is going to be run on each image for galaxies with the same redshift and magnitude. The magnitude of the galaxies injected in each iteration will be drawn based on the specified LF_shape (default = 20).
- *n_inject_max* Optional parameter to specify the maximum number of galaxies to inject at once. If the number of galaxies to be injected in an iteration is larger than this number, the code will inject the galaxies in batches to avoid overcrowding the image with the injected sources (default = *n_galaxies*). 
- *mag_bins:* The numbers of magnitude bins wanted. For a simulation run from m1 = 24.0 to m2 = 30.0 in steps of 0.5 magnitudes, there will be 13 bins (default = 13).
- *min_mag:* Brightest observed magnitude of the simulated galaxies (default = 24.0).  see *extended_mag_bins_low/high*
- *max_mag:* Faintest observed magnitude of the simulated galaxies (default = 30.0).
- *z_bins:* The numbers of redshift bins wanted. For a simulation run from z1 = 9.5 to m2 = 10.5 in steps of 0.2 magnitudes, there will be 6 bins (default = 15).
- *z_bins* Number of redshift bins (default = 16)
- *min_z:* Minimum redshift of the simulated galaxies (default = 7.5).
- *max_z:* Maximum redshift of the simulated galaxies (default = 9.0).
- *ref_uv_wl:* The wavelength in angstrom at which M<sub>UV</sub> is determined (default = 1600). 
- *n_bands:* How many filters the images have been observed in. If not specified, it will raise an error.
- *detection_band:* This is the band in which the objects are identified. If the detection image is a coadd of multiple bands, put 'det'. If not specified, it will raise an error.
- *bands:* Name of the bands. If the *detection_band* is not 'det', the detection band has to go first. If not specified, it will raise an error.
- *detection_band_combination:* Required if detection_band is 'det'. List of bands used in coadding to create the detection image.
- *coadd_type:* Required if detection_band is 'det'. 1 for a simple coadd, 2 for a noise-equalized coadd (e.g. Whitaker2019). 
- *zeropoints:* Zeropoint value corresponding to each band. The default value is 25 for each band. They are used to assign pixel values of the simulated galaxies and run SExtractor (Default = 25.0).
- *gain_values:* Gain values for each band. If not specified, it will raise an error. The gain values will be used in the SExtractor configuration file.
- *size_pix:* Pixel scale in arcsecond for the images (default = 0.08).
- *margin:* margin in arcsecond where a box of size 2*margin around the injected position will be searched for a source (sources) in the Sextracted file. The selected source will be the source in this search box whose position is closest to the injected position.
- *list_of_fields:* Text file containing the name of the fields where the simulation is going to be run. If not specified, it will raise an error.
- *path_to_images:* Directory where the science images are located.
- *path_to_results:* Directory where outputs will be placed. If not specified, it will raise an error.
- *image_name:* Heading name of the images. The naming of all science and rms images should be as follows: *image_name*_{field name}_{band name}+*imfits_end* (or *rmsfits_end* (required).
- *imfits_end:* See above for the naming of the files.
- *rmsfits_end':* See above for the naming of the files.
- *fixed_psf:* If the images are psf-matched, put name of the common psf file here. The psf image should be put in folder Files. Leave as blank if the images are not PSF-match. In that case, each image will be convolved with ``psf_{band name}.fits`` instead.
- *R_eff:* Effective radii of the galaxies at redshift 6 in kpc (default = 1.075).
- *beta_mean:* Mean of the UV slope. The injected galaxies will have spectra with slopes drawn from this mean.
- *beta_sd:* Standard deviation of the UV slopes. 
- *types_galaxies*: Number of different galaxy light profiles to be created per injection (default = 2).
- *sersic_indices:* Corresponding Sérsic indices for each *types_galaxies* (default = [1,4]).
- *fraction_type_of_galaxies:* Fraction of galaxies corresponding the the Sérsic indices given (default = [0.5,0.5]). The sum must be 1.
- *de_Vacouleur:* True of False. If true, when sersic indices == 4, the de Vacouleur profile will be used instead. (default = True).
- *ibins:* number of inclinations. The galaxies will be created for ibins inclinations with values [0,1,..,ibins-1]*0.5pi/ibins
- *ebins:* number of eccentricities. The galaxies will be created for ebins eccentricities ranging from 0 to 1.
- *min_sn:* Minimum (isophotal) signal to noise ratio in the detection band (or the first band listed in detection_band_combination) for a galaxy to be considered detected. All detected galaxies that are not blended will be run through dropout selection. So if you have S/N criteria in the dropout selection, make sure that this min_sn is safely smaller than those in the dropout selection. (default = 3.0) 
- *dropouts:* True or False. Boolean that indicates whether the user wants to run a dropout selection (default = False).
- *droptype:* Type of dropout. See dropouts.py You can customize the dropout selection criteria.
- *lin_slope:* Will be used if LF_shape is *linear*
- *exp_base:* Will be used if LF_shape is *exp*

Files
----------

Files required to run ``GLACiAR2``.

- Science images: Files with the observed images of the survey. Include all the fields and filters. See parameter *image_name* for the naming format. If the observed field of view does not cover all pixels in the image, the pixels outside of the observing area must have zero values. Otherwise, the injected galaxies may be placed in those area.

- List: Text file with the names of the fields to compute the completeness. The file name of the list should match the parameter *list_of_fields*. The file should contain at least one field.

- SExtractor parameters: ``SExtractor`` parameter file. It should be the same file used to detect the sources in the survey but the user should leave specific fields to each image/filter blank i.e. *CATALOG_NAME*, *WEIGHT_IMAGE*, *MAG_ZEROPOINT*, *GAIN*, *CHECKIMAGE_TYPE*, and *CHECKIMAGE_NAME*. We provide an example in ``SExtractor_files``.

- RMS maps: They describe the noise intensity at each pixel in relation to the science image. See parameter *image_name* for the naming format.

- PSF: Image of the PSF that we'll be used to convolve with the simulated galaxy. It depends on the instrument, telescope, and filter. If the images are not PSF-matched, provide one PSF for each band with in Files directory with the name ``psf_{band name}.fits``. If the images are PSF-matched, provide one PSF image and specify the filename under *fixed_psf* in the parameter file.

- Filter throughputs: throughput of each band. They are needed to be placed in Files directory with the name ``{band name}.tab``. Most of the hst bands are provided. See these example for format.

Modules
----------

The code is made of several modules, with the main mudule ``completeness.py``. Their descriptions are as follow.

- ``run_sextractor.py``: This module runs ``SExtractor`` in dual mode with the command ``source-extractor``. If your computer runs ``SExtractor`` with other command, you should run ``GLACiAR2`` with the command ``python completeness.py -s other_command parameters.yaml``.

- ``write_conf_files.py``: This module is called by ``run_sextractor.py``. It reads the provided SExtractor parameters file ``SExtractor_files/parameters.sex`` amd edits the parameters to be specific to each image and band (such as RMS image name, zeropoint, and gains). It then writes a new temporary configuration file to be used by ``SExtractor``.

- ``creation_of_galaxy.py``: This module performs most of the mathematics processes involved in the code. It calculates the number of galaxies to be generated according to the *LF_shape*, creates light profiles of the injected galaxies, and generates the random positions for the simulated sources.

- ``blending.py``: It identifies the detection and blending status of a source by comparing the segmentation map of the original science images to the segmentation map of the simulated images. It runs for each iteration, save the injected and recovered information of all the galaxies in that iteration.

- ``dropouts.py``: If the parameters *dropouts* is set to "True" by the user, ``blending.py`` will call this module to check if the injected galaxies pass the dropout selection criteria. The recovered magnitudes, signal-to-noise ratios, and status of all recovered galaxies will be passed to this module. The user can modify this module to the desired selection criteria.

- ``plot_completeness.py``: This module produces the plots. It always produces the plot of the completeness *C(M)* (without the dropout selection criteria) as a function of redshift and intrinsic magnitude. If the parameters *dropouts* is set to "True" by the user, it will also produce the plot of the completeness with the dropout selection.

- ``completeness.py``: This is the main module. It reads the input parameter file and calls ``creation_of_galaxy.py`` to calculate the number and flux for the simulated galaxies. It then creates the stamps with the calculated fluxes, places the galaxies in the given positions, and adds them to the science images. It later calls ``run_sextractor.py`` to run ``SExtractor`` on the images and ``blending.py`` to identify the recovered sources.

How it works
----------

This process is repeated a given number of iterations for each magnitude bin, redshift bin, and number of fields.

``GLACiAR2`` first runs ``SExtractor`` on the original science image in dual mode and creates the original catalogue. It then creates artificial galaxies stamps according to the user specified magnitude bins, redshift bins, radii, UV slope, and light profiles. In each iteration, the number of injected galaxies and the sampled magnitudes can be specified to follow either a uniform, a Schechter, a linear or an exponential probability distribution. All galaxies in an iteration will have the same magnitude, spectrum, radii, and redshift but different light profiles. These artificial sources and then added to the science image of each band at random positions. It then runs ``SExtractor`` again on the new images. 

The programme then compares the new catalogue to the original catalouge by checking a grid centered in the input position of the galaxies. The recovering status of the injected galaxies are as follows:

- Detected and isolated (status = 0)
- Detected but blended with a fainter object by less than 25% of the injected pixel area (status = 1).
- Detected but blended with a brighter object by less than 25% of the injected pixel area and that the recovered flux is within 25% of the input flux (status = 2).
- Detected but blended with a brighter object by more than 25% of the injected pixel area (status = -1).
- Detected but blended with a fainter object by more than 25% of the injected pixel area (status = -2).
- Detected below the S/N threshold (status = -3).
- Not detected by *SExtractor* (status=-4).

If required by the user, the programme can run a dropout selection criteria (dropout.py).

Outputs
----------

The code produces a set of files, images and tables. All are located in the specified *path_to_results* folder but some of them will be deleted on the fly to save memory space. We outline them below.

- **New images**.

The new images are the original images with the simulated galaxies added. Each iteration in the simulation produces one image per filter. They are saved and then deleted immediately after SExtractor runs on them. Location = *Results/images/sersic_sources\*.fits*

- **SExtractor-related files**

These are the files generated by *SExtractor*, which are the catalouges in fits format and the segmentation maps:
1. The catalouges for the original images are in the *SciImages* folder.
2. The catalogues for the images with artificial images are in the *Results/Catalogs/* folder. We only include the artificial sources that are detected. Other sources are deleted to save the memory space. For running *SExtractor*, we recommend using the provided default3.param file that detailed the catalog contents. The included fields are the following.
   * NUMBER   Running object number  
   * FLUX_ISO   Isophotal flux   [count]
   * FLUXERR_ISO   RMS error for isophotal flux   [count] 
   * MAG_ISO   Isophotal magnitude   [mag] 
   * MAGERR_ISO   RMS error for isophotal magnitude   [mag] 
   * FLUX_APER   Flux vector within fixed circular aperture(s)   [count] 
   * FLUXERR_APER   RMS error vector for aperture flux(es)   [count] 
   * MAG_APER   Fixed aperture magnitude vector   [mag] 
   * MAGERR_APER   RMS error vector for fixed aperture mag.   [mag] 
   * FLUX_AUTO    Flux within a Kron-like elliptical aperture   [count] 
   * FLUXERR_AUTO   RMS error for AUTO flux   [count] 
   * MAG_AUTO   Kron-like elliptical aperture magnitude   [mag] 
   * MAGERR_AUTO   RMS error for AUTO magnitude   [mag] 
   * KRON_RADIUS   Kron apertures in units of A or B                      
   * BACKGROUND   Background at centroid position   [count] 
   * X_IMAGE   Object position along x   [pixel] 
   * Y_IMAGE    Object position along y   [pixel] 
   * ALPHA_J2000   Right ascension of barycenter (J2000)   [deg] 
   * DELTA_J2000   Declination of barycenter (J2000)   [deg] 
   * A_IMAGE   Profile RMS along major axis   [pixel] 
   * B_IMAGE   Profile RMS along minor axis   [pixel] 
   * THETA_IMAGE   Position angle (CCW/x)   [deg] 
   * FWHM_IMAGE   FWHM assuming a gaussian core   [pixel] 
   * FWHM_WORLD   FWHM assuming a gaussian core   [deg] 
   * FLAGS   Extraction flags 
   * CLASS_STAR   S/G classifier output 
   * FLUX_RADIUS   Fraction-of-light radii   [pixel]
3. The segmentation maps (fits images) indicating the source boundaries and thier indices in the catalouges. Unless indicated otherwise, they will be deleted after each iteration. To keep these maps, run *GLACiAR2* with an option -m False e.g. ``python completeness.py -m False parameters.yaml``.

- **GLACiAR2 catalogues**

``GLACiAR2`` outputs three main tables. Their description is the following.

The first table ``{LF name}_RecoveredStats_cat{field name}.cat`` contains the summary statistics on the counts and fractions of the recovered sources and the dropouts (if specified). For each magnitude and redshift bin, tt shows the number of injected sources, and the number of recovered sources at each status. It also show the recoverd fraction and the dropout fraction (i.e. the completeness as a function of intrinsic magnitude).

Below are the columns in this table.
   
1. z: Input redshift of the simulated galaxy.
2. m: Central value of the magnitude bin.
3. N_Obj: Number of objects injected in the corresponding redshift and magnitude bin in all the iterations.
4. St=xx: Number of artificial sources recovered with status xx (7 columns).
5. N_Dropouts: Number of recovered sources with status >= 0 that passed the dropout selection criteria.
6. N_Recovered: Number of recovered sources with status >= 0.
7. Dropout_frac: Fraction of artificial sources that passed the selection criteria: N_Drop/N_Obj.
8. Rec: Fraction of not recovered artificial sources : N_Rec/N_Obj.

The second type of tables, ``{LF name}_RecoveredStats_cat{field name}_nout.cat`` and ``{LF name}_RecoveredStats_cat{field name}_dout.cat``, are the output matrices produced by the code. The \*_nout.cat shows the number of galaxies with injected magnitude M whose recovered with status >=0 and apparent magnitude m in the detection band -- N(M,m)-- at each redshift. The \*_dout.cat shows similar matrix but for the galaxies that also pass the selection criteria. For ease of use, we also provide the two tables in a Python-specific binary format (pickle file) as ``GLACiAR_output_{field name}.pickle``.

To easily check or visualize the output, the code also produces one ds9 region file for each iteration in /Results/SegmentationMaps. The region file indicates locations of the injected sources, their ``SExtractor`` ID, recovered magnitude, and status. 

- **Plots**
``GLACiAR2`` also produces a plot of the completeness, and a plot of the completeness with the selection criteria if the boolean dropouts parameters is set to True. These are the completeness as a function of intrinsic magnitude C(M).

**Note on pysysp package***
The current version ``pip install pysysp`` is for python2. To make the code runs with python3, a few lines in the installed pysysp need to be changed.
1. In pysysp.py, change "import pyfits" to "import astropy.io.fits as pyfits"
2. In pysysp.py, change "import extinction as extlaws" to "from . import extinction as extlaws"
3. In __init__.py, chage "from pysysp import StarSpectrum, BandPass, showfilters, listlaws" to "from .pysysp import StarSpectrum, BandPass, showfilters, listlaws"

