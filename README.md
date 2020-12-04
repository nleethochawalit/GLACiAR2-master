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
- *lin_slope:* Will be used if LF_shape is 'linear'
- *exp_base:* Will be used if LF_shape is 'exp'

Files
----------

Files required to run ``GLACiAR2``.

- Science images: Files with the observed images of the survey including all the fields and filters. 

- List: Text file with the names of the fields from the survey. This list is given as an input parameter. Its minimum length is one, i.e. the name of one field.

- SExtractor parameters: As previously explained, one of the steps of the code involves running ``SExtractor`` on the original images and the images with simulated galaxies. To run the software, a file defining the parameters is required. There is an example provided under ``SExtractor_files``, but we recommend for the user to change it according to their data.

- RMS maps: They describe the noise intensity at each pixel in relation to the science image. Although they are necessary only if required for the SExtractor parameters, it is strongly recommended that one of these is used in order to improve the source detection.

- PSF: Image of the PSF that we'll be used to convolve with the simulated galaxy. It depends on the instrument, telescope, and filter. If the images are not PSF-matched, provide one PSF for each band with in Files directory with the name ``psf_{band name}.fits``. If the images are PSF matched, provide one PSF image and specify the filename under *fixed_psf* in the parameter file.

- Filter throughputs: throughput of each band. They are needed to be placed in Files directory with the name ``{band name}.tab``. Most of the hst bands are provided. See these example for format.

Modules
----------

The code is made of several modules which are called by a main module, ``completeness.py``. Their descriptions are as follow.

- ``write_conf_files.py``: This module is called by ``run_sextractor.py``. It reads the provided SExtractor parameters file ``SExtractor_files/parameters.sex``. It then edits the parameters specific to each band (such as RMS image name, zeropoint, and gains), and writes new temporary configuration files. These are used by ``SExtractor``.

- ``run_sextractor.py``: This module is called by ``completeness.py``. It runs ``SExtractor`` in dual mode. This means that the photometry is performed on the location of the sources found in the detection band. The current code run SExtractor with the command ``source-extractor``. If your computer runs ``SExtractor`` with the command ``sex``, you should run ``python completeness.py -s sex parameters.yaml``.

- ``creation_of_galaxy.py``: This module performs most of the mathematics processes involved in the code. It calculates the number of galaxies to be generated according to the LF_shape, creates light profiles of the injected galaxies, and generates the random positions for the simulated sources.

- ``blending.py``: It identifies the detection and blending status of a source. It does this by comparing segmentation maps of the original science images to the segmentation maps of the simulated images. It runs for each iteration, save the injected and recovered information of all the galaxies in that iteration.

- ``dropouts.py``: If the parameters *dropouts* is set to "True" by the user, ``blending.py`` will call this module to check if the injected galaxies pass the dropout selection criteria. The recovered magnitudes, signal-to-noise ratios, and status of all recovered galaxies will be passed to this module. The user can modified this module to the desired selection criteria.

- ``plot_completeness.py``: This module produces the plots. It always produces the plot of the completeness *C(M)* (without the dropout selection criteria) as a function of redshift and intrinsic magnitude. If the parameters *dropouts* is set to "True" by the user, it will also produce the plot of the completeness (with the dropout selection).

- ``completeness.py``: This is the main module. It manages the files and calls ``creation_of_galaxy.py`` to perform the mathematical operations in order to calculate the flux for the simulated galaxies and expected magnitudes. It also calls ``run_sextractor.py`` to run the source identification software on the images as well as ``plot_completeness.py`` in order to produce plots. On its own, this module produces the features of the artificial galaxies according to the input parameters. It opens the images, create the stamps with the calculated fluxes, places the galaxies in the given positions, and adds them to the science images. Then calls the module to run SExtractor, and records the statistics regarding the recovery of sources and also the individual status and extracted properties of the simulated objects. It then produces the final tables and plots.



How it works
----------

In order to estimate the completeness of a survey, it is necessary to quantify the fraction of galaxies that are not being detected. 
``GLACiAR`` simulates artificial galaxies with particular features (inputted by the user through a parameters file) and adds them to the images of the survey. Afterwards, a source identification software is run over the original science images and the images with the simulated sources. Both catalogues are compared and the fraction of recovered artificial galaxies is measured. 

This process is repeated a given number of iterations for each magnitude bin, redshift bin, and number of fields.
For each image, there is also a number of galaxies placed in it.

An artificial galaxy stamp is created following the parameters given by the user, such as magnitude, redshift, and radius. This source is then added to the science image in random positions and the new image is saved as a new file. This file is created for each band for which the magnitude is calculated as well. To calculate the expected magnitude in each band, the program generates a synthetic spectrum, and then measures the expected flux in each band. All galaxies in one iteration have the same UV slope, which is done in order to save time.

Now, the new image information is ready to be extracted. In order to do so, \texttt{SExtractor} runs on the original science image first, creates a catalogue, and then runs over all the new images with the same SExtractor parameters. 

SExtractor runs in dual mode, this means that there is an identification band. All the sources are identified in this band first and then the information of the source (or absence of) in said position is calculated from the images on each band.

The way in which the catalogues are compared is by checking a grid centered in the input position of the galaxies. Depending on whether there is a source detected in that grid in either the segmentation map of the original science image or the segmentation map of the science image with the added artificial galaxies, a status will be assigned to the artificial galaxies. The possible statuses are:

- Detected and isolated.
- Detected and blended with a fainter object.
- Detected and blended with a brighter object.
- Detected below the S/N threshold.
- Not detected.

If required by the user, the program can run a redshift selection criteria which we set here as an example, and it is from Bernard et al (2016.)

Finally, plots are generated for C(m), S(z,m), and C(m)S(z,m).

Outputs
----------

The code produces a set of files, images and tables. Some of them are deleted for space reasons, while others are kept as a final result. We outline them below in order of appearance.

- **New images**.

The new images are the created images that contain the original observed data from the survey with the simulated galaxies added. Each iteration in the simulation produces one image. These are saved and then deleted immediately after SExtractor runs on them.

- **SExtractor Catalogues**

The source identification process generates two sets of files. One of them is the catalogues (described here), and the other one is segmentation maps.

The catalogues include a list with all the sources identified in the images and their parameters. Due to the structure of the code, the first file created for a field of the survey contains the information of the real sources as it is run on the original science image. Besides this, for each iteration of the simulation a catalogue with the information of the sources is created. The header of the catalogue is shown below:

1. NUMBER   Running object number   
2.   FLUX_ISO   Isophotal flux   [count]
3.   FLUXERR_ISO   RMS error for isophotal flux   [count] 
4.   MAG_ISO   Isophotal magnitude   [mag] 
5.   MAGERR_ISO   RMS error for isophotal magnitude   [mag] 
6.   FLUX_APER   Flux vector within fixed circular aperture(s)   [count] 
11.   FLUXERR_APER   RMS error vector for aperture flux(es)   [count] 
16.   MAG_APER   Fixed aperture magnitude vector   [mag] 
21.   MAGERR_APER   RMS error vector for fixed aperture mag.   [mag] 
26.   FLUX_AUTO    Flux within a Kron-like elliptical aperture   [count] 
27.   FLUXERR_AUTO   RMS error for AUTO flux   [count] 
28.   MAG_AUTO   Kron-like elliptical aperture magnitude   [mag] 
29.   MAGERR_AUTO   RMS error for AUTO magnitude   [mag] 
30.   KRON_RADIUS   Kron apertures in units of A or B                      
31.   BACKGROUND   Background at centroid position   [count] 
32.   X_IMAGE   Object position along x   [pixel] 
33.   Y_IMAGE    Object position along y   [pixel] 
34.   ALPHA_J2000   Right ascension of barycenter (J2000)   [deg] 
35.   DELTA_J2000   Declination of barycenter (J2000)   [deg] 
36.   A_IMAGE   Profile RMS along major axis   [pixel] 
37.   B_IMAGE   Profile RMS along minor axis   [pixel] 
38.   THETA_IMAGE   Position angle (CCW/x)   [deg] 
39.   FWHM_IMAGE   FWHM assuming a gaussian core   [pixel] 
40.   FWHM_WORLD   FWHM assuming a gaussian core   [deg] 
41.   FLAGS   Extraction flags 
42.   CLASS_STAR   S/G classifier output 
43.   FLUX_RADIUS   Fraction-of-light radii   [pixel]


- **Segmentation maps**

A segmentation map is a map with the definition of the location of a source and its borders found by a source identification software. In this case, they are produced by SExtractor. For ``GLACiAR``, they can be classified in two groups: the segmentation maps from the original science images of the survey, and the ones from the images that include the simulated galaxies.

Together with the old catalogues of the sources a segmentation map is created. This is a .fits file and it is kept.

The images with the simulated galaxies will produce new segmentation maps with the same characteristics. The only difference with the catalogue described above is the inclusion of the new detected sources.
This file is discarded in order to save space.

- **GLACiAR catalogues**

The main results produced by ``GLACiAR`` can be summarised in three tables. These contain the information of the inserted galaxies, including their input properties and their output status and properties as well. They are described below.

The first table contains information on the statistics of the results in terms of the recovered sources and the dropouts (if specified). It counts the amount of sources inserted per redshift bin and magnitude bin and it keeps track of the amount that corresponds to each detection status. It also calculates the amount of total recoveries over all the inserted simulated sources, C(m). If applicable, it counts the number of dropouts for each bin, and calculates the fraction of these over the total amount of inserted galaxies S(z,m), and over the recovered dropouts in the required redshift range over the number of recovered simulated galaxies C(m)S(z,m). Below there is an example of its structure with a brief description of the columns.

   
1. z: Input redshift of the simulated galaxy.
2. m: Median value of the magnitude bin.
3. N_Obj: Number of objects inserted for the corresponding redshift and magnitude bin in all the iterations.
4. S(0): Number of artificial sources recovered by \texttt{SExtractor} that were isolated.
5. S(2,1): Number of artificial sources recovered that were blended with a fainter object.
6. S(-1): Number of artificial sources recovered that were blended with a brighter object.
7. S(-2): Number of artificial sources that were detected by \texttt{SExtractor} with a $S/N$ under the required threshold.
8. S(-3): Number of artificial sources that were not detected by \texttt{SExtractor}.
9. N_Rec: Number of recovered artificial sources, S(0)+S(2,1).
10. N_Drop: Number of artificial sources that passed the dropout selection criteria.
11. Rec: Fraction of not recovered artificial sources : N_Rec/N_Obj.
12. Drop: Fraction of artificial sources that passed the selection criteria: N_Drop/N_Obj.


Complementary to the previous table, the algorithm produces a table with the each one of the inserted sources, their positions, the input magnitude, the blending status, and the detection status. Several tables (one for each redshift step) are produced with all the galaxies that were placed in the simulations at that redshift. It provides an insight to understand the characteristics or reasons to detect or miss an object. It also yields their detected magnitude in the detection band, and their size. In summary, it gives us information about how this source is detected instead of the input information about it. 

1. Initial Mag: Magnitude corresponding to the input flux for the star. This is not the same as Input Magnitude since the input magnitude changes depending on the beta value and size of the object.
2. Iteration: Iteration number.
3. ID Number: Identification number given by SExtractor after it runs on the image with the simulated galaxies. This number is unique for every iteration for a given magnitude and redshift.
4. Input Magnitude: Magnitude corresponding to the added flux inside all the pixels that the source includes.
5. Output Magnitude: Magnitude of the source found with SExtractor after it runs on the image with the simulated galaxies.
6. Id Status: Integer number that indicates whether a source has been recovered and/or is blended.


One last table, which is useful for redshift selection, is produced. Given that the number of bands is variable  this table is released in a Python-specific compact binary representation (using the pickle module). It does contain the ID of the object, input magnitude, status, magnitudes in all bands, and S/N for each band as well. This is an important file for redshift estimations/selection techniques. Similar to the previous tables, it contains information about the output/measured characteristics of the detected objects. It contains all the magnitudes in different bands, so it is especially useful for photometric estimations. Unlike the other two tables, this is not an ASCII file as that is not efficient and it requires too many resources.


- **Plots**

``GLACiAR`` also produces a plot of the completeness and two extra plots if the boolean dropouts parameters is set to True. The first plot corresponds to the completeness function C(m) as a function of the magnitude and the redshift. The second and third plot are the S(z,m) and S(z,m)C(m). This is only produced in case the dropout technique is applied, but given the tables produced by ``GLACiAR``, it is easily calculable with the final catalogues.

**Note on pysysp package***
The current version from pip install pysysp is for python2. To make the code runs, a few lines in the installed pysysp need to be changed.
1. In pysysp.py, change "import pyfits" to "import astropy.io.fits as pyfits"
2. In pysysp.py, change "import extinction as extlaws" to "from . import extinction as extlaws"
3. In __init__.py, chage "from pysysp import StarSpectrum, BandPass, showfilters, listlaws" to "from .pysysp import StarSpectrum, BandPass, showfilters, listlaws"

