"""
Runs SExtractor on the images with the simulated galaxies.
"""

import write_conf_files
from subprocess import check_call
from subprocess import PIPE

def science_image(name_band, detection_band, zp, g, path_to_im, 
                  path_to_results, image_name, cat, imfits_end ='_v1_drz.fits', 
                  rmsfits_end='_v1_rms_scl.fits'):
    """
    Runs SExtractor on the science image with the parameters given by the user
    and creates a catalog for the sources found in the image.
    Args:
        name_band (string) = name of the band in which the image is taken.
        detection_band (string) = name of the detection band given in the
                                  parameters file.
        zp (float) = Zeropoint value for the respective band from the
                     parameters file.
        g (float) = Gain value for the respective band. It is given
                     in the parameters file.
        path_to_im (string) = Path to the folder with the science images.
        path_to_results (string) = Path to the folder with the GLACiAR results.
        niter (integer) = Number of the current iteration.
        image_name (string) = Name of the science image (preceding the name
                              of the band) as given in the parameters file.
        cat (string) = Name of the field for which the simulation is being run.
    """
    # New parameters files are created for SExtractor.
    for i in range(len(name_band)):
        #write sextractor parameter file for the band
        write_conf_files.main(name_band[i], detection_band, zp[i], g[i],
        	path_to_im, path_to_results, 0, image_name, cat, rmsfits_end)
        # If the band SExtractor is running on is not the detection band,
        # run in dual mode with the detection band as reference.
        # NL changed to dual mode in either case cause SExtractor single mode
        # gives different number of detections from dual mode!
        if name_band[i] == detection_band: 
            p = check_call(['source-extractor ' + path_to_im + image_name + cat +'_' +
                        name_band[i] + imfits_end+',' 
                        + path_to_im + image_name + cat +'_' +
                        name_band[i] + imfits_end+' -c SExtractor_files/' +
                        'parameters_' + name_band[i] +
                        '.sex -CHECKIMAGE_NAME ' + path_to_results +
                        'SciImages/sources_' + cat + '_' + name_band[i] +
                        '_segm.fits -CATALOG_NAME '+ path_to_results +
                        'SciImages/sources_' + cat + '_' + name_band[i] +
                        '.cat'], bufsize=4096, stdin=PIPE, stdout=PIPE,
                        close_fds=True, shell=True)
        else:
            p = check_call(['source-extractor ' + path_to_im + image_name + cat+'_' +
                        detection_band + imfits_end+',' + path_to_im + 
                        image_name + cat + '_' + name_band[i] + imfits_end +
                        ' -c SExtractor_files/parameters_' +
                        name_band[i] + '.sex -CHECKIMAGE_NAME ' +
                        path_to_results + 'SciImages/sources_' + cat + '_' +
                        name_band[i] + '_segm.fits -CATALOG_NAME ' +
                        path_to_results + 'SciImages/sources_' + cat + '_' +
                        name_band[i]+'.cat -VERBOSE_TYPE QUIET'], 
                        bufsize=4096, stdin=PIPE, stdout=PIPE,
                        close_fds=True, shell=True)

def main(name_band, detection_band, zp, g, path_to_im, path_to_results, 
         niter, image_name, cat, m1, redshift, 
         rmsfits_end = '_v1_rms_scl.fits', roundnum = 0):
    """
    Runs SExtractor on the new image containing the simulated galaxies
    with the parameters given by the user and creates a new catalog for
    the objects in the image.
    Args:
        name_band (string) = name of the band in which the image is taken.
        detection_band (string) = name of the detection band given in the
                                  parameters file.
        zp (float) = Zeropoint value for the respective band from the
                     parameters file.
        g (float) = Gain value for the respective band. It is given
                     in the parameters file.
        path_to_cat (string) = Path to the folder with the science images.
                               Given in the parameters file.        
        image_name (string) = Name of the science image (preceding the name
                              of the band) as given in the parameters file.        
        The following is used to name the results
        cat (string) = Name of the field for which the simulation is being run.
        niter (integer) = Number of the current iteration.
        m1 (float) = Initial input magnitude bin for the simulated galaxy in the
                     detection band. This is used for the name of the catalog only
        redshift (float) = Redshift for the simulated galaxy.
    """
    # New parameters files are created for each band every time SExtractor
    # is run for each iteration. They are replaced in each iteration.
    write_conf_files.main(name_band, detection_band, zp, g, path_to_im, 
                          path_to_results, niter, image_name, cat, 
                          rmsfits_end,roundnum=roundnum)
    # If the band SExtractor is running on is the detection band, identify
    # identify the sources.
    # NL changed to dual mode cause SExtractor single and dual mode give 
    # different number of detections!
    if name_band == detection_band:
        p = check_call(['source-extractor ' + path_to_results + 'Results/images/'
                        'sersic_sources_' + cat + detection_band + '.fits, '
                        + path_to_results + 'Results/images/sersic_sources_' 
                        + cat + detection_band + '.fits -c '
                        'SExtractor_files/parameters_' + detection_band +
                        '.sex -CATALOG_NAME ' + path_to_results + 
                        'Results/Dropouts/source_' + cat + '_mag%.1f'%(m1) + 
                        '_z%.1f'%(redshift) +'_i' + str(niter) + '.' + 
                        str(roundnum) + '_' + detection_band + 
                        '.cat -VERBOSE_TYPE QUIET'], 
                       bufsize=4096, stdin=PIPE, stdout=PIPE, close_fds=True,
                       shell=True)
    # If the band SExtractor is running on is not the detection band, run in
    # dual mode with the detection band as reference.
    else:
        p = check_call(['source-extractor ' + path_to_results + 'Results/images/'
                        'sersic_sources_' + cat + detection_band + '.fits,' +
                        path_to_results + '/Results/images/sersic_sources_' + 
                        cat + name_band + '.fits -c SExtractor_files/'
                        'parameters_' + name_band + '.sex -CATALOG_NAME ' + 
                        path_to_results + 'Results/Dropouts/source_' + cat + 
                        '_mag%.1f'%(m1) +'_z%.1f'%(redshift) + '_i' + str(niter) + 
                        '.' + str(roundnum) + '_' + name_band + 
                        '.cat -VERBOSE_TYPE QUIET'],
                        bufsize=4096, stdin=PIPE, stdout=PIPE,
                        close_fds=True, shell=True)
