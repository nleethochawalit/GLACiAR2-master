"""Main"""
import pdb
import numpy as np
import random
import os
import yaml
import run_sextractor
import creation_of_galaxy
import blending
import plot_completeness
import pickle
from astropy.io import fits
from astropy.convolution import convolve
from astropy.cosmology import Planck15 as cosmo
from astropy import units as u
from shutil import copyfile
import warnings
import datetime
import argparse

print(datetime.datetime.now())

#Parsing arguments
msg = ""
parser = argparse.ArgumentParser(description=msg)
parser.add_argument('parameter_file', type=str, 
                    help='Name of the yaml parameter file')
parser.add_argument("-m", "--minimal", action = "store_true", 
                    default = True, help = "Only Keep minimum number of files (default: True)")
parser.add_argument("-s", "--SExtractor_command", action = "store", 
                    default = 'source-extractor',dest = 'SExtractor_command', 
                    help = "SExtractor command (default: source-extractor)")
args = parser.parse_args()
parameter_file = args.parameter_file
minimal_file = args.minimal
SExtractor_command = args.SExtractor_command
    
print("parameter file: %s"%parameter_file)
print("SExtractor command: %s"%SExtractor_command)
print("Minimal File: %s"%minimal_file)

# Read input parameters from the file 'parameters.yaml'.
# Assign a default value when possible if they are not given by the user.
stream = open(parameter_file, 'r')
parameters = yaml.load(stream,Loader=yaml.FullLoader)
if parameters['n_galaxies'] is None:
    parameters['n_galaxies'] = 100
if parameters['n_iterations'] is None:
    parameters['n_iterations'] = 20
if parameters['mag_bins'] is None:
    parameters['mag_bins'] = 13
if parameters['min_mag'] is None:
    parameters['min_mag'] = 24.0
if parameters['max_mag'] is None:
    parameters['max_mag'] = 30.0
if parameters['z_bins'] is None:
    parameters['z_bins'] = 16
if parameters['min_z'] is None:
    parameters['min_z'] = 7.5
if parameters['max_z'] is None:
    parameters['max_z'] = 9.0
if parameters['ref_uv_wl'] is None:
    parameters['ref_uv_wl'] = 1600
if parameters['n_bands'] is None:
    raise ValueError('Number of bands needed.')
if parameters['detection_band'] is None:
    raise ValueError('Input detection band.')
if parameters['bands'] is None:
    raise ValueError('Input name of the bands.')
if parameters['det_combination'] is None:
    parameters['det_combination'] = None
if parameters['coadd_type'] is None:
    parameters['coadd_type'] = 1 #1 = simple coadd, 2 = noise-:equalized coadd
if parameters['list_of_fields'] is None:
    raise ValueError('Input file with the fields.')
if (parameters['zeropoints'] is None or
        len(parameters['zeropoints']) < parameters['n_bands']):
    parameters['zeropoints'] = np.zeros(int(parameters['n_bands']))+25
if (parameters['gain_values'] is None or
        len(parameters['gain_values']) < parameters['n_bands']):
    raise ValueError('Input gain values for each band.')
if parameters['size_pix'] is None:
    parameters['size_pix'] = 0.08
if parameters['path_to_images'] is None:
    raise ValueError('Input the directory path.')
if parameters['path_to_results'] is None:
    raise ValueError('Output the directory path.')
if parameters['image_name'] is None:
    raise ValueError('Input the name of the images.')
if parameters['imfits_end'] is None:
    parameters['imfits_end'] = '_drz.fits'
if parameters['rmsfits_end'] is None:
    parameters['rmsfits_end'] = '_rms.fits'
if parameters['fixed_psf'] is None:
    parameters['fixed_psf'] = None
if parameters['R_eff'] is None:
    parameters['R_eff'] = 1.075 #kpc at redshift 6
if parameters['beta_mean'] is None:
    parameters['beta_mean'] = -2.2
if parameters['beta_sd'] is None:
    parameters['beta_sd'] = 0.4
if parameters['types_galaxies'] is None:
    parameters['types_galaxies'] = 2
if parameters['ibins'] is None:
    parameters['ibins'] = 9
if parameters['ebins'] is None:
    parameters['ebins'] = 5
if parameters['sersic_indices'][1] is None:
    parameters['sersic_indices'][1] = 4
if parameters['fraction_type'] is None:
    parameters['fraction_type'] = [0.5, 0.5]
if parameters['margin'] is None:
    parameters['margin'] = 0.3 #arcsecond. Margin for injection/recovery
if parameters['min_sn'] is None:
    parameters['min_sn'] = 0
if parameters['dropouts'] is None:
    parameters['dropouts'] = False
if parameters['de_Vacouleur'] is None:
    parameters['de_Vacouleur'] = False
if parameters['LF_shape'] is None:
    parameters['LF_shape'] = ['schechter_flat']
if parameters['lin_slope'] is None:
    parameters['lin_slope'] = 2
if parameters['exp_base'] is None:
    parameters['exp_base'] = 2
if parameters['n_inject_max'] is None:
    parameters['n_inject_max'] = parameters['n_galaxies']
if ((parameters['coadd_type'] <= 0) or (parameters['coadd_type'] > 2)):
    raise ValueError('Input coadd_type will not be recognized.')
def delete_id_fits(original_file, ids, id_keep=False):
    """Delete sextractor entries in fits format by id number"""
    
    f1 = fits.open(original_file,ignore_missing_end=True)
    if id_keep is True:
        indices = np.where(np.in1d(f1[1].data['NUMBER'], ids, 
                                   assume_unique=True))[0]
    else:
        indices = np.where(np.in1d(f1[1].data['NUMBER'], ids, 
                                   assume_unique=True, invert=True))[0]
    
    f1[1].data = f1[1].data[indices] 
    f1.writeto(original_file,overwrite=True)
        
def delete_multiple_lines(original_file, line_numbers, line_keep=False):
    """In a file, delete the lines at line number in given list"""
    if line_keep is True:
        num_lines = sum(1 for line in open(original_file))
        line_numbers = np.setdiff1d(np.arange(num_lines),line_numbers)
        
    is_skipped = False
    counter = 0
    # Create name of dummy / temporary file
    dummy_file = original_file + '.bak'
    # Open original file in read only mode and dummy file in write mode
    with open(original_file, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
        # Line by line copy data from original file to dummy file
        for line in read_obj:
            # If current line number exist in list then skip copying that line
            if counter not in line_numbers:
                write_obj.write(line)
            else:
                is_skipped = True
            counter += 1
 
    # If any line is skipped then rename dummy file as original file
    if is_skipped:
        os.remove(original_file)
        os.rename(dummy_file, original_file)
    else:
        os.remove(dummy_file)

def create_stamps(n0, size_galaxy0, Re0, types_galaxies0, ebins0, ibins0):
    """
    Creates a galaxy following a Sersic profile. (flux is normalized to 1)

    Args:
        n0 (int array) = Array with Sersic indexes.
        size_galaxy0 = Diameter of the galaxy stamp in pixels.
        Re0 (float) = Effective radius in pixels.
        types_galaxies0 (int) = Number of Sersic indexes required.
        ibins0 (int) = Number of possible inclinations for the
                         simulated galaxy.
        ebins0 (int) = Number of possible eccentricities for the
                         simulated galaxy.
    Returns:
        galaxy_grid (float array) = Stamp of a galaxy with the
                                        corresponding flux for each
                                        pixel.
        galaxy_grid_n4 (float array) = Stamp of galaxy with n=4 if
                                      "de_Vacouleur" is True. Otherwise, set
                                      to 0.
    """
    galaxy_grid = np.zeros((types_galaxies0, ebins0, ibins0, size_galaxy0,
                            size_galaxy0),dtype=np.float32)
    galaxy_grid_n4 = np.zeros((size_galaxy0, size_galaxy0),dtype=np.float32)
    ivalues = np.arange(0, 0.5, 0.5/(ibins0)) * np.pi
    evalues = np.arange(0, 1., 1. / (ebins0))
    for i in range(types_galaxies0):
        print('Creating stamps for sersic index %d'%n0[i])
        bn = creation_of_galaxy.get_bn(n0[i])
        if (parameters['de_Vacouleur'] is True and n0[i] == 4):
            galaxy_grid_n4[:, :] = creation_of_galaxy.makeSersic(n0[i], bn,
                                                                 Re0, 0, 0,
                                                                 size_galaxy0)
        else:
            print('sersic index = %d'%n0[i])
            print('eccentricity:', evalues)
            print('inclination (radian):', ivalues)
            for j in range(ebins0):
                for k in range(ibins0):
                    #print(n0[i], evalues[j], ivalues[k])
                    print(".", end = '')
                    galaxy_grid[i, j, k, :, :] = creation_of_galaxy.makeSersic(
                                                  n0[i], bn, Re0, evalues[j],
                                                  ivalues[k], size_galaxy0)
            print('\n')
    return galaxy_grid, galaxy_grid_n4


def open_images(image_name, name_band, imfits_end):
    """
    Opens science images.

    Args:
        image_name (string) = Name of the science image.
        name_band (string) = Name of the band in which the science image is
                            taken.
    Returns:
        obs_data (float array) = Data from the science image. Each cell
                                contains the flux for that pixel.
        head_data (string array) = Array with the header from the science image
                                  opened.
    """
    hdu_list = fits.open(image_name + name_band+ imfits_end,
                         ignore_missing_end=True)
    obs_data = hdu_list[0].data
    head_data = hdu_list[0].header
    hdu_list.close()
    return obs_data, head_data


def place_gal(n0, ngal, frac_n, e0, i0, flux0, frame0, x0, y0, s0, gal_g,
              gal_g_n4, psf0):
    """
    Args:
        n0 (int array) = Sersic indices
        ngal (int) = Number of galaxies to be placed
        frac_n (float) = Fraction of galaxies with n0
        e0 (int) = Eccentricity bin. Varies between 0 and ebins-1
        i0 (int) = Inclination angle bin.
                     Varies between 0 and ibins-1.
        flux0 (float) = Total corresponding flux for the artificial galaxy.
        frame0 (array) = Empty frame with the same shape as the science image.
        x0 (int array) = array with position along the x axis for new sources.
        y0 (int array) = array with position along the y axis for new sources.
        s0 (int) = radius of the galaxy stamp.
        gal_g (array) = Stamp for an artificial galaxy with n!=4.
        gal_g_n4 (array) = Stamp for artificial galaxy with n=4.
        psf0 (array) = PSF to convolve the artificial galaxy with.
    Returns:
        frame0: (array) = Frame with the same shape as science image. 
                          It has the artificial galaxy stamps in the
                          corresponding positions.
    """

    gn = 0
    # Loop for the amount of Sersic indices.
    for i in range(len(n0)):
        # Number of galaxies with given index.
        n_fraction = int(ngal * frac_n[i])
        # Run a loop for each galaxy with the given index.
        for j in range(n_fraction):
            # Scale the galaxy flux with the real flux.
            if parameters['de_Vacouleur'] is True and n0[i] == 4:
                galaxy = gal_g_n4[:, :] * flux0
            else:
                galaxy = gal_g[i, e0[gn], i0[gn], :, :] * flux0
                
            # Convolve galaxy with the PSF.
            gconv = convolve(galaxy, psf0, normalize_kernel=True)
            
            # Add galaxy stamp to the empty frame with the centre in (x0,y0).
            large_frame = np.zeros_like(frame0)
            large_frame[int(x0[gn]-s0):int(x0[gn]+s0),
                        int(y0[gn]-s0):int(y0[gn]+s0)] = gconv
            frame0 += large_frame
            gn = gn + 1
    # Convolve image with the PSF.
    return frame0

def coadd_images(path_to_results, cat, det_combination,
                 detection_band, path_to_im, image_name,
                 rmsfits_end, coadd_type):
    """
    Coadd images in det_combination
    If coadd_type = 1: simple average
    If coadd_type = 2: noise-equalized average 
    *** If coadd_type =2, assume that the *rmsfits_end images are rms image!!
    """
    #creating template frame
    hdu_list = fits.open('%sResults/images/sersic_sources_%s_%s.fits'%
                             (path_to_results,cat,det_combination[0]),
                             ignore_missing_end=True)
    head_data = hdu_list[0].header
    shape = hdu_list[0].data.shape
    hdu_list.close()
    
    frame = np.zeros(shape,dtype=np.float32)
    
    #adding bands to frame
    if coadd_type == 1:
        for ib in range(len(det_combination)):
            imhdu = fits.open('%sResults/images/sersic_sources_%s_%s.fits'%
                             (path_to_results,cat,det_combination[ib]),
                             ignore_missing_end=True)
            frame += +imhdu[0].data   
            imhdu.close()
        
    elif coadd_type == 2:
        for ib in range(len(det_combination)):
            rmshdu = fits.open('%s%s%s_%s%s'%(path_to_im, image_name, cat, 
                                 det_combination[ib], rmsfits_end),
                               ignore_missing_end=True)
            divider = np.divide(1.,rmshdu[0].data,
                                out=np.zeros(shape,dtype=np.float32), 
                                where = rmshdu[0].data!= 30000) #assume that the region outside of observation is set at 300000
            rmshdu.close()
            imhdu = fits.open('%sResults/images/sersic_sources_%s_%s.fits'%
                             (path_to_results,cat,det_combination[ib]),
                             ignore_missing_end=True)
            frame += imhdu[0].data*divider
            imhdu.close()
            
    frame /= len(det_combination)
    outfile_conv = '%sResults/images/sersic_sources_%s_%s.fits'%(
        path_to_results,cat,detection_band)     
    # Save the new fits file with the simulated galaxies.
    fits.writeto(outfile_conv, frame,head_data, overwrite=True)
                                
def create_synthetic_fits(name_hdu_list1,band,imfits_end,m,sersic_indices,
                          ninject,fraction_type_galaxies,e_rand,i_rand,flux,
                          xpos,ypos,stamp_radius,galaxy_grid,galaxy_grid_n4,
                          path_to_results,cat,return_framegal = False,
                          fixed_psf = None):
    """
    This function place fake galaxies in the images and save as fits
    in path_to_results+'Results/images/sersic_sources_'+ cat+'_'+band+'.fits'
    
    If return_framegal = True, then it will also return the fake galaxies frame
    """
    # Empty array of the size of the science image
    # where the simulated galaxies will be placed.                               
    obs, head_obs = open_images(name_hdu_list1,band,imfits_end) 
    frame = np.zeros(obs.shape,dtype=np.float32)  
    
    # if the artificial galaxy is brighter than
    # 50 mangitudes (some flux is expected, but very
    # conservative margin), place it in the image.
    if m < 50:
        # Open PSF file to convolve with the galaxy.
        if fixed_psf is None:
            psf_file = 'Files/psf_' +band + '.fits'
        else:
            psf_file = 'Files/'+fixed_psf
        ipsf=0
        hdu_psf = fits.open(psf_file)
        psf = hdu_psf[ipsf].data
        hdu_psf.close()
        #according to astropy.convolve, Kernel size must be odd in all axes
        if psf.shape[0]%2==0:
            psf = psf[0:psf.shape[0]-1,:]
        if psf.shape[1]%2==0:
            psf = psf[:,0:psf.shape[1]-1]
            
        #placing all galaxies on array of zeros
        frame_gal = place_gal(sersic_indices,ninject,fraction_type_galaxies,
                              e_rand, i_rand, flux, frame,
                              xpos, ypos, stamp_radius, galaxy_grid,
                              galaxy_grid_n4, psf)
    # if the galaxy is fainter than 50 magnitudes
    # (undetectable), save the original image.
    else:
        frame_gal = frame
        
    # new synthesized image = observed image + simulated galaxies
    new_data2 = obs + frame_gal
    outfile_conv = '%sResults/images/sersic_sources_%s_%s.fits'%(path_to_results,cat,band)
        
    # Save the new fits file with the simulated galaxies.
    fits.writeto(outfile_conv, new_data2,head_obs, overwrite=True)
    
    if return_framegal is True:
        return frame_gal

def main(minimal_file=True,SExtractor_command='source-extractor'):
    # Open the file with the name of the fields that the simulation
    # is going to be run for.
    f = open(parameters['list_of_fields'])
    k = f.readlines()
    f.close()
    cat = [str(line.split()[0]) for line in k]

    m_output = np.linspace(parameters['min_mag'], parameters['max_mag'],
                          parameters['mag_bins']) #or observed magnitudes
    mbinsize = (parameters['max_mag']-parameters['min_mag'])/\
        (parameters['mag_bins']-1.)
    nmagbins_out = parameters['mag_bins']
    
    z_total = np.linspace(parameters['min_z'], parameters['max_z'],
                          parameters['z_bins'])  #redshifts
    nzbins = parameters['z_bins']

    #intrinsic/input absolute magnitudes
    #calculated by roughly converting the output magnitude back to intrinsic
    #and round to single digit. Assuming that the observed magnitude (m) of 
    #the detection band is m at the UV reference wavelength.

    if parameters['Minput_min'] is None:
        maxlumdis = cosmo.luminosity_distance(parameters['max_z']).to(u.parsec).value
        M_input_min = np.around(parameters['min_mag']-2.5*
                            np.log10((maxlumdis/10.)**2/(1.+parameters['max_z'])),
                            decimals=1) #brightest absolute magnitude
        M_input_min = M_input_min-mbinsize #make it one magbin lower
    else:
        M_input_min = parameters['Minput_min']
        
    if parameters['Minput_max'] is None:
        minlumdis = cosmo.luminosity_distance(parameters['min_z']).to(u.parsec).value
        M_input_max = np.around(parameters['max_mag']
                                -2.5*np.log10((minlumdis/10.)**2/(1.+parameters['min_z'])),
                                decimals=1) #faintest
        M_input_max = M_input_max+mbinsize
    else:
        M_input_max = parameters['Minput_max']
        
    nmagbins_in = np.ceil((M_input_max-M_input_min)/mbinsize).astype(int)
    M_input = np.arange(nmagbins_in)*mbinsize+M_input_min
            
    nLF = len(parameters['LF_shape'])
    
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    print('\n')
    print('LF:', parameters['LF_shape'])
    print('redshift:', z_total)
    print('m_out:',m_output)
    print('M_input:',M_input)
    print('\n')
    print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    
    margin = np.round(parameters['margin']/parameters['size_pix']).astype(int)
    
    # Creates the directories that will contain the results if they
    # don't already exist.
    if not os.path.exists(parameters['path_to_results']):
        os.makedirs(parameters['path_to_results'])
    if not os.path.exists(parameters['path_to_results']+'Results'):
        os.makedirs(parameters['path_to_results']+'Results')
    if not os.path.exists(parameters['path_to_results']+'SciImages'):
        os.makedirs(parameters['path_to_results']+'SciImages')
    if not os.path.exists(parameters['path_to_results']+'Results'
                          '/SegmentationMaps'):
        os.makedirs(parameters['path_to_results']+'Results/SegmentationMaps')
    if not os.path.exists(parameters['path_to_results']+'Results/Plots'):
        os.makedirs(parameters['path_to_results']+'Results/Plots')
    if not os.path.exists(parameters['path_to_results']+'Results/images/'):
        os.makedirs(parameters['path_to_results']+'Results/images/')
    if not os.path.exists(parameters['path_to_results']+'/Results/Catalogs'):
        os.makedirs(parameters['path_to_results']+'/Results/Catalogs')
        
    #copy parameter and Sextractor files to the Results folder
    copyfile(parameter_file,parameters['path_to_results']+parameter_file)
    copyfile('SExtractor_files/parameters.sex',
             parameters['path_to_results']+'parameters.sex')
    print('Results will be in %s'%parameters['path_to_results'])
    #copy Schechter parameter if LF shape is specified as Schechter    
    for ilf in range(nLF):
        parameters['LF_shape'][ilf] = parameters['LF_shape'][ilf].lower()
        if ((parameters['LF_shape'][ilf] == 'schechter') or 
            (parameters['LF_shape'][ilf] == 'schechter_full')): 
            copyfile('Files/LF_Schechter_params.txt', 
                     parameters['path_to_results']+'LF_Schechter_params.txt')

    #loop for each field
    for ic in range(len(cat)):    
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('running field %s \n'%cat[ic])
        
        # Empty array to save the simulation results (for pickles)
        nout_completeness = np.zeros((nmagbins_out,nmagbins_in,nzbins,nLF))
        nout_dropouts = np.zeros((nmagbins_out,nmagbins_in,nzbins,nLF))
        Nin = np.zeros((nmagbins_in,nzbins,nLF))
        total_completeness = np.zeros((nmagbins_in,nzbins,nLF))
        total_dropouts = np.zeros((nmagbins_in,nzbins,nLF))
        
        #science image names (begining)
        name_hdu_list1 = (parameters['path_to_images'] +
                          parameters['image_name']+cat[ic]+'_')
        #run SExtractor on science images
        run_sextractor.science_image(parameters['bands'],
                        parameters['detection_band'], 
                        parameters['zeropoints'],
                        parameters['gain_values'], 
                        parameters['path_to_images'],
                        parameters['path_to_results'],
                        parameters['image_name'], cat[ic],
                        imfits_end =parameters['imfits_end'], 
                        rmsfits_end=parameters['rmsfits_end'],
                        SExtractor_command=SExtractor_command)
        
        # Define segmentation map from the original science image.
        segm_maps_old_file = parameters['path_to_results'] +\
            'SciImages/sources_'+cat[ic]+'_' +parameters['detection_band'] +\
                '_segm.fits'
        # Open SExtractor catalogue with the identified sources from the
        # original science image and save the following information
        # about the sources on the detection band, or the chosen main detection band
        
        if parameters['detection_band']=='det':
            ibmain = parameters['bands'].index(parameters['det_combination'][0])
            f_old_hdu = fits.open('%sSciImages/sources_%s_%s_cat.fits'%
                             (parameters['path_to_results'],cat[ic],
                              parameters['det_combination'][0]),
                             ignore_missing_end=True)
        else:
            f_old_hdu = fits.open('%sSciImages/sources_%s_%s_cat.fits'%
                             (parameters['path_to_results'],cat[ic],
                              parameters['detection_band']),
                             ignore_missing_end=True)
            ibmain = 0
                
        id_old_cat = f_old_hdu[1].data['NUMBER']
        m_old_cat = f_old_hdu[1].data['MAG_AUTO']
        f_old_cat = f_old_hdu[1].data['FLUX_ISO']
        xpos_old_cat = f_old_hdu[1].data['X_IMAGE']
        ypos_old_cat = f_old_hdu[1].data['Y_IMAGE']
        f_old_hdu.close()
        
        # Create main overall output files and write headers to save 
        # the detection stats and detection/dropout array for each field
        for ilf in range(nLF):
            curLF = parameters['LF_shape'][ilf]               
            fws = open(parameters['path_to_results']+'Results/'+curLF+
                       '_RecoveredStats_cat'+cat[ic]+'.cat', 'w')
            header = ('  z   # M_in # N(M_in) # St = 0 # St = 1  # st = 2'
                      ' # St = -1 # St = -2 # St = -3 #St = -4 # N_Dropouts'
                      ' # N_Recovered  # Dropout_frac  # Recov_frac\n')
            fws.writelines(header)
            fws.close()

            fws2 = open(parameters['path_to_results']+'Results/'+curLF+
                       '_RecoveredStats_cat'+cat[ic]+'_nout.cat', 'w')
            header = ('  z   #M_in #N_objects #n_out@m_out  #'+
                      ' #'.join('%.1f'%x for x in m_output)+'\n')
            fws2.writelines(header)  
            fws2.close()
            
            if parameters['dropouts'] is True:
                fws3 = open(parameters['path_to_results']+'Results/'+curLF+
                           '_RecoveredStats_cat'+cat[ic]+'_dout.cat', 'w')
                header = ('  z   #M_in #N_objects #d_out@m_out  #'+
                          ' #'.join('%.1f'%x for x in m_output)+'\n')
                fws3.writelines(header)  
                fws3.close()
        
        # Start the mega loops                        
        # Run for all required redshift bins.
        for iz in range(parameters['z_bins']):
            print(datetime.datetime.now())
            redshift = z_total[iz]  # Redshift of the galaxies
            print('****************************')
            print('Redshift z = %.1f'%redshift)     
            # Effective radius in pixels. Equal to 'R_eff' kpc at z = 6 
            # rescaled to the new redshift
            Re = (parameters['R_eff']*(6+1)/(redshift+1.0)) * \
                cosmo.arcsec_per_kpc_proper(redshift).value / \
                    parameters['size_pix'] 
            
            #Set the stamp diameter to 1.8 arcsecs or at least 5 Re   
            stampsize = int(np.round(np.ceil(1.8/parameters['size_pix']/2))*2) #in pixels, and is even number
            if stampsize < int(Re)*5: stampsize=int(Re)*5
            stamp_radius = stampsize/2  # Radius of the galaxy stamp.
        
            #Make galaxy stamps with specified Re, and stamp size, and various
            #ellipticity and inclinations. 
            galaxy_grid, galaxy_grid_n4 = create_stamps(
                parameters['sersic_indices'],stampsize, Re,
                int(parameters['types_galaxies']), parameters['ebins'],
                parameters['ibins'])

            #Run for all LF shapes
            for ilf in range(nLF):
                curLF = parameters['LF_shape'][ilf]   

                #Open the main overall output files to append
                fws = open(parameters['path_to_results']+'Results/'+curLF+
                           '_RecoveredStats_cat'+cat[ic]+'.cat', 'a')
                
                fws2 = open(parameters['path_to_results']+'Results/'+curLF+
                           '_RecoveredStats_cat'+cat[ic]+'_nout.cat', 'a')
               
                if parameters['dropouts'] is True:
                    fws3 = open(parameters['path_to_results']+'Results/'+curLF+
                           '_RecoveredStats_cat'+cat[ic]+'_dout.cat', 'a')
                    
                #Calculate number of galaxies to simulate according to curLF
                ngals, M_iter, refwl = creation_of_galaxy.calphi_marr(
                    M_input, curLF, redshift, 
                    msample = parameters['n_iterations'], 
                    minngal = parameters['n_galaxies'], 
                    slope = parameters['lin_slope'], 
                    expbase = parameters['exp_base'])
                
                #Note that refwl is only returned for Schechter function only
                #(specified in Schechter param file). For other LF, use the
                #specified wavelength from param.yaml
                if refwl is None: 
                    refwl = parameters['ref_uv_wl'] 
                    
                
                print('%s LF'%curLF)
                print('Number of galaxies to be created '
                      'per iteration (injecting maximum %d galaxies each time):'%
                      (parameters['n_inject_max']))
                print(ngals)
                print('\n')
                Nin[:,iz,ilf] = ngals
    
                #File to save the recovery info of the simulated galaxies.
                fwg = open(parameters['path_to_results']+'Results/'+curLF+
                           'RecoveredGalaxies_'+cat[ic]+'_z%.2f'%redshift+
                           '.cat', 'w')
                headline = ('#Mbin, Minput, niter, round_number, id, beta, '
                            'input %s mag, '%(parameters['bands'][ibmain]) +
                            'detection, dropout, mag_auto and magerr_auto pairs, ' 
                            'mag_iso and sn_iso pairs, mag_aper1 and sn_aper1 pairs, '
                            'mag_aper2 and sn_aper2 pairs,'
                            'mag_aper3 and sn_aper3 pairs in the following bands '+
                            ''.join([' %s'%x for x in parameters['bands']]) + 
                            ', Kron radius (pixels), r90, xpos, ypos \n')
                fwg.writelines(headline)    
                fwg.close()        
                
                # Empty arrays to save the stats.
                total = np.zeros(nmagbins_in)
                identified = np.zeros(nmagbins_in)
                blended_f1 = np.zeros(nmagbins_in)
                blended_f2 = np.zeros(nmagbins_in)
                blended_b = np.zeros(nmagbins_in)
                blended_l = np.zeros(nmagbins_in)
                not_indentified_sn = np.zeros(nmagbins_in)
                not_indentified     = np.zeros(nmagbins_in)
                drops = np.zeros(nmagbins_in)
                identified_total = np.zeros(nmagbins_in)
                identified_fraction = np.zeros(nmagbins_in)
                dropouts_fraction = np.zeros(nmagbins_in)
                
                nout = np.zeros([nmagbins_out,nmagbins_in])
                dout = np.zeros([nmagbins_out,nmagbins_in]) 
                
                # Run for all the required magnitudes.            
                for iM in range(nmagbins_in):
                    print('redshift %d/%d, LF %d/%d, M_input %d/%d'%
                          (iz+1,nzbins,ilf+1,nLF,iM+1,nmagbins_in))
                    print('   Input Magnitude bin M = %.1f mag:'%M_input[iM])
                    print('       simulate galaxies with M =',M_iter[iM])
                    for iit in range(parameters['n_iterations']):
                        niter = iit + 1  # Iteration number

                        Magnitude = M_iter[iM][iit] #absolute intrinsic mag for the current iteration
                        
                        # All galaxies in each iteration will have the same spectrum
                        # Assign a random beta 
                        m = np.zeros(parameters['n_bands'])
                        flux = np.zeros(parameters['n_bands'])
                        beta = random.gauss(parameters['beta_mean'],
                                            parameters['beta_sd'])
                        # Create spectrum
                        creation_of_galaxy.write_spectrum(Magnitude, beta, 
                                         redshift, parameters['path_to_results'], 
                                         absmag=True, refwl = refwl)
                        
                        # Calculate input (theoretical observed) AB mag in each band
                        for ib in range(parameters['n_bands']):      
                            m[ib] = creation_of_galaxy.mag_band(parameters['bands'][ib],
                                                                redshift,parameters['path_to_results']) 
                            flux[ib] = (10**((parameters['zeropoints'][ib] -
                                        m[ib])/2.5))
                            
                        print('   iteration %d/%d:'%(iit+1,parameters['n_iterations']))
                        print('       theoretical m =',m, end='')
                        
                        # Calculate how many rounds of injections for this 
                        # iteration. Inject total of ngals[iM] galaxies with
                        # maximum of parameters['n_inject_max'] galaxies per injection
                        nrounds = np.ceil(ngals[iM]/parameters['n_inject_max']).astype(int)
                        totinject = 0 #tally
                        
                        # Start injecting
                        for ir in range(nrounds):
                            if ir == nrounds-1: 
                                ninject = ngals[iM]-totinject
                            else:
                                ninject = parameters['n_inject_max']   
                                
                            totinject += ninject
                            
                            # Assign (xpos, ypos) = the position of the center of the
                            # galaxy in the image in pixels.        
                            xpos, ypos = creation_of_galaxy.galaxies_positions(
                                    name_hdu_list1+parameters['detection_band']+parameters['imfits_end'], 
                                    ninject, stampsize, Re)
                            
                            # Assign a random beta, incl, and elip.
                            i_rand = np.random.randint(0, parameters['ibins'],
                                                       size=ninject)
                            e_rand = np.random.randint(0, parameters['ebins'],
                                                       size=ninject)
                            # Place galaxies and make fits files in all bands.
                            for ib in range(parameters['n_bands']): 
                                create_synthetic_fits(name_hdu_list1,
                                                      parameters['bands'][ib],
                                                      parameters['imfits_end'],
                                                      m[ib],
                                                      parameters['sersic_indices'],
                                                      ninject,
                                                      parameters['fraction_type'],
                                                      e_rand,i_rand,flux[ib],
                                                      xpos,ypos,stamp_radius,
                                                      galaxy_grid,galaxy_grid_n4,
                                                      parameters['path_to_results'],
                                                      cat[ic],fixed_psf = parameters['fixed_psf'])
                            # Check if there is a separate detection band,
                            # then create the detection image
                            if parameters['detection_band']=='det':                                                                     
                                coadd_images(parameters['path_to_results'],
                                             cat[ic],
                                             parameters['det_combination'],
                                             parameters['detection_band'],
                                             parameters['path_to_images'], 
                                             parameters['image_name'],
                                             parameters['rmsfits_end'],
                                             parameters['coadd_type'])
                                
                            # Run SExtractor on the new images and generate
                            # catalogues of the detected sources.
                            for ib in range(parameters['n_bands']): 
                                do_segmentation = True if (ib == 0) else False
                                run_sextractor.main(parameters['bands'][ib],
                                                    parameters['detection_band'],
                                                    parameters['zeropoints'][ib],
                                                    parameters['gain_values'][ib],
                                                    parameters['path_to_images'],
                                                    parameters['path_to_results'],
                                                    niter, parameters['image_name'],
                                                    cat[ic], M_input[iM], redshift,
                                                    rmsfits_end = parameters['rmsfits_end'],
                                                    roundnum = ir, 
                                                    do_segmentation=do_segmentation,
                                                    SExtractor_command=SExtractor_command)
                                print('.',end='')
                            print('.')
                            
                            #remove the newly created images to save space (except detection images)
                            if minimal_file is True:
                                for ib in range(parameters['n_bands']): 
                                    os.remove(parameters['path_to_results'] + 
                                        'Results/images/sersic_sources_' +
                                        '%s_%s.fits'%(cat[ic],parameters['bands'][ib]))
                            
                            # Find the number of sources for the different categories
                            # of detection status.
                            fwg = open(parameters['path_to_results']+'Results/'+curLF+
                                       'RecoveredGalaxies_'+cat[ic]+'_z%.2f'%redshift+
                                       '.cat', 'a')
                            identified_aux, blended_f1_aux, blended_f2_aux,\
                                blended_b_aux, blended_l_aux,\
                                not_indentified_sn_aux, not_indentified_aux,\
                                nout_c, nout_d, id_nmbr = blending.main(
                                    parameters['path_to_results'],
                                    niter, ir, parameters['detection_band'],
                                    cat[ic], M_input[iM], Magnitude, beta, redshift,
                                    xpos, ypos, xpos_old_cat, ypos_old_cat, 
                                    segm_maps_old_file, m_old_cat, f_old_cat,
                                    id_old_cat, m[ibmain],
                                    parameters['zeropoints'], parameters['bands'],
                                    parameters['min_sn'], parameters['dropouts'],
                                    margin, m_output, 
                                    fwg, parameters['det_combination'])
                            fwg.close()
                            # Find the input magnitude bin of m_in and add the
                            # number of sources in each detection status for each bin.
                            # This is to combine the results from different iterations.
                            for v in range(nmagbins_in):
                                # Check where the calculated input magnitude lies and
                                # add them to the respective counter. 
                                if ((Magnitude > (M_input[v]- 0.5 * mbinsize)) and
                                    (Magnitude <= (M_input[v] + 0.5 * mbinsize))):
                                    if v != iM: warnings.warn("Something is wrong, v and iM should be equal")
                                    total[v] = total[v] + ninject 
                                    identified[v] = identified[v] + identified_aux
                                    blended_f1[v] = blended_f1[v] + blended_f1_aux
                                    blended_f2[v] = blended_f2[v] + blended_f2_aux
                                    blended_b[v] = blended_b[v] + blended_b_aux
                                    blended_l[v] = blended_l[v] + blended_l_aux
                                    not_indentified_sn[v] = (not_indentified_sn[v] +
                                                             not_indentified_sn_aux)
        
                                    not_indentified[v] = (not_indentified[v] +
                                                          not_indentified_aux)
                                    drops[v] = drops[v] + np.sum(nout_d)
                                    nout[:,v] = nout[:,v] + nout_c
                                    dout[:,v] = dout[:,v] + nout_d
                            print('      round %d/%d, inject %d, detect %d, dropout %d (in magbins)'%
                                  (ir+1,nrounds, ninject, np.sum(nout_c),np.sum(nout_d)))
                            
                            #delete big files
                            if minimal_file is True:
                                os.remove(parameters['path_to_results'] + 
                                    'Results/SegmentationMaps/Segmentation_maps_i' +
                                    '%d.%d'%(niter,ir)+ '_' + 
                                    parameters['detection_band'] + '.fits') 
                            # Only keep sextracted catalog of injected galaxies. 
                            # Delete real sources that are originally there
                            for ib in range(parameters['n_bands']):
                                delete_id_fits('%sResults/Catalogs/source_%s_mag%.1f_z%.2f_i%d.%d_%s_cat.fits'%
                                    (parameters['path_to_results'],cat[ic],
                                     M_input[iM],redshift,niter,ir,
                                     parameters['bands'][ib]), 
                                    id_nmbr, id_keep = True)
                        #finish looping over rounds
                        if totinject != ngals[iM]:
                            print('Number of injected galaxies is not equal to what expected')
                            breakpoint()
                    #finish looping over iterations
                #finish looping over input Magnitude bins
                
                #write results to the overall output files
                for s in range(nmagbins_in):
                     # Calculate fractions of different statuses.
                    identified_total[s] = (np.array(identified[s]) + 
                                           np.array(blended_f1[s]) +
                                           np.array(blended_f2[s]))
                    identified_fraction[s] = float(identified_total[s]) /\
                        np.array(total[s])
                    dropouts_fraction[s] = float(drops[s])/np.array(total[s])
                    
                    total_completeness[s, iz, ilf] = identified_fraction[s]
                    total_dropouts[s, iz , ilf] = dropouts_fraction[s]
                    
                    line = ('%.2f\t%.2f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d'
                            '\t%d\t%f\t%f\n'%
                        (redshift,M_input[s],total[s],identified[s],
                         blended_f1[s],blended_f2[s],blended_b[s],
                         blended_l[s],not_indentified_sn[s],not_indentified[s],
                         drops[s],identified_total[s],dropouts_fraction[s],
                         identified_fraction[s]))
                    fws.writelines(line) 
                    
                    line = ('%.2f\t%.2f\t%d'%(redshift, M_input[s],total[s]) +
                            ''.join('\t%d'%x for x in nout[:,s])+'\n')
                    fws2.writelines(line)
                          
                    if parameters['dropouts'] is True:
                        line = ('%.2f\t%.2f\t%d'%(redshift, M_input[s],total[s])+
                            ''.join('\t%d'%x for x in dout[:,s])+'\n')
                        fws3.writelines(line)
                    
                fws.close()
                fws2.close()
                if parameters['dropouts'] is True: fws3.close()
                
                nout_completeness[:,:,iz,ilf] = nout
                nout_dropouts[:,:,iz,ilf] = dout
            #finish looping over LF            
        #finish looping over redshift   
        #save overall output results with pickles  
        output_dict = {'Nin':Nin,'Nout_detect':nout_completeness,
                       'Nout_dropout':nout_dropouts,'magout':m_output,
                       'Magin':M_input,'redshift':z_total,
                       'LF_shape':parameters['LF_shape'],
                       'Niteration':parameters['n_iterations']} 
        outpicklefile = open(parameters['path_to_results']+
                         'Results/GLACiAR_output_'+cat[ic]+'.pickle','wb')
        pickle.dump(output_dict, outpicklefile)
        outpicklefile.close()
        
        #Generate plots (Completeness as a function of input magnitude)
        plot_completeness.main(parameters['path_to_results'],
                                parameters['LF_shape'],
                                M_input, z_total, cat[ic], total_completeness, 
                                total_dropouts)
    #finish looping over field
if __name__ == "__main__":
    main(minimal_file=minimal_file,SExtractor_command=SExtractor_command)
