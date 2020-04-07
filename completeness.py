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
# Read input parameters from the file 'parameters.yaml'.
# Assign a default value when possible if they are not given by the user.

stream = open('parameters.yaml', 'r')
parameters = yaml.load(stream,Loader=yaml.FullLoader)
if parameters['n_galaxies'] is None:
    parameters['n_galaxies'] = 100
if parameters['n_iterations'] is None:
    parameters['n_iterations'] = 100
if parameters['mag_bins'] is None:
    parameters['mag_bins'] = 20
if parameters['min_mag'] is None:
    parameters['min_mag'] = 24.1
if parameters['max_mag'] is None:
    parameters['max_mag'] = 27.9
if parameters['z_bins'] is None:
    parameters['z_bins'] = 16
if parameters['min_z'] is None:
    parameters['min_z'] = 7.0
if parameters['max_z'] is None:
    parameters['max_z'] = 9.2
if parameters['ref_uv_wavelength'] is None:
    raise ValueError('Value of lambda required.')
if parameters['n_bands'] is None:
    raise ValueError('Number of bands needed.')
if parameters['detection_band'] is None:
    raise ValueError('Input detection band.')
if parameters['bands'] is None:
    raise ValueError('Input name of the bands.')
if parameters['list_of_fields'] is None:
    raise ValueError('Input file with the fields.')
if parameters['R_eff'] is None:
    parameters['R_eff'] = 1.075 #kpc at redshift 6
if parameters['beta_mean'] is None:
    parameters['beta_mean'] = -2.2
if parameters['beta_sd'] is None:
    parameters['beta_sd'] = 0.4
if parameters['size_pix'] is None:
    parameters['size_pix'] = 0.08
if parameters['margin'] is None:
    parameters['margin'] = 0.3 #arcsecond. Margin for injection/recovery
if parameters['types_galaxies'] is None:
    parameters['types_galaxies'] = 2
if parameters['ibins'] is None:
    parameters['ibins'] = 10
if parameters['ebins'] is None:
    parameters['ebins'] = 10
if parameters['path_to_images'] is None:
    raise ValueError('Input the directory path.')
if parameters['path_to_results'] is None:
    raise ValueError('Output the directory path.')
if parameters['image_name'] is None:
    raise ValueError('Input the name of the images.')
if parameters['imfits_end'] is None:
    parameters['imfits_end'] = '_v1_drz.fits'
if parameters['rmsfits_end'] is None:
    parameters['rmsfits_end'] = '_v1_rms_scl.fits'
if parameters['sersic_indices'][1] is None:
    parameters['sersic_indices'][1] = 4
if parameters['fraction_type_galaxies'] is None:
    parameters['fraction_type_galaxies'] = [0.5, 0.5]
if (parameters['zeropoints'] is None or
        len(parameters['zeropoints']) < parameters['n_bands']):
    parameters['zeropoints'] = np.zeros(int(parameters['n_bands']))+25
if (parameters['gain_values'] is None or
        len(parameters['gain_values']) < parameters['n_bands']):
    raise ValueError('Input gain values for each band.')
if parameters['dropouts'] is None:
    parameters['dropouts'] = False
if parameters['de_Vacouleur'] is None:
    parameters['de_Vacouleur'] = False
if parameters['LF_shape'] is None:
    parameters['LF_shape'] = ['Schechter','linear','flat']
if parameters['extended_mag_bins_low'] is None:
    parameters['extended_mag_bins_low'] = 1 #should look at the pm_out and see if it's appropriate
if parameters['extended_mag_bins_high'] is None:
    parameters['extended_mag_bins_high'] = 1 
if parameters['lin_slope'] is None:
    parameters['lin_slope'] = 2
if parameters['exp_base'] is None:
    parameters['lin_slope'] = 2
def create_stamps(n0, size_galaxy0, Re0, types_galaxies0, ebins0, ibins0):
    """
    Creates a galaxy following a Sersic profile.

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
                            size_galaxy0))
    galaxy_grid_n4 = np.zeros((size_galaxy0, size_galaxy0))
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
            print('index eccentr incl')
            for j in range(ebins0):
                for k in range(ibins0):
                    print(n0[i], evalues[j], ivalues[k])
                    galaxy_grid[i, j, k, :, :] = creation_of_galaxy.makeSersic(
                                                  n0[i], bn, Re0, evalues[j],
                                                  ivalues[k], size_galaxy0)
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
            frame0[int(x0[gn]-s0):int(x0[gn]+s0),
                   int(y0[gn]-s0):int(y0[gn]+s0)] = gconv
            gn = gn + 1
    return frame0


def main(minimal_file=True):
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
    maxlumdis = cosmo.luminosity_distance(parameters['max_z']).to(u.parsec).value
    minlumdis = cosmo.luminosity_distance(parameters['min_z']).to(u.parsec).value
    M_input_min = np.around(parameters['min_mag']-2.5*
                            np.log10((maxlumdis/10.)**2/(1.+parameters['min_z'])),
                            decimals=1)
    
    M_input_max = np.around(parameters['max_mag']-2.5*
                            np.log10((minlumdis/10.)**2/\
                                     (1.+parameters['min_z'])),decimals=1)
    nmagbins_in = (np.ceil((M_input_max-M_input_min)/mbinsize)+\
        parameters['extended_mag_bins_low']+
        parameters['extended_mag_bins_high']).astype(int)

    M_input = (np.arange(nmagbins_in)-parameters['extended_mag_bins_low'])*\
        mbinsize+M_input_min
            
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
    if not os.path.exists(parameters['path_to_results']+'/Results/Dropouts'):
        os.makedirs(parameters['path_to_results']+'/Results/Dropouts')
        
    #copy parameter and Sextractor files to the Results folder
    copyfile('parameters.yaml',parameters['path_to_results']+'parameters.yaml')
    copyfile('SExtractor_files/parameters.sex',
             parameters['path_to_results']+'parameters.sex')
    #copy Schechter parameter if LF shape is specified as Schechter    
    for ilf in range(nLF):
        parameters['LF_shape'][ilf] = parameters['LF_shape'][ilf].lower()
        if parameters['LF_shape'][ilf] == 'schechter': 
            copyfile('Files/LF_Schechter_params.txt', 
                     parameters['path_to_results']+'LF_Schechter_params.txt')        
    #loop for each field
    for ic in range(len(cat)):    
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('running field '+cat[ic])
        print('\n')
        
        # Empty array to save the simulation results (for pickles)
        nout_completeness = np.zeros((nmagbins_out,nmagbins_in,nzbins,nLF))
        nout_dropouts = np.zeros((nmagbins_out,nmagbins_in,nzbins,nLF))
        Nin = np.zeros((nmagbins_in,nzbins,nLF))
        total_completeness = np.zeros((nmagbins_in,nzbins))
        total_dropouts = np.zeros((nmagbins_in,nzbins))
        
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
                       rmsfits_end=parameters['rmsfits_end'])

        # Open science image.
        obs_data_db, header_db = open_images(name_hdu_list1,
                                             parameters['detection_band'],
                                             parameters['imfits_end'])
        
        # Open segmentation maps from original science image.
        segm_science_old = fits.open(parameters['path_to_results'] +
                                     '/SciImages/sources_'+cat[ic]+'_' +
                                     parameters['detection_band'] +
                                     '_segm.fits', ignore_missing_end=True)
        segm_maps_old = segm_science_old[0].data
        segm_science_old.close()
        # Open SExtractor catalogue with the identified sources from the
        # original science image and save the following information
        # about the sources.
        # id_old_cat = SExtractor ID.
        # m_old_cat = AUTO magnitude.
        # f_old_cat = ISO flux.
        # xpos_old_cat = Centre position in the x axis.
        # ypos_old_cat = Centre position in the y axis.
        f_old_cat = open(parameters['path_to_results']+'/SciImages/sources_' +
                         cat[ic]+'_'+parameters['detection_band']+'.cat', 'r')
        k_oc = f_old_cat.readlines()
        f_old_cat.close()
        id_old_cat = [int(line.split()[0]) for line in k_oc[27:]]
        m_old_cat = [float(line.split()[27]) for line in k_oc[27:]]
        f_old_cat = [float(line.split()[25]) for line in k_oc[27:]]
        xpos_old_cat = [float(line.split()[31]) for line in k_oc[27:]]
        ypos_old_cat = [float(line.split()[32]) for line in k_oc[27:]]
        
        # Create files and write headers to save the stats for each field
        # and LF
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
        
        # Start the real loops                        
        # Run for all the required redshifts.
        for iz in range(parameters['z_bins']):
            redshift = z_total[iz]  # Redshift of the galaxies
                 
            # Effective radius in pixels. Equal to 'R_eff' kpc at z = 6 
            # rescaled to the new redshift
            Re = (parameters['R_eff']*(6+1)/(redshift+1.0)) * \
                cosmo.arcsec_per_kpc_proper(redshift).value / \
                    parameters['size_pix'] 
            
            #Set the stamp diameter to 3.5 arcsecs or at least 10 Re   
            stampsize = int(np.round(3.5/parameters['size_pix'])) #in pixels
            if stampsize < int(Re)*10: stampsize=int(Re)*10
            stamp_radius = stampsize/2  # Radius of the galaxy stamp.
        
            #Make galaxy stamps with specified Re, and stamp size, and various
            #ellipticity and inclinations
            galaxy_grid, galaxy_grid_n4 = create_stamps(
                parameters['sersic_indices'],stampsize, Re,
                int(parameters['types_galaxies']), parameters['ebins'],
                parameters['ibins'])
            
            #Run for all LF shapes
            for ilf in range(nLF):
                curLF = parameters['LF_shape'][ilf]   

                #Open the files to append
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
                    refwl = parameters['ref_uv_wavelength'] 
                    
                
                print('%s LF, redshift z = %.1f'%(curLF, redshift))
                print('Number of galaxies to be created '
                      'per iteration (injecting %d galaxies each time):'%
                      (parameters['n_galaxies']))
                print(ngals)
                print('\n')
                Nin[:,iz,ilf] = ngals
    
                #File to save the recovery info of the simulated galaxies.
                fwg = open(parameters['path_to_results']+'Results/'+curLF+
                           'RecoveredGalaxies_'+cat[ic]+'_z%.1f'%redshift+
                           '.cat', 'w')

                
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
                        
                        # assign a random beta and create spectrum
                        m = np.zeros(parameters['n_bands'])
                        flux = np.zeros(parameters['n_bands'])
                        beta = random.gauss(parameters['beta_mean'],
                                            parameters['beta_sd'])
                        creation_of_galaxy.write_spectrum(parameters['bands'][0], 
                                            Magnitude, beta, redshift, 
                                            absmag=True, refwl = refwl)
                        for ib in range(parameters['n_bands']):      
                            #calculate input (theoretical observed) AB mag in the band
                            m[ib] = creation_of_galaxy.mag_band(parameters['bands'][ib]) 
                            flux[ib] = (10**((parameters['zeropoints'][ib] -
                                        m[ib])/2.5))
                        print('   iteration %d/%d:'%(iit+1,parameters['n_iterations']))
                        print('       theoretical m =',m)
                        #print('       for bands ',parameters['bands'])
                        #calculate how many rounds of injections for this 
                        #iteration. Inject total of ngals[iM] galaxies with
                        #parameters['n_galaxies'] galaxies per injection
                        nrounds = np.ceil(ngals[iM]/parameters['n_galaxies']).astype(int)
                        
                        totinject = 0
                        for ir in range(nrounds):
                            if ir == nrounds-1: 
                                ninject = ngals[iM]-totinject
                            else:
                                ninject = parameters['n_galaxies']   
                                
                            totinject += ninject
                            
                            # (xpos, ypos) is the position of the center of the
                            # galaxy in the image in pixels.
                            xpos, ypos = creation_of_galaxy.galaxies_positions(
                                obs_data_db, ninject, stampsize, Re)
                            
                            # assign a random beta, incl, and elip.
                            i_rand = np.random.randint(0, parameters['ibins'],
                                                       size=ninject)
                            e_rand = np.random.randint(0, parameters['ebins'],
                                                       size=ninject)
                        
                            # Repeat the following process for all the bands.
                            # placing galaxies in the images
                            for ib in range(parameters['n_bands']):      
                                # Empty array of the size of the science image
                                # where the simulated galaxies will be placed.
                                obs, head_obs = open_images(name_hdu_list1,
                                                    parameters['bands'][ib],
                                                    parameters['imfits_end']) 
                                frame = np.zeros(obs.shape)  
                                
                                # if the artificial galaxy is brighter than
                                # 50 mangitudes (some flux is expected, but very
                                # conservative margin), place it in the image.
                                if m[ib] < 50:
                                    # Open PSF file to convolve with the galaxy.
                                    hdu_psf = fits.open('Files/psf_' +
                                                        parameters['bands'][ib] + 
                                                        '.fits')
                                    psf = hdu_psf[0].data
                                    hdu_psf.close()
                            
                                    frame_gal = place_gal(parameters['sersic_indices'],
                                                          ninject,parameters[
                                                          'fraction_type_galaxies'],
                                                          e_rand, i_rand, flux[ib], frame,
                                                          xpos, ypos, stamp_radius, galaxy_grid,
                                                          galaxy_grid_n4, psf)
                                # if the galaxy is fainter than 50 magnitudes
                                # (undetectable), save the original image.
                                else:
                                    frame_gal = frame
                                    
                                # observed image + simulated galaxies
                                new_data2 = obs + frame_gal
                                outfile_conv = parameters['path_to_results']+\
                                    'Results/images/sersic_sources_'+ cat[ic]+\
                                    parameters['bands'][ib]+'.fits'
                                    
                                # Save the new fits file with the simulated galaxies.
                                fits.writeto(outfile_conv, new_data2,
                                             head_obs, overwrite=True)
                                
                                # Run SExtractor on the new images. This generates
                                # catalogues with information about the
                                # detected sources.
                                run_sextractor.main(parameters['bands'][ib],
                                                    parameters['detection_band'],
                                                    parameters['zeropoints'][ib],
                                                    parameters['gain_values'][ib],
                                                    parameters['path_to_images'],
                                                    parameters['path_to_results'],
                                                    niter, parameters['image_name'],
                                                    cat[ic], M_input[iM], redshift,
                                                    rmsfits_end = parameters['rmsfits_end'],
                                                    roundnum = ir)
                            # finished looping over bands (created galaxies/sextracted)
                            
                            # Find the number of sources for the different categories
                            # of detection statuses.
                            identified_aux, blended_f1_aux, blended_f2_aux,\
                                blended_b_aux, blended_l_aux,\
                                not_indentified_sn_aux, not_indentified_aux,\
                                nout_c, nout_d = blending.main(parameters['path_to_results'],
                                                          niter, ir, 
                                                          parameters['detection_band'],
                                                          cat[ic], M_input[iM], 
                                                          Magnitude, redshift,
                                                          xpos, ypos, xpos_old_cat,
                                                          ypos_old_cat, segm_maps_old,
                                                          m_old_cat, f_old_cat,
                                                          id_old_cat, m[0],
                                                          parameters['zeropoints'],
                                                          parameters['bands'],
                                                          parameters['min_sn'],
                                                          parameters['dropouts'],
                                                          margin,
                                                          m_output, fwg)

                            # Find the input magnitude bin of m_in and add the
                            # number of sources in each detection status for each bin.
                            # This is to combine the results from different iterations.
                            for v in range(nmagbins_in):
                                # Check where the calculated input magnitude lies and
                                # add them to the respective counter. 
                                if ((Magnitude > (M_input[v]- 0.5 * mbinsize)) and
                                    (Magnitude <= (M_input[v] + 0.5 * mbinsize))):
                                    if v != iM: warnings.warn("Hmmmmmmmmm????????????")
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
                            print('      round %d/%d, inject %d, detect %d, dropout %d'%
                                  (ir+1,nrounds, ninject, np.sum(nout_c),np.sum(nout_d)))
                            #pdb.set_trace()
                            if minimal_file is True:
                                #pdb.set_trace()
                                os.remove(parameters['path_to_results'] + 
                                    'Results/SegmentationMaps/Segmentation_maps_i' +
                                    '%d.%d'%(niter,ir)+ '_' + 
                                    parameters['detection_band'] + '.fits')                   
                                for ib in range(parameters['n_bands']):
                                    os.remove(parameters['path_to_results'] +
                                        'Results/Dropouts/source_' +cat[ic] + 
                                        '_mag%.1f'%M_input[iM] + 
                                        '_z%.1f'%redshift +'_i' + 
                                        str(niter) + '.' + str(ir) + '_' + 
                                        parameters['bands'][ib] + '.cat')
                                                                          
                        #finished looping over rounds
                        if totinject != ngals[iM]: pdb.set_trace()
                    #finished looping over iterations
                #finished looping over input Magnitude bins
                fwg.close()
                
                for s in range(nmagbins_in):
                     # Calculate fractions of different statuses.
                    identified_total[s] = (np.array(identified[s]) + 
                                           np.array(blended_f1[s]) +
                                           np.array(blended_f2[s]))
                    identified_fraction[s] = float(identified_total[s]) /\
                        np.array(total[s])
                    dropouts_fraction[s] = float(drops[s])/np.array(total[s])
                    
                    total_completeness[s, iz] = identified_fraction[s]
                    total_dropouts[s, iz] = dropouts_fraction[s]
                    
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
            #finished looping over LF            
        #finished looping over redshift   
        #saving results with pickles  
        output_dict = {'Nin':Nin,'Nout_detect':nout_completeness,
                       'Nout_dropout':nout_dropouts,'magout':m_output,
                       'Magin':M_input,'redshift':z_total,
                       'LF_shape':parameters['LF_shape'],
                       'Niteration':parameters['n_iterations']} 
        outpicklefile = open(parameters['path_to_results']+
                         'Results/GLACiAR_output_'+cat[ic]+'.pickle','wb')
        pickle.dump(output_dict, outpicklefile)
        outpicklefile.close()
        
        # Generate plots.
        plot_completeness.main(parameters['path_to_results'],
                               parameters['LF_shape'],
                               M_input, parameters['min_z'],
                               parameters['max_z'], parameters['z_bins'],
                               cat[ic], total_completeness, 
                               total_dropouts)

    #finished looping over field
if __name__ == "__main__":
    main()
