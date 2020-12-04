import numpy as np
import dropouts
from astropy.io import fits
import pdb
def main(path_to_results, niter, roundnum, detection_band, cat, Mbin, Min, 
         beta, redshift, xpos, ypos, xpos_oc, ypos_oc, segm_science_fits, m_oc, f_oc, 
         id_oc, input_mag, zp, bands, min_sn, dp, margin, magbins, 
         fwg,det_combination):
    """
    Uses the information from the new and old catalogs and segmentation
    maps to find the new sources and label them according to their
    identification and blending statuses.
    Args:
        path_to_results (string) = Path to the folder with the results images.
                               Given in the parameters file.
        niter (integer) = iteration number.
        detection_band (string) = Name of the detection band given in the
                                  parameters file.
        cat (string) = Name of the field for which the simulation is run.
        Mbin (float) = Input absolute magnitude bin
        Min (float) = Initial input absolute magnitude for the simulated galaxy 
                     (at UV wavelenght).
        beta (float) = Initial input beta for the simulatted galaxy.        
        redshift (float) = Redshift for the simulated galaxy.
        xpos (int array) = Position of the simulated galaxy in the x axis.
        ypos (int array) = Position of the simulated galaxy in the y axis.
        xpos_oc (float array) = Array with the position in the x axis of the
                                centre for all the sources identified in the
                                original science image.
        ypos_oc (float array) = Array with the position in the y axis of the
                                centre for all the sources identified in the
                                original science image.
        segm_science_fits (str) = fits file for segmentation map produced by SExtractor for the
                               identified sources from the science image.
        m_oc (float array) = Array with the AB AUTO magnitude for the sources
                             identified by SExtractor in the science image.
        f_oc (float array) = Array with the ISO flux for the sources
                             identified by SExtractor in the science image.
        id_oc (int array) = Array with the ID assigned by SExtractor for each
                            one of the sources identified in the science image.
        input_mag (float) = Expected apparent magnitude of the artificial 
                            galaxies in the detection band.
        zp (float array) = Zeropoint values for all bands as given in the
                           input parameters file.
        bands (string array) = Name of the bands in which the artificial
                               galaxies will be simulated. Given in the input
                               parameters file.
        min_sn (float) = Minimum S/N ratio in the detection band for an object
                         to be considered detected by SExtractor. Given in the
                         input parameters file.
        dp (boolean) = Boolean that indicates whether the user requires to run
                       a dropout selection. Given in the input parameters file.
                       If True, the dropouts.py module will be used.
        droptype (string) = type of dropout. Choices are 'borg', 'hlf', 'candels'
        margin (int) = Number of pixels from centre in which search is performed.
        magbins (float array) = Array of the output magnitude bins.
                         The sizes of the bins must be equal. This is used to
                         report nmout
        fwg (text file) = File in which the information about the artificial
                         sources will be saved. 'RecoveredGalaxies_cat_z#.cat'

    Returns:
        identified (int) = Number of artificial galaxies from the current
                           iteration that are detected by SExtractor and that
                           are isolated.
        blended_b (int) = Number of artificial galaxies from the current
                          iteration that are detected by SExtractor and are
                          blended with previously detected brighter sources.
        blended_f (int) = Number of artificial galaxies from the current
                          iteration that are detected by SExtractor and are
                          blended with previously detected fainter sources.
        not_indentified_sn (int) = Number of artificial galaxies from the
                                   current iteration that are detected by
                                   SExtractor but are considered not identified
                                   because their S/N is below min_sn.
        not_indentified (int) = Number of artificial galaxies from the current
                                iteration that are not detected by SExtractor.
        drops (int) = Number of artificial galaxies from the current iteration
                      that passed the redshift selection criteria from
                      'dropouts.py'. If drops is set to False the value is 0.
        nmout (int arr) = array of size magbins number of galaxies
                              observed in the magbins at the input_mag                            
    """
    #open segmentation map from science image
    segm_science_old = fits.open(segm_science_fits, ignore_missing_end=True)
    segm_science = segm_science_old[0].data
    segm_science_old.close()
    
    # Open segmentation maps from simulated images, save the data
    segm_new_cat = fits.open('%sResults/SegmentationMaps/Segmentation_maps_i%d.%d_%s.fits'%(
            path_to_results,niter,roundnum,detection_band) , 
            ignore_missing_end=True)
        
    segm_sim = segm_new_cat[0].data  # Data from new images.
    segm_new_cat.close()

    # Main catalog with the identified sources from the simulated images.
    if detection_band =='det':
        selband = det_combination[0]
        locband =  bands.index(selband)
    else:
        selband = detection_band
        locband =  bands.index(selband)
        if locband != 0:breakpoint()
        
    # Information from SExtractor for the new sources (science image +
    # simulated galaxies).
    f_hdu = fits.open('%sResults/Catalogs/source_%s_mag%.1f_z%.2f_i%d.%d_%s'
                      '_cat.fits'%(path_to_results,cat,Mbin,redshift,
                                   niter,roundnum,selband),
                      ignore_missing_end=True)
                
    id_mgal = f_hdu[1].data['NUMBER']  # ID
    fauto_gal = f_hdu[1].data['FLUX_AUTO'] # auto flux
    m_gal = f_hdu[1].data['MAG_AUTO']  # AUTO magnitude
    f_gal = f_hdu[1].data['FLUX_AUTO']  # AUTO flux
    ef_gal = f_hdu[1].data['FLUXERR_AUTO'] # RMS error for flux
    sn_mgal = np.divide(f_gal, ef_gal, out=np.zeros_like(f_gal), where=ef_gal!=0) # RMS error for mag
    #radius50 = f_hdu[1].data['FLUX_RADIUS'][:,0]  # Radii
    kron_radius = f_hdu[1].data['KRON_RADIUS']*f_hdu[1].data['B_IMAGE'] #pixels
    radius90 = f_hdu[1].data['FLUX_RADIUS'][:,1]  # Radii
    f_hdu.close()
    
    # Convert the previous arrays into np.arrays
    xpos = np.array(xpos).astype(int)
    ypos = np.array(ypos).astype(int)

    # Array of the status code, initialized with -999 
    status = np.zeros(len(xpos),dtype=int) -99

    # Array of the ID of the artificial galaxy, initialized with -999 
    id_nmbr = np.zeros(len(xpos),dtype=int) - 99
  
    # Open a text file (.reg) with the formatting to be used as region file on
    # DS9. It shows location where the artificial source was originally placed.
    g = open(path_to_results + 'Results/SegmentationMaps/region_%s_mag%.1f_z%.2f_i'
             '%d.%d.reg'%(cat,Mbin,redshift,niter,roundnum),'w')

    # Loop for each artificial source.
    for i in range(len(xpos)):
        # The IF below searches for any value diferent than zero in the
        # segmentation maps of the new images (science + artificial sources)
        # within an square centered in the input position for the artificial
        # galaxy and a side of 2*margin. This is done by evaluating the sum of
        # the pixels' values. It enters the IF if the value is != 0.        
        if np.sum(segm_sim[xpos[i]-margin:xpos[i]+margin,
                           ypos[i]-margin:ypos[i]+margin]) != 0:
            
            # array with square where search is performed.
            id_mgi_aux2 = segm_sim[xpos[i]-margin:xpos[i]+margin,
                                   ypos[i]-margin:ypos[i]+margin]

            # ID of the source in the search region is recorded.
            # Choose the source that is the closest to the original xpos,ypos
            pix_index = np.nonzero(id_mgi_aux2)
            wmin = np.argmin((pix_index[0] - margin)**2 + (pix_index[1] - margin)**2)
            id_nmbr[i] = id_mgi_aux2[pix_index[0][wmin],pix_index[1][wmin]]
            
            # ID value for pixels from the science image in the position of
            # the newly found source.
            id_mgi2 = segm_science[segm_sim == id_nmbr[i]]
            idsimgal = len(id_mgi2)  # Number of pixels the source encompasses.
            
            # ID of source previously in the position of newly found source.
            id_mgi = np.max(id_mgi2) #integer 
            
            # Number of pixels of the simulated galaxy that are overlapped
            # with the old source.
            w2 = np.nonzero(id_mgi2)[0]
                            
            # If the S/N of the newly identified source is above the required,
            # threshold, enter the IF below.
            if sn_mgal[id_mgal == id_nmbr[i]] >= min_sn:
                # If there wasn't any source on the pixels where the new source
                # was found, enter the IF below.
                # These are objects detected and isolated.
                if np.sum(id_mgi2) == 0:
                    status[i] = 0  # Status for detected and isolated sources.
                    # Write the region file for the artificial sources with
                    # status=0 in colour green. The position is the input one.
                    g.write("circle %s %s 11 #color=green width=2 text={%s,%s}\n" %
                            (ypos[i], xpos[i], id_nmbr[i], 
                             m_gal[id_mgal == id_nmbr[i]]))
                # If there was a source previously on the pixels where the new
                # source was found, enter the IF below.
                # These are objects detected and blended.
                else:
                    id_blended = id_mgi  # ID of the source previously there.
                    # Find the magnitude of the old source and compare it to
                    # the input magnitude. 
                    # If the existing source is brighter, enter the IF below.
                    if m_oc[id_blended == id_oc] <= input_mag:
                        # Enter the IF below If we can recover the flux within
                        # 25% of theoretical input flux 
                        # AND if the number of pixels of the newly
                        # identified source that were previously occupied by
                        # another source is 25% or less the total number of
                        # pixels of said source (the overlap of the new source
                        # is 25% its original size). We consider this as if the
                        # artificial source was blended with a fainter object.
                        if ((np.abs((fauto_gal[id_mgal == id_nmbr[i]]/ \
                                    (10**((zp[0]-input_mag)/2.5)))-1.)  < 0.25)\
                            and (len(w2) <= 0.25*idsimgal)):
                            status[i] = 2 # Status code for blended with fainter 
                            # object but intrinsic mag is brighter
                            # Write the region file for the artificial sources
                            # with status=2 in colour blue.
                            g.write("circle %s %s 11 #color=blue width=2 text={%s,%s}\n" %
                                    (ypos[i], xpos[i], id_nmbr[i], 
                                     m_gal[id_mgal == id_nmbr[i]]))

                        # If the flux or overlap conditions aren't true, do the
                        # following
                        else:
                            status[i] = -1 # Status code for blended with brighter object.
                            # Write the region file for the artificial sources
                            # with status=-1 in colour red.
                            g.write("circle %s %s 11 #color=red width=2 text={%s,%s}\n" %
                                    (ypos[i], xpos[i], id_nmbr[i], 
                                     m_gal[id_mgal == id_nmbr[i]]))

                    # if the magnitude of the old source is fainter than the
                    # input magnitude, enter the IF below.
                    else:
                        # Make sure that it doesn't blend with large object st
                        # it is the same object as that in the science catalog.                        
                        
                        # if the (the overlap of the new source
                        # is >25% its original size). We consider this as if the
                        # artificial source was blended with a large object.
                        if (len(w2) > 0.25*idsimgal):
                            status[i] = -2  #blended with a large object.
                            # Write the region file for the artificial sources
                            # with status=-2 in colour red.
                            g.write("circle %s %s 11 #color=pink width=2 text={%s,%s}\n" %
                                    (ypos[i], xpos[i], id_nmbr[i], 
                                     m_gal[id_mgal == id_nmbr[i]]))
                        
                        else:        
                            # Status code for blended with fainter object.
                            # Different from 2 as it purely compares magnitudes.
                            status[i] = 1 # Status code for blended with fainter object.
                            # Write the region file for the artificial sources
                            # with status=1 in colour blue.
                            g.write("circle %s %s 11 #color=cyan width=2 text={%s,%s}\n" %
                                    (ypos[i], xpos[i], id_nmbr[i], 
                                     m_gal[id_mgal == id_nmbr[i]]))

            # If the S/N of the newly identified source is below the required,
            # threshold, enter the ELSE below.
            else:
                # Status for detected sources with S/N below required min_sn.
                status[i] = -3
                # Write the region file for the artificial sources with
                # status=-3 in colour Magenta.
                g.write("circle %s %s 11 #color=magenta width=2 text={%s,%s}\n" %
                                (ypos[i], xpos[i], id_nmbr[i], 
                                 m_gal[id_mgal == id_nmbr[i]]))

        # If all values of the new segmentation map within the search grid are
        # zero, the object has not been detected by SExtractor.
        else:
            status[i] = -4  # Status for sources not detected by SExtractor.
            # Write the region file for the artificial sources with status=-4
            # in colour red.
            g.write("box %s %s 11 11 0 #color=red width=2\n" %
                    (ypos[i], xpos[i]))          
    # Close the .reg file.
    g.close()

    # Initialise array for ISO magnitudes, AUTO magnitudes, and S/N in all
    # bands measured by SExtractor for the new sources.
    mag_iso = np.zeros((len(id_mgal), len(bands)))
    mag_auto = np.zeros((len(id_mgal), len(bands)))
    mag_aper1 = np.zeros((len(id_mgal), len(bands)))
    mag_aper2 = np.zeros((len(id_mgal), len(bands))) 
    mag_aper3 = np.zeros((len(id_mgal), len(bands)))
    
    magerr_auto = np.zeros((len(id_mgal), len(bands)))
    
    sn_iso = np.zeros((len(id_mgal), len(bands)))
    sn_aper1 = np.zeros((len(id_mgal), len(bands)))
    sn_aper2 = np.zeros((len(id_mgal), len(bands)))
    sn_aper3 = np.zeros((len(id_mgal), len(bands)))
      
    # Loop for the number of bands used in the simulation.
    for j in range(len(bands)):
        # Open the catalog with the identified sources from the simulated
        # images (science + artificial sources) for each band.
    
        f1 = fits.open('%sResults/Catalogs/source_%s_mag%.1f_z%.2f_i%d.%d''_%s_cat.fits'%
                       (path_to_results,cat,Mbin,redshift,niter,
                        roundnum,bands[j]),ignore_missing_end=True)

        # Save the information on the MAG and S/N for each band 
        # to be used for dropout selection.
        mag_iso[:, j] = f1[1].data['MAG_ISO'] 
        mag_auto[:,j] = f1[1].data['MAG_AUTO'] 
        mag_aper1[:,j] = f1[1].data['MAG_APER'][:,0]
        mag_aper2[:,j] = f1[1].data['MAG_APER'][:,1]
        mag_aper3[:,j] = f1[1].data['MAG_APER'][:,2]
        magerr_auto[:,j] = f1[1].data['MAGERR_AUTO'] 
        sn_iso[:,j] = np.divide(f1[1].data['FLUX_ISO'],
                                f1[1].data['FLUXERR_ISO'],
                                out = np.zeros_like(f1[1].data['FLUX_ISO']),
                                where = f1[1].data['FLUXERR_ISO']!=0)
        sn_aper1[:,j] = np.divide(f1[1].data['FLUX_APER'][:,0],
                                f1[1].data['FLUXERR_APER'][:,0],
                                out = np.zeros_like(f1[1].data['FLUX_APER'][:,0]),
                                where = f1[1].data['FLUXERR_APER'][:,0]!=0)
        sn_aper2[:,j] = np.divide(f1[1].data['FLUX_APER'][:,1],
                                f1[1].data['FLUXERR_APER'][:,1],
                                out = np.zeros_like(f1[1].data['FLUX_APER'][:,1]),
                                where = f1[1].data['FLUXERR_APER'][:,1]!=0)  
        sn_aper3[:,j] = np.divide(f1[1].data['FLUX_APER'][:,2],
                                f1[1].data['FLUXERR_APER'][:,2],
                                out = np.zeros_like(f1[1].data['FLUX_APER'][:,2]),
                                where = f1[1].data['FLUXERR_APER'][:,2]!=0) 
        f1.close()
        
    # Save the data only for sources identified as artificial sources.
    mag_iso = mag_iso[id_nmbr.astype(int)-1, :]
    mag_auto = mag_auto[id_nmbr.astype(int)-1, :]
    mag_aper1 = mag_aper1[id_nmbr.astype(int)-1, :]
    mag_aper2 = mag_aper2[id_nmbr.astype(int)-1, :]
    mag_aper3 = mag_aper3[id_nmbr.astype(int)-1, :]
    magerr_auto = magerr_auto[id_nmbr.astype(int)-1, :]
    sn_iso = sn_iso[id_nmbr.astype(int)-1, :]
    sn_aper1 = sn_aper1[id_nmbr.astype(int)-1, :]
    sn_aper2 = sn_aper2[id_nmbr.astype(int)-1, :]
    sn_aper3 = sn_aper3[id_nmbr.astype(int)-1, :]
    kron_radius = kron_radius[id_nmbr.astype(int)-1]
    radius90 = radius90[id_nmbr.astype(int)-1]
    
    wnondetect = np.nonzero(status == 4)[0]
    if len(wnondetect) > 0 :
        mag_iso[wnondetect,:] = -99
        mag_auto[wnondetect,:] = -99
        mag_aper1[wnondetect,:] = -99
        mag_aper2[wnondetect,:] = -99
        mag_aper3[wnondetect,:] = -99
        magerr_auto[wnondetect,:] = -99
        sn_iso[wnondetect,:] = -99
        sn_aper1[wnondetect,:] = -99
        sn_aper2[wnondetect,:] = -99
        sn_aper3[wnondetect,:] = -99
        kron_radius[wnondetect] = -99
        radius90[wnondetect] = -99
    
    #determine dropouts 
    # Run dropout module if dropout parameter is set to True.
    if dp is True:
        drops = dropouts.main(mag_iso, mag_auto, mag_aper1, mag_aper2,
                              sn_iso, sn_aper1, sn_aper2, status, magbins, 
                              redshift)
    else:
        drops = np.zeros(len(id_nmbr))-99
        
    #writing output to fwg file for each simulated galaxies
    
    for i in range(len(xpos)):
         # Line with information of the artificial source.          
        line = ('%.1f\t%.2f\t%d\t%d\t%d\t%.2f\t%.2f\t%d\t%d'
                %(Mbin, Min, niter, roundnum, id_nmbr[i],beta, input_mag,
                  status[i], drops[i])+ 
                ''.join(['\t%.2E\t%.2E'%(x,y) for x,y in zip(mag_auto[i,:],magerr_auto[i,:])])+
                ''.join(['\t%.2E\t%.2E'%(x,y) for x,y in zip(mag_iso[i,:],sn_iso[i,:])]) +
                ''.join(['\t%.2E\t%.2E'%(x,y) for x,y in zip(mag_aper1[i,:],sn_aper1[i,:])]) +
                ''.join(['\t%.2E\t%.2E'%(x,y) for x,y in zip(mag_aper2[i,:],sn_aper2[i,:])]) +
                ''.join(['\t%.2E\t%.2E'%(x,y) for x,y in zip(mag_aper3[i,:],sn_aper3[i,:])]) +
                '\t%.3f\t%.3f\t%.1f\t%.1f\n'%(kron_radius[i],radius90[i],xpos[i],ypos[i]))
                                         
        # Write line in the file 'RecoveredGalaxies_cat_z#.cat'
        fwg.writelines(line)
                    
    mbinsize = magbins[1]-magbins[0]

    nout_c = np.zeros(len(magbins))
    nout_d = np.zeros(len(magbins))
    for i in range(len(magbins)):
        cwin = np.nonzero((mag_auto[:,locband] >= magbins[i]-mbinsize*0.5) &
                         (mag_auto[:,locband] < magbins[i]+mbinsize*0.5) &
                         (status >= 0))[0]
        dwin = np.nonzero((mag_auto[:,locband] >= magbins[i]-mbinsize*0.5) &
                         (mag_auto[:,locband] < magbins[i]+mbinsize*0.5) &
                         (status >= 0) & (drops>0))[0]
        
        nout_c[i] = len(cwin)
        nout_d[i] = len(dwin)

    identified = len(np.nonzero(status == 0)[0])
    blended_f1 = len(np.nonzero(status == 1)[0])
    blended_f2 = len(np.nonzero(status == 2)[0])
    blended_b = len(np.nonzero(status == -1)[0])
    blended_l = len(np.nonzero(status == -2)[0])
    not_indentified_sn = len(np.nonzero(status == -3)[0])
    not_indentified = len(np.nonzero(status == -4)[0])
    
    # Return number of sources for each of the following categories.
    return identified, blended_f1, blended_f2, blended_b, blended_l, \
        not_indentified_sn, not_indentified, nout_c, nout_d, id_nmbr

if __name__ == "__main__":
    main()
