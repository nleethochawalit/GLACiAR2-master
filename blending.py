import numpy as np
import pickle
import dropouts
import os
from astropy.io import fits
import pdb

def main(path_to_results, niter, roundnum, detection_band, cat, Mbin, Min, 
         redshift, xpos, ypos, xpos_oc, ypos_oc, segm_science, m_oc, f_oc, 
         id_oc, input_mag, zp, bands, min_sn, dp, margin, magbins, fwg):
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
        redshift (float) = Redshift for the simulated galaxy.
        xpos (int array) = Position of the simulated galaxy in the x axis.
        ypos (int array) = Position of the simulated galaxy in the y axis.
        xpos_oc (float array) = Array with the position in the x axis of the
                                centre for all the sources identified in the
                                original science image.
        ypos_oc (float array) = Array with the position in the y axis of the
                                centre for all the sources identified in the
                                original science image.
        segm_science (array) = Segmentation map produced by SExtractor for the
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
    # Open segmentation maps from simulated images, save the data
    segm_new_cat = fits.open('%sResults/SegmentationMaps/Segmentation_maps_i%d.%d_%s.fits'%(
            path_to_results,niter,roundnum,detection_band) , ignore_missing_end=True)
    segm_sim = segm_new_cat[0].data  # Data from new images.
    segm_new_cat.close()


    # Catalog with the identified sources from the simulated images.
    f = open('%sResults/Dropouts/source_%s_mag%.1f_z%.1f_i%d.%d_%s'
             '.cat'%(path_to_results,cat,Mbin,redshift,niter,roundnum,detection_band))
    k = f.readlines()
    f.close()

    # Information from SExtractor for the new sources (science image +
    # simulated galaxies).
    id_mgal = [int(line.split()[0]) for line in k[27:]]  # ID
    fauto_gal = [float(line.split()[25]) for line in k[27:]] # auto flux
    f_gal = [float(line.split()[1]) for line in k[27:]]  # Isophotal flux
    ef_gal = [float(line.split()[2]) for line in k[27:]]  # RMS error for flux
    m_gal = [float(line.split()[27]) for line in k[27:]]  # AUTO magnitude
    merr_gal = [float(line.split()[28]) for line in k[27:]]  # AUTO magnitude
    sn_mgal = np.array(np.array(f_gal)/np.array(ef_gal))  # RMS error for mag
    xpos_nc = [float(line.split()[32]) for line in k[27:]]  # Position in x
    ypos_nc = [float(line.split()[31]) for line in k[27:]]  # Position in y
    #radius = [float(line.split()[42]) for line in k[27:]]  # Radii

    # Convert the previous arrays into np.arrays
    xpos = np.array(xpos).astype(int)
    ypos = np.array(ypos).astype(int)
    xpos_oc = np.array(xpos_oc)
    ypos_oc = np.array(ypos_oc)
    xpos_nc = np.array(xpos_nc)
    ypos_nc = np.array(ypos_nc)
    m_oc = np.array(m_oc)
    f_oc = np.array(f_oc)
    id_oc = np.array(id_oc)
    id_mgal = np.array(id_mgal)
    m_gal = np.array(m_gal)
    merr_gal = np.array(merr_gal)
    fauto_gal = np.array(fauto_gal)
    

    # Array of the status code, initialized with -999 
    status = np.zeros(len(xpos),dtype=int) -999

    # Array of the ID of the artificial galaxy, initialized with -999 
    id_nmbr = np.zeros(len(xpos),dtype=int) - 999
  
    # Open a text file (.reg) with the formatting to be used as region file on
    # DS9. It shows location where the artificial source was originally placed.
    g = open(path_to_results + 'Results/SegmentationMaps/region_%s_mag%.1f_z%.1f_i'
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
                    # If the source is brighter, enter the IF below.
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
                            status[i] = 2 # Status code for blended with fainter object but intrinsic mag is birghter
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
                        if ((len(w2) > 0.25*idsimgal) or 
                            (np.abs((fauto_gal[id_mgal == id_nmbr[i]]/ \
                                    (10**((zp[0]-input_mag)/2.5)))-1.) > 0.25)):
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
                            g.write("circle %s %s 11 #color=blue width=2 text={%s,%s}\n" %
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
    magerr_auto = np.zeros((len(id_mgal), len(bands)))
    sn = np.zeros((len(id_mgal), len(bands)))
    # Loop for the number of bands used in the simulation.
    for j in range(len(bands)):
        # Open the catalog with the identified sources from the simulated
        # images (science + artificial sources) for each band.
        f1 = open(path_to_results + 'Results/Dropouts/source_%s_mag%.1f_'
                  'z%.1f_i%d.%d''_%s.cat'%(cat,Mbin,redshift,niter,
                                           roundnum,bands[j]))
        k1 = f1.readlines()
        f1.close()

        # Save the information on the ISO mag, AUTO mag, and S/N for each band.
        mag_iso[:, j] = [float(line.split()[3]) for line in k1[27:]]
        mag_auto[:, j] = [float(line.split()[27]) for line in k1[27:]]
        magerr_auto[:,j] = [float(line.split()[28]) for line in k1[27:]]
        sn[:, j] = [float(line.split()[1])/float(line.split()[2]) for
                    line in k1[27:]]
        
    # Save the data only for sources identified as artificial sources.
    mag_iso = mag_iso[id_nmbr.astype(int)-1, :]
    mag_auto = mag_auto[id_nmbr.astype(int)-1, :]
    magerr_auto = magerr_auto[id_nmbr.astype(int)-1, :]
    sn = sn[id_nmbr.astype(int)-1, :]
    
    wnondetect = np.nonzero(status == 4)[0]
    if len(wnondetect) > 0 :
        mag_iso[wnondetect,:] = -999
        mag_auto[wnondetect,:] = -999
        magerr_auto[wnondetect,:] = -999
        sn[wnondetect,:] = -999
        
    #determine dropouts 
    # Run dropout module if dropout parameter is set to True.
    if dp is True:
        drops = dropouts.main(mag_iso, mag_auto, sn, status, magbins, redshift)
    else:
        drops = np.zeros(len(id_nmbr))-99
        
    #writing output to fwg file for each simulated galaxies
    for i in range(len(xpos)):
         # Line with information of the artificial source.
         if status[i] > -4:
             if m_gal[id_mgal == id_nmbr[i]] != mag_auto[i,0]:
                 print('something is wrong')
                 pdb.set_trace()
                 
         line = ('%.1f\t%.2f\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d'
                 %(Mbin, Min, niter, roundnum, id_nmbr[i],input_mag,
                   mag_auto[i,0],magerr_auto[i,0],status[i], drops[i]) + 
                 ''.join(['\t%.2f\t%.2f'%(x,y) for x,y in zip(mag_iso[i,:],sn[i,:])]) +
                 '\n')
                                          
         # Write line in the file 'RecoveredGalaxies_cat_z#.cat'
         fwg.writelines(line)
                    
    mbinsize = magbins[1]-magbins[0]

    nout_c = np.zeros(len(magbins))
    nout_d = np.zeros(len(magbins))
    for i in range(len(magbins)):
        cwin = np.nonzero((mag_auto[:,0] >= magbins[i]-mbinsize*0.5) &
                         (mag_auto[:,0] < magbins[i]+mbinsize*0.5) &
                         (status >= 0))[0]
        dwin = np.nonzero((mag_auto[:,0] >= magbins[i]-mbinsize*0.5) &
                         (mag_auto[:,0] < magbins[i]+mbinsize*0.5) &
                         (status >= 0) & (drops==1))[0]
        
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
        not_indentified_sn, not_indentified, nout_c, nout_d

if __name__ == "__main__":
    main()
