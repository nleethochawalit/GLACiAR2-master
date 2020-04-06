"""
Writes the configuration files with the parameters to run SExtractor.
"""


def main(name_band, detection_band, zp, g, path_to_im, path_to_results, 
         niter, image_name, cat, rmsfits_end,roundnum = 0):
    """
    Writes and saves the configuration files needed to run SExtractor
    according to the parameters given by the user. It saves new files
    for each band and  they are replaced for each iteration.
    Args:
        name_band (string) = name of the band in which the image is taken.
        detection_band (string) = name of the detection band given in the
                                parameters file.
        zp (float) = Zeropoint value for the respective band. It is given
                         in the parameters file.
        g (float) = Gain value for the respective band. It is given
                         in the parameters file.
        path_to_cat (string) = Path to the folder with the science images.
                               Given in the parameters file.
        niter (integer) = iteration number.
        image_name (string) = Name of the science image (before the band) as
                              given in the parameters file.
        cat (string) = Name of the field for which the simulation is run.
        rmsfits_end (string) = the ending of the rms file names
    """
    # Open the example sextractor file from where we will read the lines.
    f = open('SExtractor_files/parameters.sex', 'r')
    k = f.readlines()
    f.close()
    # k_new will have the same written lines than the default file
    # (lines saved in k), and then the ones that should be changed
    # according to the input parameters, are overwritten. A new file is
    # saved for each band with k_new as the text.
    k_new = k
    if name_band == detection_band:
        for i in range(len(k)):
            if k[i] == 'MAG_ZEROPOINT\n':
                k_new[i] = 'MAG_ZEROPOINT    '+str(zp) + '\n'
            if k[i] == 'GAIN\n':
                k_new[i] = 'GAIN             '+str(g) + '\n'
            if k[i] == 'CHECKIMAGE_TYPE\n':
                k_new[i] = 'CHECKIMAGE_TYPE   SEGMENTATION\n'
            if k[i] == 'CHECKIMAGE_NAME\n':
                k_new[i] = 'CHECKIMAGE_NAME'+\
                ' %sResults/SegmentationMaps/Segmentation_maps_i%d.%d_%s'%(
                        path_to_results,niter,roundnum,name_band)+'.fits \n'
            if k[i] == 'CATALOG_NAME\n':
                k_new[i] = 'CATALOG_NAME test \n'
            if k[i] == 'WEIGHT_IMAGE\n':
                k_new[i] = 'WEIGHT_IMAGE     '+\
                '%s%s%s_%s%s, %s%s%s_%s%s \n'%(path_to_im, image_name, cat, 
                                               detection_band, rmsfits_end, 
                                               path_to_im, image_name, cat, 
                                               detection_band, rmsfits_end)
                
        nf = open('SExtractor_files/parameters_' + name_band + '.sex', 'w')
        nf.writelines(k_new)
        nf.close()
    else:
        for i in range(len(k)):
            if k[i] == 'MAG_ZEROPOINT\n':
                k_new[i] = 'MAG_ZEROPOINT    '+str(zp)+'\n'
            if k[i] == 'GAIN\n':
                k_new[i] = 'GAIN             '+str(g)+'\n'
            if k[i] == 'CHECKIMAGE_TYPE\n':
                k_new[i] = 'CHECKIMAGE_TYPE\n'
            if k[i] == 'CHECKIMAGE_NAME\n':
                k_new[i] = 'CHECKIMAGE_NAME \n'
            if k[i] == 'CATALOG_NAME\n':
                k_new[i] = 'CATALOG_NAME test \n'
            if k[i] == 'WEIGHT_IMAGE\n':
                k_new[i] = 'WEIGHT_IMAGE     %s%s%s_%s%s, %s%s%s_%s%s \n'%(
                        path_to_im, image_name, cat, detection_band, 
                        rmsfits_end, path_to_im, image_name, cat, 
                        name_band, rmsfits_end)
                
        nf = open(('SExtractor_files/parameters_' + name_band + '.sex'), 'w')
        nf.writelines(k_new)
        nf.close()


if __name__ == '__main__':
    main()
