"""
Dropout selection
It applies a dropout redshift selection for galaxies from the BoRG
survey at z~10.
Mangitudes for this case:
    # mag f160w = mag[0][]
    # mag f098m = mag[1][]
    # mag f125w = mag[2][]
    # mag f606w = mag[3][]

This file can be modified by the user.
"""

import numpy as np
import pdb

def drop_borg25(mag_aper, sn_aper1, st, magbins, redshift):
    drops = np.zeros(mag_aper.shape[0],dtype=int)  # Initialise dropouts status array.    
    #mag bands are: f350,f098w,f105w,f125w,f140w, f160w
    for i in range(mag_aper.shape[0]):
        if (st[i] >= 0 and sn_aper1[i,5]>4):    
            #Check for Yband dropout z=8
            if ((sn_aper1[i,3] >6) and
                (mag_aper[i,2] - mag_aper[i,3] > 0.45) and
                (mag_aper[i,3] - mag_aper[i,5] < 0.5) and
                ((mag_aper[i,2] - mag_aper[i,3]) > 
                 1.5*(mag_aper[i,3] - mag_aper[i,5])+0.45) and
                (mag_aper[i,1] - mag_aper[i,3] > 1.75) and
                ((mag_aper[i,3] - mag_aper[i,5]) < 
                 0.15*(mag_aper[i,1] - mag_aper[i,3] - 1.75)+0.02)):
                    if sn_aper1[i,0]<1.5: 
                        drops[i] = -1
                    if sn_aper1[i,0]<1: 
                        drops[i] = 1
          
            #Check for Jband dropout
            if ((drops[i] == 0) and
                (sn_aper1[i,3] >6) and
                 (mag_aper[i,2] - mag_aper[i,4] > 1.5) and
                 (mag_aper[i,4] - mag_aper[i,5] < 0.3) and
                 ((mag_aper[i,2] - mag_aper[i,4]) > 
                  5.33*(mag_aper[i,4] - mag_aper[i,5])+0.7) and
                 (mag_aper[i,1] - mag_aper[i,4] > 1.5) and
                 ((mag_aper[i,1] - mag_aper[i,4]) > 
                  5.33*(mag_aper[i,4] - mag_aper[i,5])+0.7)):
                    if sn_aper1[i,0]<1.5: 
                        drops[i] = -2
                    if sn_aper1[i,0]<1: 
                        drops[i] = 2
            #Check for Hband dropout
            if ((drops[i] == 0) and
                (sn_aper1[i,3] >6) and
                (mag_aper[i,3] - mag_aper[i,5] > 1.3)):
                    if ((sn_aper1[i,0]<1.5) and (sn_aper1[i,1]<1.5) 
                        and (sn_aper1[i,3]<1.5)): 
                        drops[i] = -3
                    if ((sn_aper1[i,0]<1) and (sn_aper1[i,1]<1) 
                        and (sn_aper1[i,3]<1)): 
                        drops[i] = 3
        
    return drops

def drop_borg(mag_iso, mag_auto, sn, st, magbins, redshift):
    drops = np.zeros(sn.shape[0],dtype=int)  # Initialise dropouts status array.    
    for i in range(sn.shape[0]):
        #Yband dropout
        if (redshift < 9 and st[i] >= 0 and 
            (mag_iso[i, 1] - mag_iso[i, 2]) > 1.75 and 
            (mag_iso[i,2] - mag_iso[i,0] < 0.02 + 0.15*(mag_iso[i,1]-mag_iso[i,2]-1.75)) and 
            sn[i, 3] < 1.5):
            drops[i] = 1
        #Jband dropout
        if (redshift >= 9 and st[i] >= 0 and 
            (mag_iso[i, 2] - mag_iso[i, 0]) > 1.5 and 
            sn[i, 3] < 1.5 and sn[i, 1] < 1.5):
            drops[i] = 1
    return drops

def drop_candels(mag_iso, mag_auto, sn, st, magbins, redshift):
    #from Bouwens2011 (ignore the non-detection criterion.)
    drops = np.zeros(sn.shape[0],dtype=int)  # Initialise dropouts status array.    
    for i in range(sn.shape[0]):
        if (redshift < 7.5 and redshift > 6.5 and st[i] >= 0 and 
            (mag_iso[i, 1] - mag_iso[i, 2]) > 0.7 and 
            (mag_iso[i,2] - mag_iso[i,3]) < 0.45 and
            (mag_iso[i,1] - mag_iso[i,2] > 0.42 + 0.14*(mag_iso[i,2]-mag_iso[i,3]))):
            drops[i] = 1
        #Yband dropout
        if (redshift >= 7.5 and redshift < 8.5 and st[i] >= 0 and 
            (mag_iso[i, 2] - mag_iso[i, 3]) > 0.45 and 
            (mag_iso[i, 3] - mag_iso[i, 0]) < 0.5 and
            sn[i, 3] > 3.5 and sn[i, 0] > 3 and sn[i, 1] < 3):
            drops[i] = 1
    return drops

def drop_hlf_archived(mag_iso, mag_auto, sn, sn_aper1, sn_aper2,st, magbins, redshift):
    #from Bouwens2015 
    #  f160w f105w f125w f435w f606w f775w f814w
    drops = np.zeros(sn.shape[0],dtype=int)  # Initialise dropouts status array.   
    #calculate chisq optical for non-detection bands
    chisq_optY1 = np.sum(np.sign(sn_aper1[:,3:7])*sn_aper1[:,3:7]**2, axis=1)
    chisq_optY2 = np.sum(np.sign(sn_aper2[:,3:7])*sn_aper2[:,3:7]**2, axis=1)
    #Make sure that at least 1 band is available. Flux will be zero if it's 
    #not in the observed field.
    chisq_optY1[chisq_optY1==0] = 99
    chisq_optY2[chisq_optY2==0] = 99
    chisq_ir = np.sum(np.sign(sn_aper1[:,0:3])*sn_aper1[:,0:3]**2, axis=1)
    
    drops = np.where((st>=0) & ((mag_iso[:, 1] - mag_iso[:, 2]) > 0.45) & 
                     ((mag_iso[:, 2] - mag_iso[:, 0]) < 0.5) &
                     ((mag_iso[:, 1] - mag_iso[:, 2]) > 
                      0.75*(mag_iso[:, 2]- mag_iso[:, 0])+0.525) & 
                     (sn[:,0]>5.5) & (sn[:,2]>5.5) & (chisq_ir>25) & 
                     (sn[:,3]<2) & (sn[:,4]<2) & (sn[:,5]<2) & 
                     (sn[:,6]<2) & (chisq_optY1<4) &(chisq_optY2<3),1,0)

    return drops

def drop_hlf(mag_iso, mag_auto, mag_aper1, mag_aper2, 
             sn_iso, sn_aper1, sn_aper2, st):
    """
    From Bouwens2015
    List of bands: 0=f435w, 1=f606w, 2=f775w, 3=f814w, 4=f850lp
    5=f105w, 6=f125w, 7=f140w, 8=f160w
    aper1 = 0.2" diameter
    aper2 = 0.35" diameter
    
    Will use colors and S/N from aperture2
    Detected at 5sigma redward of the break
    Output: 
        drops array with values: 0 if not dropouts
        4 if redshift 4, 
        1 if redshift 5, 
        6 if redshift 6, 
        7 if redshift7
    """
    drops = np.zeros(st.shape[0],dtype=int)  # Initialise dropouts status array. 
    chisq_red = np.sum(np.sign(sn_aper2[:,5:9])*sn_aper2[:,5:9]**2, axis=1)

    #redshift 4
    drops = np.where((st >= 0) & 
                     (mag_aper2[:,0]-mag_aper2[:,1] > 1.0) &
                     (mag_aper2[:,2]-mag_aper2[:,6] < 1.0) &
                     (mag_aper2[:,0]-mag_aper2[:,1] > 1.6*(mag_aper2[:,2]-mag_aper2[:,6])+1.) &
                     (chisq_red > 25),
                     4, drops)
    #redshift 5
    drops = np.where((st >= 0) & 
                     (mag_aper2[:,1]-mag_aper2[:,2] > 1.2) &
                     (mag_aper2[:,4]-mag_aper2[:,8] < 1.3) &
                     (mag_aper2[:,1]-mag_aper2[:,2] > 0.8*(mag_aper2[:,4]-mag_aper2[:,8])+1.2) &
                     (sn_aper2[:,0] < 2) & (chisq_red > 25),
                     1, drops)
    #redshift 6                 
    drops = np.where((st >= 0) & 
                     (mag_aper2[:,2]-mag_aper2[:,4] > 1) &
                     (mag_aper2[:,5]-mag_aper2[:,8] < 1) &
                     (mag_aper2[:,2]-mag_aper2[:,4] > 0.78*(mag_aper2[:,5]-mag_aper2[:,8])+1.0) &
                     (sn_aper2[:,0] < 2) &
                     ((mag_aper2[:,1]-mag_aper2[:,4] > 2.7)|(sn_aper2[:,1] < 2)) &
                     (chisq_red > 25),
                     6, drops)
    # #redshift 7
    # chisq_optY1 = np.sum(np.sign(sn_aper1[:,0:3])*sn_aper1[:,0:3]**2, axis=1)
    # chisq_optY2 = np.sum(np.sign(sn_aper2[:,0:3])*sn_aper2[:,0:3]**2, axis=1)
    # drops = np.where((mag_aper2[:,4]-mag_aper2[:,5] > 0.7) &
    #                  (mag_aper2[:,6]-mag_aper2[:,8] < 0.45) &
    #                  (mag_aper2[:,4]-mag_aper2[:,5] > 0.8*(mag_aper2[:,6]-mag_aper2[:,8])+0.7) &
    #                  ((mag_aper2[:,3]-mag_aper2[:,6] > 1.0)|(sn_aper2[:,3] < 1.5))
    #                  (sn_aper2[:,0] < 2) & (sn_aper2[:,1] < 2) & 
    #                  (sn_aper2[:,2] < 2) & 
    #                  (chisq_optY1 < 4) & (chisq_optY2 < 3) & (chisq_red > 25),
    #                  7,drops)
    return drops

def main(mag_iso, mag_auto, mag_aper1, mag_aper2, sn_iso, sn_aper1, sn_aper2, 
         st, magbins, redshift, droptype):
    """
    Dropout selection for galaxies from BoRG at z~10
    Args:
        mag_iso (float array) = ISO magnitudes for each band and object.
        mag_auto (float array) = AUTO magnitudes for each band and object.
        mag_aper1,2 (float array) = APER magnitudes for each band and object.
                 Correspond to first,second aperture defined in sextractor parameters.                       
        sn_iso (float array) = Signal to noise (ISO) for each band and object.
        sn_aper1,2 (float array) = Signal to noise (APER) for each band and object.
        st  (int array) = Indicates the identification status of the object
                (non detected, blended, detected).
        magbins (float array) = bins of observed magnitude bins 
        redshift (float) = redshift of the galaxies so that the dropout
                            criterion can be determined accordingly
       droptype (string) = type of dropout. Choices are 'borg', 'hlf', 'candels'      
       
       **Note that each droptype doesn't require all data**
    Returns:
        drops (int array) = array of size st. 1 if it's a dropout. 0 if it's not.
    """
    if droptype.lower() == 'borg':
        drops = drop_borg(mag_iso, mag_auto, sn_iso, st, magbins, redshift)
    if droptype.lower() == 'candels':
        drops = drop_candels(mag_iso, mag_auto, sn_iso, st, magbins, redshift)
    if droptype.lower() == 'hlf':
        drops = drop_hlf(mag_iso, mag_auto, mag_aper1, mag_aper2, 
                         sn_iso, sn_aper1, sn_aper2, st)
    if droptype.lower() == 'borg25':
        drops = drop_borg25(mag_aper1, sn_aper1, st, magbins, redshift)    
    return drops

if __name__ == "__main__":
    main()
