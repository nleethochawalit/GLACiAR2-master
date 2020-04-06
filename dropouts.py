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
    #from Bouwens2011
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

def main(mag_iso, mag_auto, sn, st, magbins, redshift):
    """
    Dropout selection for galaxies from BoRG at z~10
    Args:
        mag_iso (float array) = ISO magnitudes for each band and object.
        mag_auto (float array) = AUTO magnitudes for each band and object.
        sn (float array) = Signal to noise for each band and object.
        st  (int array) = Indicates the identification status of the object
                (non detected, blended, detected).
        magbins (float array) = bins of observed magnitude bins 
        redshift (float) = redshift of the galaxies so that the dropout
                            criterion can be determined accordingly
       
    Returns:
        drops (int array) = array of size st. 1 if it's a dropout. 0 if it's not.
    """
    #drops = drop_borg(mag_iso, mag_auto, sn, st, magbins, redshift)
    drops = drop_candels(mag_iso, mag_auto, sn, st, magbins, redshift)
    return drops

if __name__ == "__main__":
    main()
