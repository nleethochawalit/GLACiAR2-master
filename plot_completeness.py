import numpy as np
import matplotlib.pyplot as plt


def main(path_to_cat,LF_shape, xarr, ymin, ymax, zbins, cat, 
         cs, ds, plotline = None):
    """
    Generates plots for C(m), S(z,m), and S(z,m)C(m)
    Args:
        path_to_cat (string) = Path to the folder with the science images.
                               Given in the parameters file.
        LF_shape (string arr) = name of lum functions
                                dimension must match the last dim of cso, dso
        xarr (float array) = input magnitude array assigned by the user.
        ymin (float) = Minimum input redshift assigned by the user.
        ymax (float) = Maximum input redshift assigned by the user.
        zbins (int) = Number of redshift bins dessignated by user.
        cat (string) = Name of the field for which the simulation is run.
        cs (int array) = Array with fraction of recovered sources for each
                         input magnitude and redshift.
        ds (int array) = Array with fraction of dropouts for each 
                         input magnitude and redshift.
        plotline (array) = Array of [npoint,2] for x,y pairs to be plot on 
                            the image
    """

    x = xarr
    y = np.linspace(ymin, ymax, zbins)
    xmin = np.min(x)
    xmax = np.max(x)
    spacing_x = round(5*(np.amax(x) - np.amin(x))/len(x), 1)
    spacing_y = round(5*(ymax - ymin)/(zbins), 1)

    with np.errstate(divide='ignore', invalid='ignore'):
        csm = np.true_divide(ds, cs)  # cover in case of division by 0.

    for lf in range(len(LF_shape)):
        plt.imshow(cs[:,:,lf], extent=[min(y), max(y), max(x), min(x)], cmap='RdPu',
                   origin="upper")
        plt.xlabel('$z$', fontsize=16)
        plt.ylabel('Input Absolute Magnitude', fontsize=16)
        plt.yticks(np.arange(xmin, xmax, spacing_x), fontsize=16)
        plt.xticks(np.arange(ymin, ymax, spacing_y), fontsize=16)
        plt.colorbar().set_label(label='$C(M)$', size=16)
        plt.savefig(path_to_cat+'Results/Plots/Completeness_Field'+
                    cat+'_'+LF_shape[lf]+'.pdf')
        plt.close()
    
    
        # Enter the IF below any dropout is different than 0 ()
        curds = ds[:,:,lf]
        if np.sum(curds) != 0.0:
            plt.imshow(curds, extent=[min(y), max(y), max(x), min(x)],
                       cmap='GnBu', origin="upper")
            plt.xlabel('$z$', fontsize=16)
            plt.ylabel('Input Absolute MAgnitude', fontsize=16)
            plt.yticks(np.arange(xmin, xmax, spacing_x), fontsize=16)
            plt.xticks(np.arange(ymin, ymax, spacing_y), fontsize=16)
            plt.colorbar().set_label(label='$C(m)S(z,m)$', size=16)
            plt.savefig(path_to_cat+'Results/Plots/Dropouts_Field'+cat+'_'+LF_shape[lf]+'.pdf')
            plt.close()
    
            plt.imshow(csm[:,:,lf], extent=[min(y), max(y), max(x), min(x)],
                       cmap='YlGn', origin="upper")
            plt.xlabel('$z$', fontsize=16)
            plt.ylabel('Input Absolute Magnitude', fontsize=16)
            plt.yticks(np.arange(xmin, xmax, spacing_x), fontsize=16)
            plt.xticks(np.arange(ymin, ymax, spacing_y), fontsize=16)
            plt.colorbar().set_label(label='$S(z,m)$', size=16)
            plt.savefig(path_to_cat+'Results/Plots/CS_Field'+cat+'_'+LF_shape[lf]+'.pdf')
            plt.close()
        