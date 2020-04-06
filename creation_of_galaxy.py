"""
Contains the main equations and calculations needed for the program.
"""

import numpy as np
import random
import pysysp
import scipy.special
import scipy.optimize
import scipy.integrate
from astropy.io import fits
from astropy.table import Table
from astropy.cosmology import Planck15 as cosmo
from astropy import units as u
def gammainc(s, x):
    """
    Define and return the value of the gamma incomplete function.
    """
    def integrand(t, s): return t**(s-1) * np.exp(-1*t)
    gi = scipy.integrate.quad(integrand, 0, x, args=s)[0]
    return gi


def gamma(s):
    """
    Define and return the value of the gamma function.
    """
    def integrand(t, s): return t**(s-1) * np.exp(-1*t)
    gi = scipy.integrate.quad(integrand, 0, np.inf, args=s)[0]
    return gi


def get_bn(n0):
    """
    Calculates the parameter bn from the Sersic profile.

    Args:
        n0 (int) = Sersic index.

    Returns:
        bn (float) = Value of bn.
    """
    def errfunc(bn, n):
        return abs(scipy.special.gamma(2*n0) -
                   2*scipy.special.gammainc(2*n0, bn) *
                   scipy.special.gamma(2*n0))
    bn = scipy.optimize.fmin(errfunc, 1., args=(n0,), disp=False)[0]
    return bn


def get_Ie_n1(flux0, minor_axis, major_axis):
    """
    Calculates de parameter Ie for the Sersic profile with n = 1,
    which corresponds to the intensity at the radius that encloses
    half of the total light of the galaxy,  the effective radius.

    Args:
        flux0 (float) = Total flux for the galaxy.
        minor_axis (float) = Minor axis of the ellipse in pixels.
        major_axis (float) = Major axis of the ellipse in pixels.

    Returns:
        ie (float) = Value of Ie
    """
    ie = 1 / (2 * minor_axis * major_axis * np.pi)
    return ie


def get_Ie(bn0, n0, flux0, re):
    """
    Calculates de parameter Ie for the Sersic profile with n = 4,
    which corresponds to the intensity at the radius that encloses
    half of the total light of the galaxy, the effective radius.

    Args:
        bn0 (float) = bn parameter.
        flux0 (float) = Total flux for the galaxy.
        re (float) = effective radius in pixels.
    Returns:
        ie (float) = Value of Ie.
    """
    ie = ((bn0)**(2 * n0)) / (re**2 * 2 * np.pi * n0 * gamma(2*n0))
    return ie


def calculate_distance(x1, y1, x2, y2):
    """
    Calculates the distance between two pairs of coordinates.
    Args:

    Returns:
       dis (float) = distance in the same units as the coordinates.
    """
    dis = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    return dis


def makeSersic(n0, bn0, re, ell, inc_angle, size_galaxy):
    """
    Calculates the flux for each pixel following a Sersic profile.

    Args:
        n0 (int) = Sersic index.
        bn0 (float) = bn parameter.
        re (float) = Effective radius in pixels.
        ell (float) = Eccentricity. Varies between 0 and 1. 0 is circle
        inc_angle (float) = Inclination angle in radians. Varies
                            between 0 and Pi/2.
        size_galaxy (int) = Diameter of the galaxy stamp in pixels.
    Returns:
        fl (float) = Flux for each pixel.
    """
    stamp = np.zeros((size_galaxy,size_galaxy))
    s2 = size_galaxy / 2
    major_axis = re
    minor_axis = re * (1.-ell)
    I_e = ((bn0)**(2*n0)) / (2*np.pi*n0*major_axis*minor_axis*gamma(2*n0))
    def f(x, y):
        x_aux = (x-s2)*np.cos(inc_angle) + (y-s2)*np.sin(inc_angle)
        y_aux = -(x-s2)*np.sin(inc_angle) + (y-s2)*np.cos(inc_angle)
        radius = np.sqrt((x_aux/major_axis)**2 + (y_aux/minor_axis)**2)
        return I_e * np.exp(-bn0*((radius)**(1./n0)))

    for i in range(size_galaxy):
        def g(x):
            return i - 1./2.
        def h(x):
            return i + 1./2.
        for j in range(size_galaxy):
            fl = scipy.integrate.dblquad(f, j-1./2., j+1./2., g, h, 
                                 epsabs=1.49e-08, epsrel=1.49e-08)[0]
            stamp[i,j] = fl
    return stamp


def galaxies_positions(image_data, nsources, size, re):
    """
    Provides the position in the image for the simulated galaxies.

    Args:
        image_data (float array) = Science image. It corresponds to
                                   an array with the counts value for
                                   each pixel.
        nsources (int) = Number of simulated galaxies per iteration.
        size (int) = Diameter of the galaxy stamp in pixels.
        re (float) = Effective radius (in pixels).
    Returns:
        xpos (int) = Position of the simulated galaxy in the x axis.
        ypos (int) = Position of the simulated galaxy in the y axis.
    """
    s2 = size / 2
    xpos, ypos = np.zeros(nsources), np.zeros(nsources)
#    for i in range(nsources):
#        xr, yr = 1141, 1156
#        while ((xr == 0 and yr == 0) or image_data[xr, yr] == 0):
#            xr = random.randrange(s2 + 1, image_data.shape[0] - s2 - 1, 1)
#            yr = random.randrange(s2 + 1, image_data.shape[1] - s2 - 1, 1)
#            if (image_data[xr, yr] != 0):
#                xpos[i] = int(xr)
#                ypos[i] = int(yr)
    for j in range(nsources):
        d2 = calculate_distance(xpos[j], ypos[j], xpos, ypos)
        w1 = np.where(np.logical_and(np.array(d2) <= re,
                      np.array(d2) > 0.0))[0]
        while ((xpos[j] == 0 and ypos[j] == 0) or
               (image_data[int(xpos[j]), int(ypos[j])] == 0) or
               (len(w1) >= 1)):
            xr = random.randrange(s2 + 1, image_data.shape[0] - s2 - 1, 1)
            yr = random.randrange(s2 + 1, image_data.shape[1] - s2 - 1, 1)
            xpos[j] = int(xr)
            ypos[j] = int(yr)
            d2 = calculate_distance(xpos[j], ypos[j], xpos, ypos)
            w1 = np.where(np.logical_and(np.array(d2) <= 5 * re,
                          np.array(d2) > 0.0))[0]
    #xpos[0] = int(1138)
    #ypos[0] = int(1152)
    return xpos, ypos


def spectrum(a, x, b, redshift):
    """
    Creates observed spectrum following a Lyman break galaxy model.
    Args:
        a (float) = Normalisation factor for the spectrum.
        x (float array) = Array with points for the observed wavelength axis.
        b (float) = Value of the UV spectral slope.
        redshift (float) = Input redshift of the artificial source.
    Returns:
        f (float array )= flux for each wavelength for the spextrum
    """
    x = x / (1.+redshift)
    f = np.piecewise(x, [x < 1216, x >= 1216], [0, lambda x: a*(x**b)/
        (1.+redshift)])
    return f


def write_spectrum(detectionband, mag, beta, redshift,absmag=False,refwl=1600):
    """
    Saves the spectrum (spec.fits) of the simulated Lyman break galaxy
    so the magnitudes expected in each filter can be calculated.
    If absmag = False: take (mag, beta, redshift), returned observed flambda
    If absmag = True: take (Mag, beta, redshift), returned observed flambda
    Observed flambda has zero point magnitude = 0
    Args:
        detectionband (str) = name of the detection band  
        mag (float) = Input mangitude of the artificial source in the 
                      detection band. If absmag = True, it should be absolute
                      magtitude. Otherwise, it's apparent magnitude.
        beta (float) = Value of the UV spectral slope.  
        redshift (float) = Input redshift of the artificial source.
        absmag (boolean) = Whether the magnitude is absolute magnitude
        refwl (float) = Must provide if absmag is True. Wavelength for absmag 
                        reference in angstrom.
    """
    wavelength = np.arange(0, 30000, 1.)
    if absmag:
        c = 2.99792458e18  # In angstrom
        Lumdis = cosmo.luminosity_distance(redshift).to(u.parsec).value
        #define f_lambdae as L_lambdae/(4pi D_L^2) where f_lambdao = f_lambdae/(1+z)
        f_lambdae = ((c * 10**(-mag/2.5)) / (refwl**2))*(10./Lumdis)**2
        #this flux has zero_point magnitude = 0
        norm = f_lambdae/refwl**beta
    else:        
        ftilde = pysysp.StarSpectrum()
        ftilde.wavelength = wavelength
        ftilde.flux = spectrum(1., wavelength, beta, redshift)
        bandpass = pysysp.BandPass('Files/' + detectionband + '.tab')
        mtilde = ftilde.apmag(bandpass,mag='AB')+48.6
        #NL: sorry for making this a bit confusing. But +48.6 is to keep the flux
        #flambda to have zero_point magnitude = 0. This is to keep it close to 
        #the original code.
        norm = 10**((mtilde-mag)/2.5)
        
    spec = spectrum(norm, wavelength, beta, redshift)
    data = [[wavelength],[spec]]
    hdu = fits.PrimaryHDU()
    hdu.header['CTYPE1'] = 'Wavelength'
    hdu.header['CTYPE2'] = 'Flux'
    t = Table(data, names=('Wavelength', 'flux'))
    data = np.array(t)
    fits.writeto('Files/spec.fits', data,header=None, overwrite=True)


def mag_band(name_band):
    """
    Args
        name_band (string) = name of the band for which the magnitude of the
                             simulated galaxy is calculated.
    Returns
        mag (float): Expected AB magnitude in "name_band".
    """
    synth_spec = pysysp.StarSpectrum('Files/spec.fits')
    flux = pysysp.BandPass('Files/' + name_band + '.tab')
    mag = synth_spec.apmag(flux, mag='AB') + 48.6
    return mag

def absmag_refwl(redshift,refwl):
    #Calculate absolute magnitude at refwl (rest frame) from spec.fits
    
    c = 2.99792458e18 #in angstrom/s
    spec = pysysp.StarSpectrum('Files/spec.fits')
    flambda_e = (1+redshift)*np.interp(refwl*(1+redshift),
                 spec.wavelength[0],spec.flux[0])
    flambda_10pc = flambda_e*(cosmo.luminosity_distance(redshift).\
                              to(u.parsec).value/10)**2 
    absmag = -2.5*np.log10(refwl**2*(flambda_10pc)/c)
    return absmag

def schechter_func(m, mstar, alpha):
    """
    not a normalized schechter function
    m is Absolute magnitude
    """
    f = 10**(0.4*(alpha+1.)*(mstar-m))*np.exp(-10**(0.4*(mstar-m)))
    return f

def calphi(shapefn,redshift, slope=2):
    """
    calculate relative number of galaxies whose spectra is 'Files/spec.fits'
    for a given LF, redshift 
    --input--
    shapefn (string) : current choices are 
        'flat'     - return array of 1s
        'shechter' - must provide Files/LF_Schechter_params.txt and redshift
        'linear'   - must provide slope, Default = 2
    """
    if shapefn.lower() == 'flat': phi = 1.
    if shapefn.lower() == 'linear':
        constant = 1.
        refwl  = 1450.
        absmag = absmag_refwl(redshift,refwl)
        phi = (absmag+27.)*slope+constant  #use M=-27 as to make sure that phi is positive
        if phi < 0: phi = 1
        
    if shapefn.lower() == 'schechter':
        f = open('Files/LF_Schechter_params.txt')
        k = f.readlines()
        f.close()
        zmed   = np.array([float(line.split()[0]) for line in k[1:]])
        mstar  = np.array([float(line.split()[3]) for line in k[1:]])
        alpha  = np.array([float(line.split()[4]) for line in k[1:]])
        refwl  = np.array([float(line.split()[6]) for line in k[1:]]) 
        
        #picking the value closest to the redshift
        wmin = np.argmin(np.abs(zmed-redshift))
        mstar = mstar[wmin]
        alpha = alpha[wmin]
        refwl = refwl[wmin]
        
        absmag = absmag_refwl(redshift,refwl) 
        
        #calculate phi
        phi = schechter_func(absmag,mstar,alpha)
        
    return phi

def calphi_marr(mbins, shapefn, redshift, msample = 5, minngal = 100, 
                slope = 2, expbase = 2):
    """
    calculate relative number of galaxies with input magnitude marr 
    (apparent magnitude). 
    Normalized so that each bin has at least minngal galaxy.
    Also returned a set of random magnitudes for each bins. The magnitudes are
    according to the LF
    
    --input--
    mbins : intrinsic absolute magnitude bins 
    shapefn (string) : current choices are 
        'flat'     - return array of 1s
        'shechter' - must provide Files/LF_Schechter_params.txt and redshift
        'linear'   - must provide slope, Default = 2
        'exp'      - must provide expbase, Default = 2
    minngal (int) : minimum number of galaxies in each bin, Default = 100
    msample (int) : number of sampling in each magnitude bin. This is usually
                    the number of iterations. (For each iteration j and each 
                    mbins(i), the simulated galaxies will have one
                    magnitude m_j (within the bin of mbins(i), then 
                    phiarr(mbins(i)) galaxies will be injected into the image)
    """
    if shapefn.lower() == 'flat': 
        phiarr, marr = flat(mbins, minngal, msample)
        refwl = None
        
    if shapefn.lower() == 'linear':
        phiarr, marr = linear(mbins, minngal, msample, slope)
        refwl = None
        
    if shapefn.lower() == 'schechter':
        phiarr, marr, refwl = schechter(mbins, redshift, minngal, msample)
        
    if shapefn.lower() == 'exp':
        phiarr, marr = exponential(mbins, minngal, msample, expbase)
        refwl = None
        
    return phiarr, marr, refwl

def schechter(marr,redshift,min_ngal,msample):
    """
    The Schechter function
    Calculate number of galaxies in each bin, according to Schechter function
    at the specified redshift. The function is specified in 
    Files/LF_Schechter_params.txt
    
    The function is normalized such that the magnitude bin with the smallest 
    number of galaxies has min_ngal galaxies.
    
    --- INPUT ---
    marr (float array): array of ABSOLUTE magnitude bins. Each bin is centered at 
                        specified M with a binsize of M[1]-M[0]
    redshift (float): the code will take the parameters in LF_Schechter_params
                      where the redshift is the closest to this input z
    min_ngal (float): minimum number of galaxies in each mag bin
    msample (int)   : number of sampling in each magnitude bin.
    --- OUTPUT ---
    ngals (array) = array of number of galaxies in each magnitude bin
    mlist_out (list) = list of np.array. Each array contains msample absolute 
                        magnitude of galaxies in the magnitude bin
    refwl (float) = reference wavelength for the magnitude
    """
    
    # parameters of the Schechter function
    f = open('Files/LF_Schechter_params.txt')
    k = f.readlines()
    f.close()
    zmed   = np.array([float(line.split()[0]) for line in k[1:]])
    mstar  = np.array([float(line.split()[3]) for line in k[1:]])
    alpha  = np.array([float(line.split()[4]) for line in k[1:]])
    refwl  = np.array([float(line.split()[6]) for line in k[1:]]) 
        
    #picking the value closest to the redshift
    wmin = np.argmin(np.abs(zmed-redshift))
    mstar = mstar[wmin]
    alpha = alpha[wmin]
    refwl = refwl[wmin]
    if alpha > -1: mzero = mstar-2.5*np.log10(alpha+1)
    
    #calculate number of galaxies
    halfmbin = 0.5*(marr[1]-marr[0])    
    integral = np.zeros(len(marr))
    for i,m in enumerate(marr):    
        integral[i] = scipy.integrate.quad(schechter_func, m-halfmbin, m+halfmbin,
             args=(mstar,alpha))[0]
    norm = min_ngal/np.amin(integral)
    ngals = (norm*integral).astype(int)
    print('Schechter function:')
    for i,m in enumerate(marr): print('generating %d galaxies at magnitude %1.1f'%(ngals[i],m))
    
    #randomly draw magnitude values from Schechter function using the 
    #rejection method    
    mlist_out = []
    for i,m in enumerate(marr):
        mlow = m-halfmbin
        mhigh = m+halfmbin
        
        #find the maximum prob in that magnitude bin
        #if alpha is less than -1, then schechter is strictly increasing fn
        #if alpha is greater than -1, then schechter has negative slope if mag
        #nitude is greater(fainter) than mstar-2.5*np.log10(alpha+1)
        pdfvals = schechter_func(np.array([mhigh,mlow]),mstar,alpha)
        maxpdf = np.amax(pdfvals)
        minpdf = np.amin(pdfvals)
        if (alpha > -1):
            if (mhigh > mzero) and (mlow < mzero):
                maxpdf = schechter_func(mzero,mstar,alpha)
                
        submarr = np.zeros(msample)
        for j in range(msample):
            while submarr[j] == 0:
                rand_xy = np.random.uniform(low=[mlow,minpdf],high=[mhigh,maxpdf])
                fval = schechter_func(rand_xy[0],mstar,alpha)
                if rand_xy[1] < fval:
                    submarr[j] = rand_xy[0]
        mlist_out.append(submarr)
    return ngals,mlist_out,refwl


def flat(marr,min_ngal, msample):
    """
    randomly assign magnitudes of galaxies with a uniform distribution.
    Each magnitude bin has min_ngal galaxies.
    
    --- INPUT ---
    marr (float array): array of magnitude bins. Each bin is centered at 
                        specified M with a binsize of M[1]-M[0]
    min_ngal (float): minimum number of galaxies in each mag bin
    --- OUTPUT ---
    ngals (array) = array of number of galaxies in each magnitude bin
    mlist_out (list) = list of np.array. Each array contains msample 
                        magnitudes of galaxies in the magnitude bin
    """
    halfmbin = 0.5*(marr[1]-marr[0])
    mlist_out = []
    ngals = np.zeros(len(marr),dtype=int)
    ngals[:] = min_ngal
    for m in marr:
        submarr = np.random.uniform(m-halfmbin,m+halfmbin,size=msample)
        mlist_out.append(submarr)
    return ngals,mlist_out

def linear(marr, min_ngal, msample, slope):
    """
    randomly assign magnitudes of galaxies with a linear distribution.
    slope is number per magnitude
    """    
    def lin_eq(x,x0,m,y0):
        y = (x-x0)*m+y0
        return y    

    ngals = lin_eq(marr,marr[0],slope,min_ngal)  
    
    #randomly draw magnitude values from linear PDF using the 
    #rejection method    
    mlist_out = []
    halfmbin = 0.5*(marr[1]-marr[0])
    for i,m in enumerate(marr):
        mlow = m-halfmbin
        mhigh = m+halfmbin
        
        pdfvals = lin_eq(np.array([mhigh,mlow]),marr[0],slope,min_ngal)
        
        maxpdf = np.amax(pdfvals)
        minpdf = np.amin(pdfvals)
                
        submarr = np.zeros(msample)
        for j in range(msample):
            while submarr[j] == 0:
                rand_xy = np.random.uniform(low=[mlow,minpdf],high=[mhigh,maxpdf])
                fval = lin_eq(rand_xy[0],marr[0],slope,min_ngal)
                if rand_xy[1] < fval:
                    submarr[j] = rand_xy[0]
        mlist_out.append(submarr)
    return ngals,mlist_out
    

def exponential(marr, min_ngal, msample, expbase):
    """
    randomly assign magnitudes of galaxies with an exponential distribution.
    """    
    def exp_eq(x,a,x0,y0):
        y = y0*(a**(x-x0))
        return y  

    ngals = exp_eq(marr,expbase,marr[0],min_ngal).astype(int)
    
    #randomly draw magnitude values from linear PDF using the 
    #rejection method    
    mlist_out = []
    halfmbin = 0.5*(marr[1]-marr[0])
    for i,m in enumerate(marr):
        mlow = m-halfmbin
        mhigh = m+halfmbin
        
        pdfvals = exp_eq(np.array([mhigh,mlow]),expbase,marr[0],min_ngal)  
        
        maxpdf = np.amax(pdfvals)
        minpdf = np.amin(pdfvals)
                
        submarr = np.zeros(msample)
        for j in range(msample):
            while submarr[j] == 0:
                rand_xy = np.random.uniform(low=[mlow,minpdf],high=[mhigh,maxpdf])
                fval = exp_eq(rand_xy[0],expbase,marr[0],min_ngal)  
                if rand_xy[1] < fval:
                    submarr[j] = rand_xy[0]
        mlist_out.append(submarr)
    return ngals,mlist_out
    