#for ASTR 136, author: Olivier Hervet

# crude calculator of the number of counts per pixel expected for a V magnitude of a star in the AB system for the Lick Nickel instrument
#atmospheric extinction is considered
# the (unknown?) effective area of the Nickel is not taken into account. This mean a 100% photon transmission through the 
# optical system. So the result is more an upper limit than an accurate estimation.

#Note: The number of counts need to be well above 100 and below 60k for each exposure


#import numpy as np
from pylab import*
from astropy.io import ascii
from scipy.interpolate import interp1d

c = 2.998e8 #speed fo light [m.s-1]
h = 6.6260755e-27 # Planck Constant	 [erg s]

#extinction part
#atm extinction at Lick is availabe here: https://mthamilton.ucolick.org/techdocs/standards/lick_mean_extinct.html
font = 'atm_extinction_lick.txt'
table_ext = ascii.read(font)
lambda_ext = table_ext['wavelength[nm]']
magn_ext = table_ext['extinction[magnitude/airmass]']
f_ext = interp1d(lambda_ext,magn_ext)
#plot(lambda_ext,magn_ext)



def nexpose(mV,t, El):
    #inputs: V magnitude of a star, t:exposure time [sec], El: elevation of the star [deg]
    #output: number of counts in image pixel per exposure
    
    #central wavelength of the filter in nanometers. For Blue 450, for Visible  600, and for Red  675
    #https://old.aip.de/en/research/facilities/stella/instruments/data/johnson-ubvri-filter-curves
    #V: Effective central wavelength 547.7nm, FWHM 99.1nm
    
    #I considered here a usual B filter, I didn't check with the actual filter on the Nickel
    
    lambda_V = 547.7e-9 #wavelength [m]
    Dlambda_V = 99.1e-9
    nu_V = c/lambda_V   #frequency [Hz]
    
    #observed magnitude, correct from the atm extinction
    mV += f_ext(lambda_V*1e9)/sin(El*pi/180)
    
    
    #energy flux
    fv = 10**(-0.4*(mV + 48.6)) * nu_V # erg.cm-2.s-1
    
    #photon flux (divide by energy of photons)
    ph_flux = fv/(h*nu_V) # ph.cm-2.s-1
    
    #mirror surface Nickel
    A = pi*(50)**2  #[cm2]
    
    #total photon flux on mirror ~ photon flux on focal plane
    ph_flux *= A #ph.s-1
    
    #Nickel fov
    #N_fov = 6.3**2 #arcmin2
    
    #pixel fov (0.184 arcsec/pixel)
    pix_fov = (0.184)**2 #arcsec2
    #considering a binning of 2x2 (4 camera pixel per image pixel)
    pix_fov *= 4 #arcsec2
    
    #PSF on nickel apparently ~2arcsec FWHM
    #use ~3 arcsec from previous observations
    PSF = pi*(3/2.)**2 #arcsec2
    
    #photon flux per image pixel for a point source
    ph_flux_impix = ph_flux *pix_fov/PSF
    
    
    #quantum efficiency V
    Queff = 0.828 #ph.elec-1 (approx)
    
    #Gain Nickel
    G = 1.70 #e-.cts
    
    #counts flux per image pix
    c_flux_impix = ph_flux_impix*Queff/G
    
    print('You should expect a maximum of', round(c_flux_impix*t,2),'count/sec per image pixel')
