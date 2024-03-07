___________________________________________________________________________________________________________________________________________________________________
___________________________________________________________________________________________________________________________________________________________________
___________________________________________________________________________________________________________________________________________________________________
# rovibAbsSpec_IDL

___________________________________________________________________________________________________________________________________________________________________
GENERAL DESCRIPTION:
This IDL script models a multi-component absorption spectrum in the form of cross-sectional ground-based telescope data.

___________________________________________________________________________________________________________________________________________________________________
DATA DESCRIPTION:
The data (gvtauSpec_JWK.dat) contains the neccessary variables (i.e. partition functions, Q(T); transition wavenumber, (e.g. f12_10); etc.). The data were
processed from raw telescope data, wavelength calibrated, and normalized. As such, the data is an absorption spectrum from ~1950 cm^-1 to 2200 cm^-1. 

___________________________________________________________________________________________________________________________________________________________________
CODE DESCRIPTION:
The code first defines the constants and other molecular properties. Then these molecular properties are used to calculate the optical depth of each transition as 
a function of wavenumber. These optical depths are the input for radiative transfer equations (i.e. e^-tau) from which the spectrum of each component is generated. 
Then for physical reasons, the lines are weighted with a veiling function, and then convolved with their respective doppler broadening. Then the total spectrum is 
constructed and is convolved with the instrument resolution. The code then plots the results in four .png files for each of the 12CO v=1-0, 12CO v=2-1, 13CO v=1-0, 
and C18O v=1-0 transitions. 

___________________________________________________________________________________________________________________________________________________________________
RUNNING THE CODE:
1) Download the IDL script (rovibAbsSpec_JWK.pro), the data file containing the spectrum (gvtauSpec_JWK.dat), and the remaining .csv files which contain the CO
   molecular properties data
 
2) In a terminal, cd into the directory that now contains the script and the data

3) Run the script by typing the following into the command line:

           restore,'/home/jwkern/Research/.../gvtauSpec_JWK.dat'
           .r rovibAbsSpec_JWK.pro
___________________________________________________________________________________________________________________________________________________________________
___________________________________________________________________________________________________________________________________________________________________
___________________________________________________________________________________________________________________________________________________________________
