from __future__ import print_function
# -*- coding: iso-8859-1 -*-
'''
Function: meridional_eigenmode
author:   Mike Bell   
Date:     May 2018
Description: calculates a normalised meridional eigenmode  
''' 
#          version 1:  15 May 2018 is first version
#                      5 April 2019 corrected exponent to use -y*y/4 in place of -y*y/2
#                      22 July 2019 -  dy_p introduced so that amplitudes of eigenfunctions are of order 1.  
# 
#-------------------------------------------------------------------------------------------------

def meridional_eigenmode( j_merid, c_m, lats, file_out):

# inputs 

# j_merid             - integer number of mode to calculate (0 <= j_merid <= 4)   
# c_m                 - the separation constant for the vertical eigenmode
# lats[jy]            - the latitudes in radians at which to calculate the modes
# file_out            - the output unit for standard diagnostics

# returns

# merid_mode[jy]      - the vertical projection factor for each of the modes
# dy_p                - the north-south grid-spacing in the y_p (y tilde) coordinate (a scalar - independent of latitude) 

    import math
    import numpy as np
    import sys
    from constants import earth_radius, OMEGA  

    lev_diag_user = 3    # >5 is used as test for diagnostic print-outs below. The level can be set by the user in each function independently  

#----------------------------------------------------------------------------------------------------------
#  1. Calculate normalised latitudes
#----------------------------------------------------------------------------------------------------------

    beta = 2.0 * OMEGA / earth_radius

    norm_fac_y = math.sqrt( 2.0 * beta / c_m ) 

    y_p = norm_fac_y * earth_radius * lats    # y_p stands for y' ; it has the same dimensions as lats

#----------------------------------------------------------------------------------------------------------
#  2. Calculate Hermite polynomials
#----------------------------------------------------------------------------------------------------------
    Hermite = np.zeros_like(lats)  
    if j_merid == 0 : 
      Hermite[:] = 1.0 
    elif j_merid == 1 :
      Hermite = y_p 
    elif j_merid == 2 :
      Hermite = y_p*y_p - 1.0 
    elif j_merid == 3 :
      Hermite = y_p*y_p*y_p - 3. * y_p  
    elif j_merid == 4 :
      Hermite = y_p*y_p*y_p*y_p - 6. * y_p*y_p + 3
    elif j_merid == 5 :
      Hermite = y_p*y_p*y_p*y_p*y_p  - 10. * y_p*y_p*y_p + 15 * y_p  
    else : 
      print ( "error 1 in meridional_eigenmode: you need to include more Hermite polynomials") 
      sys.exit()

#----------------------------------------------------------------------------------------------------------
#  3. Calculate the eigenmode pre-normalisation 
#----------------------------------------------------------------------------------------------------------

#    exponent_y = - 0.5 * y_p * y_p    line corrected 5 April 2019 

    exponent_y = - 0.25 * y_p * y_p 
    merid_mode = Hermite * np.exp( exponent_y ) 

#----------------------------------------------------------------------------------------------------------
#  3. Normalise the eigenmode  
#----------------------------------------------------------------------------------------------------------
    dy_p = y_p[1] - y_p[0]      #  the latitude spacing is uniform so dy_p can be taken to be a scalar 

    norm_sq = dy_p * np.sum ( merid_mode * merid_mode )  
    norm_amp = 1.0 / math.sqrt ( norm_sq ) 
    merid_mode = norm_amp * merid_mode
 
#----------------------------------------------------------------------------------------------------------
#  3. Diagnostic print-out of the eigenmode  
#----------------------------------------------------------------------------------------------------------

#  \int He_m^2 dy = \sqrt{2 pi} m! .   So the normalisation factor should be proportional to 1 / ( sqrt(m!) ). 
# Multiplying norm_fac by sqrt(m!) should give a constant. Note that 0! = 1. 

    j_factorial = math.factorial(j_merid)
    exact_integral = math.sqrt( j_factorial )      
    check_const = exact_integral * norm_amp 

    if lev_diag_user > 2:
       print(' meridional mode: j_merid, norm_amp, check_const = ', j_merid, norm_amp, check_const, file=file_out) 

    if lev_diag_user > 5:
       print('jy, merid_mode[jy] for c_m, j_merid = ', c_m, j_merid, file=file_out)
       if lev_diag_user > 5:
          for jy in range ( len( merid_mode ) ) :
             print(jy, merid_mode[jy], file=file_out)

#----------------------------------------------------------------------------------------------------------
    return merid_mode, dy_p 
#----------------------------------------------------------------------------------------------------------
