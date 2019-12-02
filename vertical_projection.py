from __future__ import print_function
# -*- coding: iso-8859-1 -*-
'''
Function: vertical_projection
author:   Mike Bell   
Date:     May 2018
Description: calculates the vertical projection factors of the wind stresses on the eigenmodes  
''' 
#          version 1:  15 May 2018 is first version
#          version 2:  22 July 2019 -  normalisation of depth introduced so that amplitudes are properly dimensional 
#
#-------------------------------------------------------------------------------------------------

def vertical_projection( p_mode, dep_rho_levels, MLD, H_bottom, file_out):

# note that z is height so is negative; dep is depth so is positive 

# inputs 

# p_mode[m, k]        - normalised pressure modes on the rho grid  
# dep_rho_levels[k]   - the depths of the density levels
# MLD                 - depth of penetration of winds (in metres) - the mixed layer depth
# H_bottom            - depth of bottom of domain used for eigenfunctions 
# file_out            - the output unit for standard diagnostics

# returns

# fac_mode[m] - the vertical projection factor for each of the modes

    import math
    import numpy as np
    from constants import rho_0  

#----------------------------------------------------------------------------------------------------------
#  1. Calculate the vertical grid and grid spacings

# the kth rho level lies between the kth and (k+1)th rho levels. So rho levels are called _ph in places 
#----------------------------------------------------------------------------------------------------------
    n_modes, n_dz_half = p_mode.shape 

    fac_mode = np.zeros ( n_modes )

    dz_ph_rho = np.zeros ( n_dz_half )
    dep_w_levels = np.zeros ( n_dz_half+1 ) 

    dep_w_levels[0] = 0.

    for k in range(n_dz_half) :
       dep_w_levels[k+1] = dep_w_levels[k] + 2.0 * ( dep_rho_levels[k] - dep_w_levels[k] )

    for k in range(n_dz_half) :
       dz_ph_rho[k] = - dep_w_levels[k+1] + dep_w_levels[k]    # note dz is negative

    fac_mode = np.zeros ( n_modes ) 
    for k in range(n_dz_half) :

       if dep_w_levels[k+1] < MLD : 
          fac_mode[:] = fac_mode[:] - p_mode[:, k] * ( dz_ph_rho[k] / H_bottom )   #  note these contributions are +ve 

       else :

# find the depth covered by the ML within this model layer 
          dep_remain = MLD - dep_w_levels[k]

# interpolate to the middle of this layer
          dep_target = dep_w_levels[k] + 0.5 * dep_remain 
          dep_0 = dep_rho_levels[k-1]
          dep_1 = dep_rho_levels[k] 
          alpha = (dep_1 - dep_target) / (dep_1 - dep_0)
          for n in range (n_modes) : 
             pmode_targ = alpha * p_mode[n,k-1] + (1.0 - alpha) * p_mode[n,k]   
#             fac_mode[n] = fac_mode[n] - p_mode[n, k] * dep_remain     # corrected 28 Sep 2018
             fac_mode[n] = fac_mode[n] + pmode_targ * ( dep_remain / H_bottom )        # this contribution is +ve 
          break 

    fac_mode = fac_mode / (rho_0 * MLD ) 
 
    print('fac_mode for vertical projections (should be just less than 9.8 times surface value) = ', fac_mode, file=file_out)

#----------------------------------------------------------------------------------------------------------
    return fac_mode 
#----------------------------------------------------------------------------------------------------------
