from __future__ import print_function
# -*- coding: iso-8859-1 -*-
'''
Function: dep_w_levels
author:   Mike Bell   
Date:     May 2018
Description: calculate the depths of w levels from rho levels  
''' 
#          version 1:  15 May 2018 is first version
#
#-------------------------------------------------------------------------------------------------

def dep_w_levels(dep_rho_levels, file_out):

# inputs 

# dep_rho_levels[k]   - the depths of the density levels
# file_out            - the output unit for standard diagnostics

# returns

# dep_w_levels[k]     - the depth of w (vertical velocity) levels

    import numpy as np 

    lev_diag_user = 5

#----------------------------------------------------------------------------------------------------------
#  1. Calculate the vertical grid and grid spacings

# the kth rho level lies between the kth and (k+1)th rho levels. So rho levels are called _ph in places 
#----------------------------------------------------------------------------------------------------------
    n_dz_half = len(dep_rho_levels) 
    dep_w_levels = np.zeros ( n_dz_half+1 ) 

    dep_w_levels[0] = 0.

    for k in range(n_dz_half) :
       dep_w_levels[k+1] = dep_w_levels[k] + 2.0 * ( dep_rho_levels[k] - dep_w_levels[k] )

    if lev_diag_user > 10 : 
       print(' k,  dep_rho_levels[k], dep_levels[k] ', file=file_out)
       print([ (k,  dep_rho_levels[k], dep_levels[k]) for k in range ( n_dz_half ) ], file=file_out)
     
#----------------------------------------------------------------------------------------------------------
    return dep_w_levels
#----------------------------------------------------------------------------------------------------------
