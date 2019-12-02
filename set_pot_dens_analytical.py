from __future__ import print_function
# -*- coding: iso-8859-1 -*-
'''
Function: set_pot_dens_analytical
author:   Mike Bell   
Date:     May 2018
Description: calculates the potential density at each model level as an analytical function of depth specified by the user 
''' 
#          version 1:  15 May 2018 is first version
#
#-------------------------------------------------------------------------------------------------

def set_pot_dens_analytical ( k_bottom_w, dep_rho_levels, dens_type, dens_diff, exp_scale_depth, file_out):

# note that z is height so is negative; dep is depth so is positive 

# inputs 

#    k_bottom_w			 # level of the bottom (the bottom is at w-level k ; k = 0 is the top level) ; the number of rho-levels is k_bottom_w   
#    dens_type                   # options 'linear', 'exponential' 
#    dens_diff                   # difference in density between top and bottom kg m^{-3}; 6.0 is sensible  
#    exp_scale_depth             # depth scale in metres for exponential case; 500. is sensible  

# returns

# pot_dens [k]                   # potential density profile ; uses k_bottom_w as dimension to avoid arrays with spurious numbers in them 

    import math
    import numpy as np
    from dep_w_levels import dep_w_levels

    lev_diag_user = 3 

#----------------------------------------------------------------------------------------------------------

    pot_dens = np.zeros( k_bottom_w )

    dep_w_levels = dep_w_levels( dep_rho_levels, file_out ) 

    if dens_type == 'linear' : 

       factor = dens_diff / dep_w_levels[k_bottom_w] 

       for k in range ( k_bottom_w ) : 
          pot_dens[k] =  - factor * dep_rho_levels[k]
 
    else :    #    dens_type = 'exponential'

      for k in range ( k_bottom_w ) : 
         scaled_height =  - dep_rho_levels[k] / exp_scale_depth
         pot_dens[k] = dens_diff * math.exp( scaled_height ) 

    if lev_diag_user > 10: 
       print(' k, pot_dens[k] ', file=file_out)
       for k in range ( k_bottom_w ) :
          print(k, pot_dens[k], file=file_out)

#----------------------------------------------------------------------------------------------------------
    return pot_dens
#----------------------------------------------------------------------------------------------------------
