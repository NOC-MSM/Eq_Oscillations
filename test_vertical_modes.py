from __future__ import print_function
# -*- coding: iso-8859-1 -*-
'''
Function: test_vertical_eigenmodes
author:   Mike Bell   
Date:     May 2018
Description: checks the calculation of the vertical eigenmodes for given density field as a function of depth; 

This function is not used at the moment 

''' 
#          version 1:  15 May 2018 is first version
#
#-------------------------------------------------------------------------------------------------

def test_vertical_eigenmodes():

# note that z is height so is negative; dep is depth so is positive 

# inputs 

# dep_rho_levels[k]   - the depths of the density levels
# pot_rho[k]          - the potential densities at these levels
# n_modes             - the number of modes to return
# file_out            - the output unit for standard diagnostics

# returns

# eigen_out           - the eigenvalues of the first n_modes
# p_mode              - normalised pressure modes on the rho grid  

   import numpy as np
   from read_fields import read_one_field
   from set_pot_dens_analytical import set_pot_dens_analytical
   from vertical_eigenmodes import vertical_eigenmodes
   from vertical_eigenmodes_shooting import vertical_eigenmodes_shooting

#----------------------------------------------------------------------------------------------------------
#  0. Specify the user inputs 
#----------------------------------------------------------------------------------------------------------

   region_lc = 'pac_e' 
   j_lat_row_min = 50      #  latitude of row of densities chosen between 0 and 99 
   j_lat_row_max = 51      #  latitude of row of densities chosen between 0 and 99 

# Adam chose k_bottom_w = 74 which gives different results 
   k_bottom_w = 65   # w-level of the bottom (k = 0 is at z = 0; the bottom is at w-level k)  level 65 is at 3995 metres in the 75 level grid
                      # the number of rho depth levels is k_bottom_w; the number of w-levels is k_bottom_w + 1

   dir_data_in = '/data/local/frmk/AdamBlaker201810/'
   vertical_density_profile_input_file = dir_data_in + 'VN206HF_4h_200604_sig0_' + region_lc + '.nc'

   ll_vertical_density_profile_analytical = False    # if False uses a data input file to set the density profile

   if ll_vertical_density_profile_analytical :
       dens_type = 'exponential'     # options 'linear', 'exponential' (use lower case)
       dens_diff = -6.0         #   density at top minus density at the bottom (level k_bottom_w)
       exp_scale_depth = 500.   # depth scale in metres for exponential case

   ll_vertical_modes_shooting = True     # True => use shooting method; False => use matrix solution method

   n_vertical_modes = 3  

#----------------------------------------------------------------------------------------------------------
#  1. Open the output file and output the user inputs 
#----------------------------------------------------------------------------------------------------------

   if ll_vertical_density_profile_analytical :
      expt_name = 'analytical_'+dens_type+str(dens_diff)+str(exp_scale_depth) + '_' 

   else: 
      expt_name = region_lc + str(j_lat_row_min) + '_'+ str(j_lat_row_max) + '_'

   if ll_vertical_modes_shooting :
     expt_name = expt_name + 'shooting_'
   else : 
     expt_name = expt_name + 'matrix_'

   expt_name = expt_name + str(k_bottom_w) 

   filestem = './std_out/v_modes/'+expt_name
   file_out = open(filestem+'.txt','w')

   print(' region_lc, j_lat_row_min, j_lat_row_max = ', region_lc, j_lat_row_min, j_lat_row_max, file=file_out)

   print(' k_bottom_w  = ', k_bottom_w, file=file_out)

   print(' vertical_density_profile_input_file = ', vertical_density_profile_input_file, file=file_out)

   print(' ll_vertical_density_profile_analytical = ', ll_vertical_density_profile_analytical, file=file_out)
   if ll_vertical_density_profile_analytical :
       print(' dens_type, dens_diff, exp_scale_depth = ', dens_type, dens_diff, exp_scale_depth, file=file_out)

   print(' ll_vertical_modes_shooting = ', ll_vertical_modes_shooting, file=file_out)

#----------------------------------------------------------------------------------------------------------
#  2. Set the densities 
#----------------------------------------------------------------------------------------------------------

   deps = read_one_field ( vertical_density_profile_input_file, 'deptht', file_out )
   print ( 'deps = ', deps , file = file_out ) 

# just a check that the depth levels agree: time_series_input_file is not used hereafter 
   time_series_input_file = dir_data_in + 'VN206HF_4h_pacific_458_538_bdy_pressures.nc'
   dep_rho_levels  = read_one_field(time_series_input_file, 'gdept_0', file_out )
   print ( 'dep_rho_levels = ', dep_rho_levels , file = file_out ) 

   if ll_vertical_density_profile_analytical :
      pot_rho = set_pot_dens_analytical( k_bottom_w, deps, dens_type, dens_diff, exp_scale_depth, file_out )

   else:
      pot_rho_in = read_one_field ( vertical_density_profile_input_file, 'vosigmai', file_out )
      print ( ' pot_rho_in.shape = ', pot_rho_in.shape, file = file_out ) 

      pot_rho = np.zeros ( k_bottom_w ) 
      for k in range ( k_bottom_w ) : 
         for j in range ( j_lat_row_min, j_lat_row_max ) : 
            pot_rho[k] = pot_rho[k] + pot_rho_in[0,k,j]    # these are zonally averaged density fields  
         pot_rho[k] = pot_rho[k] / float ( j_lat_row_max - j_lat_row_min ) 

      if len(dep_rho_levels) != len(deps) :
        print(' input data set vertical dimensions do not match. Exiting ', file=file_out)
        sys.exit()

   print ( ' k   pot_rho[k] ' , file=file_out)
   for k in range( k_bottom_w ) : 
      print ( k, pot_rho[k] , file=file_out)

#----------------------------------------------------------------------------------------------------------
#  3. Calculate the eigenvalues  
#----------------------------------------------------------------------------------------------------------

# note that dz_layers values are all negative 
   if ll_vertical_modes_shooting:
       c_mode, p_mode, dz_layers = vertical_eigenmodes_shooting( k_bottom_w, deps, pot_rho, n_vertical_modes, file_out )
   else:
       c_mode, p_mode, dz_layers = vertical_eigenmodes( k_bottom_w, deps, pot_rho, n_vertical_modes, file_out )

#----------------------------------------------------------------------------------------------------------
test_vertical_eigenmodes()