from __future__ import print_function
# -*- coding: iso-8859-1 -*-
'''
Function: vertical_eigenmodes_shooting, shoot
author:   Mike Bell   
Date:     May 2018
Description: calculate the vertical eigenmodes for given density field as a function of depth using a shooting method
''' 
#          version 1:  15 May 2018 is first version
#          version 2:  22 July 2019 -  normalisation of depth introduced so that amplitudes are properly dimensional 
#
#-------------------------------------------------------------------------------------------------

def shoot (k_bottom_w, dz_ph_rho, N_sq_dz_w, lambda_test ): 

# inputs 

# k_bottom_w          - w-level of the bottom (the bottom is at w-level k ; k = 0 is the top level)    
# dz_ph_rho           - depths of the layers containing the densities (at which rho, u and v are held) 
# N_sq_dz_w           - N^2 times dz_w evaluated at the w levels 
# lambda_test         - test value of lambda (the eigenvalue)  

# outputs

# n_zeros             - number of times the solution changes sign 
# h_bdy               - value of the solution at the far boundary i.e. h[k_bottom_w]
# h[k]                - the solution at all levels 

# integrates d^2 h / dz^2  = - lambda_test N^2 h from k=1 to k=k_bottom_w-1 setting h[0] = 0 and h[1] = 1.0 and using 
# (h[k+1] - h[k])/dz_rho[k] - (h[k] - h[k-1])/dz_rho[k-1]  = - lambda N_sq_w[k] dz_w[k] h[k]  

   import numpy as np 

   h = np.zeros( k_bottom_w + 1) 
   h[1] = 1.0 
   n_zeros = 0

   for k in range (  1, k_bottom_w ) : 
      temp = - lambda_test*N_sq_dz_w[k]*h[k] + (h[k] - h[k-1])/dz_ph_rho[k-1]
      h[k+1] = h[k] + dz_ph_rho[k]*temp 
      if k <  k_bottom_w - 1 and h[k+1]*h[k] < 0.0 :      # corrected 28 Sep 2018. Only count changes in sign above the lower boundary 
        n_zeros = n_zeros + 1        

   h_bdy = h[k_bottom_w] 

   return n_zeros, h_bdy, h

def vertical_eigenmodes_shooting(k_bottom_w, dep_rho_levels, pot_rho, n_modes, file_out):

# note that z is height so is negative; dep is depth so is positive 

# inputs 

# k_bottom_w          - w-level of the bottom (the bottom is at w-level k ; k = 0 is the top level)    
# dep_rho_levels[k]   - the depths of the density levels
# pot_rho[k]          - the potential densities at these levels
# n_modes             - the number of modes to return
# file_out            - the output unit for standard diagnostics

# returns

# c_mode[m]           - the eigenvalues (as equivalent phase speed c_e) of the first n_modes
# p_mode[m, k]        - normalised pressure modes on the rho grid  
# dz_ph_rho           - depths of the layers containing the densities (at which rho, u and v are held) 
# H_bottom            - depth of bottom of domain used for eigenfunctions 

    import math
    import numpy as np 
    import sys as sys
    from scipy.linalg import eig 
    from dep_w_levels import dep_w_levels
    from constants import gdef, rho_0

    lev_diag_user = 5    # >2 and >5 are used as tests for diagnostic print-outs below. The level can be set by the user in each function independently  


#----------------------------------------------------------------------------------------------------------
#  1. Calculate the vertical grid and grid spacings

# the kth rho level lies between the kth and (k+1)th rho levels. So rho levels are called _ph in places 
#----------------------------------------------------------------------------------------------------------

    dz_ph_rho = np.zeros ( k_bottom_w )
    dz_w      = np.zeros ( k_bottom_w + 1 )

    dep_w_levels = dep_w_levels( dep_rho_levels, file_out ) 

    for k in range(k_bottom_w) :
       dz_ph_rho[k] = - dep_w_levels[k+1] + dep_w_levels[k]    # note dz is negative

    dz_w[0] = - 2. *  dep_rho_levels[0]
    for k in range(1, k_bottom_w) : 
       dz_w[k] = - dep_rho_levels[k] + dep_rho_levels[k-1] 
    dz_w[k_bottom_w] = - 2. *  ( dep_w_levels[k_bottom_w] - dep_rho_levels[k_bottom_w-1] )

#----------------------------------------------------------------------------------------------------------
#  2. Calculate the Brunt Vaisala frequency (it should be positive !) 
#----------------------------------------------------------------------------------------------------------
    N_sq_dz_w = np.zeros ( k_bottom_w )   # values are not required (or calculated) at k=0 and k=k_bottom_w

    for k in range(k_bottom_w - 1) :
       N_sq_dz_w[k+1] =  - ( gdef / rho_0 ) * ( pot_rho[k+1] - pot_rho[k] )    # N^2 dz_w  at w levels 

    if lev_diag_user > 10:
       print('k, pot_rho[k], dep_rho_levels[k], dep_levels[k+1], dz_ph_rho[k], dz_w[k], N_sq_dz_w[k] ', file=file_out)
       for k in range(k_bottom_w) :
          print(k, pot_rho[k], dep_rho_levels[k], dep_levels[k+1], dz_ph_rho[k], dz_w[k], N_sq_dz_w[k], file=file_out)

#----------------------------------------------------------------------------------------------------------
#  3. Find eigenvalues lambda_mode and eigenfunctions h
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
# 3.0 initialise with reasonable first guesses for lambda_mode[0] = c_e^{-2}; c_0 is about 3 m/s and start loop over modes   
#----------------------------------------------------------------------------------------------------------

    iter_max = 100   # this should be a large enough value
    max_tol = 1.E-6  # this should be small enough 
    lambda_0 = 0.01   
    lambda_1 = 0.1      

    p_mode = np.zeros ( (n_modes, k_bottom_w) ) 
    c_mode = np.zeros ( n_modes ) 

    for j_mode in range(n_modes): 
 
#----------------------------------------------------------------------------------------------------------
# 3.1 find a value of lambda_0 smaller than the eigenvalue for this mode 
#----------------------------------------------------------------------------------------------------------
       n_zeros_0, h_bdy_0, h_0 = shoot (k_bottom_w, dz_ph_rho, N_sq_dz_w, lambda_0 ) 

       if lev_diag_user > 5:
          print (' n_zeros_0, h_bdy_0, lambda_0 = ', n_zeros_0, h_bdy_0, lambda_0, file=file_out)

       for iter in range(iter_max) : 
          if n_zeros_0 >  j_mode : 
             lambda_1 = lambda_0 
             lambda_0 = 0.5 * lambda_0 
             n_zeros_0, h_bdy_0, h_0 = shoot (k_bottom_w, dz_ph_rho, N_sq_dz_w, lambda_0 ) 
             print ( 'n_zeros_0, h_bdy_0, lambda_0 = ', n_zeros_0, h_bdy_0, lambda_0, file=file_out)

          else : 
             break 

       if n_zeros_0 > j_mode :
         print(' failed to find a sufficiently small value for lambda_0: exiting ', lambda_0, file=file_out)
         sys.exit()

#----------------------------------------------------------------------------------------------------------
# 3.2 find a value of lambda_1 larger than the eigenvalue for this mode 
#----------------------------------------------------------------------------------------------------------

       n_zeros_1, h_bdy_1, h_1 = shoot (k_bottom_w, dz_ph_rho, N_sq_dz_w, lambda_1 )
       if lev_diag_user > 5:
          print ( 'n_zeros_1, h_bdy_1, lambda_1 = ', n_zeros_1, h_bdy_1, lambda_1, file=file_out)

       for iter in range(iter_max) : 
          if n_zeros_1 <  j_mode+1 : 
             lambda_1 = 2.0 * lambda_1 
             n_zeros_1, h_bdy_1, h_1 = shoot (k_bottom_w, dz_ph_rho, N_sq_dz_w, lambda_1 ) 
             if lev_diag_user > 5:
                print ( 'n_zeros_1, h_bdy_1, lambda_1 = ', n_zeros_1, h_bdy_1, lambda_1, file=file_out)

          else : 
             break 

       if n_zeros_1 < j_mode+1 :
         print ( ' failed to find a sufficiently small value for lambda_0: exiting ', lambda_0, file=file_out)    
         sys.exit()
        
#----------------------------------------------------------------------------------------------------------
#  3.3 Halve the interval between lambda_0 and lambda_1 until the solution is accurate enough  
#----------------------------------------------------------------------------------------------------------
       for iter in range(iter_max) :

          lambda_new = 0.5 * (lambda_0 + lambda_1)  
 
          n_zeros, h_bdy, h_mode = shoot (k_bottom_w, dz_ph_rho, N_sq_dz_w, lambda_new ) 
          if lev_diag_user > 5:
             print ( 'n_zeros, h_bdy, lambda_new = ', n_zeros, h_bdy, lambda_new, file=file_out)

          if abs(h_bdy) < max_tol and n_zeros_0 + 1 == n_zeros_1 :    

# we have found the solution ! set initial values for next mode before "exiting" iterations
             lambda_0 = lambda_new
             lambda_1 = 2. * lambda_new

             break     

          if n_zeros < j_mode+1 : 
            lambda_0  = lambda_new 
            n_zeros_0 = n_zeros
          else :
            lambda_1  = lambda_new 
            n_zeros_1 = n_zeros
           
#----------------------------------------------------------------------------------------------------------
#  4. Find the p_mode solution and normalise it   
#----------------------------------------------------------------------------------------------------------
       H_bottom = dep_w_levels[k_bottom_w+1]
       p_sum = 0.0 
       for k in range( k_bottom_w ) : 
          p_mode[j_mode, k] = - ( h_mode[k+1] - h_mode[k] )  / dz_ph_rho[k] 
          p_sum = p_sum - p_mode[j_mode,k] * p_mode[j_mode,k] * ( dz_ph_rho[k]  / H_bottom ) 

       recip_amp = 1.0 / math.sqrt(p_sum) 

       p_mode[j_mode] = recip_amp * p_mode[j_mode]  

       if lambda_new < 0.0 : 
          print ( ' lambda_new < 0.0. This should not happen. Exiting ', lambda_new, file=file_out)
          sys.exit()

       c_mode[j_mode] = 1.0 / math.sqrt( lambda_new ) 
       if lev_diag_user > 3:
          print ( ' lambda_new = ', lambda_new, file=file_out)

    if lev_diag_user > 2:
       print ( ' H_bottom = ', H_bottom, file=file_out)
       print ( ' c_mode = ', c_mode, file=file_out)
       print ( ' p_mode[:,0] = ', p_mode[:,0], file=file_out)

# used this to look for resonantly interacting modes - but this can't explain missing peaks in power spectrum 
       for j_merid in range (6): 
          print ( ' j_merid, (2j_merid+1)*c_mode  = ', j_merid, (2*j_merid+1)*c_mode, file = file_out)   
       if lev_diag_user > 5:

          for jmode in range ( n_modes ):
             print ( ' jmode, p_mode[jmode,k]', jmode, [ p_mode[jmode, k] for k in range (k_bottom_w ) ], file=file_out)  

#----------------------------------------------------------------------------------------------------------
    return c_mode, p_mode, dz_ph_rho, H_bottom 
#----------------------------------------------------------------------------------------------------------
