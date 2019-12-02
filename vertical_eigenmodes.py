from __future__ import print_function
# -*- coding: iso-8859-1 -*-
'''
Function: vertical_eigenmodes
author:   Mike Bell   
Date:     May 2018
Description: calculate the vertical eigenmodes for given density field as a function of depth 
''' 
#          version 1:  15 May 2018 is first version
#
#-------------------------------------------------------------------------------------------------

def vertical_eigenmodes(k_bottom, dep_rho_levels, pot_rho, n_modes, file_out):

# note that z is height so is negative; dep is depth so is positive 

# inputs 

# k_bottom            - level of the bottom (the bottom is at w-level k ; k = 0 is the top level)    
# dep_rho_levels[k]   - the depths of the density levels
# pot_rho[k]          - the potential densities at these levels
# n_modes             - the number of modes to return
# file_out            - the output unit for standard diagnostics

# returns

# c_mode[m]           - the eigenvalues (as equivalent phase speed c_e) of the first n_modes
# p_mode[m, k]        - normalised pressure modes on the rho grid  
# dz_ph_rho           - depths of the layers containing the densities (at which rho, u and v are held) 

    import math
    import numpy as np 
    from scipy.linalg import eig 
    from dep_w_levels import dep_w_levels
    from constants import gdef, rho_0

#----------------------------------------------------------------------------------------------------------
#  1. Calculate the vertical grid and grid spacings

# the kth rho level lies between the kth and (k+1)th rho levels. So rho levels are called _ph in places 
#----------------------------------------------------------------------------------------------------------
    dz_ph_rho = np.zeros ( k_bottom )
    dz_w      = np.zeros ( k_bottom + 1 )

    dep_w_levels = dep_w_levels( dep_rho_levels, file_out ) 

    for k in range(k_bottom) :
       dz_ph_rho[k] = - dep_w_levels[k+1] + dep_w_levels[k]    # note dz is negative

    dz_w[0] = - 2. *  dep_rho_levels[0]
    for k in range(1, k_bottom) : 
       dz_w[k] = - dep_rho_levels[k] + dep_rho_levels[k-1] 
    dz_w[k_bottom] = - 2. *  ( dep_w_levels[k_bottom] - dep_rho_levels[k_bottom-1] )

#----------------------------------------------------------------------------------------------------------
#  2. Calculate the Brunt Vaisala frequency (it should be positive !) 
#----------------------------------------------------------------------------------------------------------
    N_sq_dz_w = np.zeros ( k_bottom + 1 ) 

    for k in range(k_bottom - 1) :
       N_sq_dz_w[k+1] =  - ( gdef / rho_0 ) * ( pot_rho[k+1] - pot_rho[k] )    # N^2 dz_w  at w levels 

    print('k, pot_rho[k], dep_rho_levels[k], dep_levels[k+1], dz_ph_rho[k], dz_w[k], N_sq_dz_w[k] ', file=file_out)
    for k in range(k_bottom) :
       print(k, pot_rho[k], dep_rho_levels[k], dep_levels[k+1], dz_ph_rho[k], dz_w[k], N_sq_dz_w[k], file=file_out)

#----------------------------------------------------------------------------------------------------------
#  3. Calculate the full matrix A defined at all dep_w_levels (including the boundaries) down to k_bottom. 
#     This is NOT the matrix whose eigenvalues we calculate
#----------------------------------------------------------------------------------------------------------
    n_A = k_bottom + 1 
    A = np.zeros( (n_A,n_A) )  

#----------------------------------------------------------------------------------------------------------
#  3.1 Calculate the elements of A; it is dz_w[k] * d^2 h / dz^2, which is a symmetric matrix
#      We do not need to calculate A[0,1] or A[n_A-1,n_A-2]
 
#----------------------------------------------------------------------------------------------------------
    for k in range (n_A-1) : 
       A[k, k+1]   = 1.0 / dz_w[k]    
       A[k+1, k]   = A[k, k+1] 

    for k in range (1, n_A - 1):  
       A[k, k] = - A[k+1, k] - A[k, k+1]  
 
#----------------------------------------------------------------------------------------------------------
#  3.2 Calculate the arrays whose eigenvalues will be calculated 
#----------------------------------------------------------------------------------------------------------

    B = A[1:n_A-1,1:n_A-1]
    print(' B.shape = ', B.shape, file=file_out)

    C = np.zeros( (n_A-2,n_A-2) )
    for k in range(n_A-2) :
       C[k,k] = - N_sq_dz_w[k+1]
       print(' k, B[k,k], C[k,k] = ', k, B[k,k], C[k,k], file=file_out)

    for k in range(n_A - 3) :
       print(' k, B[k,k+1], B[k+1,k] = ', k, B[k,k+1], B[k+1,k], file=file_out)

    print(' A[1,0:3] = ', A[1,0:3], file=file_out)
    print(' B[0,0:2] = ', B[0,0:2], file=file_out)

    print(' A[n_A-2,n_A-3:n_A-1] = ', A[n_A-2,n_A-3:n_A-1], file=file_out)
    print(' B[n_A-3,n_A-4:n_A-2] = ', B[n_A-3,n_A-4:n_A-2], file=file_out)

    print(' C[0,0], C[1,1] = ', C[0,0], C[1,1], file=file_out)
    print(' C[n_A-3,n_A-3], C[n_A-2,n_A-2] = ', C[n_A-4,n_A-4], C[n_A-3,n_A-3], file=file_out)

#----------------------------------------------------------------------------------------------------------
#  5. Calculate and sort the eigenvalues; the eigenvalues are 1.0 / c_e^2 
#----------------------------------------------------------------------------------------------------------

    eigenvalues, vertical_modes = eig (B, C)

    print(' eigenvalues = ', eigenvalues, file=file_out)

    indices = np.argsort(eigenvalues)

#----------------------------------------------------------------------------------------------------------
#  6. Gather the eigenvalues and calculate the pressure modes 
#----------------------------------------------------------------------------------------------------------
    c_mode = np.zeros ( n_modes ) 
    p_mode = np.zeros ( (n_modes, k_bottom ) )

    h_mode = np.zeros ( k_bottom+1 )

    for n in range ( n_modes ): 
       ind = indices[n]
       eigen = eigenvalues[ind]
       print('ind, eigen = ', ind, eigen, file=file_out)
       c_mode[n] = 1.0 / math.sqrt(eigen)

       h_mode[1:k_bottom] = vertical_modes[:,ind]

       p_sum = 0.0 
       for k in range( k_bottom ) : 
          p_mode[n, k] = ( h_mode[k+1] - h_mode[k] ) / dz_ph_rho[k]
          p_sum = p_sum - dz_ph_rho[k] * p_mode[n,k] * p_mode[n,k]

       recip_amp = 1.0 / math.sqrt(p_sum) 

       p_mode[n] = - recip_amp * p_mode[n]  

    print(' c_mode = ', c_mode, file=file_out)
    print(' k, p_mode[:,k]', file=file_out)
    for k in range( k_bottom ) :
       print(k, p_mode[:,k], file=file_out)

#----------------------------------------------------------------------------------------------------------
    return c_mode, p_mode, dz_ph_rho
#----------------------------------------------------------------------------------------------------------
