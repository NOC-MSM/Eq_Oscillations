# -*- coding: iso-8859-1 -*-
'''
Function: test_vertical_eigenmodes
author:   Mike Bell   
Date:     May 2018
Description: checks the calculation of the vertical eigenmodes for given density field as a function of depth 
''' 
#          version 1:  15 May 2018 is first version
#
#-------------------------------------------------------------------------------------------------


def test_function_calls():

# note that z is height so is negative; dep is depth so is positive 

# inputs 

# dep_rho_levels[k]   - the depths of the density levels
# pot_rho[k]          - the potential densities at these levels
# n_modes             - the number of modes to return
# file_out            - the output unit for standard diagnostics

# returns

# eigen_out           - the eigenvalues of the first n_modes
# p_mode              - normalised pressure modes on the rho grid  

   def G(t): 
      return t*t

   c = 5.0 

   def function (t, x) : 
      y = [G(t) + x[0],  - c * x[1] ] 
      return y   

   for i in range (3) : 

      t = float(i) 
      x = [0, 1]

      v = function (t, x) 
      print x, v 

#----------------------------------------------------------------------------------------------------------
test_function_calls()