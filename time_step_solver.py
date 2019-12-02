# -*- coding: iso-8859-1 -*-
from __future__ import print_function
'''
Function: time_step_solver
author:   Mike Bell
Date:     October 2018
Description: Calculate the wind-driven simulation of a particular meridional and vertical mode
             for user specified segments starting from "correct" values of v and dv_dt
             ***updated October 2018 to include least squares solution***
             ***Version that includes boundary pressure effect***
'''
#          version 1:  15 May 2018 is first version
#
#-------------------------------------------------------------------------------------------------

def fun_ode ( vector, t, tck_tauy, tck_f_taux, omega_sq, Rayleigh_damp_coef, const ) :
   import scipy.interpolate as interpol
   v = vector[0]
   q = vector[1]
   tauy_temp   = interpol.splev(t, tck_tauy, der = 0)
   f_taux_temp = interpol.splev(t, tck_f_taux, der = 0)
   dv_dt =  q + tauy_temp
   dq_dt =  - 2.0*Rayleigh_damp_coef*dv_dt - omega_sq*v - f_taux_temp  + const
   return [ dv_dt, dq_dt ]

#-------------------------------------------------------------------------------------------------

def best_fit_ics ( t_spline, v_true_seg, t_seg, tck_tauy, tck_f_forcex, k_spline, vector_initial,
                   omega_sq, Rayleigh_damp_fac, const, file_out ) : 

# Purpose: find the best fit initial conditions for a given choice of parameters and forcing 
# Mike Bell  6 February 2019  

   import numpy as np
   import numpy.linalg as linalg
   import scipy.integrate as integ
   import scipy.interpolate as interpol

   lev_diag_user = 3     #    > 5 and > 10 are used as tests for diagnostic print-outs below

# 0. Dimension the least squares matrix and vector 

   len_time_segments = len (t_seg) 
   A = np.zeros ( (len_time_segments, 2) )      
   b = np.zeros ( len_time_segments )

# 1. Find the particular solution, sol_p 

   sol_p = integ.odeint( fun_ode, vector_initial, t_seg, args=(tck_tauy, tck_f_forcex, omega_sq, Rayleigh_damp_fac, const) )  

   b = v_true_seg - sol_p[:,0] 

# 2. Set up a spline solution for zero forcing and a zero value for the constant forcing 

   zero_force = np.zeros_like (t_spline) 
   tck_zero   = interpol.splrep(t_spline, zero_force, k = k_spline)
   zero = 0.0 

# 3. Calculate the A and B homogeneous solutions and load them into matrix A ; these solutions have no forcing  

   vector_A = [ 1., 0. ]
   sol_A = integ.odeint( fun_ode, vector_A, t_seg, args=(tck_zero, tck_zero, omega_sq, Rayleigh_damp_fac, zero) )  
   A[:,0] = sol_A[:,0]

   vector_B = [ 0., 1. ]
   sol_B = integ.odeint( fun_ode, vector_B, t_seg, args=(tck_zero, tck_zero, omega_sq, Rayleigh_damp_fac, zero) )  
   A[:,1] = sol_B[:,0] 

   if lev_diag_user > 10 : 
      print ( ' omega_sq, Rayleigh_damp_fac, const = ', omega_sq, Rayleigh_damp_fac, const, file = file_out ) 
      print ( ' ipt, sol_p[ipt, 0], sol_A[ipt,0], sol_B[ipt,0] ', file = file_out ) 
      for ipt in range ( 0, len_time_segments, 10 ) : 
         print ( ipt, sol_p[ipt, 0], sol_A[ipt,0], sol_B[ipt,0], file = file_out ) 

# 4. Calculate the least squares solution and return it 

   x, sum_res, rank, sing_vals = linalg.lstsq( A, b, rcond = None )      

   if lev_diag_user > 5 : print ( 'solution x, sum_res, rank, sing_vals = ', x, sum_res, rank, sing_vals, file = file_out )

   best_fit_sol = sol_p[:,0]  + x[0] * sol_A[:,0]  + x[1] * sol_B[:,0] 

   return best_fit_sol

#-------------------------------------------------------------------------------------------------

def mean_std( a_seg ) : 

# return the mean and standard deviation of the input numpy array (a_seg) 

   import numpy as np
   import math
   import sys 

   mean_a = np.mean (a_seg ) 
   mean_a_sq = np.mean (a_seg*a_seg )

   var_a = mean_a_sq - mean_a * mean_a 

   if var_a > 0.0 : 
      std_a = math.sqrt(var_a) 
   elif var_a < -1.0E-10 :  
      print ( ' mean_std is taking square root of a significant negative value ', var_a ) 
      sys.exit()
   else : 
      std_a = 0.0 

   return mean_a, std_a 

#-------------------------------------------------------------------------------------------------

def time_step_solver(u_ts, v_ts, f_taux_ts, tauy_ts, f_press_w_e_ts,  
                     input_file_time_step, omega_sq, Rayleigh_damp_fac, 
                     len_time_segments, first_segment, inc_segments, max_no_segments, 
                     len_halo_spline, k_spline, l_Rayleigh, l_wind_x, l_wind_y, l_press, l_const, l_force_fix,
                     l_best_fit_initial_cond, file_out):

# Inputs
#
#    time-series (_ts) of zonally integrated projections onto this meridional and vertical mode
# u_ts                      - actual  zonal velocities (only used to output du_dt_sgmts
# v_ts                      - actual  meridional velocities  
# f_taux_ts                 - zonal wind stresses (multiplied by f) 
# tauy_ts                   - meridional wind stresses
# f_press_w_e_ts            - pressure on western boundaries minus pressure on eastern boundaries (multiplied by f)
#
# input_file_time_step      - 
# omega_sq                  - natural frequency of this mode 
# Rayleigh_damp_fac         - Rayleigh friction coefficient 

#                            User controlled definition of time segments
# len_time_segments         - number of time-inputs/outputs in a segment  
# first_segment             - start of first segment (integer) 
# inc_segments              - increment between segments
# max_no_segments           - (maximum) number of segments to generate for each meridional wavenumber
# len_halo_spline           - halo round time period used to make spline calculation accurate near the start and end of the period
# l_Rayleigh                - True implies include Rayleigh damping in the least squares solution
# l_wind_x                   -       "             zonal wind forcing                     "
# l_wind_y                   -       "             merid wind forcing                     "
# l_press                    -       "             pressure forces on E & W bdys          " 
# l_const                    -       "              a constant value                       " 
# l_force_fix                -       "             use fixed forcing in the "fit" solution depending on l_wind_x, l_wind_y, l_press
# l_best_fit_initial_cond   - True implies finds best fit initial condition; False uses the NEMO model state as initial condition 
# file_out                  - diagnostics output file 

# Returns
# 
# v_sim_sgmts                 - simulated     sets of v time-series
# v_fit_sgmts                 - simulated     sets of v time-series using best fit solution
# v_sgmts                     - corresponding sets of time-series generated from the actual meridional velocities
# t_sgmts                     - corresponding time for plotting
# du_dt_sgmts                 - du_dt calculated from the actual zonal velocities 
# dtauy_dt_ts                 - time series of dtauy_dt_ts (padded with zeroes at start and end of time-series) 

   import math
   import numpy as np
   import numpy.linalg as linalg
   import scipy.integrate as integ
   import scipy.interpolate as interpol

   lev_diag_user = 2     #    > 5 and > 10 are used as tests for diagnostic print-outs below

#----------------------------------------------------------------------------------------------------------
# 0. initialise outputs and start loop over segments
#----------------------------------------------------------------------------------------------------------

   v_sim_sgmts     = np.zeros ( ( max_no_segments,  len_time_segments ) )  # simulated segments of v time-series
   v_fit_sgmts     = np.zeros ( ( max_no_segments,  len_time_segments ) )  # simulated segments of v time-series using best fit solution 
   v_sgmts         = np.zeros ( ( max_no_segments,  len_time_segments ) )  # corresponding segments of NEMO v time-series
   t_sgmts         = np.zeros ( ( max_no_segments,  len_time_segments ) )  # corresponding time for plotting
   du_dt_sgmts     = np.zeros ( ( max_no_segments,  len_time_segments ) )  # du/dt segments calculated from NEMO u time-series 

   mean_b_seg     = np.zeros ( max_no_segments )
   mean_b_sq_seg  = np.zeros ( max_no_segments )
   res_fit_seg    = np.zeros ( max_no_segments )   # sum of square residuals for least squares fit to d^2 v / dt^2
   res_sim_seg    = np.zeros ( max_no_segments )   # sum of square residuals for original dynamics based simulation

   om_sq_fit_seg = np.zeros ( max_no_segments )
   om_sq_sim_seg = np.zeros ( max_no_segments )

   Ray_fit_seg = np.zeros ( max_no_segments )
   Ray_sim_seg = np.zeros ( max_no_segments )

   tau_x_fit_seg = np.zeros ( max_no_segments )
   tau_x_sim_seg = np.zeros ( max_no_segments )

   tau_y_fit_seg = np.zeros ( max_no_segments )
   tau_y_sim_seg = np.zeros ( max_no_segments )

   press_w_e_fit_seg = np.zeros ( max_no_segments )
   press_w_e_sim_seg = np.zeros ( max_no_segments )

   const_fit_seg = np.zeros ( max_no_segments )

   mean_orig_seg    = np.zeros( max_no_segments )
   mean_sq_orig_seg = np.zeros( max_no_segments )
   mse_sim_seg      = np.zeros( max_no_segments )
   mse_fit_seg      = np.zeros( max_no_segments )
   
   dtauy_dt_ts = np.zeros_like (v_ts)

   for iseg in range ( max_no_segments ) :

#----------------------------------------------------------------------------------------------------------
# 1. Specify segments of the time-series and calculate spline fits to segments extended by halos
#----------------------------------------------------------------------------------------------------------

      it_seg_start = first_segment + iseg * inc_segments
      it_seg_end   = it_seg_start + len_time_segments     # the segment end is it_seg_end-1 (the last point is not included)

      it_spline_start = it_seg_start - len_halo_spline
      it_spline_end   = it_seg_end   + len_halo_spline

      it_seg = np.arange( it_seg_start, it_seg_end )
      t_seg = input_file_time_step * it_seg     # "actual" times
      it_spline = np.arange( it_spline_start, it_spline_end )

# calculate the segments for the spline fits 
      t_spline      = input_file_time_step * it_spline     # "actual" times
      u_spline      = u_ts[it_spline_start:it_spline_end]
      v_spline      = v_ts[it_spline_start:it_spline_end]
      f_taux_spline = f_taux_ts[it_spline_start:it_spline_end]
      f_pr_w_e_spline = f_press_w_e_ts[it_spline_start:it_spline_end]
      tauy_spline   = tauy_ts[it_spline_start:it_spline_end]
        
# debugging
      if lev_diag_user > 10 :
         print(' len(t_spline), len(v_spline) = ', len(t_spline), len(v_spline))#, file=file_out)
         print('t_spline = ', t_spline)#, file=file_out)
         print('v_spline = ', v_spline)#, file=file_out)
         print('f_taux_spline = ', f_taux_spline)#, file=file_out)
         print("tauy_ts=",tauy_ts)

# calculate the spline fits themselves
      tck_u      = interpol.splrep(t_spline, u_spline, k = k_spline)
      tck_v      = interpol.splrep(t_spline, v_spline, k = k_spline)
      tck_f_taux = interpol.splrep(t_spline, f_taux_spline, k = k_spline)
      tck_f_pr_w_e = interpol.splrep(t_spline, f_pr_w_e_spline, k = k_spline)
      tck_tauy   = interpol.splrep(t_spline, tauy_spline, k = k_spline)

# calculate the values on the time segments 
      du_dt_seg    = interpol.splev(t_seg, tck_u, der = 1)
      v_seg        = interpol.splev(t_seg, tck_v, der = 0)
      dv_dt_seg    = interpol.splev(t_seg, tck_v, der = 1)
      f_taux_seg   = interpol.splev(t_seg, tck_f_taux, der = 0)
      f_pr_w_e_seg = interpol.splev(t_seg, tck_f_pr_w_e, der = 0)
      tauy_seg     = interpol.splev(t_seg, tck_tauy, der = 0)
      dtauy_dt_seg = interpol.splev(t_seg, tck_tauy, der = 1)
      
      dtauy_dt_ts[it_seg_start: it_seg_end] = dtauy_dt_seg

      if lev_diag_user > 10 :
         print('it, time, v_time, f_taux_time, tauy_time', file=file_out)
         for it in range ( 0, len(it_seg), 10  ):
            print(it, t_seg[it], v_seg[it], f_taux_seg[it], tauy_seg[it] , file=file_out)

#----------------------------------------------------------------------------------------------------------
# 2. Calculate a least squares solution for the segment  minimise |(A x - b)|^2
#----------------------------------------------------------------------------------------------------------

# 2a. Calculate A and b 

      dim_A = 1
      if l_Rayleigh : 
        dim_A = dim_A + 1
      if not(l_force_fix) : 
        if l_wind_x : 
          dim_A = dim_A + 1
        if l_wind_y : 
          dim_A = dim_A + 1
        if l_press : 
          dim_A = dim_A + 1
      if l_const : 
        dim_A = dim_A + 1

        
      A = np.zeros ( (len_time_segments, dim_A) )      
      b = np.zeros ( len_time_segments )

      b = interpol.splev(t_seg, tck_v, der = 2)          # b is d2v_dt2_seg

      indx = 0 
      if l_Rayleigh : 
        A[:,indx]   = - interpol.splev(t_seg, tck_v, der = 1)      #  for Rayleigh friction 
        indx = indx+1

      A[:,indx] = - interpol.splev(t_seg, tck_v, der = 0)      #  - v_seg (this is the - omega_sq v term)  
      indx = indx+1

      if l_force_fix : 
         if l_wind_x : 
            b = b + interpol.splev(t_seg, tck_f_taux, der = 0)          # + f taux       
      else : 
         if l_wind_x : 
            A[:,indx] = - interpol.splev(t_seg, tck_f_taux, der = 0)          # - f taux 
            indx = indx+1

      if l_force_fix : 
         if l_wind_y : 
            b = b - interpol.splev(t_seg, tck_tauy, der = 1)              # - d tauy / dt 
      else : 
         if l_wind_y : 
            A[:,indx] = interpol.splev(t_seg, tck_tauy, der = 1)              # d tauy / dt 
            indx = indx+1

      if l_force_fix : 
         if l_press : 
            b = b + interpol.splev(t_seg, tck_f_pr_w_e, der = 0)        # pressure on west minus east; +ve sign  
      else : 
         if l_press : 
            A[:,indx] = - interpol.splev(t_seg, tck_f_pr_w_e, der = 0)        # pressure on west minus east; -ve sign  
            indx = indx+1

      if l_const : 
         A[:,indx] = 1.0       
         indx = indx+1 
        
      if lev_diag_user > 5 : 
         print( ' indx, [ A[k_seg,indx] for k_seg in range (len_time_segments) ] ' , file=file_out)
         for indx in range (dim_A) : 
            print(indx, [ A[k_seg,indx] for k_seg in range (len_time_segments) ]  , file=file_out)
         print( 'b = ', b, file=file_out)

# 2b. Find x such that  | A x - b |^2 is a minimum   

      x, sum_res, rank, sing_vals = linalg.lstsq( A, b, rcond = None )      

# 2c. Calculate sum of square residuals for 0) d^2 v / dt^2 1) least squares fit and 2) original dynamics based simulation  

# 0) d^2 v / dt^2  ; the naming of mean_b_seg and mean_b_sq_seg is a bit out-of-date; mean_b_orig_seg would be better

      b_orig = interpol.splev(t_seg, tck_v, der = 2)          # b_orig is d2v_dt2_seg
      mean_b_seg[iseg]   = np.mean(b_orig) 
      mean_b_sq_seg[iseg] = np.mean (b_orig*b_orig)   

# 1) least squares fit and 2) original dynamics based simulation  

      b_fit = np.zeros_like(b) 
      for indx in range ( dim_A ) :  
         b_fit[:] = b_fit[:] + A[:,indx] * x[indx] 

      diffs = (b - b_fit)*(b - b_fit)
      res_fit_seg[iseg] =  np.mean ( diffs )

# 2) original dynamics based simulation

      b_sim = - omega_sq * v_seg - 2.0 * Rayleigh_damp_fac * dv_dt_seg

      if l_wind_x :
         b_sim = b_sim - f_taux_seg 

      if l_wind_y :
         b_sim = b_sim + dtauy_dt_seg

      if l_press : 
         b_sim = b_sim - f_pr_w_e_seg 

      diffs = (b_orig - b_sim)*(b_orig - b_sim)
      res_sim_seg[iseg] =  np.mean ( diffs )
  
# 2d. Unpack x into more easily understood  terms  
      if lev_diag_user > 5 :
        print( 'iseg = ', iseg, file=file_out)
        print( 'x, sum_res = ', x, sum_res, file=file_out)
        print( 'av_res = ', math.sqrt( sum_res / (len_time_segments - dim_A) )  , file=file_out)

      indx = 0 
      if l_Rayleigh : 
         Rayleigh_fit = 0.5 * x[indx] 
         if lev_diag_user > 3 : print( 'Rayleigh_fit (units: per sec) = ', Rayleigh_fit , file=file_out)
         indx = indx+1
      else : 
         Rayleigh_fit = 0.0 

      Ray_fit_seg[iseg] = Rayleigh_fit

      om_sq_fit = x[indx] 
      if lev_diag_user > 3 : print( ' omega_sq, om_sq_fit = ', omega_sq, om_sq_fit, file=file_out)
      indx = indx+1
      om_sq_fit_seg[iseg] = om_sq_fit

      if l_force_fix : 
         if l_wind_x : 
            tau_x_fit = 1.0 
         else : 
            tau_x_fit = 0.0	 
      else :
         if l_wind_x : 
            tau_x_fit = x[indx]
            indx = indx+1
         else :
            tau_x_fit = 0.0
      if lev_diag_user > 3 : print( 'factor for zonal winds = ', tau_x_fit, file=file_out)
      tau_x_fit_seg[iseg] = tau_x_fit

      if l_force_fix : 
         if l_wind_y : 
            tau_y_fit = 1.0 
         else : 
            tau_y_fit = 0.0	 
      else :
         if l_wind_y : 
            tau_y_fit = x[indx] 
            indx = indx+1
         else :
            tau_y_fit = 0.0 
      if lev_diag_user > 3 : print( 'factor for meridional winds = ', tau_y_fit , file=file_out)
      tau_y_fit_seg[iseg] = tau_y_fit

      if l_force_fix : 
         if l_press : 
            press_w_e_fit = 1.0 
         else : 
            press_w_e_fit = 0.0
      else : 
         if l_press : 
            press_w_e_fit = x[indx] 
            indx = indx+1
         else: 
            press_w_e_fit = 0.0
      if lev_diag_user > 3 : print( ' factor for pressure on western boundary = ', press_w_fit, file=file_out)
      press_w_e_fit_seg[iseg] = press_w_e_fit 

      if l_const : 
         const_fit = x[indx]        
         if lev_diag_user > 3 : print( ' constant offset  = ', const_fit, file=file_out)
         indx = indx+1
      else: 
         const_fit = 0.0
      const_fit_seg[iseg] = const_fit

#----------------------------------------------------------------------------------------------------------
# 3. Integrate the pair of first order ODEs over the time segment and save the resulting time-series 
#----------------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------------
# 3a. Calculate the initial conditions for the segment
#----------------------------------------------------------------------------------------------------------

      t_initial = t_seg[0]
      v_initial = interpol.splev(t_initial, tck_v, der = 0)
      dv_dt_initial = interpol.splev(t_initial, tck_v, der = 1)

      q_initial = dv_dt_initial
      if l_wind_y : 
         q_initial = q_initial - interpol.splev(t_initial, tck_tauy, der = 0)

      vector_initial = [ v_initial, q_initial ]

# diagnostics
      tauy_initial   = interpol.splev(t_initial, tck_tauy, der = 0)
      f_taux_initial = interpol.splev(t_initial, tck_f_taux, der = 0)

      if lev_diag_user > 5 :
         print('iseg = ', iseg, file=file_out)
         print('t_initial, v_initial = ', t_initial, v_initial, file=file_out)
         print('dv_dt_initial, q_initial, tauy_initial = ', dv_dt_initial, q_initial, tauy_initial, file=file_out)
         print('- omega_sq*v_initial, f_taux_initial = ', - omega_sq*v_initial, f_taux_initial, file=file_out)


      t_sgmts[iseg]     = t_seg # *4 because time increments are 4 hours
      du_dt_sgmts[iseg] = du_dt_seg
      v_sgmts[iseg]     = v_seg    # this is the NEMO model solution for this mode in this segment 

      mean_orig_seg[iseg]    = np.mean( v_seg )     
      mean_sq_orig_seg[iseg] = np.mean( v_seg*v_seg )   
#----------------------------------------------------------------------------------------------------------
# 3b. Generate time-series using the original equations  - pressure forces are NO LONGER commented out
#----------------------------------------------------------------------------------------------------------
      const = 0.0
 
      om_sq_sim_seg[iseg] = omega_sq
      Ray_sim_seg[iseg] = Rayleigh_damp_fac

      f_forcex_spline = np.zeros_like(f_taux_spline)
      tau_x_sim_seg[iseg] = 0.0
      press_w_e_sim_seg[iseg] = 0.0

      if l_wind_x : 
         f_forcex_spline = f_forcex_spline + f_taux_spline
         tau_x_sim_seg[iseg] = 1.0
      if l_press : 
         f_forcex_spline = f_forcex_spline + f_pr_w_e_spline 
         press_w_e_sim_seg[iseg] = 1.0
      tck_f_forcex = interpol.splrep(t_spline, f_forcex_spline, k = k_spline)

      forcey_spline = np.zeros_like(tauy_spline)
      tau_y_sim_seg[iseg] = 0.0
      if l_wind_y : 
         forcey_spline = forcey_spline + tauy_spline
         tau_y_sim_seg[iseg] = 1.0
      tck_forcey = interpol.splrep(t_spline, forcey_spline, k = k_spline)

      if l_best_fit_initial_cond : 
#         print ( 'sim_sgmts , omega_sq, iseg = ', omega_sq, iseg, file = file_out ) 
         v_sim_seg = best_fit_ics ( t_spline, v_sgmts[iseg], t_seg, tck_forcey, tck_f_forcex, k_spline, vector_initial, 
                                    omega_sq, Rayleigh_damp_fac, const, file_out)
      else : 
         sol = integ.odeint( fun_ode, vector_initial, t_seg, args=(tck_forcey, tck_f_forcex, omega_sq, Rayleigh_damp_fac, const) )  
         v_sim_seg = sol[:,0]

      v_sim_sgmts[iseg] = v_sim_seg
      v_res_seg = v_sim_seg - v_seg
#      mse_sim_seg = np.mean(v_res_seg*v_res_seg)      # corrected 5 April 2019 
      mse_sim_seg[iseg] = np.mean(v_res_seg*v_res_seg) 

         
#----------------------------------------------------------------------------------------------------------
# 3c. Generate time-series using the best fit solution
#----------------------------------------------------------------------------------------------------------

      f_forcex_spline_fit = tau_x_fit * f_taux_spline 
      if l_press : 
         f_forcex_spline_fit = f_forcex_spline_fit + press_w_e_fit * f_pr_w_e_spline

      tck_f_forcex_fit = interpol.splrep(t_spline, f_forcex_spline_fit, k = k_spline)

      tck_y_fit   = interpol.splrep(t_spline, tau_y_fit * tauy_spline, k = k_spline)

      if l_best_fit_initial_cond : 
#         print ( 'sim_fit , omega_sq, iseg = ', omega_sq, iseg, file = file_out ) 
         v_fit_seg = best_fit_ics ( t_spline, v_sgmts[iseg], t_seg, tck_y_fit, tck_f_forcex_fit, k_spline, vector_initial,
                                    om_sq_fit, Rayleigh_fit, const_fit, file_out ) 
      else : 
         sol_fit = integ.odeint( fun_ode, vector_initial, t_seg, args=(tck_y_fit, tck_f_forcex_fit, om_sq_fit, Rayleigh_fit, const_fit) )
         v_fit_seg = sol_fit[:,0]
 
      v_fit_sgmts[iseg] = v_fit_seg       
      v_res_seg = v_fit_seg - v_seg
#      mse_fit_seg = np.mean(v_res_seg*v_res_seg)    # corrected 5 April 2019
      mse_fit_seg[iseg] = np.mean(v_res_seg*v_res_seg) 

#----------------------------------------------------------------------------------------------------------
# 4. Print out results summary 
#----------------------------------------------------------------------------------------------------------
      if lev_diag_user > 3 :
         print('v_sim_sgmts = ', [ v_sim_sgmts[iseg,it] for it in range (len_time_segments) ], file=file_out)
         print('v_sgmts =     ', [ v_sgmts[iseg,it] for it in range (len_time_segments) ], file=file_out)
         print(' error =    ', [ (v_sim_sgmts[iseg,it] - v_sgmts[iseg,it]) for it in range (len_time_segments) ], file=file_out)

   if lev_diag_user > 1 :

      mean_b = np.mean ( mean_b_seg ) 
      mean_b_sq = np.mean ( mean_b_sq_seg ) 
      var_b  = mean_b_sq - mean_b * mean_b
      res_sim_norm = np.mean ( res_sim_seg ) / var_b 
      res_fit_norm = np.mean ( res_fit_seg ) / var_b 

      mean_orig    =  np.mean ( mean_orig_seg ) 
      mean_sq_orig =  np.mean ( mean_sq_orig_seg ) 
      var_orig = mean_sq_orig - mean_orig * mean_orig
      mse_sim_norm = np.mean ( mse_sim_seg ) / var_orig
      mse_fit_norm = np.mean ( mse_fit_seg ) / var_orig
       
      mean_om_sq_sim, std_om_sq_sim = mean_std( om_sq_sim_seg )
      mean_om_sq_fit, std_om_sq_fit = mean_std( om_sq_fit_seg )

      mean_Ray_sim, std_Ray_sim = mean_std( Ray_sim_seg )
      mean_Ray_fit, std_Ray_fit = mean_std( Ray_fit_seg )

      mean_tau_x_sim, std_tau_x_sim = mean_std( tau_x_sim_seg )
      mean_tau_x_fit, std_tau_x_fit = mean_std( tau_x_fit_seg )

      mean_tau_y_sim, std_tau_y_sim = mean_std( tau_y_sim_seg )
      mean_tau_y_fit, std_tau_y_fit = mean_std( tau_y_fit_seg )

      mean_press_w_e_sim, std_press_w_e_sim = mean_std( press_w_e_sim_seg )
      mean_press_w_e_fit, std_press_w_e_fit = mean_std( press_w_e_fit_seg )

      mean_const_fit, std_const_fit = mean_std( const_fit_seg ) 
      mean_const_norm = mean_const_fit / var_b
      std_const_norm  = std_const_fit / var_b

# print changed to include var_orig on 10 April 2019 
#      print( 'res_sim_norm, mse_sim_norm, mean_om_sq_sim, std_om_sq_sim, mean_tau_x_sim, std_tau_x_sim, mean_tau_y_sim, std_tau_y_sim, mean_Ray_sim, std_Ray_sim, var_b ', file=file_out)
#      print( res_sim_norm, mse_sim_norm, mean_om_sq_sim, std_om_sq_sim, mean_tau_x_sim, std_tau_x_sim, mean_tau_y_sim, std_tau_y_sim, mean_Ray_sim, std_Ray_sim, var_b, file=file_out)

      print( 'res_sim_norm, mse_sim_norm, mean_om_sq_sim, std_om_sq_sim, mean_tau_x_sim, std_tau_x_sim, mean_tau_y_sim, std_tau_y_sim, mean_Ray_sim, std_Ray_sim, var_b, var_orig ', file=file_out)
      print( res_sim_norm, mse_sim_norm, mean_om_sq_sim, std_om_sq_sim, mean_tau_x_sim, std_tau_x_sim, mean_tau_y_sim, std_tau_y_sim, mean_Ray_sim, std_Ray_sim, var_b, var_orig, file=file_out)
      print( 'mean_press_w_e_sim, std_press_w_e_sim', file=file_out)
      print( mean_press_w_e_sim, std_press_w_e_sim, file=file_out)

      print( 'res_fit_norm, mse_fit_norm, mean_om_sq_fit, std_om_sq_fit, mean_tau_x_fit, std_tau_x_fit, mean_tau_y_fit, std_tau_y_fit, mean_Ray_fit, std_Ray_fit, mean_const_norm, std_const_norm ', file=file_out)
      print( res_fit_norm, mse_fit_norm, mean_om_sq_fit, std_om_sq_fit, mean_tau_x_fit, std_tau_x_fit, mean_tau_y_fit, std_tau_y_fit, mean_Ray_fit, std_Ray_fit, mean_const_norm, std_const_norm, file=file_out)
      print( 'mean_press_w_e_fit, std_press_w_e_fit', file=file_out)
      print( mean_press_w_e_fit, std_press_w_e_fit, file=file_out)

   return v_sim_sgmts, v_fit_sgmts, v_sgmts, t_sgmts, du_dt_sgmts, dtauy_dt_ts
#--------------------------------------------------------------------------------------------------------
