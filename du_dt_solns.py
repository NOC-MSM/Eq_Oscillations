# -*- coding: iso-8859-1 -*-
'''
Functions: Slow_du_dt_calc
author:   Mike Bell
Date:     June 2019 
Description: Build matrix of contributions to du/dt (useful for least squares solution)   
'''
#          version 1:  17 September 2018
#
#-------------------------------------------------------------------------------------------------

def plot_du_dt( mode_ts, c_mode, j_merid, plotfilestem, input_file_time_step, period_in_days, ttl, file_out ): 

    from Slow_plots import Slow_plots
    
    filename = plotfilestem+'_' + c_mode + '_' + str(j_merid) + '.pdf'    
    mean_sq_mode_ts = Slow_plots(filename, file_out, input_file_time_step, period_in_days, mode_ts, c_mode, ttl)

#-------------------------------------------------------------------------------------------------

def du_dt_calc(du_dt_cmpt, j_merid, beta_y_fac, v_mode, u_mode, force_x_mode, 
               plotfilestem, input_file_time_step, period_in_days, ttl, file_out):

#   Inputs
# du_dt_cmpt  contributions to du_dt for all meridional modes: [n_meridional_modes-1, 4, n_times_soln]
#             both input and output (because it is only slowly built up) 
# j_merid     meridional mode number of inputs below  
# beta_y_fac  beta y = beta_y_fac * y_tilde 
# v_mode      solution for v [n_times_soln] 
# force_x_mode  forcing in x direction [n_times_soln]
# 
# Following are for diagnostic outputs only
#
# plotfilestem           - stem of filename of plots
# input_file_time_step  
# period_in_days         - period of this mode to use for smoothing the time-series before plotting
#                          Note this can be set to zero in this routine to avoid doing this smoothing
# ttl                    - title  
# file_out
#
#    Returns
# u_merid_modes 

    import math
    import numpy as np 
    
    lev_diag = 10 
    
    n_merid_1, n_cmpts, n_times_soln = du_dt_cmpt.shape

# calculate contribution from force_x_mode
    if j_merid < n_merid_1 : 
       du_dt_cmpt[j_merid,0,:] = - u_mode[:]
       du_dt_cmpt[j_merid,1,:] = force_x_mode[:]
    
       if lev_diag > 5: 
          plot_du_dt( force_x_mode, 'force_x', j_merid, plotfilestem, input_file_time_step, period_in_days, ttl, file_out  )   
      

#  normalisation factor for mode j_merid is 1. / sqrt( math.factorial(j_merid) ) 
    norm_fac_j = 1.0 / math.sqrt( math.factorial(j_merid) ) 

# calculate contribution from  beta y v  into mode j_merid + 1    

# calculate contribution from  beta y v  into mode j_merid - 1    
    if j_merid > 0 : 
       fac_jm1 = math.sqrt( math.factorial(j_merid-1) ) * norm_fac_j
       v_force = float(j_merid) * beta_y_fac *  fac_jm1 * v_mode
       du_dt_cmpt[j_merid-1,2,:] =  v_force[:]    

       if lev_diag > 5: 
          plot_du_dt( v_force, 'v_minus', j_merid+1, plotfilestem, input_file_time_step, period_in_days, ttl, file_out )   

    if j_merid + 1 < n_merid_1 : 
       fac_jp1 = math.sqrt( math.factorial(j_merid+1) ) * norm_fac_j
       v_force = beta_y_fac * fac_jp1 * v_mode
       du_dt_cmpt[j_merid+1,3,:] =  v_force[:]

       if lev_diag > 5: 
          plot_du_dt( v_force, 'v_plus', j_merid-1, plotfilestem, input_file_time_step, period_in_days, ttl, file_out )   
    
    return du_dt_cmpt

#-------------------------------------------------------------------------------------------------

def du_dt_integrate(du_dt_ts, u_initial, Rayleigh_damp_fac, input_file_time_step, file_out):

#   Inputs
# du_dt_ts                  (du_dt + ru) for one mode: [n_times_soln]
# u_initial                 initial value for u  
# Rayleigh_damp_fac         Rayleigh damping factor  
# input_file_time_step      length of each time-step (in seconds) 
# file_out
#
#    Returns
# u_ts         integral 

   import numpy as np 
   
   n_times_soln = len(du_dt_ts)
   time = input_file_time_step * np.arange( n_times_soln ) 
   exp_time = np.exp( Rayleigh_damp_fac * time )
   source = du_dt_ts * exp_time * input_file_time_step   # this includes the time-step dt 
   
# use trapezoidal rule to solve [ exp(rt') u(t') ]_0^t = \int_0^t S(t') exp(rt') dt'

   u_integral = np.zeros ( n_times_soln )    
   u_integral[0] = u_initial
   for itime in range ( n_times_soln-1 ): 
      u_integral[itime+1] = u_integral[itime] + 0.5 * ( source[itime] + source[itime+1] ) 
   
   u_ts =  u_integral / exp_time
   
   return u_ts
   
#-------------------------------------------------------------------------------------------------

def u_plot(plotfilename, time_axis, u_sim_ts, u_sim_ts_0, u_true_ts, file_out):

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(15,4))
    plt.plot(time_axis, u_true_ts, color='k', linewidth=2, linestyle='-',  label = 'NEMO')
    plt.plot(time_axis, u_sim_ts,  color='r', linewidth=2, linestyle=':', label = 'sim')
    plt.plot(time_axis, u_sim_ts_0,color='b', linewidth=2, linestyle='--', label = 'lsf')
    
    plt.xlabel("Time (days)", fontsize=12)
    plt.ylabel(r'$ \tilde{u}_{m,n} \; (\mathrm{m}^2 \mathrm{s}^{-1})$', fontsize=12)
    plt.legend()

    plt.savefig(plotfilename, transparent=True, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')

    return
#-------------------------------------------------------------------------------------------------
