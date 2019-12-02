# -*- coding: iso-8859-1 -*-
'''
Function: test_time_step_solver
author:   Mike Bell   
Date:     May 2018
Description: Test the time-step solver for the wind-driven simulation of a particular meridional and vertical mode 
             and user specified segments starting from "correct" values of v and dv_dt      

''' 
#          version 1:  15 May 2018 is first version
#
#-------------------------------------------------------------------------------------------------

def test_time_step_solver():

    import sys
    import math 
    import numpy as np 
    from time_step_solver import time_step_solver 

#----------------------------------------------------------------------------------------------------------
# 0. Set user choices  
#----------------------------------------------------------------------------------------------------------

    date = '2018_05_15'
    expt_name  = 'test_tss'

    filestem = './std_out/'+date+'/'+expt_name
    file_out = open(filestem+'.txt','w')


# 0.1 Specify the time-series to generate for testing 
    n_times = 30          # number of points in the time-series
    time_step       = 4. * 3600.   # 4 hour time step 

    print >> file_out, ' n_times, time_step = ' , n_times, time_step

# 0.2 Specify the winds as Y = Y_c cos alpha t + Y_s sin alpha t
#                          F = F_c cos alpha t + F_s sin alpha t

    alpha = 0.02 / time_step 
    Y_c =   3.0 / alpha 
    Y_s =  -1.0 / alpha
    F_c =   2.0
    F_s =   3.0

    print >> file_out, ' alpha, Y_c, Y_s, F_c, F_s = ', alpha, Y_c, Y_s, F_c, F_s

# 0.3 Specify the response frequency and damping and the amplitudes of the general solutions
    omega   = 0.03 / time_step
    denom = omega * omega - alpha * alpha
  
    r_damp  = 0.0               # damping is not implemented 
    v_c_gen =  5.0 /  denom 
    v_s_gen =  3.0 /  denom 

    print >> file_out, ' omega, r_damp, v_c_gen, v_s_gen =  ', omega, r_damp, v_c_gen, v_s_gen

# 0.4 set the time segments to simulate  

    len_time_segments =  10      # (number of time-outputs in a segment) 6*10 is 10 days for 4 hourly time-outputs 
    first_segment     =   3      # start of first segment
    inc_segments      =   5      # increment between segments  
    max_no_segments   =   2      # (maximum) number of segments to generate for each meridional wavenumber
    len_halo_spline   =   3      # halo round time period used to make spline calculation accurate near the start and end of the period

    print >> file_out, ' len_time_segments = ', len_time_segments
    print >> file_out, ' first_segment     = ', first_segment
    print >> file_out, ' inc_segments     = ', inc_segments
    print >> file_out, ' max_no_segments   = ', max_no_segments
    print >> file_out, ' len_halo_spline = ', len_halo_spline  

    last_tstep_reqd =  first_segment + (max_no_segments-1)*inc_segments + len_time_segments + len_halo_spline
    if last_tstep_reqd > n_times : 
       print >> file_out, '  last_tstep_reqd > n_times; please avoid this choice ', last_tstep_reqd, n_times
       sys.exit()


#----------------------------------------------------------------------------------------------------------
# 1. Calculate the wind forcing and solutions   
#----------------------------------------------------------------------------------------------------------

    Y = np.zeros( n_times ) 
    F = np.zeros_like (Y)   

    for it in range ( n_times ) : 
       alpha_t =  alpha * time_step * it 
       Y [it] = Y_c * math.cos(alpha_t) + Y_s * math.sin(alpha_t)  
       F [it] = F_c * math.cos(alpha_t) + F_s * math.sin(alpha_t)  

    print >> file_out, ' Y = ', Y
    print >> file_out, ' F = ', F

#----------------------------------------------------------------------------------------------------------
# 1.1 Calculate the meridional velocity response
#----------------------------------------------------------------------------------------------------------
    v_part = np.zeros_like (Y)     # this is the particular solution 
    v_gen  = np.zeros_like (Y)     # this is the particular solution 


    v_c_part =   (alpha * Y_s - F_c ) / denom  
    v_s_part = - (alpha * Y_c + F_s ) / denom  

    for it in range ( n_times ) : 
       alpha_t = alpha * time_step * it 
       v_part[it] = v_c_part * math.cos(alpha_t) + v_s_part * math.sin(alpha_t)  
       omega_t = omega * time_step * it
       v_gen[it]  = v_c_gen * math.cos(omega_t)  + v_s_gen * math.sin(omega_t)  

    print >> file_out, ' v_part = ', v_part
    print >> file_out, ' v_gen  = ', v_gen

    v = v_part + v_gen
    print >> file_out, ' v  = ', v

#----------------------------------------------------------------------------------------------------------
# 2. Test whether time_step_solver produces matching results    
#----------------------------------------------------------------------------------------------------------

    v_sim_seq, v_seq = time_step_solver(v, F, Y, time_step, omega*omega, len_time_segments, first_segment, inc_segments, max_no_segments, len_halo_spline, file_out)

#--------------------------------------------------------------------------------------------------------
test_time_step_solver()