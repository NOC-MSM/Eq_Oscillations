# -*- coding: iso-8859-1 -*-
from __future__ import print_function
'''
Function: Eq_Osc_Wind
author:   Mike Bell
Date:     May 2018
Description: calculates time-series of zonal mean equatorial waves ( inertia-gravity, Yanai, Kelvin and Rossby) in the low-order baroclinic
             vertical modes (modes 1, 2 and 3) driven by the local wind forcing
             ***Version that includes boundary pressure effect***
'''
#          version 1:  15 May 2018 is first version
#          version 2:  22 July 2019 -  normalisation of depth introduced so that amplitudes are properly dimensional 
#
#-------------------------------------------------------------------------------------------------

def find_nearest(array,value):
    idx,val = min(enumerate(array), key=lambda x: abs(x[1]-value))
    return idx

def Eq_Osc_Wind(region):

# Inputs
#    region     string: e.g. 'PAC'  options are: 'PAC', 'PAC_E', 'PAC_W', 'ATL', 'IND', 'IND_E', 'IND_W'
#
# Returns: null

    import sys
    import math as math
    import numpy as np
    import numpy.fft as fft
    import numpy.linalg as linalg
    import scipy.integrate as integ
    import scipy.interpolate as interpol
    import matplotlib.pyplot as plt
    from read_fields import read_one_field, read_three_fields
    from constants import deg2rad, rad2deg, earth_radius, gdef, OMEGA, rho_0
    from set_pot_dens_analytical import set_pot_dens_analytical
    from vertical_eigenmodes_shooting import vertical_eigenmodes_shooting
    from vertical_projection import vertical_projection
    from mode_plots import vmode_plots, mmode_plots
    from Osc_plots import Osc_plots, Osc_plots2
    from Slow_plots import Slow_plots
    from Force_plots import Force_plots 
    from du_dt_solns import du_dt_calc, du_dt_integrate, u_plot   
    from meridional_eigenmode import meridional_eigenmode
    from time_step_solver import time_step_solver
    from Calc_moc import calc_moc, rmodes 
    from Save_moc import save_moc
    from Save_modes import save_modes
    from evaluate import evaluate

    region_lc = region.lower()   # converts to lower case

#----------------------------------------------------------------------------------------------------------
# 0. Set user choices
#----------------------------------------------------------------------------------------------------------

    lev_diag_user = 6    # >5 and >10 are used as tests for diagnostic print-outs below. 
                          # The level can be set by the user in each function independently

    # Control path for output files
    date = '2020_08_05/'+ region_lc
    sig_type          = '04_sig2'
    k_bottom_w        = 65   # w-level of the bottom (k = 0 is at z = 0; the bottom is at w-level k)  level 65 is at 3995 metres in the 75 level grid
                             # the number of rho depth levels is k_bottom_w; the number of w-levels is k_bottom_w + 1
                             # 65 is the usual choice ; 73 is a sensible alternative; 75 is not great because some values are missing   

    # Set expt_name to be used for output files
    expt_name  = region_lc + '_' + sig_type + '_' + str(k_bottom_w) + '_Rf_Xt_Yt_Pt_Of_Bt_Ft_Ct_Sf_Apr_Eq_strat_ls_ics_22_100_segs'

    # Set the options for the least squares solver
    l_Rayleigh        = False    #  True implies include Rayleigh damping        in the least squares solution      (R) 
    l_wind_x          = True     #        "              zonal wind forcing                     "                   (X)
    l_wind_y          = True     #        "              meridional wind forcing                "                   (Y)   
    l_press           = True     #        "              pressure forcing (on east & west boundaries)               (P)
    l_press_orig      = False    #        "      use original (approximate) calculation of pressure differences     (O)
    l_press_full      = True     #        "      use full pressure (with barotropic) in new formulation             (B) 
    l_force_fix       = True     #        "      forces are fixed in the least squares fit solution (depending on l_wind_x, l_wind_y, l_press) (F)
    l_const           = True     #        "      include a constant value in the least squares solution             (C)
    l_surface_mode    = False    #  True gives surface modes; False uses normal modes with w=0 at the lower boundary (S)

# 0.1 specify the time-series data sets to input: pressures on e & w boundaries, zonally integrated horizontal velocities and zonally integrated surface wind stresses
    # The fields are assumed to have the shapes:  press_w, press_e, u, v: [it, kz, jy] and taux, tauy[it, jy]     
    # Specify the directory for the input files
    dir_data_in = '/noc/msm/scratch/meso-clip/atb299/'
#   dir_data_in = '/data/users/frmk/AdamBlaker201810/'
    time_series_input_file = dir_data_in + 'VN206HF_' + region_lc + '_458_538_bdy_press_tot.nc'   # standard case
    input_file_time_step =  4. * 3600.    # 4 hourly timesteps ; this information could be read from each of the files
 
# 0.1.1 specify the vertical density profile to use to calculate the vertical modes
    # This is a single vertical density profile rho[kz]
    ll_vertical_density_profile_analytical = False    # if False uses a data input file to set the density profile

    if ll_vertical_density_profile_analytical :
       dens_type       = 'exponential'    # options 'linear', 'exponential' (use lower case)
       dens_diff       = -6.0        #   density at top minus density at the bottom (level k_bottom_w)
       exp_scale_depth = 500.   # depth scale in metres for exponential case

    else :
       # 200610 or 200604
       vertical_density_profile_input_file = dir_data_in + 'NEMO_HF2/PDens/VN206HF_4h_2006' + sig_type +'_' + region_lc + '.nc'
       j_lat_row_min     = 30        #  min latitude of row of densities (between 0 and 99) 
       j_lat_row_max     = 71        #  max latitude of row of densities (between 0 and 99) 

# 0.2 specify the number of vertical modes to calculate and the depth of the wind penetration to use as inputs;
    n_vertical_modes     = 6   
    MLD                  = 50.       # depth of wind penetration in metres (this is described as the MLD, mixed layer depth, in the documentation)  
                                     # results should be relatively insensitive to MLD if MLD <= 50.

# 0.3 set the list of meridional wavenumbers to calculate
    n_meridional_modes   = 6    

# 0.4 set the Rayleigh damping coefficient 
    if l_Rayleigh : 
      Rayleigh_damp_fac = 1.E-7      # units are per second. So exp(-t/tau) with tau = 10 days uses 1.0/864000 which is about 1.e-6.   
    else :
       Rayleigh_damp_fac = 0.0

# 0.5 set the time segments to simulate
    len_time_segments    = 100       # (number of time-outputs in a segment) 6*10 is 10 days for 4 hourly time-outputs 
    first_segment        =   3       # start of first segment
    inc_segments         = 100       # increment between segments (hours)
    max_no_segments      =  22       # (maximum) number of segments to generate for each meridional wavenumber
    len_halo_spline      =   3       # halo round time period used to make spline calculation accurate near the start and end of the period
    k_spline             =   3       # use 3 or 5 

# 0.6 set whether to use the model state for the initial conditions or find a best fit initial state
    l_best_fit_initial_cond_v = True #  initial conditions for v solutions (soln of d^2 v / dt^2)
    l_best_fit_initial_cond_u = True #  initial conditions for u solutions (soln of du/dt)

# 0.7 control of calculation of power spectra and time series plots ; the code checks that inc_segments <= len_time_segments for this option 
    l_compute_spectra    = True 
    l_plot_slow          = True 
    l_plot_forces        = True 
    l_calc_du_dt         = False
    
    if l_calc_du_dt or l_plot_forces :
       if len_time_segments != inc_segments : 
          print (' if l_calc_du_dt = True choose len_time_segments and inc_segments to be the same' )
          sys.exit()
      
    # n_times_soln is the number of points in the spline solutions 
    n_times_soln = len_time_segments * max_no_segments

# 0.8 control of movie loops of reconstituted velocities ; the code checks that inc_segments <= len_time_segments for this option
    l_reconstruct     = True
    ch_recon          = 'orig'   #   'orig'  or  'sim' or 'fit' 

    # set up numpy array of logicals to control which vertical and meridional modes are used in the reconstruction
    #   np.zeros sets all to False; np.ones set all to True  
    l_recon_modes     = np.ones( ( n_vertical_modes, n_meridional_modes ), dtype = bool )   
#   l_recon_modes     = np.array( [ [1,1,0], [1,0,0], [1,0,0] ], dtype = bool )   # first list is the meridional modes for "first" vertical mode  

#----------------------------------------------------------------------------------------------------------
# 1. Print out user choices
#----------------------------------------------------------------------------------------------------------

    filestem = './std_out/'+ date + '/' + expt_name
    file_out = open(filestem+'.txt','w')

    print(' region, region_lc = ', region, region_lc, file=file_out)
    
    print(' time_series_input_file = ', time_series_input_file, file=file_out)

    print(' k_bottom_w  = ', k_bottom_w, file=file_out)
    if ll_vertical_density_profile_analytical :
       print(' dens_type, dens_diff, exp_scale_depth = ', dens_type, dens_diff, exp_scale_depth, file=file_out)
    else :
       print(' vertical_density_profile_input_file = ', vertical_density_profile_input_file, file=file_out)
       print(' j_lat_row_min, j_lat_row_max = ', j_lat_row_min, j_lat_row_max, file=file_out)

    print(' n_vertical_modes = ', n_vertical_modes, file=file_out)
    print(' MLD (depth_wind_penetration) = ', MLD, file=file_out)
    print(' l_surface_mode = ', l_surface_mode, file=file_out)
    print(' n_meridional_modes = ', n_meridional_modes, file=file_out)
    print(' Rayleigh_damp_fac = ', Rayleigh_damp_fac,  file=file_out)

    print(' len_time_segments = ', len_time_segments, file=file_out)
    print(' first_segment     = ', first_segment, file=file_out)
    print(' inc_segments     = ', inc_segments, file=file_out)
    print(' max_no_segments   = ', max_no_segments, file=file_out)
    print(' len_halo_spline = ', len_halo_spline, file=file_out)

    if len_halo_spline > first_segment :
       print('  len_halo_spline > first_segment; please avoid this choice ', len_halo_spline, first_segment, file=file_out)
       sys.exit()

    print(' l_Rayleigh, l_wind_x, l_wind_y, l_const = ', l_Rayleigh, l_wind_x, l_wind_y, l_const, file=file_out)
    print(' l_press, l_press_orig, l_press_full = ', l_press, l_press_orig, l_press_full, file=file_out)
    print(' l_force_fix = ', l_force_fix, file=file_out )

    print(' l_best_fit_initial_cond_u, l_best_fit_initial_cond_v = ', l_best_fit_initial_cond_u, l_best_fit_initial_cond_v, file=file_out ) 

    print(' l_compute_spectra, l_plot_slow, l_plot_forces, l_calc_du_dt = ', l_compute_spectra, l_plot_slow, l_plot_forces, l_calc_du_dt, file=file_out) 

    if l_reconstruct or l_compute_spectra : 
       if inc_segments > len_time_segments :
          print ( ' spectra and movie loops have gaps for this case. So setting these logicals to false ', file = file_out ) 
          l_reconstruct = False 
          l_compute_spectra = False

    if l_reconstruct : 
       print(' ch_recon = ', ch_recon, '; reconstruction includes the following modes', file=file_out)

       for m_vert in range ( n_vertical_modes ) :     
          for j_merid in range ( n_meridional_modes ) :     
              if l_recon_modes[ m_vert, j_merid ] : 
                 print( ' meridional mode', j_merid, '; vertical mode ', m_vert, file=file_out)

    print('-----------------------------', file=file_out)


#----------------------------------------------------------------------------------------------------------
# 2. Read in the time series input files
#----------------------------------------------------------------------------------------------------------

    e3t_0                    = read_one_field("/noc/msm/scratch/climate/atb299/NEMO_domains/75_lev_gridinfo/mesh_zgr.nc", 'e3t_0', file_out )
    deps                     = read_one_field(time_series_input_file, 'gdept_0', file_out )
    times                    = read_one_field(time_series_input_file, 'time_counter', file_out )
    v                        = read_one_field(time_series_input_file, 'vvel_int', file_out )    
    u_own_grid               = read_one_field(time_series_input_file, 'uvel_int', file_out )

    if l_press_full: 
       press_new = 'press_tot_w_e_dz'
    else : 
       press_new = 'press_w_e_dz'


    press_w_dz, press_e_dz, press_w_e_dz  = read_three_fields(time_series_input_file, 'press_w_dz', 'press_e_dz', press_new, file_out )
    lats_u, taux, tauy                    = read_three_fields(time_series_input_file, 'latitude', 'taux', 'tauy', file_out )
   
    n_times, n_lats = taux.shape
    n_times_y, n_lats_y = tauy.shape
    n_times_v, n_deps, n_lats_v = v.shape

    print(' n_times, n_times_y, n_times_v = ', n_times, n_times_y, n_times_v, file=file_out)
    print(' n_lats, n_lats_y, n_lats_v = ', n_lats, n_lats_y, n_lats_v, file=file_out)
    print(' n_deps = ', n_deps, file=file_out)

    last_tstep_reqd =  first_segment + (max_no_segments-1)*inc_segments + len_time_segments + len_halo_spline
    if last_tstep_reqd > n_times :
       print('  last_tstep_reqd > n_times; please avoid this choice ', last_tstep_reqd, n_times, file=file_out)
       sys.exit()

    if n_times != n_times_y or n_times != n_times_v  :
       print(' input data set time dimensions do not match. Exiting ', file=file_out)
       sys.exit()

    if n_lats_y != n_lats_v or n_lats != n_lats_y  :
       print(' input data set meridional dimensions do not match. Exiting ', file=file_out)
       sys.exit()

# 2.1 Set press_w_e_dz depending on user input (l_press_orig) 

    if l_press_orig : 
       press_w_e_dz = press_w_dz - press_e_dz

# 2.2 Adjust the latitudes if necessary 
#   The tracer and u fields have grids that coincide with the equator. There can be problems with the 0.0 values for lats.u being masked.   
#   So we change lats_u to a simple np array
#   Commenting out as this seems to give problems ??  
#    lats_u = lats_u.data      
 
    print('len(lats_u), lats_u[0], lats_u[2] =',len(lats_u), lats_u[0], lats_u[2], file=file_out)

    if lev_diag_user > 5:
       print(' lats_u[0:3] = ', lats_u[0:3], file=file_out)
       print(' lats_u[35:45] = ', lats_u[35:45], file=file_out)
    
# the latitudes for the tauy points are missing so we construct them  
    lats_v = np.zeros( n_lats ) 
    for jy in range (0, n_lats-1) :                                  
       lats_v[jy] = 0.5 * ( lats_u[jy] + lats_u[jy+1] )    
    lats_v[n_lats-1] = lats_u[n_lats-1] + 0.5 * ( lats_u[n_lats-1] - lats_u[n_lats-2] )
        
    if lev_diag_user > 5:
       print(' lats_v[0:3] = ', lats_v[0:3], file=file_out)

# convert latitudes in degrees to radians
    lats_u = deg2rad * lats_u
    lats_v = deg2rad * lats_v
    
#----------------------------------------------------------------------------------------------------------
# 3. Calculate taux and press_w_dz averaged onto the grid on which v and tauy are held
#----------------------------------------------------------------------------------------------------------

    beta = 2.0 * OMEGA / earth_radius
    y_u = earth_radius * lats_u
    y_v = earth_radius * lats_v

    taux_y_v = np.zeros( (n_times, n_lats) )
    press_w_e_dz_y_v = np.zeros( (n_times, n_deps, n_lats) )

    for jy in range ( n_lats-1 ) :
       taux_y_v[:, jy] = 0.5 * ( taux[:, jy] + taux[:, jy+1] )
       press_w_e_dz_y_v[:,:,jy] = 0.5 * ( press_w_e_dz[:,:,jy] + press_w_e_dz[:,:,jy+1] ) 

# the next line uses taux at the wrong latitude but it should not make much difference to results because meridional modes should be small there
    taux_y_v[:, n_lats-1] = taux[:, n_lats-1]
    press_w_e_dz_y_v[:,:, n_lats-1] = press_w_e_dz[:,:, n_lats-1] 

    u_y_v = np.zeros( (n_times, n_deps, n_lats) )
    for jy in range ( n_lats-1 ) :
       u_y_v[:,:,jy] = 0.5 * ( u_own_grid[:,:, jy] + u_own_grid[:,:, jy+1] )
    u_y_v[:,:, n_lats-1] = u_own_grid[:,:, n_lats-1]

#----------------------------------------------------------------------------------------------------------
# 4. Calculate the normalised vertical eigenmodes and the equivalent phase speeds
#----------------------------------------------------------------------------------------------------------
    if ll_vertical_density_profile_analytical :

      pot_rho = set_pot_dens_analytical( k_bottom_w, deps, dens_type, dens_diff, exp_scale_depth, file_out )

    else:
      dep_rho_levels = read_one_field ( vertical_density_profile_input_file, 'deptht', file_out )
      pot_rho_in = read_one_field ( vertical_density_profile_input_file, 'vosigmai', file_out )
      print ( ' pot_rho_in.shape = ', pot_rho_in.shape, file = file_out ) 

# 2019/07/17 start: avoid using missing data in calculation of pot_rho 

      pot_rho_min = 0.1    # missing values are actually 0.0 

      pot_rho = np.zeros ( k_bottom_w ) 
      for k in range ( k_bottom_w ) :
         icount = 0  
         for j in range ( j_lat_row_min, j_lat_row_max ) : 
            if pot_rho_in[0,k,j] > pot_rho_min: 
               icount = icount + 1
               pot_rho[k] = pot_rho[k] + pot_rho_in[0,k,j]    # pot_rho_in is a zonally averaged density field 
         if icount > 0 : 
            pot_rho[k] = pot_rho[k] / float(icount) 
         else : 
            print( ' icount = 0 in calculation of pot_rho. You need to reduce k_bottom_w: j = ', j, file = file_out )

# 2019/07/17 end

      if len(dep_rho_levels) != n_deps :
        print(' input data set vertical dimensions do not match. Exiting ', file=file_out)
        sys.exit()

    if lev_diag_user > 10:
       print ( ' k   pot_rho[k] ' , file=file_out)
       for k in range( k_bottom_w ) : 
          print ( k, pot_rho[k] , file=file_out)

# note that dz_layers values are all negative 
    c_mode, p_mode, dz_layers, H_bottom = vertical_eigenmodes_shooting( k_bottom_w, deps, pot_rho, n_vertical_modes, l_surface_mode, file_out )

    vert_proj_fac_mode = vertical_projection( p_mode, deps, MLD, H_bottom, file_out )    # note that vert_proj_fac_mode values are all positive 

    if lev_diag_user > 5:
        plotfilename = filestem + '_vert_modes.pdf' 
        vmode_plots(plotfilename, np.cumsum(dz_layers),p_mode)  # ATB Added ability to plot vertical mode structure
        
#----------------------------------------------------------------------------------------------------------
# 5. Begin loop over vertical mode projections
#----------------------------------------------------------------------------------------------------------

    if l_reconstruct : 
       v_recon_orig = np.zeros_like (v)      #   v_recon will contain velocity time-series reconstructed from selected fields and modes 
       v_recon_sim  = np.zeros_like (v)      #   v_recon will contain velocity time-series reconstructed from selected fields and modes 
       moc   = np.zeros_like (v)             #   Compute the model MOC from the velocities read in
       mocR  = np.zeros_like (v)             #   Compute the MOC from the velocities in v_recon_orig (mode projection onto v)
       mocRt = np.zeros_like (v)             #   Compute the MOC from the velocities in v_recon_sim  (simulation from surface wind and pressure fields)

    if l_plot_slow: 
       mean_sq_v_slow = np.zeros( (n_vertical_modes, n_meridional_modes) )
       mean_sq_v_diff = np.zeros_like(mean_sq_v_slow) 

    for m_vert in range(n_vertical_modes) :

       mvp1 = m_vert + 1   # this is the baroclinic vertical mode number 

       c_m = c_mode[m_vert]
       vert_proj_fac = vert_proj_fac_mode[m_vert]
       vert_mode = p_mode[m_vert]

       v_mode = np.zeros ( (n_times, n_lats) )
       press_w_e_mode = np.zeros ( (n_times, n_lats) )

       u_mode = np.zeros ( (n_times, n_lats) )  
       
       u_ts_modes = np.zeros( (n_meridional_modes, n_times) )   
       du_dt_sol_modes = np.zeros( (n_meridional_modes, n_times_soln) )   

#  du_dt_cmpt stores du/dt for up to 4 components in each meridional mode. It is calculated by projecting - r u + X^p + f v onto the meridional modes 
#             using a ladder operator expression. Only the first n_meridional_modes-1 modes can be calculated correctly.  
#             For j_merid = 0 only two components are non-zero; see du_dt_calc
       if l_calc_du_dt: 
          du_dt_cmpt    = np.zeros( (n_meridional_modes-1, 4, n_times_soln) )     

       for k in range ( k_bottom_w ) :     # changed to n_deps by MJB on 2 Nov 2018 

#  24/07/2019 removed gdef * rho_0 from expressions for v_mode, u_mode and press_w_e_mode
          v_mode[:,:] = v_mode[:,:] - vert_mode[k] * v[:,k,:] * ( dz_layers[k] / H_bottom )   
          u_mode[:,:] = u_mode[:,:] - vert_mode[k] * u_y_v[:,k,:] * ( dz_layers[k] / H_bottom )  

# next expression is based on (23) from version 1.1 of documentation      
# correction 4 March 2019: f_press_w[k] is p_w[k] * dz_layers[k]. So it is incorrect to multiply p_w[k] by dz_layers[k]   
#          f_press_w_mode[:,:] = f_press_w_mode[:,:] - gdef * vert_mode[k] * f_press_w[:,k,:] * dz_layers[k]        

          press_w_e_mode[:,:] = press_w_e_mode[:,:] + vert_mode[k] * ( press_w_e_dz_y_v[:,k,:] / H_bottom ) 
       press_w_e_mode = press_w_e_mode / rho_0     
            
#----------------------------------------------------------------------------------------------------------
# 5B. calculate and plot out all the meridional modes for this vertical wavenumber  
#----------------------------------------------------------------------------------------------------------

       if lev_diag_user > 5 or l_calc_du_dt :               
          merid_mode_set     = np.zeros( (n_meridional_modes, n_lats) ) 
          ttl = r'$\mathrm{mode} \; m = $' + str(mvp1)
          for j_merid in range (n_meridional_modes) :
             merid_mode_set[j_merid], dy_tilde = meridional_eigenmode( j_merid, c_m, lats_v, file_out )
             omega_sq = beta * c_m * (1.0 + 2.0*j_merid)    
         
          plotfilename = filestem + '_merid_modes_mvp1_' + str(mvp1) + '.pdf' 
          mmode_plots(plotfilename, lats_v*rad2deg, merid_mode_set, title=ttl)  # ATB Added ability to plot meridional mode structure

#----------------------------------------------------------------------------------------------------------
# 6. Begin loop over the meridional modes
#----------------------------------------------------------------------------------------------------------

       for j_merid in range (n_meridional_modes) :
          print("")
          print("Vertical mode", m_vert, ", Meridional mode", j_merid, "== Farrar & Durland (2012):","bc",m_vert+1,"m",j_merid)

          merid_mode, dy_tilde = meridional_eigenmode( j_merid, c_m, lats_v, file_out )      
          omega_sq = beta * c_m * (1.0 + 2.0*j_merid)   # Corrected: 26/09/2018. square of natural frequency of this mode (Eq 32) 
          period_in_days =  2.0 * math.pi / ( 86400. * math.sqrt(omega_sq)  ) 

          print('m_vert+1, j_merid, c_m, omega_sq, period in days = ', m_vert+1, j_merid, c_m, omega_sq, period_in_days, file=file_out)
                
# calculate the projection of the winds onto this meridional mode
          f_taux_ts = np.zeros( n_times )
          tauy_ts   = np.zeros( n_times )

          for jy in range ( n_lats ) :
             f_taux_ts[:] = f_taux_ts[:] + beta * y_v[jy] * taux_y_v[:, jy] * merid_mode[jy] * dy_tilde
             tauy_ts[:]   = tauy_ts[:]   + tauy[:, jy] * merid_mode[jy] * dy_tilde

# calculate the projection on this vertical mode
          f_taux_ts = vert_proj_fac * f_taux_ts
          tauy_ts   = vert_proj_fac * tauy_ts
      
          if l_calc_du_dt: 
             taux_ts   = np.zeros( n_times )
             for jy in range ( n_lats ) :
                taux_ts[:] = taux_ts[:] + taux_y_v[:, jy] * merid_mode[jy] * dy_tilde
             taux_ts = vert_proj_fac * taux_ts

          if lev_diag_user > 10:        
                print('f_taux_ts,tauy_ts=', f_taux_ts[:2], tauy_ts[:2])#, file=file_out )
                print('f_taux_ts.shape,tauy_ts.shape=', f_taux_ts.shape, tauy_ts.shape)#, file=file_out )


#----------------------------------------------------------------------------------------------------------
# 7. Calculate the projections of the meridional velocities onto this meridional mode
#----------------------------------------------------------------------------------------------------------
          u_ts = np.zeros( n_times )
          v_ts = np.zeros( n_times )
          f_press_w_e_ts = np.zeros( n_times ) 

          for jy in range (n_lats ) :
              u_ts[:] = u_ts[:] + u_mode[:, jy] * merid_mode[jy] * dy_tilde
              v_ts[:] = v_ts[:] + v_mode[:, jy] * merid_mode[jy] * dy_tilde
              f_press_w_e_ts[:] = f_press_w_e_ts[:] + beta * y_v[jy] * press_w_e_mode[:, jy] * merid_mode[jy] * dy_tilde
          u_ts_modes[j_merid] = u_ts
  
          if lev_diag_user > 10:        
              print('v_ts=', v_ts[:2])#, file=file_out )
          
          if j_merid < n_meridional_modes -1 : 
             if l_calc_du_dt: 
                press_w_e_ts = np.zeros( n_times )
                for jy in range (n_lats ) :
                   press_w_e_ts[:] = press_w_e_ts[:] + press_w_e_mode[:, jy] * merid_mode[jy] * dy_tilde
                
#----------------------------------------------------------------------------------------------------------
# 8. Calculate the wind-driven simulation of this meridional mode; starting from values of v_ts for user specified segments
#----------------------------------------------------------------------------------------------------------
          v_sim_seq, v_fit_seq, v_seq, t_seq, du_dt_seq, dtauy_dt_ts = time_step_solver(u_ts, v_ts, f_taux_ts,    
                                                         tauy_ts, f_press_w_e_ts, input_file_time_step, omega_sq, Rayleigh_damp_fac,
                                                         len_time_segments, first_segment, inc_segments, max_no_segments, 
                                                         len_halo_spline, k_spline, l_Rayleigh, l_wind_x, l_wind_y, l_press, 
                                                         l_const, l_force_fix, l_best_fit_initial_cond_v, file_out)

          du_dt_sol_modes[j_merid] = np.reshape(du_dt_seq[:,:inc_segments], n_times_soln)
#----------------------------------------------------------------------------------------------------------
# 9. Plot out results for this mode: specific values of m_vert and j_merid  
#----------------------------------------------------------------------------------------------------------

          print(v_sim_seq.shape,v_seq.shape,t_seq.shape, file=file_out)

          fileend = 'm_vert_' + str(m_vert) + '_j_merid_'+ str(j_merid) + '.pdf' 

# Compute power spectra using norm = "ortho" to obtain a unitary transform 
          if l_compute_spectra:  
          #print(inc_segments, max_no_segments)
            tsfull=np.reshape(v_seq[:,:inc_segments],max_no_segments*inc_segments)      
            v_seq_ps=np.abs(np.fft.rfft(tsfull,axis=0,norm="ortho"))
            samples_per_day=4./24.                                      #    
            freqy = np.fft.rfftfreq(len(tsfull), d=samples_per_day)
            #print(freqy.shape, v_seq_ps.shape)
            ulim=find_nearest(freqy,0.5)
            llim=1

# If simulated splines can form a continuous time series, compute power spectra for these too...
            ttl = "for mode (m,n)= (" + str(m_vert+1) + ", " + str(j_merid) + ")"
            if inc_segments == len_time_segments:
             print("Computing PS for splines also...")
             tsspline=np.reshape(v_sim_seq[:,:inc_segments],max_no_segments*inc_segments)        
             v_sim_seq_ps=np.abs(np.fft.rfft(tsspline,axis=0,norm="ortho"))
#             m_vertfd = m_vert+1
#             print(type(m_vertfd),type(j_merid), type(str(j_merid)))

             days_per_step  = input_file_time_step / 86400.
             end_segments = first_segment + max_no_segments*inc_segments
             t_seg_starts = days_per_step * np.arange(first_segment, end_segments, inc_segments)

             Osc_plots2(filestem, fileend, t_seg_starts, t_seq.T/86400., v_seq.T, v_sim_seq.T, v_fit_seq.T, 
                    freqy[llim:ulim], v_seq_ps[llim:ulim], v_sim_seq_ps[llim:ulim], title=ttl) # time (s) converted to days
            else:
             Osc_plots(filestem, fileend, t_seg_starts, t_seq.T/86400.,v_seq.T,v_sim_seq.T, v_fit_seq.T,
                    freqy[llim:ulim], v_seq_ps[llim:ulim]) # time (s) converted to days
            

#----------------------------------------------------------------------------------------------------------
# 9B. Calculate "slow" running mean time-series and an estimate of them from the forcing   
#----------------------------------------------------------------------------------------------------------
          if l_plot_slow:
             plotfilename = filestem + '_v_slow_m_vert_' + str(m_vert) + '_j_merid_'+ str(j_merid) + '.pdf' 
#             forced_ts =  (dtauy_dt_ts - f_taux_ts - f_press_w_e_ts ) / omega_sq
             forced_ts =  - (f_taux_ts + f_press_w_e_ts ) / omega_sq
             y_label = r'$\tilde{v}_{m,n} \; (\mathrm{m}^2 s^{-1})$'
             mean_sq_v_slow[m_vert,j_merid], mean_sq_v_diff[m_vert,j_merid] = Slow_plots(plotfilename, file_out, input_file_time_step, 
                                                                                     period_in_days, v_ts, y_label, ttl, add_ts = forced_ts ) 
         
          if l_plot_forces:
             plotfilename = filestem + '_force_m_vert_' + str(m_vert) + '_j_merid_'+ str(j_merid) + '.pdf'         
             end_segments = first_segment + max_no_segments*inc_segments 
             Force_plots(plotfilename, input_file_time_step/86400., first_segment, end_segments, dtauy_dt_ts, f_taux_ts, f_press_w_e_ts, ttl, file_out)   

#----------------------------------------------------------------------------------------------------------
# 9C. Calculate source terms for du/dt; these are stored in du_dt_cmpt
#----------------------------------------------------------------------------------------------------------
# the u modes are calculated from fv and X^* at the timesteps of the spline solutions
          if l_calc_du_dt: 

             if ch_recon == 'orig' : 
                v_sol = v_seq
             elif ch_recon == 'sim' :
                v_sol = v_sim_seq
             elif ch_recon == 'fit' :
                v_sol = v_fit_seq

             v_sol     = np.reshape(v_sol[:,:inc_segments], n_times_soln)      #  this is the time-series of v for this mode 
           
             force_x_ts = taux_ts + press_w_e_ts  
             istart = first_segment
             istop  = istart + n_times_soln
#  beta y = sqrt(0.5 * beta * c_m) * tilde_y           
             beta_y_fac = math.sqrt( 0.5 * beta * c_m )
             plotfilestem = filestem + '_m_vert_' + str(m_vert) 
             du_dt_cmpt = du_dt_calc(du_dt_cmpt, j_merid, beta_y_fac, v_sol, u_ts[istart:istop], force_x_ts[istart:istop], 
                             plotfilestem, input_file_time_step, period_in_days, ttl, file_out )
        
#----------------------------------------------------------------------------------------------------------
# 10. Calculate the contribution to v_recon from this mode (specific values of m_vert and j_merid)  
#----------------------------------------------------------------------------------------------------------

          if l_reconstruct : 
             if l_recon_modes[ m_vert, j_merid ] :
                v_recon_orig,v_recon_sim = rmodes(v_recon_orig,v_recon_sim,max_no_segments,first_segment,inc_segments,len_time_segments,k_bottom_w,v_sim_seq,v_fit_seq,vert_mode,merid_mode)

#----------------------------------------------------------------------------------------------------------
# 10B. Compare du_dt solutions with NEMO u time-series for these modes (still within for m_vert "loop")   
#----------------------------------------------------------------------------------------------------------

       if l_calc_du_dt: 
       
          time_axis = np.reshape(t_seq[:,:inc_segments], n_times_soln)   
          time_axis = time_axis / 86400.                 

          for j_merid in range (n_meridional_modes-1) :

             if j_merid == 0 : 
                n_cmpts = 3    #  no contribution from mode j_merid - 1   
             else : 
                n_cmpts = 4    #   contributions from j_merid - 1, j_merid, j_merid + 1 
        
             du_dt_matrix = np.zeros ( (n_times_soln, n_cmpts) )    #  matrix input to the least squares solve
             for icmpt in range ( n_cmpts ) : 
                du_dt_matrix[:,icmpt] =  du_dt_cmpt[j_merid, icmpt,:]

             du_dt_sol_ref = du_dt_sol_modes[j_merid] 
          
             x_lsq, sum_res_du_dt, rank, sing_vals = linalg.lstsq( du_dt_matrix, du_dt_sol_ref, rcond = None )
             Ray = x_lsq[0]      
             sum_du_dt = np.sum( du_dt_sol_ref*du_dt_sol_ref ) 
             rat_res_du_dt = sum_res_du_dt / sum_du_dt    

             u_ts = u_ts_modes[j_merid]      #  this is the time-series of u for this mode 
             u_sol_ref = u_ts[istart:istop]

# find particular solutions for both sets  
             u_part_full = np.zeros_like(u_sol_ref) 
             u_part_lsq  = np.zeros_like(u_sol_ref)
             for icmpt in range (1, n_cmpts ) :    #    - r u is not included in this loop! 
                u_part_icmpt = du_dt_integrate(du_dt_cmpt[j_merid, icmpt], 0.0, Ray, input_file_time_step, file_out )             
                u_part_full = u_part_full + u_part_icmpt 
                u_part_lsq  = u_part_lsq  + x_lsq[icmpt] * u_part_icmpt

# find homogeneous solution (the same one for both sets)  
             u_initial = u_sol_ref[0]
             u_zero = np.zeros( n_times_soln ) 
             u_sol_homog = du_dt_integrate(u_zero, u_initial, Ray, input_file_time_step, file_out )       

             if ( l_best_fit_initial_cond_u == False ) :
                u_sol_full = u_part_full + u_sol_homog
                u_sol_lsq  = u_part_lsq  + u_sol_homog
# moved these two lines on 30 July 2019
#                sum_res_full = np.mean( (u_sol_full - u_sol_ref)**2 )              
#                sum_res_lsq  = np.mean( (u_sol_lsq  - u_sol_ref)**2 )              
             else : 
                error_part_full = u_sol_ref - u_part_full
                error_part_lsq  = u_sol_ref - u_part_lsq
                
                matrix = np.zeros ( (n_times_soln, 1) )
                matrix[:,0] = u_sol_homog
        
                x_full_ics, sum_res_full, rank, sing_vals = linalg.lstsq( matrix, error_part_full, rcond = None ) 
                x_lsq_ics,  sum_res_lsq,  rank, sing_vals = linalg.lstsq( matrix, error_part_lsq,  rcond = None ) 

                u_sol_full = x_full_ics * u_sol_homog + u_part_full        
                u_sol_lsq  = x_lsq_ics  * u_sol_homog + u_part_lsq
# new position of these two lines on 30 July 2019
             mean_res_full = np.mean( (u_sol_full - u_sol_ref)**2 )              
             mean_res_lsq  = np.mean( (u_sol_lsq  - u_sol_ref)**2 )              
             mean_u_sol_ref = np.mean ( u_sol_ref ) 
             mean_ref_sq   = np.mean( (u_sol_ref - mean_u_sol_ref)**2 )    
             rat_res_full_ref =  mean_res_full / mean_ref_sq
             rat_res_lsq_ref  =  mean_res_lsq / mean_ref_sq
         
             print ( ' u solutions: j_merid, m_vert, Ray =  ', j_merid, m_vert, Ray, file = file_out ) 
             print ( ' rat_res_du_dt, sum_res_du_dt, sum_du_dt = ', rat_res_du_dt, sum_res_du_dt, sum_du_dt, file = file_out )  
             print ( ' mean_res_full, mean_res_lsq, mean_ref_sq, x_lsq = ', mean_res_full, mean_res_lsq, mean_ref_sq, x_lsq, file = file_out ) 
             print ( ' rat_res_full_ref, rat_res_lsq_ref = ', rat_res_full_ref, rat_res_lsq_ref, file = file_out )
        
             plotfilename = filestem + '_u_m_vert_' + str(m_vert) + '_j_merid_'+ str(j_merid) + '.pdf'  
             u_plot( plotfilename, time_axis,  u_sol_full, u_sol_lsq, u_sol_ref, file_out )


    if l_reconstruct :
       moc,mocR,mocRt = calc_moc(moc,mocR,mocRt,v,v_recon_orig,v_recon_sim,e3t_0[0,:],k_bottom_w)

       k1 = find_nearest(deps[:],1000)
       k2 = find_nearest(deps[:],2000)
       
       print("Values at          ",'{0:.0f}'.format(deps[k1]),"m ",'{0:.0f}'.format(deps[k2]),"m")
       rms1=np.mean(np.abs(moc[first_segment:inc_segments*max_no_segments,k1,40]))
       rms2=np.mean(np.abs(moc[first_segment:inc_segments*max_no_segments,k2,40]))
       print("MOC          RMS = ",'{0:.2f}'.format(rms1)," ",'{0:.2f}'.format(rms2))

#      rms1=np.mean(np.abs(mocR[first_segment:inc_segments*max_no_segments,k1,40]))
#      rms2=np.mean(np.abs(mocR[first_segment:inc_segments*max_no_segments,k2,40]))
#      print("MOCR         RMS = ",'{0:.2f}'.format(rms1)," ",'{0:.2f}'.format(rms2))
#      cc1 = np.corrcoef(moc[first_segment:inc_segments*max_no_segments,k1,40],mocR[first_segment:inc_segments*max_no_segments,k1,40])
#      cc2 = np.corrcoef(moc[first_segment:inc_segments*max_no_segments,k2,40],mocR[first_segment:inc_segments*max_no_segments,k2,40])
#      print("Corr. MOCR v MOC = ",'{0:.2f}'.format(cc1[1,0]**2)," ",'{0:.2f}'.format(cc2[1,0]**2))

       rms1=np.mean(np.abs(mocRt[first_segment:inc_segments*max_no_segments,k1,40]))
       rms2=np.mean(np.abs(mocRt[first_segment:inc_segments*max_no_segments,k2,40]))
       print("MOCRt        RMS = ",'{0:.2f}'.format(rms1)," ",'{0:.2f}'.format(rms2))
       cc1 = np.corrcoef(moc[first_segment:inc_segments*max_no_segments,k1,40],mocRt[first_segment:inc_segments*max_no_segments,k1,40])
       cc2 = np.corrcoef(moc[first_segment:inc_segments*max_no_segments,k2,40],mocRt[first_segment:inc_segments*max_no_segments,k2,40])
       print("Corr. MOCR v MOC = ",'{0:.2f}'.format(cc1[1,0]**2)," ",'{0:.2f}'.format(cc2[1,0]**2))
 
       results = evaluate(moc[first_segment:inc_segments*max_no_segments,k1,40],mocRt[first_segment:inc_segments*max_no_segments,k1,40])
       print(results)

       #fname="/noc/users/atb299/Python/"+str(expt_name)+"_moc_recon_nvert_"+str(n_vertical_modes)+"_nmerid_"+str(n_meridional_modes)+".nc"
       fname = filestem + '_MOC_m_vert_' + str(n_vertical_modes) + '_n_merid_'+ str(n_meridional_modes) + '.nc'  
       print(fname)
       save_moc(fname,deps,lats_v*rad2deg,moc,mocR,mocRt)
       
       #fname="/noc/users/atb299/Python/"+str(expt_name)+"_modes_nvert_"+str(n_vertical_modes)+"_nmerid_"+str(n_meridional_modes)+".nc"
       fname = filestem + '_MODES_m_vert_' + str(n_vertical_modes) + '_n_merid_'+ str(n_meridional_modes) + '.nc'  
       print(fname)
       save_modes(fname,deps[:k_bottom_w],lats_v*rad2deg,p_mode.T,merid_mode_set.T)

#----------------------------------------------------------------------------------------------------------
# 11. Generate movie loops of the original time series of <v'> and reconstructed time-series   
#----------------------------------------------------------------------------------------------------------
    for itime in range ( 0, n_times, 500 ) : 
       for kz in range (0, k_bottom_w, 30 ) : 
          print ( 'v[',itime,kz,']      ',  [ v[itime, kz, jy] for jy in range (0, n_lats, 20) ], file = file_out )      
          print ( 'v_recon_sim[',itime,kz,']', [ v_recon_sim[itime, kz, jy] for jy in range (0, n_lats, 20) ], file = file_out ) 
          print ( 'v_recon_orig[',itime,kz,']', [ v_recon_orig[itime, kz, jy] for jy in range (0, n_lats, 20) ], file = file_out ) 

#----------------------------------------------------------------------------------------------------------
# 12. Output stats for slow modes    
#----------------------------------------------------------------------------------------------------------
    if l_plot_slow: 
       for m_vert in range (n_vertical_modes) : 
          print( ' m_vert, mean_sq_v_slow[m_vert] = ', m_vert, mean_sq_v_slow[m_vert], file = file_out )  
          print( ' m_vert, mean_sq_v_diff[m_vert] = ', m_vert, mean_sq_v_diff[m_vert], file = file_out )  
          print( ' m_vert, mean_sq_v_diff/mean_sq_v_slow = ', m_vert, mean_sq_v_diff[m_vert]/mean_sq_v_slow[m_vert], file = file_out )  
    

#----------------------------------------------------------------------------------------------------------
    return   #      v_sim_seq, v_seq, t_seq  # MJB these do not need to be returned
#----------------------------------------------------------------------------------------------------------
