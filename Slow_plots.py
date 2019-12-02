# -*- coding: iso-8859-1 -*-
'''
Functions: Slow_plots 
author:   Mike Bell
Date:     September 2018
Description: Produce line plots of slow timeseries 
'''
#          version 1:  17 September 2018
#
#-------------------------------------------------------------------------------------------------

def Slow_plots(plotfilename, file_out, input_file_time_step, period_in_days, v_ts, y_label, ttl, add_ts=None):

# return values depend on whether add_ts=None 


    import numpy as np 
    import matplotlib.pyplot as plt

    lev_diag = 20

    len_ts = len(v_ts) 
    
    half_window = 0.5 * 86400. * period_in_days / input_file_time_step   # number of time-steps in average to be used
    
    it_off_start = int( round( -half_window ) ) # negative integer (offset to start of window)
    it_off_end   = int( round(  half_window ) ) # positive integer (offset to end of window) 
 
    len_window = it_off_end - it_off_start    
    len_slow_ts = len_ts + it_off_start - it_off_end
    
    if lev_diag > 10: 
       print( ' half_window, len_window, len_slow_ts = ', half_window, len_window, len_slow_ts, file=file_out) 
    
    time_slow   = np.zeros(len_slow_ts) 
    v_slow      = np.zeros(len_slow_ts)
    forced_slow = np.zeros(len_slow_ts)
        
    for islow in range (len_slow_ts) : 
       it = islow - it_off_start
       time_slow[islow] = ( input_file_time_step * it ) / 86400.   # time in days
       v_slow[islow] = np.mean( v_ts[islow:islow+len_window]) 
       if not ( add_ts is None ) : 
          forced_slow[islow] = np.mean( add_ts[islow:islow+len_window]) 
    
# Plot time series of v_slow and forced_slow
    fig = plt.figure(figsize=(15,8))
    
    ax1 = plt.subplot2grid((5, 4), (0, 0),colspan=4,rowspan=2)
    plt.plot(time_slow,v_slow,color='k',linewidth=2, linestyle='-', label = 'NEMO')
    if not ( add_ts is None ) : 
       plt.plot(time_slow,forced_slow,color='r',linewidth=2, linestyle=':',label = 'forced')

    plt.title("Time series of slow variations " + ttl, fontsize=16)
    plt.xlabel("Time (days)", fontsize=12)
#    plt.ylabel(r'$\tilde{v}_{m,n} \; (\mathrm{m}^2 s^{-1})$', fontsize=12)
    plt.ylabel(y_label, fontsize=12)
    plt.legend()
    
    plt.savefig(plotfilename, transparent=True, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')
    
    mean_sq_v_slow = np.mean(v_slow*v_slow)
    
    if add_ts is None :
    
       return mean_sq_v_slow
       
    else : 
       v_diff = v_slow - forced_slow
       mean_sq_diff = np.mean(v_diff*v_diff)  
       return mean_sq_v_slow, mean_sq_diff

#-------------------------------------------------------------------------------------------------

