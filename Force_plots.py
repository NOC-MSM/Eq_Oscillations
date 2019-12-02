# -*- coding: iso-8859-1 -*-
'''
Functions: Force_plots 
author:   Mike Bell
Date:     Apr 2019
Description: Produce line plots of forces as timeseries: note that tauy is not differentiated ?  
'''
#          version 1:  17 September 2018
#
#-------------------------------------------------------------------------------------------------

def Force_plots(plotfilename, input_file_time_step, start, end,
                dtauy_dt_ts, f_taux_ts, f_press_w_e_ts, ttl, file_out):

    import numpy as np 
    import matplotlib.pyplot as plt

    
    len_ts = len(dtauy_dt_ts) 
    time_ts = input_file_time_step * np.arange(len_ts) 

# Plot time series of v_slow and forced_slow
    fig = plt.figure(figsize=(15,8))
    
    ax1 = plt.subplot2grid((5, 4), (0, 0),colspan=4,rowspan=2)
    plt.plot(time_ts[start:end],dtauy_dt_ts[start:end],color='0.7',linewidth=2, linestyle='-', label=r'$ \mathrm{d} \tilde{Y} /\mathrm{d}t $')
    plt.plot(time_ts[start:end],f_taux_ts[start:end],color='r',linewidth=2, linestyle=':', label = r'$f \tilde{X}$')
    plt.plot(time_ts[start:end],f_press_w_e_ts[start:end],color='b',linewidth=2, linestyle='--', label = r'$f( \tilde{p}_w - \tilde{p}_E )$')

    plt.title("Time series of forcings " + ttl, fontsize=16)
    plt.xlabel("Time (days)", fontsize=12)
    plt.ylabel("Forcing $(\mathrm{m}^2 s^{-2})$", fontsize=12)
    plt.legend()
    
    plt.savefig(plotfilename, transparent=True, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')

    return

#-------------------------------------------------------------------------------------------------

