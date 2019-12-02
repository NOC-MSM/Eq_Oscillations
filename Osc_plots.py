# -*- coding: iso-8859-1 -*-
'''
Functions: Osc_plots, Osc_plots2
author:   Adam Blaker
Date:     September 2018
Description: Produce line plots of the output from the function "Eq_Osc_Wind"
'''
#          version 1:  17 September 2018
#                      24 April 2019   removed factor 1.E-4 mulitiplying time-series before output
#
#-------------------------------------------------------------------------------------------------

def Osc_plots(filestem, fileend, t_seg_starts, time_axis, v_seq, v_sim_seq, v_sim_fit_seq, freqy, v_seq_ps):

    import matplotlib.pyplot as plt
    import numpy as np


    # Plot time series of v_seq and v_sim_seq 
    # Removed multipliciation by 1.0E-4 in following 3 plt.plot statements: 10/06/2019 
    # (for consistency with Osc_plots2)
      
    fig = plt.figure(figsize=(15,8))
    plotfilename = filestem + '_timeseries_' + fileend

    ax1 = plt.subplot2grid((5, 4), (0, 0),colspan=4,rowspan=2)
    plt.plot(time_axis,v_seq,color='k',linewidth=2, linestyle='-', label = 'NEMO')
    plt.plot(time_axis,v_sim_seq,color='r',linewidth=2, linestyle=':', label = 'sim')
#    plt.plot(time_axis,v_sim_fit_seq,color='b',linewidth=2, linestyle='--', label = 'lsf')

    heights = np.zeros_like(t_seg_starts)
    plt.plot(t_seg_starts, heights, 'ro') 

    plt.title("Time series of mode projection and spline fits", fontsize=16)
    plt.xlabel("Time (days)", fontsize=12)
    plt.ylabel(r'$ \tilde{v}_{m,n} \; (\mathrm{m}^2 \mathrm{s}^{-1})$', fontsize=12)
    plt.legend()

    plt.savefig(plotfilename, transparent=True, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')
        
    # Plot power spectrum of data
    fig = plt.figure(figsize=(15,8))
    plotfilename = filestem + '_power_spectrum_' + fileend

    ax2 = plt.subplot2grid((5, 4), (3, 0),colspan=4,rowspan=2)
    plt.plot(freqy[1:],v_seq_ps[1:])#,freqy[llim:ulim],v_seq_ps[llim:ulim])
    plt.title("Power spectrum of v_seq (solid black above)", fontsize=16)
    plt.xlabel("Frequency (cpd)", fontsize=12)
    plt.ylabel("Power", fontsize=12)
    plt.legend()

    plt.show()

    plt.savefig(plotfilename, transparent=True, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')
    return

#-------------------------------------------------------------------------------------------------

def Osc_plots2(filestem, fileend, t_seg_starts, time_axis, v_seq, v_sim_seq, v_sim_fit_seq, freqy, v_seq_ps, v_sim_seq_ps, title=None):

    import matplotlib.pyplot as plt
    import numpy as np


    # Plot time series of v_seq and v_sim_seq
    fig = plt.figure(figsize=(15,8))
    plotfilename = filestem + '_timeseries_' + fileend

    ax1 = plt.subplot2grid((5, 4), (0, 0),colspan=4,rowspan=2)
#    plt.plot(time_axis,v_seq*1.E-4,color='k',linewidth=2, linestyle='-')   factor 1.E-4 removed from following 3 lines (24/04/2019)
    plt.plot(time_axis,v_seq,color='k',linewidth=2, linestyle='-') #  label = 'NEMO')
    plt.plot(time_axis,v_sim_seq,color='r',linewidth=2, linestyle=':') # label = 'sim')
#    plt.plot(time_axis,v_sim_fit_seq,color='b',linewidth=2, linestyle='--') #  label = 'lsf')

    heights = np.zeros_like(t_seg_starts)
    plt.plot(t_seg_starts, heights, 'ro') 

    plt.title("Time series: "+title, fontsize=16)
    plt.xlabel("Time (days)", fontsize=12)
    plt.ylabel(r'$ \tilde{v}_{m,n} \; (\mathrm{m}^2 \mathrm{s}^{-1})$', fontsize=12)
#     plt.legend()     # there are multiple legends for each line (because it is plotted for each segment so labels & legends are commented out 

    fig.savefig(plotfilename, transparent=True, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')
    
    # Plot power spectrum of data
    fig = plt.figure(figsize=(15,8))
    plotfilename = filestem + '_power_spectrum_' + fileend

    ax2 = plt.subplot2grid((5, 4), (3, 0),colspan=4,rowspan=2)
    plt.plot(freqy[1:],v_seq_ps[1:],color='k' )     #  label = 'NEMO')  #,freqy[llim:ulim],v_seq_ps[llim:ulim])
    plt.plot(freqy[1:],v_sim_seq_ps[1:],color='r', linestyle=':' ) # label = 'sim')    #,freqy[llim:ulim],v_seq_ps[llim:ulim])
    plt.title("Power spectra: "+title, fontsize=16)
    plt.xlabel("Frequency (cpd)", fontsize=12)
    plt.ylabel(r'$\mathrm{Power} \; (\mathrm{m}^2 \mathrm{s}^{-1})$', fontsize=12)
#     plt.legend()

    fig.savefig(plotfilename, transparent=True, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')

    return

#-------------------------------------------------------------------------------------------------
