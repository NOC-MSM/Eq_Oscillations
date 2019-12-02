# -*- coding: iso-8859-1 -*-
'''
Functions: vmode_plots, mmode_plots
author:   Adam Blaker
Date:     September 2018
Description: 
   vmode_plots: Plot vertical mode structure as a function of depth for all modes
   mmode_plots: plot meridional mode strucutre as function of latitude
'''
#          version 1:  17 September 2018
#                      22 July 2019 - corrected "Phase" to "Amplitude" in both vmode and mmode
#
#-------------------------------------------------------------------------------------------------

def vmode_plots(plotfilename, depths, p_mode):

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(6,12))

    ax1 = plt.subplot2grid((2, 2), (0, 0),colspan=2,rowspan=2)
    for m in range (p_mode.shape[0]) :
       if m == 0 : 
          plt.plot(p_mode[m,:],depths,linewidth=4, color='r', label=r'$m=$'+str(m+1))
       else :
          plt.plot(p_mode[m,:],depths,linewidth=2, linestyle='--', label=r'$m=$'+str(m+1))
       
    plt.xlim([-4.6,4.6])
    plt.title("Vertical mode structure", fontsize=16)
    plt.xlabel("Amplitude", fontsize=12)
    plt.ylabel("Depth (m)", fontsize=12)
    plt.legend()
        
    fig.savefig(plotfilename, transparent=True, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')
    return

def mmode_plots(plotfilename, lats_v, merid_mode_set, title=None):

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(15,4))

    ax1 = plt.subplot2grid((2, 2), (0, 0),colspan=2,rowspan=2)
    for n in range (merid_mode_set.shape[0]) :
       if n == 0 : 
          plt.plot(lats_v,merid_mode_set[n,:], linewidth=4, color='b', label=r'$n=$'+str(n))
       else : 
          plt.plot(lats_v,merid_mode_set[n,:], linewidth=2, linestyle='--', label=r'$n=$'+str(n))

    plt.title("Meridional mode structures for "+title, fontsize=16)
    plt.xlabel("Latitude", fontsize=12)
    plt.ylabel("Amplitude", fontsize=12)
    plt.legend()
        
    fig.savefig(plotfilename, transparent=True, bbox_inches='tight', pad_inches=0.1)
    plt.close('all')
    return
