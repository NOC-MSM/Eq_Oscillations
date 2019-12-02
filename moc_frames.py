########################################################
# Produce images with multiple (3) simulations 
########################################################

import numpy as np
import multiprocessing

def make_plot(tt):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    fig = plt.figure(figsize=(24,13))

    ax1 = plt.subplot2grid((10, 8), (0, 0),colspan=2,rowspan=6)
    ax2 = plt.subplot2grid((10, 8), (0, 2),colspan=2,rowspan=6)
    ax3 = plt.subplot2grid((10, 8), (0, 4),colspan=2,rowspan=6)
    ax4 = plt.subplot2grid((10, 8), (0, 6),colspan=2,rowspan=6)
    ax5 = plt.subplot2grid((10, 8), (7, 0),colspan=8,rowspan=3)
    axC = fig.add_axes([0.91, 0.43, 0.02, 0.45])

    cmap = mpl.cm.seismic
    norm = mpl.colors.Normalize(vmin=-200, vmax=200)
    cb1  = mpl.colorbar.ColorbarBase(axC, cmap=cmap, norm=norm, orientation='vertical',extend='both')
    cb1.set_label('Sv',fontsize=20)
    cb1.ax.tick_params(labelsize=16) 

    ax1.set_title("Pacific MOC", fontsize=20)
    ax1.set_xlabel("Latitude", fontsize=20)
    ax1.set_ylabel("Depth (m)", fontsize=20)
    ax1.set_xlim([-10.5,10.5])
    ax1.xaxis.set_tick_params(labelsize=16)
    ax1.yaxis.set_tick_params(labelsize=16)
    ax1.invert_yaxis()

    ax2.set_title("Mode (1,0) only", fontsize=20)
    ax2.set_xlabel("Latitude", fontsize=20)
    ax2.xaxis.set_tick_params(labelsize=16)
    ax2.set_xlim([-10.5,10.5])
    ax2.axes.get_yaxis().set_visible(False)
    ax2.invert_yaxis()

    ax3.set_title("All modes to (3,2)", fontsize=20)
    ax3.set_xlabel("Latitude", fontsize=20)
    ax3.xaxis.set_tick_params(labelsize=16)
    ax3.set_xlim([-10.5,10.5])
    ax3.axes.get_yaxis().set_visible(False)
    ax3.invert_yaxis()

    ax4.set_title("All modes to (6,5)", fontsize=20)
    ax4.set_xlabel("Latitude", fontsize=20)
    ax4.xaxis.set_tick_params(labelsize=16)
    ax4.set_xlim([-10.5,10.5])
    ax4.axes.get_yaxis().set_visible(False)
    ax4.invert_yaxis()

    ax5.set_xlim(0,len(moc)/6)
    ax5.set_xlabel('Time (days)', fontsize=20)
    ax5.set_ylabel('MOC (Sv)', fontsize=20)
    ax5.xaxis.set_tick_params(labelsize=16)
    ax5.yaxis.set_tick_params(labelsize=16)
    ax5.plot(tx,moc[:,46,40],color='k',linewidth=3,label="Pacific MOC")
    ax5.plot(tx,mocRt11[:,46,40],color='r',linewidth=2,label="Mode (1,0)")
    ax5.plot(tx,mocRt33[:,46,40],color='b',linewidth=2,label="All modes to (3,2)")
    ax5.plot(tx,mocRt66[:,46,40],color='c',linewidth=2,label="All modes to (6,5)")
    ax5.legend(bbox_to_anchor=(0.1, 1.02, .8, .102), loc=3, ncol=4, mode="expand", borderaxespad=0., fontsize=16)

    # Ensure these values are consistent with those used when running Eq_Osc_Wind_1906.py
    first_segment     =   3      # start of first segment
    inc_segments      = 100      # increment between segments (hours)
    max_no_segments   =  22      # (maximum) number of segments to generate for each meridional wavenumber
    input_file_time_step =  4. * 3600.  # 4 hourly timesteps ; this information could be read from each of the files
    days_per_step     = input_file_time_step / 86400.
    end_segments      = first_segment + max_no_segments*inc_segments
    t_seg_starts      = days_per_step * np.arange(first_segment, end_segments, inc_segments)
    heights = np.zeros_like(t_seg_starts)
    ax5.plot(t_seg_starts, heights, 'ro') 

    #PMOC plot
    cs = ax1.pcolormesh(lats,deps[:72],moc[tt,:72,:], cmap="seismic",vmin=-200,vmax=200)
    ax1.plot([0],[1000], marker="D",markersize=10,color='k')
 
    #PMOC mode (1,0) plot
    cs = ax2.pcolormesh(lats,deps[:72],mocRt11[tt,:72,:], cmap="seismic",vmin=-200,vmax=200)
    ax2.plot([0],[1000], marker="D",markersize=10,color='k')

    #PMOC mode (3,3) plot
    cs = ax3.pcolormesh(lats,deps[:72],mocRt33[tt,:72,:], cmap="seismic",vmin=-200,vmax=200)
    ax3.plot([0],[1000], marker="D",markersize=10,color='k')

    #PMOC mode (6,6) plot
    cs = ax4.pcolormesh(lats,deps[:72],mocRt66[tt,:72,:], cmap="seismic",vmin=-200,vmax=200)
    ax4.plot([0],[1000], marker="D",markersize=10,color='k')
    
    #Time series plot
    tmark = ax5.axvline(x=tx[tt],linewidth=3)

    # Draw and save >> Edit path as required
    plt.savefig(odir+file+"_{:04d}.png".format(tt), dpi=64, transparent=False, bbox_inches='tight', pad_inches=0.1)
    plt.close()
    return None


def moc_frames(moc,mocR,mocRt11,mocRt33,mocRt66,odir):
    tt = np.arange( 402 )#len(moc) )
    p = multiprocessing.Pool()
    p.map(make_plot, tt)
    return None


########################################################
########################################################
# Read MOC data
import netCDF4
import os

from numpy import pi
deg2rad, rad2deg = pi / 180.0,  180.0 / pi

# Choose number of days
ndays = 67
nt = ndays * 6

fpath = "/noc/msm/working/valor/MB_NIWs/"

#file = "pac_moc_recon_nvert_1_nmerid_1_nat100"
#file = "pac_moc_recon_nvert_1_nmerid_1_lsf100"
file = "pac_moc_recon_nvert_1_nmerid_1_lsf100_73L"
fs = netCDF4.Dataset(fpath+file+'.nc')
moc    = fs.variables['moc'][:nt,:,:]
mocR   = fs.variables['mocR'][:nt,:,:]
mocRt11  = fs.variables['mocRt'][:nt,:,:]
lats = fs.variables['latitude'][:]*rad2deg
deps = fs.variables['depth'][:]
tx   = fs.variables['time_counter'][:nt]
fs.close()

#file = "pac_moc_recon_nvert_3_nmerid_3_nat100"
#file = "pac_moc_recon_nvert_3_nmerid_3_lsf100"
file = "pac_moc_recon_nvert_3_nmerid_3_lsf100_73L"
fs = netCDF4.Dataset(fpath+file+'.nc')
mocRt33  = fs.variables['mocRt'][:nt,:,:]
fs.close()

#file = "pac_moc_recon_nvert_6_nmerid_6_nat100"
#file = "pac_moc_recon_nvert_6_nmerid_6_lsf100"
file = "pac_moc_recon_nvert_6_nmerid_6_lsf100_73L"
fs = netCDF4.Dataset(fpath+file+'.nc')
mocRt66  = fs.variables['mocRt'][:nt,:,:]
fs.close()

# Ensure that the output directory exists
odir = '/noc/users/atb299/Python/animation_pngs/'
os.system('mkdir '+odir)

moc_frames(moc,mocR,mocRt11,mocRt33,mocRt66,odir)

# Command line animation via ffmpeg
# On command line, something like: ffmpeg -framerate 24 -i pac_moc_recon_nvert_6_nmerid_6_lsf100_73L_%04d.png -vcodec mpeg4 -vb 50M ../pac_moc_recon.mp4

gen_movie = 'ffmpeg -framerate 24 -i /noc/users/atb299/Python/animation_pngs/'+file+'_%04d.png -vcodec mpeg4 -vb 50M '+file+'.mp4'
os.system(gen_movie)

