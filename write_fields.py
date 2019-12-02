# -*- coding: iso-8859-1 -*-
# Name:    write_fields
# Functions:  write_press_diff_field
# Author:  Adam Blaker
# Date:    Oct 2018
# Purpose: writes fields computed during calculation of forces

import netCDF4
import numpy as np


#  j_out = j_max + 1 - j_min
#  press_w_out  = np.zeros( (kmax, j_out) )
#  press_e_out  = np.zeros( (kmax, j_out) )
#  lat_out      = np.zeros( j_out )
#  taux_out     = np.zeros( j_out )
#  tauy_out     = np.zeros( j_out )

def write_press_diff_field(lats, deps, times, press_w, press_e, ofile ):
    ncid_out=netCDF4.Dataset(ofile,'w',format='NETCDF3_CLASSIC')
    print ("output_file = ", ofile)

    #Make the dimensions
    jmax = len(lats)
    kmax = len(deps) 
    nmax = len(times)
    y=ncid_out.createDimension('y',jmax)
    z=ncid_out.createDimension('z',kmax)
    t=ncid_out.createDimension('t',nmax)
    press_diff = press_w - press_e

    #Define the variables
    time   =ncid_out.createVariable('time_counter','d',('t'))
    gdept_0=ncid_out.createVariable('gdept_0','d',('z'))
    lat    =ncid_out.createVariable('latitude','d',('y'))
    press_w_e=ncid_out.createVariable('press_w_e','d',('t','z','y'))

    #Write the variables
    gdept_0[:]=deps
    lat[:]=lats
    press_w_e[:]=press_diff

    ncid_out.close()

    return
