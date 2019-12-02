# -*- coding: iso-8859-1 -*-
from __future__ import print_function
'''
Created on 8 Mar 2019

@author: Mike Bell 

check consistency of T and S fields with tracer mask for Adam's ORCA025
 
'''
#--------------------------------------------------------------------------------------------------------------

def check_tmask() : 


#----------------------------------------
# 0. Preliminary declarations
# ---------------------------------------

    import sys
    import numpy as np 
    import math 
    from read_fields import read_lsm, read_one_field

#----------------------------------------
# choose files to use 
# --------------------------------------- 

    dir_in = '/data/local/frmk/AdamBlaker201810'
    file_in = dir_in + '/VN206HF_4h_20060403_1_T.nc'

    file_std = open('std_out/check_tmask.txt','w') 
    print ( 'Starting ; check_tmask, dir_in, file_in = ', dir_in, file_in, file=file_std )  


#----------------------------------------
# 1. Read in the model fields that will be used 
# ---------------------------------------

    mask_T =  read_lsm (dir_in) 

    t_all = read_one_field(file_in, 'votemper', file_std )
    tn = t_all[0].data
    s_all = read_one_field(file_in, 'vosaline', file_std )
    sn = s_all[0].data
    kmax, jmax, imax = tn.shape
    print ( ' kmax, jmax, imax = ', kmax, jmax, imax, file = file_std ) 

#----------------------------------------
# 1. Compare the fields
# ---------------------------------------

    tn_mask = np.where( tn == 0.0, False, True)
    sn_mask = np.where( sn == 0.0, False, True) 
    ts_mask = np.logical_and( tn_mask, sn_mask)      # ts_mask is land where BOTH tn and sn suggest it is. 

    miss_match = np.where( ts_mask != mask_T, 1, 0) 

    icount = np.sum(miss_match)

    print ( ' icount = ', icount , file = file_std ) 

check_tmask()
#----------------------------------------------------------------------------------------------------------
