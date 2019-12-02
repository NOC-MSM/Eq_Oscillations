# -*- coding: iso-8859-1 -*-
'''
Function: Run_Eq_Osc_Wind
author:   Mike Bell   
Date:     May 2018
Description: calculates time-series of zonal mean equatorial waves ( inertia-gravity, Yanai) in the low-order baroclinic 
             vertical modes (modes 1, 2 and 3) driven by the local wind forcing  
''' 
#          version 1:  15 May 2018 is first version
#
#-------------------------------------------------------------------------------------------------

def Run_Eq_Osc_Wind():

# Inputs     
#    region     string: e.g. 'PAC'  options are: 'PAC', 'PAC_E', 'PAC_W', 'ATL', 'IND', 'IND_E', 'IND_W' 
#
# Returns: null 

    from Eq_Osc_Wind import Eq_Osc_Wind 

#----------------------------------------------------------------------------------------------------------
# 0. Set list of regions to loop over. Other user options are set in section 0 of Eq_Osc_Wind   
#----------------------------------------------------------------------------------------------------------

    list_regions_full = ['PAC', 'PAC_E', 'PAC_W', 'ATL', 'IND', 'IND_E', 'IND_W']

    list_regions = ['PAC', 'ATL', 'IND']
    list_regions = ['PAC']

    for region in list_regions: 

       Eq_Osc_Wind(region)

#----------------------------------------------------------------------------------------------------------
Run_Eq_Osc_Wind()
#----------------------------------------------------------------------------------------------------------
 
