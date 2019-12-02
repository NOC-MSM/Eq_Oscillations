# -*- coding: iso-8859-1 -*-
from __future__ import print_function
'''
Functions: read_one_field, read_three_fields
author:   Mike Bell   
Date:     May 2018
Description: reads input fields  
''' 
#          version 1:  15 May 2018 is first version
#
#-------------------------------------------------------------------------------------------------

def read_lsm(dir_in):

   import numpy as np 
   import netCDF4

   f = netCDF4.Dataset(dir_in+'/mask.nc','r')
   mask_orca = f.variables["tmask"][:]
   mask_config = mask_orca[0]
   f.close()

   print ( 'shape(mask_config) = ', np.shape(mask_config) )
   print ( 'type(mask_config)  = ', type(mask_config) )

   return mask_config 

def read_one_field(input_file, field_id, file_out):

    import numpy as np 
    import netCDF4

    g = netCDF4.Dataset(input_file,'r')

    field = g.variables[field_id][:]
    print(' input_file, field_id, field.shape = ', input_file, field_id, field.shape, file=file_out)
    print(' type(field) = ', type(field), file=file_out)

    if np.ndim(field) == 1: 
       print(' field[0:2] = ', field[0:2], file=file_out)

    g.close()

    return field 

def read_three_fields(input_file, id_0, id_1, id_2, file_out):

    import numpy as np 
    import netCDF4

    field_0 = read_one_field(input_file, id_0, file_out)
    field_1 = read_one_field(input_file, id_1, file_out)
    field_2 = read_one_field(input_file, id_2, file_out)

    return field_0, field_1, field_2 

#-------------------------------------------------------------------------------------------------