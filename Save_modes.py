# Create and save netcdf output
def save_modes(fname,deps,lats,vmodes,mmodes):
    import netCDF4
    import numpy as np
    
    ncid_out=netCDF4.Dataset(fname,'w',format='NETCDF4',clobber=True)
    #Make the dimensions
    shape1=vmodes.shape
    shape2=mmodes.shape
    z=ncid_out.createDimension('depth',shape1[0])
    y=ncid_out.createDimension('latitude',shape2[0])
    m=ncid_out.createDimension('vmode',shape1[1])
    n=ncid_out.createDimension('mmode',shape2[1])
    
    #Define the variables
    #========
    vmvar=ncid_out.createVariable('vmodes','f4',('depth','vmode'), fill_value=0., zlib=True)
    vmvar.coordinates="depth vmode"
    
    mmvar=ncid_out.createVariable('mmodes','f4',('latitude','mmode'), fill_value=0., zlib=True)
    mmvar.coordinates="latitude mmode"
    
    #========
    gphitvar=ncid_out.createVariable('latitude','f4',('latitude'), zlib=True)
    gphitvar.units="degrees_north"
    gphitvar.standard_name = "latitude"
    gphitvar.long_name = "Latitude"
    #========
    gdepvar=ncid_out.createVariable('depth','f4',('depth'), zlib=True)
    gdepvar.units="m"
    gdepvar.standard_name = "depth"
    gdepvar.long_name = "depth"
    
    #Write to the variables
    vmvar[:]=vmodes
    mmvar[:]=mmodes
    gphitvar[:]=lats
    gdepvar[:]=deps
    
    ncid_out.close()
    return None
