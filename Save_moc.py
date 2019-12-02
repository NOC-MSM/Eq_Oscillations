# Create and save netcdf output
def save_moc(fname,deps,lats,moc,mocR,mocRt):
    import netCDF4
    import numpy as np
    
    ncid_out=netCDF4.Dataset(fname,'w',format='NETCDF4',clobber=True)
    #Make the dimensions
    shape=moc.shape
    y=ncid_out.createDimension('latitude',shape[2])
    z=ncid_out.createDimension('depth',shape[1])
    t=ncid_out.createDimension('time_counter',None)
    
    tx = np.arange(0,(2220-1/6)/6,1/6)  # Timeseries
    
    #Define the variables
    #========
    MOCvar=ncid_out.createVariable('moc','f4',('time_counter','depth','latitude'), fill_value=0., zlib=True)
    MOCvar.coordinates="time_counter depth latitude"
    
    MOCRvar=ncid_out.createVariable('mocR','f4',('time_counter','depth','latitude'), fill_value=0., zlib=True)
    MOCRvar.coordinates="time_counter depth latitude"
    
    MOCRTvar=ncid_out.createVariable('mocRt','f4',('time_counter','depth','latitude'), fill_value=0., zlib=True)
    MOCRTvar.coordinates="time_counter depth latitude"
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
    #========
    tcntvar=ncid_out.createVariable('time_counter','f4',('time_counter'),zlib=True)
    tcntvar.units="days since 2006-04-01 00:00:00"
    tcntvar.axis="T"
    tcntvar.standard_name="time"
    tcntvar.long_name="Time axis"
    tcntvar.calendar="noleap"
    tcntvar.time_origin="2006-04-01 00:00:00"
    
    #Write to the variables
    MOCvar[:]=moc
    MOCRvar[:]=mocR
    MOCRTvar[:]=mocRt
    
    gphitvar[:]=lats
    gdepvar[:]=deps
    tcntvar[:]=tx
    
    ncid_out.close()
    return None
