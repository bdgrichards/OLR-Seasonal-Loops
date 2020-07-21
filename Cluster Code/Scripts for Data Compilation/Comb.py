import xarray as xr
import numpy as np

outp = xr.open_mfdataset('./Output_Data_temp_case/*.nc', combine='by_coords')

outp.to_netcdf('Combined_data_temp_case.nc')

