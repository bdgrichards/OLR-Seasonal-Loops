import xarray as xr
import numpy as np

outp = xr.open_mfdataset('./Output_Data_ceres_const_rh/*.nc', combine='by_coords')

outp.to_netcdf('Combined_data_ceres_const_rh.nc')
