import numpy as np
import sys
import climlab
import xarray as xr

def sat_vap_pres(T):
    '''
    T must be in K, returns pressure in Pa - checked against online calculators to be correct
    '''
    return 100*6.112*np.exp(17.67*(T-273.15)/((T-273.15)+243.5))

def q_from_rh(rh, T, p): 
    # Taking gas constant for dry air as 287.05 J kg^-1 K^-1
    # Taking gas constant for water vapor as 461.52 J kg^-1 K^-1 
    
    # This is the equation at the end of the cell above
    q = (rh/100)*(287.05/461.52) * (sat_vap_pres(T)/(100*p))
    
    return q

def truncate(temps, rhums, levels, pres):
    
    for i in range(len(levels)):
        if pres >= levels[i]:
            arr_index = i
            break
            
    levels_arr = np.array(levels)
    rhums_arr = np.array(rhums)
    temps_arr = np.array(temps)
    
    levels_new = levels_arr[arr_index:]
    temps_new = temps_arr[:,arr_index:]
    rhums_new = rhums_arr[:,arr_index:]
    
    return levels_new, temps_new, rhums_new 

def lat_lon_index(latval, lonval, latsvals, lonsvals):
    '''
    A lat - lon indexing system
    starts at lat = 90, lon = 0 and numbers left to right, one row at a time
    '''
   
    lat_ind = np.where(latsvals == latval)[0][0]
    lon_ind = np.where(lonsvals == lonval)[0][0]
    
    ind = lat_ind*len(lonsvals) + lon_ind
    return ind

#******************************************************************************

def compute():
    
    # Values on a 1x1 degree grid
    lats = np.flip([float(i) for i in np.arange(-90,91,1)])
    lons = np.array([float(i) for i in np.arange(0,360,1)])
    
    # Import Data, ignoring 2000 because it starts in March
    data_t = xr.open_dataset('./Data/ERA5/era5_t_2000_2020.nc').sel(latitude=lats, longitude=lons).sel(time=slice("2001-01-01", "2020-12-30"))
    data_ts = xr.open_dataset('./Data/ERA5/era5_t2m_2000_2020.nc').sel(latitude=lats, longitude=lons, expver=1).sel(time=slice("2001-01-01", "2020-12-30"))
    data_q = xr.open_dataset('./Data/ERA5/era5_q_2000_2020.nc').sel(latitude=lats, longitude=lons).sel(time=slice("2001-01-01", "2020-12-30"))
    data_sp = xr.open_dataset('./Data/ERA5/era5_sp_2000_2020.nc').sel(latitude=lats, longitude=lons, expver=1).sel(time=slice("2001-01-01", "2020-12-30"))
    
    levels = data_t.level.values
    month_list = np.arange(1,13)
    
    olr_monav = [] 
    ts_monav = []
    
    for i_mon in range(len(month_list)):

        times_list = data_t.sel(time=data_t['time.month'] == month_list[i_mon]).time.values
        
        olr_list = []

        for j in range(len(times_list)):
            
            # Get the correct pressure levels, when surface pressure is < 1000 mBar
            surf_p = data_sp.sel(latitude=latval, longitude=lonval, time=times_list[j]).sp.values[()]/100
            levels_correct = data_t.sel(level=slice(0,surf_p)).level.values
            
            # Create a column state with the correct pressure levels
            state = climlab.column_state(lev=levels_correct)
            
            # Correct lowest pressure level from 1000 mbar set by default
            state.Tatm.domain.lev.bounds[-1] = surf_p
            state.Tatm.domain.lev.delta[-1] = state.Tatm.domain.lev.bounds[-1] - state.Tatm.domain.lev.bounds[-2]

            # Temperatures
            state.Tatm[:] = data_t.sel(latitude=latval, longitude=lonval, time=times_list[j], level=slice(0,surf_p)).t.values
            state.Ts[:] = data_ts.sel(latitude=latval, longitude=lonval, time=times_list[j]).t2m.values
            
            # Humidities
            h2o = climlab.radiation.water_vapor.ManabeWaterVapor(state=state);
            h2o.q = data_q.sel(latitude=latval, longitude=lonval, time=times_list[j], level=slice(0,surf_p)).q.values

            # Couple water vapor to radiation
            rad = climlab.radiation.RRTMG(state=state, specific_humidity=h2o.q)

            # Run the model
            rad.compute();
            
            # Save output
            olr_list.append(rad.OLR[0])

        olr_monav.append(np.mean(olr_list))
        ts_monav.append(np.mean(data_ts.sel(latitude=latval, longitude=lonval, time=data_ts['time.month'] == month_list[i_mon]).t2m.values))

    # Empty arrays to put data in
    ts_input = np.zeros((1,1,12))
    olr_input = np.zeros((1,1,12))
    
    for i in range(len(ts_input[0,0,:])):
        ts_input[0,0,i] = ts_monav[i]
        olr_input[0,0,i] = olr_monav[i]

    ts_data = xr.DataArray(ts_input, dims=('lat', 'lon', 'month'), coords={'lat': [latval], 'lon': [lonval], 'month':month_list})
    olr_data = xr.DataArray(olr_input, dims=('lat', 'lon', 'month'), coords={'lat': [latval], 'lon': [lonval], 'month':month_list})

    h = xr.Dataset({'ts': ts_data, 'olr_calc':olr_data})
    
    # A unique identifier for the output file, for later concatenation
    h_index = lat_lon_index(latval, lonval, lats, lons)
    h.to_netcdf('./Output_Data_ceres/'+str(h_index)+'.nc')
    
    return

#*******************************************************************************

## allow compute() to be called from command line:  #
##    USAGE:
##    python compute_stuff.py 1. 2.

if __name__ == "__main__":   # this bit will only get executed if I run the file from command line
    latval = float(sys.argv[1])
    lonval = float(sys.argv[2])
    compute()
