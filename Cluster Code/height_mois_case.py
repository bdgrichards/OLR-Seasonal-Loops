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
    # Taking gas constant for dry air as 287.05 J kg^-1 K^-1 from Marshall and Plumb
    # Taking gas constant for water vapor as 461.52 J kg^-1 K^-1 from Wikipedia
    
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

def PolyArea(x,y):
    xvals = x[:]
    yvals = y[:]
        
    if xvals[len(xvals)-1] != xvals[0] and yvals[len(yvals)-1] != yvals[0]:
        xvals.append(xvals[0])
        yvals.append(yvals[0])
        
    sorting_list = []
    simple_polygons = []
    total_area = 0.
    
    sorting_list.append([xvals,yvals])
    
    while len(sorting_list) > 0:
        x1 = sorting_list[0][0]
        y1 = sorting_list[0][1]
        inters, intloc, ind1, ind2 = test_selfint(x1,y1)
        
        if inters:
            x2, x3 = break_poly(x1,intloc[0],ind1,ind2)
            y2, y3 = break_poly(y1,intloc[1],ind1,ind2)
            sorting_list.append([x2,y2])
            sorting_list.append([x3,y3])
            sorting_list.remove([x1,y1])
        else:
            simple_polygons.append([x1,y1])
            sorting_list.remove([x1,y1])
    
    for i in simple_polygons:
        total_area += SimplePolyArea(i[0],i[1])
    
    return total_area

def test_selfint(x,y):
    xvalst = x[:]
    yvalst = y[:]
    
    if xvalst[len(xvalst)-1] != xvalst[0] and yvalst[len(yvalst)-1] != yvalst[0]:
        xvalst.append(xvalst[0])
        yvalst.append(yvalst[0])
    
    for i in range(len(xvalst)-1):
        for j in range(len(xvalst)-1):
            if j > i:              
                xsec = [xvalst[i],xvalst[i+1],xvalst[j],xvalst[j+1]]
                ysec = [yvalst[i],yvalst[i+1],yvalst[j],yvalst[j+1]]
                intersect, intloc = isIntersect(xsec,ysec)
                if intersect:
                    return True, intloc, i, j
    else:
        return False, 0, 0, 0
    
def isIntersect(x,y):
    x1,x2,x3,x4 = x[0],x[1],x[2],x[3]
    y1,y2,y3,y4 = y[0],y[1],y[2],y[3]
    
    # Sets gradient as very large number to stop divide by 0 errors
    if x1 == x2:
        A1 = 1E100
    else:
        A1 = (y1-y2)/(x1-x2)
        
    if x3 == x4:
        A2 = 1E100
    else:
        A2 = (y3-y4)/(x3-x4)
    
    b1 = y1-A1*x1
    b2 = y3-A2*x3
    
    xaxis1 = np.linspace(x1,x2,100)
    xaxis2 = np.linspace(x3,x4,100)
    yaxis1 = A1*xaxis1 + b1
    yaxis2 = A2*xaxis2 + b2

    if (max(x1,x2) <= min(x3,x4)):
        return False, [0,0]

    if (A1 == A2):
        return False, [0,0]
       
    Xa = (b2 - b1) / (A1 - A2)
    
    buffer = 0.0001 # to fix floating point arithmetic errors 
    if ( (Xa <= buffer + max( min(x1,x2), min(x3,x4) )) or
         (Xa >= (buffer*(-1) + min( max(x1,x2), max(x3,x4))) ) ):
        return False, [0,0]
    else:
        return True, [Xa,A1*Xa+b1]

def break_poly(a,aint,aind1,aind2):
    anew1 = a[:aind1+1] + [aint] + a[1+aind2:]
    anew2 = [aint] + a[aind1+1:aind2+1] + [aint]
    return anew1, anew2

def SimplePolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

def lat_lon_index(latval, lonval):
    '''
    Uses the air temperature lat lon indexing system
    starts at lat = 90, lon = 0 and numbers left to right, one row at a time
    '''
    
    lats = np.array([ 90. ,  87.5,  85. ,  82.5,  80. ,  77.5,  75. ,  72.5,  70. ,
        67.5,  65. ,  62.5,  60. ,  57.5,  55. ,  52.5,  50. ,  47.5,
        45. ,  42.5,  40. ,  37.5,  35. ,  32.5,  30. ,  27.5,  25. ,
        22.5,  20. ,  17.5,  15. ,  12.5,  10. ,   7.5,   5. ,   2.5,
         0. ,  -2.5,  -5. ,  -7.5, -10. , -12.5, -15. , -17.5, -20. ,
       -22.5, -25. , -27.5, -30. , -32.5, -35. , -37.5, -40. , -42.5,
       -45. , -47.5, -50. , -52.5, -55. , -57.5, -60. , -62.5, -65. ,
       -67.5, -70. , -72.5, -75. , -77.5, -80. , -82.5, -85. , -87.5,
       -90. ])
    lons = np.array([  0. ,   2.5,   5. ,   7.5,  10. ,  12.5,  15. ,  17.5,  20. ,
        22.5,  25. ,  27.5,  30. ,  32.5,  35. ,  37.5,  40. ,  42.5,
        45. ,  47.5,  50. ,  52.5,  55. ,  57.5,  60. ,  62.5,  65. ,
        67.5,  70. ,  72.5,  75. ,  77.5,  80. ,  82.5,  85. ,  87.5,
        90. ,  92.5,  95. ,  97.5, 100. , 102.5, 105. , 107.5, 110. ,
       112.5, 115. , 117.5, 120. , 122.5, 125. , 127.5, 130. , 132.5,
       135. , 137.5, 140. , 142.5, 145. , 147.5, 150. , 152.5, 155. ,
       157.5, 160. , 162.5, 165. , 167.5, 170. , 172.5, 175. , 177.5,
       180. , 182.5, 185. , 187.5, 190. , 192.5, 195. , 197.5, 200. ,
       202.5, 205. , 207.5, 210. , 212.5, 215. , 217.5, 220. , 222.5,
       225. , 227.5, 230. , 232.5, 235. , 237.5, 240. , 242.5, 245. ,
       247.5, 250. , 252.5, 255. , 257.5, 260. , 262.5, 265. , 267.5,
       270. , 272.5, 275. , 277.5, 280. , 282.5, 285. , 287.5, 290. ,
       292.5, 295. , 297.5, 300. , 302.5, 305. , 307.5, 310. , 312.5,
       315. , 317.5, 320. , 322.5, 325. , 327.5, 330. , 332.5, 335. ,
       337.5, 340. , 342.5, 345. , 347.5, 350. , 352.5, 355. , 357.5])                
    lat_ind = np.where(lats == latval)[0][0]
    lon_ind = np.where(lons == lonval)[0][0]
    
    ind = lat_ind*len(lons) + lon_ind
    return ind



#******************************************************************************



def compute_height():
    
    #=====================================================
    # Import Data
    #=====================================================
    
    path = './Data/'
    
    temp = xr.open_dataset(path + 'temp.nc', decode_times=False)
    rhum = xr.open_dataset(path + 'rhum.nc', decode_times=False)
    pres = xr.open_dataset(path + 'pres.nc')
    ts2m = xr.open_dataset(path + 'ts_2m.nc', decode_times=False)
    
    lats = temp.lat.values
    lons = temp.lon.values
    levels = temp.level.values
    
    #=====================================================
    # Calculate the height
    #=====================================================
    
    month_list = np.arange(1,13)
    ts_monav = np.array(ts2m.sel(lat=latval, lon=lonval).air.values)
    temp_profile_monav_raw = np.array(temp.sel(lat=latval, lon=lonval).air.values)
    rhum_profile_monav_raw = np.array(rhum.sel(lat=latval, lon=lonval).rhum.values)
    pres_in = pres.sel(lat=latval, lon=lonval).pres.values[()]
    olr_monav = []

    levels_new, temp_profile_monav, rhum_profile_monav = truncate(temp_profile_monav_raw, rhum_profile_monav_raw, levels, pres_in)
    
    shum_profile_monav = []
    
    for k in range(len(rhum_profile_monav[:,0])):
        shum_section = []
        for l in range(len(rhum_profile_monav[0,:])):
            shum_section.append(q_from_rh(rh=rhum_profile_monav[k,l], T=temp_profile_monav[k,l], p=levels[l]))
        shum_profile_monav.append(shum_section)

    shum_profile_monav = np.array(shum_profile_monav)
    
    ts_mean = np.mean(ts_monav)
    ts_diff = [tval - ts_mean for tval in ts_monav]
    
    for i_mon in range(len(month_list)):
        temp_vals_mean = np.mean(temp_profile_monav, axis=0)
        temp_vals = [dat + ts_diff[i_mon] for dat in temp_vals_mean]
        shum_vals = shum_profile_monav[i_mon]

        #  State variables (Air and surface temperature)
        state = climlab.column_state(lev=levels_new)
        state.Tatm[:] = np.flip(temp_vals)
        state.Ts[:] = ts_monav[i_mon]

        #  Humidities
        h2o = climlab.radiation.water_vapor.ManabeWaterVapor(qStrat=5e-6, state=state);
        h2o.q[len(h2o.q)-len(shum_vals):len(h2o.q)] = np.flip(shum_vals)

        #  Couple water vapor to radiation
        rad = climlab.radiation.RRTMG(state=state, specific_humidity=h2o.q)

        #  Run the model
        rad.compute();
        olr_monav.append(rad.OLR[0])

    height_val = PolyArea(list(ts_monav),olr_monav) / (max(ts_monav) - min(ts_monav))
    
    h_data = xr.DataArray(height_val, dims=('lat', 'lon'), coords={'lat': [latval], 'lon': [lonval]})
    h = xr.Dataset({'height': h_data})

    h_index = lat_lon_index(latval, lonval)

    h.to_netcdf('./Output_Data_mois_case/'+str(h_index)+'.nc')
    
    return


#*******************************************************************************



## allow compute() to be called from command line:  #
##    USAGE:
##    python compute_stuff.py 1. 2.

if __name__ == "__main__":   # this bit will only get executed if I run the file from command line
    latval = float(sys.argv[1])/2
    lonval = float(sys.argv[2])/2
    compute_height()
