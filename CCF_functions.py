def Surface_Temp_Advection(Us,Vs,Ts):
    ## In this case we use the xarray differentiate function
    ## as such all data must be xarray dataarrays
    ## U and V are in units of M/s
    ## Hence, this function ouputs in K/s
    a = 6.371e6 ## Radius of Earth in Meters
    lat = Ts.lat
    lon = Ts.lon
    ## interpolate any missing values
    Ts = Ts#.interpolate_na('lat','linear')
    Us = Us#.interpolate_na('lat','linear')
    Vs = Vs#.interpolate_na('lat','linear')
    Ts = Ts.assign_coords({'lat':lat*np.pi/180,'lon':lon*np.pi/180}) ## convert to radians
    lat_derivative = Ts.differentiate('lat')## central differences method
    lon_derivative = Ts.differentiate('lon')## Must interpolate back to original grid so compatible with next calculation
    lon = lon_derivative.lon*180/np.pi ## convert back to degrees
    lat = lat_derivative.lat*180/np.pi
    lat_derivative = lat_derivative.assign_coords({'lat':lat,'lon':lon})
    lon_derivative = lon_derivative.assign_coords({'lat':lat,'lon':lon})
    lat_radians = lat * np.pi/180 ## need to still save radians one last time
    lon_radians = lon * np.pi/180
    longrid,latgrid = np.meshgrid(lon_radians,lat_radians)
    coeff_lon = Us.interp({'lon':lon,'lat':lat}).assign_coords({'lat':lat,'lon':lon})/(a*np.cos(latgrid))
    coeff_lat = Vs.interp({'lon':lon,'lat':lat}).assign_coords({'lat':lat,'lon':lon})/a
    T_adv = -coeff_lon*lon_derivative - coeff_lat*lat_derivative
    return(T_adv)


    ## Write funcitons for the various terms in the EIS equation

def Potential_Temperature(T,P):
    ## T in Kelvin, P in Pa
    R_cp = 0.286
    P0 = 100000 ## Pa
    theta = T*(P0/P)**R_cp
    return(theta)

def LTS(Ts,Ps,T_700):
## T in Kelvin, P in Pa
    R_cp = 0.286
    P0 = 100000 ## Pa
    theta_700 = T_700*(P0/70000)**R_cp
    theta_s = Ts*(P0/Ps)**R_cp
    lts = theta_700 - theta_s
    return(lts)

def Sat_vapor_pres(T):
    ## T in Kelvin
    A = 2.53e11 #Pa
    B = 5420 #K
    es = A*np.exp(-B/T)
    return(es)

def Sat_mixing_ratio(T,P):
    ## T in Kelvin, P in Pa
    A = 2.53e11 #Pa
    B = 5420 #K
    Rd = 287.053 #J/kg/K
    Rv = 461.52 #J/kg/K
    # Calculate the saturation vapor pressure using a Clausius-Clapeyron numerical solution
    es = A*np.exp(-B/T)
    qs = (Ra/Rv)*es/(P-es)
    return(qs)

def Moist_Adiabatic_LR(T,P):
    ## T in Kelvin, P in Pa
    A = 2.53e11 #Pa
    B = 5420 #K
    Rd = 287.053 #J/kg/K
    Rv = 461.52 #J/kg/K
    g = 9.81 #m/s/s
    cp = 1005 #J/kg/K
    Lv = 2.5e6 #J/kg
    
    # Calculate the saturation vapor pressure using a Clausius-Clapeyron numerical solution
    es = A*np.exp(-B/T)
    qs = (Rd/Rv)*es/(P-es)
    ## Now calculate the moist adiabatic lapse rate
    gamma = (g/cp)*(1 - (1+Lv*qs/(Rd*T))/(1+Lv**2*qs/(cp*Rv*T**2)))
    return(gamma)

def Relative_Humidity(T,P,q):
    A = 2.53e11 #Pa
    B = 5420 #K
    Rd = 287.053 #J/kg/K
    Rv = 461.52 #J/kg/K
    # Calculate the saturation vapor pressure using a Clausius-Clapeyron numerical solution
    es = A*np.exp(-B/T)
    # Now calculate the actual vapor pressure
    e = P*(q/((Rd/Rv)+q))
    RH = e/es
    return(RH)

def z_LCL(T,P,q):
    #First have to calculate the relative humidity
    A = 2.53e11 #Pa
    B = 5420 #K
    Rd = 287.053 #J/kg/K
    Rv = 461.52 #J/kg/K
    # Calculate the saturation vapor pressure using a Clausius-Clapeyron numerical solution
    es = A*np.exp(-B/T)
    # Now calculate the actual vapor pressure
    e = P*(q/((Rd/Rv)+q))
    ## Then relative humidity is the fraction of vapor pressure over saturation
    RH = e/es
    ## Now calculate the dew point temperature using this numerical approximation
    Td = T - (100-RH)/5
    ## Now the LCL with another approximation from Lawrence et al. (2005)
    LCL = (T - Td)/8
    return LCL
    
def Estimated_Inversion_Strength(Ts, T_700, Ps, q, z700):
    ## First calculate the LTS
    R_cp = 0.286
    P0 = 100000 ## Pa
    theta_700 = T_700*(P0/70000)**R_cp
    theta_s = Ts*(P0/Ps)**R_cp
    lts = theta_700 - theta_s
    
    
    ## Now calculate the moist adiabatic LR at 850 hPa
    T_850 = (T_700+Ts)/2
    A = 2.53e11 #Pa
    B = 5420 #K
    Rd = 287.053 #J/kg/K
    Rv = 461.52 #J/kg/K
    g = 9.81 #m/s/s
    cp = 1005 #J/kg/K
    Lv = 2.5e6 #J/kg
    # Calculate the saturation vapor pressure using a Clausius-Clapeyron numerical solution
    es = A*np.exp(-B/T_850)
    qs = (Rd/Rv)*es/(85000-es)
    ## Now calculate the moist adiabatic lapse rate
    gamma_850 = (g/cp)*(1 - (1+Lv*qs/(Rd*T_850))/(1+Lv**2*qs/(cp*Rv*T_850**2)))
    
    
    ## Now calculate the LCL from surface variables
    es = A*np.exp(-B/Ts)
    # Now calculate the actual vapor pressure
    e = Ps*(q/((Rd/Rv)+qs))
    ## Then relative humidity is the fraction of vapor pressure over saturation
    RH = e/es
    ## Now calculate the dew point temperature using this numerical approximation
    Td = Ts - (100-RH)/5
    ## Now the LCL with another approximation from Lawrence et al. (2005)
    LCL = (Ts - Td)/8
    
    
    ## Now put it all together
    EIS = lts - gamma_850*(z_700 - LCL)
    return(EIS)
    
