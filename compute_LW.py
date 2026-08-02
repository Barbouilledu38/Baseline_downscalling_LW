# -*- coding: utf-8 -*-

import xarray as xr
import numpy as np
import pandas as pd
from pyproj import Transformer
import numba as nb

from pathlib import Path
import sys
import time
import rioxarray
from rasterio.enums import Resampling

############################ Functions ################################
    
sigma = 5.67e-8

@nb.njit(cache = True)
def altit_G(
    dem : np.ndarray,
    z0 : float | np.ndarray,
    T0 : float | np.ndarray,
    lapse_rate : float,
    )-> float | np.ndarray:
        
    return T0 + (dem-z0)*lapse_rate

@nb.njit(cache = True)
def LW(
    LW_AROME : np.ndarray,
    T2m_AROME : np.ndarray,
    z_AROME : np.ndarray,
    z_DEM : np.ndarray,
    SVF : np.ndarray,
    eps_a : np.ndarray,
    lapse_rate = -6.5e-3,
    )-> np.ndarray:
    
    """
    [Input]
    - LW_AROME : np.ndarray (nx_A,ny_A) of AROME downwelling LW flux
    - T2m_AROME : np.ndarray (nx_A,ny_A) of AROME temperature
    - zs : np.ndarray (nx,ny) of goal elevation
    - SVF : np.ndarray (nx,ny) of SVF
    - lapse_rate : float of standard altitudinal gradiant of temperature
    - force_LWu : option to add (True) the neighboor contribution to the downwelling LW flux
    
    [Output]
    - LW : np.ndarray (nx,ny) of downwelling LW flux
    """
    
    # Downscall the air temperature from the AROME topographie based on a linear interpolation 
    # with a standard lapse of -6.5 K/km
    T2m_lapse_rate = altit_G(T0 = T2m_AROME, z0 = z_AROME, dem =  z_DEM, lapse_rate = lapse_rate)
    
    # Atmospherical LW is as emitted by a black body of same surface and emissivity eps_a_fine
    sigma = 5.67*1e-8
    LWd_fine = eps_a*sigma*T2m_lapse_rate**4
    
    return LWd_fine*SVF
    
##################################### Baseline ############################################

def baseline_LW(
    save_name : str,
    forcing_file : str,
    topo_params_file : str,
    run_dir : str = ".") :

    path_run_dir = Path(run_dir)

    # Loading forcing file in epsg 4326, topographic parameters file in epsg 2154
    ds_forçage = xr.open_dataset(path_run_dir / forcing_file)
    ds_topo = xr.open_dataset(path_run_dir / topo_params_file)
    
    # Creating a new variable eps_a in forcage xr.Dataset before reggridding
    ds_forçage["eps_a"] = (["latitude","longitude"], ds_forçage.LWdown.values/(sigma*ds_forçage.Tair.values**4))
    
    # Changing the coordinates to epsg 2154, the epsg of the goal DEM
    ## Writing crs for the projection
    ds_forçage = ds_forçage.rio.write_crs("EPSG:4326")
    ds_topo = ds_topo.rio.write_crs("EPSG:2154")
    
    ## Projecting/interpolating the coarse AROME on the fine DEM with rio.reproject_match
    # and bilinear option. Now calculation may be made easily
    ds_forçage = ds_forçage.rio.reproject_match(ds_topo,
                                                resampling=Resampling.bilinear)
    
    # The forcing file contains only a date, no time coordinate
    LW_fine = LW(LW_AROME = ds_forçage.LWdown.values,
                        T2m_AROME = ds_forçage.Tair.values,
                        z_AROME = ds_forçage.ZS.values,
                        z_DEM = ds_topo.ZS.values,
                        SVF = ds_topo.svf.values,
                        eps_a = ds_forçage.eps_a.values)
                        
    # Sauvegarde sous .nc
    ds_LW = xr.Dataset(
        data_vars=dict(
            LW=(["y", "x"], LW_fine)
        ),
        coords=dict(
            x=("x", ds_topo.x.values),
            y=("y", ds_topo.y.values),
        )
    )
    ds_LW.LW.attrs={'units': '$W.m^{-2}$', 'standard_name': 'LWd', 'long_name': 'Downwelling Long Wave'}
    
    ds_LW.to_netcdf(path_run_dir / save_name)

########################### Arguments ######################################

if len(sys.argv) != 4:
    print("Usage: python3 compute_LW.py run_dir forcing.nc nom_experience")
    sys.exit(1)
    
start_time = time.time()

baseline_LW(
    run_dir = sys.argv[1],
    forcing_file = sys.argv[2],
    topo_params_file = f"topo_params_{sys.argv[3]}.nc",
    save_name = f"LW_{sys.argv[3]}.nc")
    
end_time = time.time()

print(f"Computing the downscalled LW flux : {round(end_time - start_time,3)} s")
