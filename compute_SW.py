# -*- coding: utf-8 -*-

############################################# Imports #########################################

import xarray as xr
import numpy as np
import numba as nb
from pathlib import Path
import sys
import time
import rioxarray
from rasterio.enums import Resampling

############################################# SW #########################################

def Baseline_SW(
    run_dir : str = ".",
    fic_forcing : str = "forcing.nc",
    fic_topo_params : str = "topo_params.nc",
    fic_shadow : str = "shadow.nc",
    fic_SW_downscalled : str = "ds_SW.nc",
    date_input : str = "2026010100",
    x_dim_AROME = "longitude",
    y_dim_AROME = "latitude"):
    
    """
    Multiply the SW forcing interpolated on the topographic parameters grid by the shadow_mask variable of the fic_shadow netcdf
    of the corresponding time. The projection grid for diffus radiation is the Sky View Factor. 
    
    The shadow mask file must contain only a date in the time coordinate in order to parallelize on the time coordinate.
    
    [Input]
    - run_dir : str = adresse of the run_dir
    - fic_forcing : str = name of forcing file in the run_dir
    - fic_topo_params : str = name of topo parameters file in the run_dir, precedently 
            calculated with compute_svf.py
    - fic_shadow : str = name of the shadow netcdf precedently calculated with compute_shadow.py
    - fic_SW_downscalled : str = name of the saved result
    
    [Output]
    - Save the downscalled SW forçing in a netcdf, indexed on the topo_params x and y coordinates 
            (Lambdert 93, i.e. epsg 2154) with variables SWdir and SWdif
    """
    
    run_dir = Path(run_dir)
    fic_forcing = run_dir / fic_forcing
    fic_topo_params = run_dir / fic_topo_params
    fic_shadow = run_dir / fic_shadow
    fic_SW_downscalled = run_dir / fic_SW_downscalled
    
    # Error on missing file in working directory
    
    try:
        xr.open_dataset(fic_forcing).close()
        xr.open_dataset(fic_topo_params).close()
        xr.open_dataset(fic_shadow).close()
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Fichier manquant dans le répertoire de run : {e}") from e
        
    # Loading forcing file in epsg 4326, topographic parameters file in epsg 2154
    ds_forcing = xr.open_dataset(fic_forcing)
    ds_topo = xr.open_dataset(fic_topo_params)
    
    # Changing the coordinates to epsg 2154, the epsg of the goal DEM
    ## Writing crs for the projection
    ds_forcing = ds_forcing.rio.write_crs("EPSG:4326")
    ds_topo = ds_topo.rio.write_crs("EPSG:2154")
    
    # Projecting/interpolating the coarse AROME on the fine DEM with rio.reproject_match
    # and bilinear option. Now calculation may be made easily
    ds_forcing = ds_forcing.rio.reproject_match(ds_topo,
                                                resampling=Resampling.bilinear)
                                                
    # Loading the shadow mask netcdf
    ds_shadow = xr.open_dataset(fic_shadow)
    
    # Selecting the date for the shadow mask 
    mm,dd,hh = date_input[4:6],date_input[6:8],date_input[8:10]
    ds_shadow_t = ds_shadow.sel(time=f"2026-{mm}-{dd}T{hh}:00:00.000000000", method = 'nearest')
                                               
    # save the result inside a netcdf        
    ds_SW = xr.Dataset(
        data_vars=dict(
            SWdir=(["y", "x"], ds_forcing.DIR_SWdown.values*ds_shadow_t.shadow_mask.values),
            SWdif=(["y", "x"], ds_forcing.SCA_SWdown.values*ds_topo.svf.values)
        ),
        coords=dict(
            x=("x", ds_topo.x.values),
            y=("y", ds_topo.y.values),
        )
    )
    
    ds_SW.SWdir.attrs={'units': '$W.m^{-2}$', 'standard_name': 'SWdir', 'long_name': 'Direct downwelling shortwave'}
    ds_SW.SWdif.attrs={'units': '$W.m^{-2}$', 'standard_name': 'SWdif', 'long_name': 'Diffuse downwelling shortwave'}
    
    ds_SW.to_netcdf(fic_SW_downscalled)
                                                
    print("SW OK")
    
############################################# Calling function #########################################

if len(sys.argv) != 5:
    print("Usage: python3 compute_SW.py run_dir forcing.nc nom_experience dateinput")
    sys.exit(1)
    
start_time = time.time()
    
# Computing SW
Baseline_SW(
    run_dir = sys.argv[1],
    fic_forcing = sys.argv[2],
    fic_topo_params = f"topo_params_{sys.argv[3]}.nc",
    fic_shadow = f"shadow_mask_{sys.argv[3]}.nc",
    fic_SW_downscalled = f"SW_{sys.argv[3]}.nc",
    date_input = sys.argv[4],
    x_dim_AROME = "longitude",
    y_dim_AROME = "latitude")

end_time = time.time()

print(f"Computing SW downscalling : {round(end_time - start_time,3)} s")









