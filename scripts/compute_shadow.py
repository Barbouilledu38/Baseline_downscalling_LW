# -*- coding: utf-8 -*-

from utils import *

from topocalc import gradient
from topocalc import viewf
from topocalc import horizon

import pvlib

############################################# Core function ###############################################

############# Position soleil et date #################

def from_dates_to_solar_angles(
    fic_solar_pos : str,
    fic_dem : str,
    time_slice,
    run_dir : str = "." ,
    ):
    
    path_run_dir = Path(run_dir)
    
    ds_topo = xr.open_dataset(path_run_dir / fic_dem)
    
    xx,yy = np.meshgrid(ds_topo.x.values,ds_topo.y.values)
    lons,lats = convert_epsg_pts(xx,yy, epsg_src=2154, epsg_tgt=4326)
    
    nt,nx,ny = len(time_slice),lons.shape[0],lons.shape[1]
    solar_pos = np.zeros((nt, nx, ny,2))
    
    for t, timestamp in enumerate(time_slice): 
        
        # création d’un tableau de dates avec le même forme que lats.ravel()
        time_array = np.full(lats.ravel().shape, timestamp, dtype='datetime64[ns]')
    
        solar = pvlib.solarposition.get_solarposition(
            latitude  = lats.ravel(),   # tableau 1D de tous les pixels
            longitude = lons.ravel(),
            time      = time_array,
        )
        solar_pos[t, :, :, 0] = solar["azimuth"].values.reshape(nx, ny)
        solar_pos[t, :, :, 1] = solar["elevation"].values.reshape(nx, ny)
                
    ds_solar = xr.Dataset(
        data_vars=dict(
            elevation=(["time","y", "x"], solar_pos[:,:,:,1]),
            azimuth=(["time","y", "x"], solar_pos[:,:,:,0]),
        ),
        coords=dict(
            time=("time", time_slice),
            x=("x", ds_topo.x.values),
            y=("y", ds_topo.y.values),
        )
    )
    ds_solar.elevation.attrs={'units': 'degre', 'standard_name': 'elevation', 'long_name': 'solar elevation'}
    ds_solar.azimuth.attrs = {'units': 'degre', 'standard_name': 'azimuth', 'long_name': 'solar azimuth'}
    
    ds_solar.to_netcdf(path_run_dir / fic_solar_pos)
    
######################## Projection et ombres ######################

@nb.njit(cache=True)
def proj(
    slope : float,
    aspect : float,
    azimuth : float,
    elevation : float,
    )-> float:
    
    fact = np.cos(np.pi/2 - elevation)*np.cos(slope) + np.sin(np.pi/2 - elevation)*np.sin(slope)*np.cos(azimuth - aspect)
    
    if fact > 0 : # fact < 0 : à l'ombre de sa propre pente
        return fact
    
    else :
        return 0

def horizon_angles_tab(
    DEM : np.ndarray,
    spacing : float,
    azimuths = np.arange(360),
    )->np.ndarray :

    a,b = DEM.shape
    cosha = np.zeros((len(azimuths),a,b))

    for i,azimuth in enumerate(azimuths):
        
        cosha[i,:,:] = horizon.horizon(azimuth-180,DEM.astype(np.float64),spacing) # Attention : où bien placer l'azimth zéro ?
        
    return cosha

@nb.njit(cache=True) 
def shadow(
    cosha : np.ndarray, # en degrés, 3D
    azimuths_t : np.ndarray, # en degrés, 3D
    elevations_t : np.ndarray, # cos(elevation) en degrés, 3D
    slope : np.ndarray,
    aspect : np.ndarray,
    )-> np.ndarray :
    
    _,a,b = cosha.shape
    res = np.ones((a,b))
    
    for i in range(a):
        for j in range(b):
            
            azi_index = int(round(azimuths_t[i,j])) % 360
                            
            if cosha[azi_index,i,j] > elevations_t[i,j] :
                res[i,j] = 0
                
            else :
                res[i,j] = proj(
                                slope = slope[i,j],
                                aspect = aspect[i,j],
                                azimuth = azimuths_t[i,j],
                                elevation = elevations_t[i,j])
    return res 

@nb.njit(cache=True)
def shadow_tab(
    cosha : np.ndarray, # np.ndarray (360,nx,ny)
    elevations : np.ndarray, # np.cos(elevation) np.ndarray (nt,nx,ny)
    azimuths : np.ndarray, # np.ndarray (nt,nx,ny)
    slope : np.ndarray, # np.ndarray (nx,ny)
    aspect : np.ndarray, # np.ndarray (nx,ny)
    )-> np.ndarray:
    
    nt,nx,ny = elevations.shape
    shad = np.ones((nt,nx,ny))
    
    for t in range(nt):
        
        shad[t,:,:] = shadow(
            cosha = cosha,
            azimuths_t = azimuths[t],
            elevations_t = elevations[t],
            slope = slope,
            aspect = aspect)
        
    return shad

def shadow_dataset(
    fic_topo_params : str = "topo_params.nc",
    fic_solar_pos : str = "ds_solar.nc",
    fic_ha : str = "horizon_angle.nc",
    fic_shadow : str = "ds_shadow.nc",
    run_dir: str  = "." ,
    force_recompute_solar = False,
    force_recompute_horizon = False
    ):
    
    # Constructing a time slice with pandas for all the year 2026
    start = pd.Timestamp('2026-01-01')
    end   = pd.Timestamp('2026-01-03')
    slice_t = slice(start, end)
    time_slice = pd.date_range(start, end, freq='h')
    
    path_run_dir = Path(run_dir)
    
    ds_topo = xr.open_dataset(path_run_dir / fic_topo_params)
    
    dx = np.median(np.diff(ds_topo.x.values)) 
    xx, yy = np.meshgrid(ds_topo.x.values,ds_topo.y.values)
    
    # Changement systeme de coordonnées de Lambert 93 (epsg 2154) à epsg 4326 pour from_dates_to_solar_angles
    lons,lats = convert_epsg_pts(xx,yy, epsg_src=2154, epsg_tgt=4326)
    
    solar_exists = Path(path_run_dir / fic_solar_pos).exists() and not force_recompute_solar
    if not solar_exists :
    
        path = Path(fic_solar_pos)
        path.unlink(missing_ok=True)
        
        from_dates_to_solar_angles(
            run_dir = run_dir,
            fic_dem = fic_topo_params,
            fic_solar_pos = fic_solar_pos,
            time_slice = time_slice)
    
    ds_solar = xr.open_dataset(fic_solar_pos)
    
    print("Solar features OK")
    
    horizon_exists = Path(fic_ha).exists() and not force_recompute_horizon
    if not horizon_exists :
        
        path = Path(fic_ha)
        path.unlink(missing_ok=True)
    
        cosha = horizon_angles_tab(
                    DEM = ds_topo.ZS.values,
                    spacing = dx)*(180/np.pi)

        # Sauvegarde sous .nc
        ds_ha = xr.Dataset(
            data_vars=dict(
                cosha=(["azimuth","y", "x"], cosha)
            ),
            coords=dict(
                azimuth=("azimuth", np.arange(360)),
                y=("y", ds_topo.y.values),
                x=("x", ds_topo.x.values),
            )
        )
        ds_ha.cosha.attrs={'units': 'degres','standard_name': 'Horizon angle', 'long_name': 'Horizon angle'}
        ds_ha.to_netcdf(fic_ha)
        
    ds_ha = xr.open_dataset(fic_ha)
    
    print("Horizon angles OK")
    
    shadows = shadow_tab(
        cosha = ds_ha.cosha.values,
        elevations = np.cos(ds_solar.elevation.values*(np.pi/180)),
        azimuths = ds_solar.azimuth.values,
        slope = ds_topo.slope.values,
        aspect = ds_topo.aspect.values)
    
    # Sauvegarde sous .nc
    ds_shadow = xr.Dataset(
        data_vars=dict(
            shadow_mask=(["time","y", "x"], shadows[:,:,:])
        ),
        coords=dict(
            time=("time", time_slice),
            x=("x", ds_topo.x.values),
            y=("y", ds_topo.y.values),
        )
    )
    ds_shadow.shadow_mask.attrs={'standard_name': 'Topo SW proj', 'long_name': 'Topographic SW projection factor'}
    ds_shadow.to_netcdf(fic_shadow)
    
    print("Shadow mask OK")

################################## Arguments ############################################

if len(sys.argv) != 8:
    print("Usage: python3 compute_shadow.py run_dir fic_topo.nc ds_solar.nc fic_ha.nc fic_shadow.nc force_recompute_solar force_recompute_horizon_angle")
    sys.exit(1)
    
shadow_dataset(
    run_dir = sys.argv[1],
    fic_topo_params = sys.argv[2],
    fic_solar_pos = sys.argv[3],
    fic_ha = sys.argv[4],
    fic_shadow = sys.argv[5],
    force_recompute_solar = sys.argv[6],
    force_recompute_horizon = sys.argv[7])
