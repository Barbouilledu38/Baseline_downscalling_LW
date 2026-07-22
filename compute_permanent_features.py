# -*- coding: utf-8 -*-

############################################# To modify #########################################

"""
Indétermination si degré ou radiant étant un obstacle à la compréhension du code
A commenter avec générosité
Le caclcul des angles d'horizon peut être utilisé pour retourer le SVF plutôt que le calculer séparément
"""

############################################# Imports #########################################

from utils import *

from topocalc import gradient
from topocalc import viewf
from topocalc import horizon

import pvlib

############################################# SVF #########################################

def compute_svf_regular_grid(
    run_dir: str = "." ,
    fic_topo : str  = "dem.nc",
    save_name : str = "topo_params.nc",
    ):
    """
    Compute SVF for the whole .nc topography in input
    """
    
    path_run_dir = Path(run_dir)
        
    ds_topo = xr.open_dataset(path_run_dir / fic_topo)
        
    # La projection Lambert 93 (ou epsg 2154) est en mètre donc conforme pour le calcul, pas besoin de changer en 32632
    #xs, ys = convert_epsg_pts(ds_topo.x.values,ds_topo.x.values, epsg_src = epsg_init, epsg_tgt=32632)
    dx = np.median(np.diff(ds_topo.x.values))
    dy = np.median(np.diff(ds_topo.y.values))

    svf = viewf.viewf(np.double(ds_topo.ZS.values), dx)[0]
    slope, aspect = gradient.gradient_d8(ds_topo.ZS.values, dx, dy)
    
    # Sauvegarde sous .nc
    ds_topo_params = xr.Dataset(
        data_vars=dict(
            ZS=(["y", "x"], ds_topo.ZS.values.astype(np.float64)),
            svf=(["y", "x"], svf.astype(np.float64)),
            slope=(["y", "x"], slope.astype(np.float64)),
            aspect=(["y", "x"], aspect.astype(np.float64)),
        ),
        coords=dict(
            x=("x", ds_topo.x.values),
            y=("y", ds_topo.y.values),
        )
    )
    
    ds_topo_params.ZS.attrs={'units': 'm', 'standard_name': 'elevation', 'long_name': 'elevation'}
    ds_topo_params.x.attrs = {'units': 'm'}
    ds_topo_params.y.attrs = {'units': 'm'}
    ds_topo_params.slope.attrs = {'units': 'rad'}
    ds_topo_params.aspect.attrs = {'units': 'rad'}
    ds_topo_params.svf.attrs = {'units': 'ratio', 'standard_name': 'svf', 'long_name': 'Sky view factor'}
    
    ds_topo_params.to_netcdf(path_run_dir / save_name)
    
    print("Topographic parameters OK")
    
############################################# Solar parameters #########################################

def from_dates_to_solar_angles(
    time_slice,
    run_dir : str,
    fic_dem : str,
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
    
    return ds_solar
    
############################################# Projection et ombres #########################################

@nb.njit(cache=True)
def proj(
    slope : float,
    aspect : float,
    azimuth : float,
    elevation : float,
    )-> float:
    
    fact = np.cos(np.pi/2 - elevation*(np.pi/180))*np.cos(slope) + np.sin(np.pi/2 - elevation*(np.pi/180))*np.sin(slope)*np.cos(azimuth*(np.pi/180) - aspect)


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
    ha = np.zeros((len(azimuths),a,b))

    for i,azimuth in enumerate(azimuths):
        
        # azimuth-180 because the function horizon.horizon must receive azimuth bewteen -180/180°. With the 0 beeing the north
        ha[i,:,:] = horizon.horizon(azimuth-180,DEM.astype(np.float64),spacing)
        
    return ha

@nb.njit(cache=True) 
def shadow(
    ha : np.ndarray, # en degrés, 3D
    azimuths_t : np.ndarray, # en degrés, 2D
    elevations_t : np.ndarray, # en degrés, 2D
    slope : np.ndarray, # radiant
    aspect : np.ndarray, # radiant
    )-> np.ndarray :
    
    _,a,b = ha.shape
    res = np.ones((a,b))
    
    for i in range(a):
        for j in range(b):
            
            azi_index = int(round(azimuths_t[i,j])+180) % 360
                            
            if ha[azi_index,i,j] > elevations_t[i,j] :
                res[i,j] = 0
                
            else :
                res[i,j] = proj(slope = slope[i,j],
                                aspect = aspect[i,j],
                                azimuth = azimuths_t[i,j],
                                elevation = elevations_t[i,j])
    return res 

@nb.njit(cache=True)
def shadow_tab(
    ha : np.ndarray, # np.ndarray (360,nx,ny)
    elevations : np.ndarray, # np.cos(elevation) np.ndarray (nt,nx,ny)
    azimuths : np.ndarray, # np.ndarray (nt,nx,ny)
    slope : np.ndarray, # np.ndarray (nx,ny)
    aspect : np.ndarray, # np.ndarray (nx,ny)
    )-> np.ndarray:
    
    nt,nx,ny = elevations.shape
    shad = np.ones((nt,nx,ny))
    
    for t in range(nt):
        
        shad[t,:,:] = shadow(
            ha = ha,
            azimuths_t = azimuths[t],
            elevations_t = elevations[t],
            slope = slope,
            aspect = aspect)
        
    return shad

def shadow_dataset(
    fic_topo_params : str = "topo_params.nc",
    fic_shadow : str = "ds_shadow.nc",
    run_dir: str  = "." ):
    
    path_run_dir = Path(run_dir)
    
    ############# Solar parameters ##################
    
    print("Computing solar parameters")
    
    # Constructing a time slice with pandas for all the year 2026
    start = pd.Timestamp('2026-01-01 00:00:00')
    end   = pd.Timestamp('2026-01-02 00:00:00')
    slice_t = slice(start, end)
    time_slice = pd.date_range(start, end, freq='h')
    
    path_run_dir = Path(run_dir)
    
    ds_topo = xr.open_dataset(path_run_dir / fic_topo_params)
    
    dx = np.median(np.diff(ds_topo.x.values)) 
    xx, yy = np.meshgrid(ds_topo.x.values,ds_topo.y.values)
    
    # Changement systeme de coordonnées de Lambert 93 (epsg 2154) à epsg 4326 pour from_dates_to_solar_angles
    lons,lats = convert_epsg_pts(xx,yy, epsg_src=2154, epsg_tgt=4326)

    ds_solar = from_dates_to_solar_angles(
        run_dir = run_dir,
        fic_dem = fic_topo_params,
        time_slice = time_slice)
        
    #ds_solar.to_netcdf(path_run_dir / "ds_solar.nc")
    
    print("Solar features OK")
    
    ############# Horizon angles ##################
    
    print("Computing horizon angles")
    
    ha = horizon_angles_tab(
                DEM = ds_topo.ZS.values,
                spacing = dx)

    # Sauvegarde sous .nc
    ds_ha = xr.Dataset(
        data_vars=dict(
            ha=(["azimuth","y", "x"], -np.arccos(ha)*(180/np.pi)+90)
        ),
        coords=dict(
            azimuth=("azimuth", np.arange(360)),
            y=("y", ds_topo.y.values),
            x=("x", ds_topo.x.values),
        )
    )
    ds_ha.ha.attrs={'units': 'degres','standard_name': 'Horizon angle', 'long_name': 'Horizon angle'}
    
    #ds_ha.to_netcdf(path_run_dir / "ds_ha.nc")
    
    print("Horizon angles OK")
    
    ############# Shadow mask ##################
    
    print("Computing shadow mask")
    
    shadows = shadow_tab(
        ha = ds_ha.ha.values, # degrees
        elevations = ds_solar.elevation.values, # degrees
        azimuths = ds_solar.azimuth.values, # degrees
        slope = ds_topo.slope.values, # radiant
        aspect = ds_topo.aspect.values) # radiant
    
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
    
    ds_shadow.to_netcdf(path_run_dir / fic_shadow)
    
    print("Shadow mask OK")
    
############################################# Calling function #########################################

if len(sys.argv) != 4:
    print("Usage: python3 compute_permanent_features.py run_dir dem.nc nom_experience")
    sys.exit(1)
    
# Computing svf
compute_svf_regular_grid(
    run_dir = sys.argv[1],
    fic_topo = sys.argv[2],
    save_name = f"topo_params_{sys.argv[3]}.nc")
    
# Computing shadow mask
shadow_dataset(
    run_dir = sys.argv[1],
    fic_topo_params = f"topo_params_{sys.argv[3]}.nc",
    fic_shadow = f"shadow_mask_{sys.argv[3]}.nc")
    

    
    
    

























