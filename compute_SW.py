# -*- coding: utf-8 -*-

############################################# Imports #########################################

from utils import *

############################################# SW #########################################

def Baseline_SW(
    run_dir : str = ".",
    fic_forcing : str = "forcing.nc",
    fic_topo_params : str = "topo_params.nc",
    fic_shadow : str = "shadow.nc",
    fic_SW_downscalled : str = "SW_downscalled.nc",
    x_dim_AROME = "longitude",
    y_dim_AROME = "latitude"):
    
    """
    Based on a spatial and temporal extent, select inside de shadow netcdf the SW projection
    map and multiply it by the SW forcing interpolated on the same grid. The projection grid 
    for diffus radiation is the Sky View Factor
    
    [Input]
    - run_dir : str = adresse of the run_dir
    - fic_forcing : str = name of forcing file in the run_dir
    - fic_topo_params : str = name of topo parameters file in the run_dir, precedently 
            calculated with compute_svf.py
    - fic_shadow : str = name of the shadow netcdf precedently calculated with compute_shadow.py
    - fic_SW_downscalled : str = name of the saved result
    
    [Output]
    - Save the downscalled SW forçing in a netcdf, indexed on the topo_params x and y coordinates 
            (Lambdert 93, i.e. epsg 2154)with variables SWdir and SWdif
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
        
    # Opening meteorological forçing file
    
    ds_forcage = xr.open_dataset(fic_forcing)   
    
    xx,yy = np.meshgrid(ds_forcage[x_dim_AROME].values,ds_forcage[y_dim_AROME].values)
    xs_, ys_ = convert_epsg_pts(xx,yy, epsg_src =4326, epsg_tgt=2154)
    xs_AROME, ys_AROME = xs_[0],ys_[:,0]
    
    xmin = min(xs_AROME)
    xmax = max(xs_AROME)
    ymin = min(ys_AROME)
    ymax = max(ys_AROME)
    
    print("Forcing file OK")
                
    # Opening topo params

    ds_topo_params = xr.open_dataset(fic_topo_params)
    
    print("Topographique params crop OK")
        
    # Opening shadow netcdf
    ds_shadow = xr.open_dataset(fic_shadow)
            
    # Les points pour les interpolations
    xx,yy = np.meshgrid(ds_topo.x.values,ds_topo.y.values)
    nx,ny = xx.shape
    points = np.column_stack((np.ravel(xx),np.ravel(yy)))
    
    ################## Calcul ###########################
    SW_fine = np.zeros((nx,ny,2))
        
    # dir : Interpolated points multiplicated by the shadow mask
    SW_fine[:,:,0] = interpolate( 
        xs = xs,
        ys = ys,
        arr = ds_forcage.SWdir.values, #.astype(np.float64),
        points = points).reshape(nx, ny)*ds_shadow_t.shadow_mask.values
        
    # dif : Interpolated points multiplicated by the sky view factor
    SW_fine[:,:,1] = interpolate( 
        xs = xs,
        ys = ys,
        arr = ds_forcage_t.SWdif.values, #.astype(np.float64),
        points = points).reshape(nx, ny)*ds_topo.svf.values
        
    # Sauvegarde sous .nc
    ds_SW = xr.Dataset(
        data_vars=dict(
            SWdir=(["time","y", "x"], SW_fine[:,:,0]),
            SWdif=(["time","y", "x"], SW_fine[:,:,1])
        ),
        coords=dict(
            x=("x", ds_topo.x.values),
            y=("y", ds_topo.y.values),
        )
    )
    
    ds_SW.SWdir.attrs={'units': '$W.m^{-2}$', 'standard_name': 'SWdir', 'long_name': 'Direct downwelling shortwave'}
    ds_SW.SWdif.attrs={'units': '$W.m^{-2}$', 'standard_name': 'SWdif', 'long_name': 'Diffuse downwelling shortwave'}
    
    ds_SW.to_netcdf(fic_res)
    
    print("SW OK")
    
############################################# Calling function #########################################

if len(sys.argv) != 4:
    print("Usage: python3 compute_SW.py run_dir forcing.nc nom_experience")
    sys.exit(1)
    
# Computing SW
Baseline_SW(
    run_dir = sys.argv[1],
    fic_forcing = sys.argv[2],
    fic_topo_params = f"topo_params_{sys.argv[3]}.nc",
    fic_shadow : f"shadow_mask_{sys.argv[3]}.nc",
    fic_SW_downscalled : str = f"SW_{sys.argv[3]}.nc",
    x_dim_AROME = "longitude",
    y_dim_AROME = "latitude")











