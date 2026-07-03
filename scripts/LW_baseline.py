# -*- coding: utf-8 -*-

from utils import *

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
def Downscale_lw_u(
    Ts : np.ndarray, 
    LWd : np.ndarray, 
    xs : np.ndarray,
    ys : np.ndarray,
    delta = 2000,
    )-> np.ndarray:
    
    """
    [Input]
    - Ts : np.ndarray (nx,ny) of the surface temperature previously downscalled to the goal resolution
    - LWd : np.ndarray (nx,ny) of the surface LW downwelling flux without neighboors 
                illumination at goal resolution
    - xs, ys : np.ndarray (nx*ny) coordinates (m) of the goal resolution
    - delta : float describing the maximal distance considering illumination from neighboors in the LW

    [Output]
    - LWu : np.ndarray (nx,ny) neighboors contribution to LW flux
    
    """
    
    nx, ny = LWd.shape
    LWu = np.zeros((nx, ny))
    
    # The surface emissivity is considered equal to 1 whatever the ground nature
    eps_surf = np.full_like(LWd, 1)
    
    for i in range(nx):
        for j in range(ny):
            mask = (xs-xs[ny*i+j])**2 + (ys-ys[ny*i+j])**2 <= delta**2
            LWu[i,j] = np.mean((1-eps_surf[mask])*LWd[mask] + eps_surf[mask]*sigma*Ts[mask]**4)
            
    return LWu

#@nb.njit(cache = True)
def LW(
    LW_AROME : np.ndarray,
    T2m_AROME : np.ndarray,
    z_AROME : np.ndarray,
    x_AROME : np.ndarray,
    y_AROME : np.ndarray,
    zs : np.ndarray,
    xs : np.ndarray,
    ys : np.ndarray,
    SVF : np.ndarray,
    lapse_rate = -6.5e-3,
    force_LWu : bool = False,
    )-> np.ndarray:
    
    """
    [Input]
    - LW_AROME : np.ndarray (nx_A,ny_A) of AROME downwelling LW flux
    - T2m_AROME : np.ndarray (nx_A,ny_A) of AROME temperature
    - z_AROME,x_AROME,y_AROME : np.ndarray (nx_A,ny_A) of AROME coordinates
    - zs : np.ndarray (nx,ny) of goal elevation
    - xs,ys : np.ndarray (nx,ny) of the goal coordinates
    - SVF : np.ndarray (nx,ny) of SVF
    - lapse_rate : float of standard altitudinal gradiant of temperature
    - force_LWu : option to add (True) the neighboor contribution to the downwelling LW flux
    
    [Output]
    - LW : np.ndarray (nx,ny) of downwelling LW flux
    """
    
    nx,ny = SVF.shape
    
    # Points of fine grid on which to interpolate
    points = np.column_stack((xs, ys))
    
    # creating coarse (then fine by bilinear interpolation) array of atmospherique emissivity
    eps_a_coarse = LW_AROME/(sigma*T2m_AROME**4)
    eps_a_fine = interpolate(xs = x_AROME, ys = y_AROME, arr = eps_a_coarse, points = points)
    
    # Downscalling temperature with altitudinal gradiant
    z_AROME_fine = interpolate(xs = x_AROME, ys = y_AROME, arr = z_AROME, points = points)
    T2m_AROME_fine = interpolate(xs = x_AROME, ys = y_AROME, arr = T2m_AROME, points = points)
    T2m_fine = altit_G(T0 = T2m_AROME_fine, z0 = z_AROME_fine, dem =  zs.ravel(), lapse_rate = lapse_rate)
    
    # Atmospherical LW is as emitted by a black body of same surface and emissivity eps_a_fine
    LWd_fine = (eps_a_fine*sigma*T2m_fine**4).reshape(nx,ny)
    if force_LWu :
        LWu = Downscale_lw_u(
            Ts = T2m_fine, # By default, surface temperature is  set equal to 2m temperature
            LWd  = LWd_fine,
            xs = xs,
            ys = ys)

        return LWd_fine*SVF + LWu*(1-SVF)
    
    else :
        return LWd_fine*SVF,T2m_fine,T2m_AROME_fine,eps_a_fine
    
##################################### Baseline ############################################

def baseline_LW(
    save_name : str,
    forcing_file : str,
    topo_params_file : str,
    run_dir : str = "." ,
    force_LWu : bool = False):

    path_run_dir = Path(run_dir)

    # Loading forcing file in epsg 4326 
    ds_forçage = xr.open_dataset(path_run_dir / forcing_file)
    
    # Changing the format of the coordinate array for calculus
    xx_AROME, yy_AROME = np.meshgrid(ds_forçage.longitude.values,ds_forçage.latitude.values)
    xs_AROME,ys_AROME = np.ravel(xx_AROME),np.ravel(yy_AROME)
    
    # Changement systeme de coordonnées de epsg 4326 à Lambert 93 (epsg 2154)
    xs_A_2154, ys_A_2154 = convert_epsg_pts(xs_AROME,ys_AROME, epsg_src =4326, epsg_tgt=2154)

    # Ouverture données topo     
    ds_topo = crop_dataset(
        fic_to_ds = path_run_dir / topo_params_file,
        xmin = np.min(xs_A_2154),
        xmax = np.max(xs_A_2154),
        ymin = np.min(ys_A_2154),
        ymax = np.max(ys_A_2154),
        x_dim = "x",
        y_dim = "y",
        epsg_i = 2154,
        epsg_ds = 2154)
    
    xx,yy = np.meshgrid(ds_topo.x.values,ds_topo.y.values)
    nx,ny = xx.shape
    xs,ys = np.ravel(xx), np.ravel(yy)
        
    LW_fine = np.zeros((len(ds_forçage.valid_time),nx,ny))
    
    for i,date in enumerate(ds_forçage.valid_time[:1]) :
        
        # Select specifique date inside xr.Dataset forcings
        ds_forçage_t = ds_forçage.sel(valid_time=date)
        LW_fine[i,:,:] = LW(LW_AROME = ds_forçage_t.LWdown.values,
                            T2m_AROME = ds_forçage_t.t2m.values,
                            z_AROME = ds_forçage_t.z.values,
                            x_AROME = xs_A_2154,
                            y_AROME = ys_A_2154,
                            zs = ds_topo.ZS.values,
                            xs = xs,
                            ys = ys,
                            SVF = ds_topo.svf.values,
                            force_LWu = force_LWu)
    # Sauvegarde sous .nc
    ds_LW = xr.Dataset(
        data_vars=dict(
            LW=(["valid_time","y", "x"], LW_fine)
        ),
        coords=dict(
            time=("valid_time", ds_forçage.valid_time.values),
            x=("x", ds_topo.x.values),
            y=("y", ds_topo.y.values),
        )
    )
    ds_LW.LW.attrs={'units': '$W.m^{-2}$', 'standard_name': 'LWd', 'long_name': 'Downwelling Long Wave'}
    """
    ds_LW.to_netcdf(path_run_dir / save_name)
    """
    
    return ds_LW
    
########################### Arguments ######################################

if len(sys.argv) != 6:
    print("Usage: python3 baseline_lw.py run_dir forcing.nc topo_params.nc save_name.nc force_LWu")
    sys.exit(1)

baseline_LW(
    run_dir = sys.argv[1],
    forcing_file = sys.argv[2],
    topo_params_file = sys.argv[3],
    save_name = sys.argv[4],
    force_LWu = sys;argv[5])
    

