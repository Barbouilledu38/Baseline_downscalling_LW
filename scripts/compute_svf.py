# -*- coding: utf-8 -*-

from utils import *

from topocalc import gradient
from topocalc import viewf
from topocalc import horizon

############################################# Core function ####################################

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
            ZS=(["y", "x"], dem_arr.astype(np.float64)),
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
    
############################### Arguments ########################################################

if len(sys.argv) != 3:
    print("Usage: python3 compute_svf.py run_dir fic_topo.nc")
    sys.exit(1)
    
Baseline_SW(
    run_dir = sys.argv[1],
    fic_topo = sys.argv[2])


