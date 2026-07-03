# -*- coding: utf-8 -*-

import xarray as xr
import numpy as np
import pandas as pd
from pyproj import Transformer
import numba as nb

from pathlib import Path

############################ Dictionnaire variables et constantes ######################

sigma = 5.67e-8

############################ Functions ####################################

######## Select area ##########

def crop_dataset(
    fic_to_ds : str,
    epsg_i : int,
    epsg_ds : int,
    x_dim,
    y_dim,
    xmin = None,
    xmax = None,
    ymin = None,
    ymax = None,
    )-> xr.Dataset :
    
    ds = xr.open_dataset(fic_to_ds)
    
    # Reindexing on increasing coordinates 
    ds = ds.sortby(x_dim)
    ds = ds.sortby(y_dim)

    if epsg_i != epsg_ds :
        xmin,ymin = convert_epsg_pts(xmin,ymin, epsg_src = epsg_i, epsg_tgt=epsg_ds)
        xmax,ymax = convert_epsg_pts(xmax,ymax, epsg_src = epsg_i, epsg_tgt=epsg_ds)
        
    ds = ds.sel({x_dim: slice(xmin, xmax), y_dim: slice(ymin, ymax)})
    
    return ds

def slice_by_time(ds: xr.Dataset,
                  start: str | pd.Timestamp,
                  stop:  str | pd.Timestamp,
                  dim: str = "valid_time") -> xr.Dataset:
    
    if dim not in ds.coords:
        raise KeyError(f"'{dim}' n’est pas une coordonnée du dataset "
                       f"(coordonnées disponibles : {list(ds.coords)})")

    start_ts = pd.to_datetime(start)
    stop_ts  = pd.to_datetime(stop)

    ds_slice = ds.sel({dim: slice(start_ts, stop_ts)})

    return ds_slice

######## Interpolation #########

@nb.njit()
def find_closest(arr: np.ndarray, 
                 val: float) -> int:
    """
    [Input]
    - array 1D rangé dans l'ordre croissant
    - valeur
    [Output]
    - index tel que arr[i] <= val < arr[i+1]
    """
    a, b = 0, arr.shape[0] - 2          # on veut l'index gauche
    while a < b:
        milieu = (a + b + 1) // 2
        if arr[milieu] <= val:
            a = milieu
        else:
            b = milieu - 1
    return a

@nb.njit()
def bilinear_interpolation(
    xs: np.ndarray, 
    ys: np.ndarray,
    arr: np.ndarray,
    points: np.ndarray,
    ) -> np.ndarray:
    """
    [Input]
    - coordonnées x de la grille (Nx,) triées croissantes
    - coordonnées y de la grille (Ny,) triées croissantes
    - array (Ny, Nx) des valeurs à interpoler pour chaque points de la grille 
    - points cibles (N, 3)  colonnes x, y, z
    
    [Output]
    - (N,) de valeurs SVF interpolées (NaN si hors grille)
    """
    N = points.shape[0]
    result = np.empty(N, dtype=np.float64)
 
    x_min, x_max = xs[0], xs[-1]
    y_min, y_max = ys[0], ys[-1]
 
    for k in range(N):
        px = points[k, 0]
        py = points[k, 1]
 
        # Hors domaine renvoit NaN
        if px < x_min or px > x_max or py < y_min or py > y_max:
            result[k] = np.nan
            continue
 
        ix = find_closest(xs, px)
        iy = find_closest(ys, py)
        
        # Renvoit un Nan si les proches voisins sont des Nan
        if (np.isnan(arr[iy,ix]) == True or arr[iy+1,ix+1] == True or
            arr[iy+1,ix] == True or arr[iy,ix+1] == True):
            
            result[k] = np.nan
            continue
            
        x0, x1 = xs[ix], xs[ix + 1]
        y0, y1 = ys[iy], ys[iy + 1]
 
        # Poids de pondération = distance relative à l'écart entre les proches voisins
        tx = (px - x0) / (x1 - x0)
        ty = (py - y0) / (y1 - y0)
 
        # Coins de la cellule  arr[iy, ix] → ligne = y, colonne = x
        res00 = arr[iy,     ix    ]
        res10 = arr[iy,     ix + 1]
        res01 = arr[iy + 1, ix    ]
        res11 = arr[iy + 1, ix + 1]
 
        result[k] = (
            res00 * (1 - tx) * (1 - ty)
            + res10 *      tx  * (1 - ty)
            + res01 * (1 - tx) *      ty
            + res11 *      tx  *      ty
        )
 
    return result

@nb.njit()
def interpolate(
    xs : np.ndarray,
    ys : np.ndarray,
    arr : np.ndarray,
    points: np.ndarray,
    ) -> np.ndarray:
    """
    [Input]
    - xs : valeurs 1D des abscisses du tableau
    - ys : valeurs 1D des ordonnées du tableau
    - arr : tableau 2D des valeurs à interpoler
    - points : ndarray des points où interpoler (Nbre_points, 2) (x,y)
    [Output]
    - np.ndarray (N,) de valeurs interpolées (Nan hors domaine)
    """
    
    # Garantit que les axes sont croissants pour si jamais sortie de xarray
    if xs[-1] < xs[0]:
        xs  = xs[::-1].copy()
        arr = arr[:, ::-1].copy()
    if ys[-1] < ys[0]:
        ys  = ys[::-1].copy()
        arr = arr[::-1, :].copy()
 
    pts = np.ascontiguousarray(points)#, dtype=np.float64)
 
    return bilinear_interpolation(xs, ys, arr, pts)

########### Altitude gradiant for temperature ###########

def altitudinal_gradiant(
    diff_ZS : float | np.ndarray, #ZS fine - ZS coarse
    Ts : float | np.ndarray, # Ts coarse interpolé sur la grille fine
    )-> float | np.ndarray: # Ts fine corrigé par un laspe rate moyen tropospherique
    
    laspe_rate = -6.5e-3 #K/m
    
    return Ts + diff_ZS*laspe_rate
    
########## Convert coordinates ########################

def convert_epsg_pts(xs,ys, epsg_src, epsg_tgt=32632):
    """
    Simple function to convert a list fo poitn from one projection to another oen using PyProj

    Args:
        xs (array): 1D array with X-coordinate expressed in the source EPSG
        ys (array): 1D array with Y-coordinate expressed in the source EPSG
        epsg_src (int): source projection EPSG code
        epsg_tgt (int): target projection EPSG code

    Returns: 
        array: Xs 1D arrays of the point coordinates expressed in the target projection
        array: Ys 1D arrays of the point coordinates expressed in the target projection
    """
    #print('Convert coordinates from EPSG:{} to EPSG:{}'.format(epsg_src, epsg_tgt))
    trans = Transformer.from_crs("epsg:{}".format(epsg_src), "epsg:{}".format(epsg_tgt), always_xy=True)
    Xs, Ys = trans.transform(xs, ys)
    return Xs, Ys

