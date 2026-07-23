# Downscalling downwelling irradiance in complexe terrain

The following code aim to **downscall solar and thermal irrandiance in complexe terrain** (i.e. mountains). It is seperated in three 
parts, as represented in the following flowchart : (i) **Computing permanent topographic parameters** independant of the meteorological variables (_compute_permanent_features.py_),
computing the downscalled (ii) **short wave** (_compute_SW.py_) and (iii) **long wave** (_compute_SW.py_) irradiance based on constant 
topographic parameters and time dependant meteorological forcing. 

<img width="1823" height="1229" alt="full_baseline_flowchart_formal" src="https://github.com/user-attachments/assets/4a5aedda-114e-4b19-a7cc-3c2a9f4df4d7" />

## _Compute_permanent_features.py_

This part of code is to be launched only one time per topographic grid. It uses the **Topocalc library** (https://github.com/USDA-ARS-NWRC/topocalc/blob/main/README.md). 
Please note that the Topocalc library shall only function with numpy version inferior to 2.0. It provides the other functions (_Compute_SW.py_ & _Compute_LW.py_) two netcdf files gridded on the input digital elevation model :

(i) **topo_params_{nom_experience}.nc** wich is the input DEM with added slope, aspect and sky view factor variables. The sky view factor may be seen as the transfert function of an atmospheric isotropic diffuse radiation.

(ii) **shadow_mask_{nom_experience}.nc**, containing a unique variable shadow_mask, the transfert function of the short wave direct irradiance. 

## _Compute_SW.py_

Takes in input a meteorological forcing - with variables SWdir, SWdif and coordinates lat, lon - to be downscalled to the topographic grid represented by the output of _Compute_permanent_features.py_, as well as
the two outputs of _Compute_permanent_features.py_. Returns **SW.nc** of downscalled direct (resp. diffuse) solar irradiance SWdir (resp. SWdif).

### Method



## _Compute_LW.py_

Takes in input a meteorological forcing - with variables LWd, T2m and ZS and coordinates lat, lon - to be downscalled to the topographic grid represented by the output of _Compute_permanent_features.py_, as well as 
**topo_params_{nom_experience}.nc** from _Compute_permanent_features.py_. Returns **LW.nc** of downscalled thermal irradiance from the atmosphere and the surrounding slopes (in option).

### Method

La contribution atmosphérique dans le flux descendant long wave sur un élément de surface est décrit par celui émis par une surface corps noir de température T2m (downscalled) d'émissivité celle d'AROME :
$$\epsilon_{atm} = \frac{LW_{AROME}}{\sigma T_{2m,AROME}^4}$$ \\
La composante atmosphérique du rayonnement LW sur une surface plane est alors :
$$LW_{atm} = \epsilon_{atm}\sigma T_{2m}^4$$
La composante surfacique à partir de la température de surface déscendue d'échèle ou interpolée à partir de l'émissivité de la surface $\epsilon_s$ ($=0.98$ pour la neige) :
$$LW_{s} = avg_\mathcal{N}[\underbrace{\epsilon_s\sigma T_s^4}_\text{émis }+ \underbrace{(1-\epsilon_s)\epsilon_{atm}\sigma T_{2m}^4}_\text{atm reflechis}]$$
D'où un bilan complet modulé par le SVF :
$$LW_d = \underbrace{\epsilon_{atm}\sigma T_{2m}^4\times SVF}_\text{atm}  + \underbrace{avg_\mathcal{N}[\epsilon_s\sigma T_s^4+(1-\epsilon_s)\epsilon_{atm}\sigma T_{2m}^4] \times (1-SVF)}_\text{surfaces voisines}$$
\begin{figure}[!h]

