::: center
\
**Downscalling downwelling irradiance in complexe terrain**\
:::

The following code aim to **downscall solar and thermal irrandiance in
complexe terrain** (i.e. mountains). It is seperated in three parts, as
represented in figure [1](#fig:flowchart){reference-type="ref"
reference="fig:flowchart"} : (i) **Computing permanent topographic
parameters** independant of the meteorological variables
(*compute_permanent_features.py*), computing the downscalled (ii)
**short wave** (*compute_SW.py*) and (iii) **long wave**
(*compute_SW.py*) irradiance based on constant topographic parameters
and time dependant meteorological forcing.

<figure id="fig:flowchart" data-latex-placement="H">
<img src="./full_baseline_flowchart_formal.png" />
<figcaption>Flowchart of the downscalling baseline.</figcaption>
</figure>

# *Compute_permanent_features.py* {#compute_permanent_features.py .unnumbered}

This part of code is to be launched only one time per topographic grid.
It uses the [**Topocalc
library**](https://github.com/USDA-ARS-NWRC/topocalc/blob/main/README.md).
Please note that the Topocalc library shall only function with numpy
version inferior to 2.0. It provides the other functions
(*Compute_SW.py* & *Compute_LW.py*) two netcdf files gridded on the
input digital elevation model :

- **topo_params\_{nom_experience}.nc** wich is the input DEM with added
  slope, aspect and sky view factor variables. The sky view factor may
  be seen as the transfert function of an atmospheric isotropic diffuse
  radiation.

- **shadow_mask\_{nom_experience}.nc**, containing a unique variable
  shadow_mask, the transfert function of the short wave direct
  irradiance.

# *Compute_SW.py* {#compute_sw.py .unnumbered}

Takes in input a meteorological forcing - with variables SWdir, SWdif
and coordinates lat, lon - to be downscalled to the topographic grid
represented by the output of *Compute_permanent_features.py*, as well as
the two outputs of *Compute_permanent_features.py*. Returns **SW.nc** of
downscalled direct (resp. diffuse) solar irradiance SWdir (resp. SWdif).

## Method {#method .unnumbered}

Le rayonnement solaire direct est multiplié par le mask d'ombrage,
prenant en compte la projection sur la pente, l'ombrage propre et
projeté par le voisinnage :
$$SW_{dir,SUB} = Shadow\_mask \times SW_{dir,GRID}$$ Le rayonnement
solaire diffus, considéré comme isotrope, est multiplié par le $SVF$ :
$$SW_{dif,SUB} = SVF \times SW_{dif,GRID}$$

# *Compute_LW.py* {#compute_lw.py .unnumbered}

Takes in input a meteorological forcing - with variables LWd, T2m and ZS
and coordinates lat, lon - to be downscalled to the topographic grid
represented by the output of *compute_permanent_features.py*, as well as
topo_params\_nom_experience.nc from Compute_permanent_features.py.
Returns \*\*LW.nc\*\* of downscalled thermal irradiance from the
atmosphere and the surrounding slopes (in option).

## Method {#method-1 .unnumbered}

La contribution atmosphérique dans le flux descendant long wave sur un
élément de surface est décrit par celui émis par une surface corps noir
de température T2m (downscalled) d'émissivité celle d'AROME. La calcul
est fait sur la grille AROME :
$$\epsilon_{atm} = \frac{LW_{AROME}}{\sigma T_{2m,AROME}^4}$$\
La température de l'atmosphère à 2 mètre est déscendue d'échelle (en
attendant une meilleur méhode) par interplation linéaire verticale
suivant un gradiant standard $\Gamma$ :
$$T_\text{2m,SUB} = T_{2m,GRID} + \Gamma \times \Delta z$$ La composante
atmosphérique du rayonnement LW sur une surface plane est alors :
$$LW_{atm,SUB} = \epsilon_{atm}\sigma T_{2m,SUB}^4$$ La composante
atmophérique du rayonnement LW sur une surface complexe est alors,
considérant le rayonnement LW comme isotrope :
$$LW_{atm,SUB} = SVF \times \epsilon_{atm}\sigma T_{2m,SUB}^4$$ La
composante surfacique provenant de l'irradiance des pentes adjacentes
sur l'élément de surface considéré est modélisé comme le rayonnement du
surface corps noirs à la température $T_s$ de l'élément de surface
considéré, d'émissivité $\epsilon_s$, dont les pentes voisines occupes
(1-SVF) d'angle solide :
$$LW_{s,SUB} = (1-SVF)\times \epsilon_s\sigma T_s^4$$ L'approximation de
considérer l'irradiance des surfaces voisines équivalente à celle
produite par un ensemble de surfaces radiativement équivalent à
l'élément de surface considéré est vérifié sur la figure
[2](#fig:voisins){reference-type="ref" reference="fig:voisins"}.

<figure id="fig:voisins" data-latex-placement="H">
<img src="./Compare_method_neighbours_EDELWEISS.jpg" />
<figcaption>Flowchart of the downscalling baseline.</figcaption>
</figure>

L'irradiance complète pour un élément de surface donné revient donc à
une somme de flux atmosphérique et surfacique pondérés par le SVF
$$LW_d = \underbrace{SVF\times \epsilon_{atm}\sigma T_{2m}^4}_\text{atmosphère}  + \underbrace{(1-SVF)\times \epsilon_s\sigma T_s^4 }_\text{surfaces voisines}$$
