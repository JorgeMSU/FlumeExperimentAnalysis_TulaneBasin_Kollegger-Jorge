# FlumeExperimentAnalysis_TulaneBasin_Kollegger-Jorge

Below we include a description of the codes used to generate Figures 1 to 4 of Kollegger et al. 

The first step was to extract the elevation information of the fluvial surface over time for our analysis. The experiment was well documented with photographs taken every 15 minutes, and a laser scanner gathered topographical information every hour. This topographical information was then converted into a Digital Elevation Model (DEM) that we used to observe how the fluvial surface changed with time (Yu et al. 2017). 

We used the script called “TopoMatrix.m” to define two 3D matrices with the fluvial surface elevations over time for each sea-level scenario. The "LMLP.mat" matrix for the Low Magnitude, Long Period (LMLP) scenario, and the "Anew.mat" matrix for the High Magnitude, Short Period (HMSP) scenario. Note that this script is not required to reproduce the plots included in the manuscript; “Anew.mat” and "LMLP.mat" are also included in the Zenodo repository. The HMSP scenario is associated with hours 50-540 of the experiment (and thus in the SEAD repository), a sea-level amplitude of 12.25mm, and a period of 24.5 hours. The LMLP scenario is associated with hours 50-540 in the SEAD repository, a sea-level amplitude of 3.06mm, and a period of 98 hours. Both sea-level scenarios had a background sea-level rise rate of 0.25mm/hr superimposed on their respective sea-level oscillations. Sediment (quartz dominated) and water input were held constant at the respective rates of 3.9 x 10-4kg/s, and 1.7 x 10-4m3/s, and entered the basin via a weir. The link to the SEAD repository can be found in the acknowledgments section of the manuscript. For more details on the experiment, see Yu et al. 2017.

After loading the matrices "Anew.map" or "LMLP.mat", we use the following scripts:
-	“RadialAverageMatrixResiduals.m” to generate figures 2D and 2E (HMSP scenario).
-	“AverageProfile.m” to generate figures 3C and 3D (HMSP scenario).
-	“LMLPRadialAverageMatrixResiduals.m” to generate figure 4A (LMLP scenario).
-	“LMLPAverageProfile.m” to generate figures 4B, 4C, and 4D (LMLP scenario).

With these scripts, we calculate the radial averages between a location 100mm downstream of the inlet (to avoid boundary effects) and the shoreline with a spacing of 5mm. We locate the shoreline using the curve intersect function (Hölz 2021). In some occasions, the shoreline extends beyond the domain covered by the topographic scans (1.3m from the inlet); in this case we fix the shoreline location until it retreats again. 
