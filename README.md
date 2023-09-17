# 3D_DensityMap
Python script to compute 3D density map for an ensemble of conformational models


## General description

Using a trajectory that comprises of ensemble members, and their corresponding population fractions (weights),
this script calculates 3D density map, by mapping the coordinates of ensemble members on a grid of uniform dimensions
(_nx_, _ny_, _nz_) = (100,100,100) and whose voxel size is _dg_ = 2 Å. Voxel occupancy is a sum of weights of those ensemble members whose atomic coordinates were found
in a respective voxel.

In this example, 3D density map is calculated using all atoms of conformational models. Alternatively it is possible to calculate 3D density map of backbone atoms only, or of particular residues of interest.
Calculation of 3D density map requires superimposition of conformers in the ensemble. In the current version of the script, superimposition is done on
the first frame, which in the case of a given ensemble represents the most populated ensemble member. 
It is optional to what reference structure the ensemble will be aligned. For example, one can align to some independent reference structure that is not part of the ensemble, e.g. crystal structure model.


3D density map for an ensemble is saved in the MRC2014<sup>2</sup> file format using _mrcfile_ python library<sup>3</sup> as "3D_density_map_g.mrc" file. 3D density map can be opened in 3D rendering software UCSF Chimera<sup>4
</sup>. 


![F_posterior](https://github.com/mpopara/3D_DensityMap/assets/40856779/e605d4d0-8a01-4893-897c-cd967f2c7231)

Illustration outlining the extent of the ensemble at 50%, 68% and 90% of the density-weighted volume is adapted from Dittrich _et al_, 2023.<sup>1</sup>

Besides calculating total 3D density map of an ensemble, this script also calculated so-called half-density maps ("3D_density_map_h1.mrc" and "3D_density_map_h2.mrc").
Half maps are computed using odd and even frames of the ensemble. Half maps can be subsequentlyy used to calulate so called Fourier Shell Correlation (FSC) curves<sup>5</sup> , which give
precision estimate for the ensemble.<sup>1</sup> 

Additionally, this script calculates density maps for different quartiles of cumulative weights, i.e. Q<sub>1</sub>(0.00-0.25), Q<sub>2</sub>(0.25-0.5), Q<sub>3</sub>(0.5-0.75) and Q<sub>4</sub>(0.75-1.00). Moreover, calculated are density maps for cumulative
sum of weights of 0.25, 0.5, 0.75 and 1. Such analysis informs whether the 3D density map changes between highly and less populated ensemble members. 
This is particularly convenient when comparing 3D density maps between two ensembles- for example, similarity between ensembles may be high when accounting the most populated ensemble members,
but low upon addition of less populated structures.


## Input file requirements

* ensemble of conformational models provided as trajectory in any of the mdtraj compatible formats (dcd, nc, xtc..)
* topology as .pdb file
* .dat file containing weights (population fractions) of ensemble members. This space-delimited file is of a size N<sub>conformers</sub> x 2, where the first column contains indices of the ensemble members,
 and the second column contains their corresponding weights. This script assumes that the order of ensemble members in the trajectory file follows the same order as in the weights file.

## Dependencies
_compute_density.py_ is a python script built on Python 3.8.8. Script was tested using provided exemplary data files under the following configuration:

* Windows 10
* Python 3.8.8
* mdtraj 1.9.4
* numpy 1.23.0
* numba 0.53.1
* mrcfile 1.3.0



## References
1. Dittrich, J.; Popara, M.; Kubiak, J.; Dimura, M.; Schepers, B.; Verma, N.; Schmitz, B.; Dollinger, . P.; Kovacic, F.; Jaeger, K. E.;
Seidel, C. A. M.; Peulen, T. O.; Gohlke, H., Resolution of Maximum Entropy Method-Derived Posterior Conformational Ensembles of a Flexible System Probed by FRET and Molecular Dynamics Simulations.
J Chem Theory Comput 2023, 19 (8), 2389-2409.

2. Cheng, A.; Henderson, R.; Mastronarde, D.; Ludtke, S. J.; Schoenmakers, R. H. M.; Short, J.; Marabini, R.; Dallakyan, S.; Agard, D.; Winn, M.  MRC2014: Extensions to the MRC format header for
electron cryo-microscopy and tomography. J. Struct. Biol. 2015, 192, 146−150.

3. Burnley, T.; Palmer, C. M.; Winn, M. Recent developments in the CCP-EM software suite. Acta Crystallogr., Sect. D: Struct. Biol. 2017, 73, 469−477.

4. Pettersen, E. F.; Goddard, T. D.; Huang, C. C.; Couch, G. S.; Greenblatt, D. M.; Meng, E. C.; Ferrin, T. E. UCSF Chimera-A visualization system for exploratory research and analysis. J. Comput. Chem. 2004, 25, 1605−1612.

5. Harauz, G. and van Heel, M. (1986) Exact Filters for General Geometry Three Dimensional Reconstruction. Optik, 73, 146-156.

## Authors

* Milana Popara
* Thomas-Otavio Peulen
