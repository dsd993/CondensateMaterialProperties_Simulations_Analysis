These set of files can be used to run the slab simulations for computing concentrations and to cut of the dense phase section after equilibration for the production runs.  
The example set of files given here are for the diblock E-K sequence (80 chains), simulated using the Martini 3 model.  See 
EKV_15.dat file for the sequence of diblock E-K.

File details:

Initial coordinate file --> PRO_SOL_IONS.gro

Final coordinate file --> final_denseEKV15.gro

Topology file --> PRO_SOL_IONS.top

GROMACS input script file --> prod_Martini_NPT.mdp

Command used to compute the density profile from slab simulations: "gmx density" within the GROMACS software package.