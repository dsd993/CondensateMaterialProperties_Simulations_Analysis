These set of files can be used to run the dense phase simulations for computing viscosity based on the Green-Kubo relation.  The example set of 
files given here are for the diblock E-K sequence (80 chains), simulated using the Martini 3 model.  See 
EKV_15.dat file for the sequence of diblock E-K.

Note than an extension of this run to longer times, as discussed in the manuscript, can be done for the purpose of 
computing diffusion coefficients at long times.

File details:

Initial coordinate file --> PRO_SOL_IONS.gro

Final coordinate file --> prod.gro

Topology file --> PRO_SOL_IONS.top

GROMACS input script file --> prod_Martini_NVT.mdp


Command used to compute mean-squared displacement within GROMACS from which diffusion coefficient was extracted as discussed in the Methods
section of the manuscript: "gmx msd" within the GROMACS software package.

For computing viscosity based on the Green-Kubo relation, we used the codes shared on a Github repository as part
of a recent Journal of Chemical Information and Modeling (JCIM) article by Prass et al., 2023 (link: https://doi.org/10.1021/acs.jcim.3c00947).