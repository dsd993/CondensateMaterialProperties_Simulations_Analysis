These set of files can be used to run the initial equilibrium P = 0 atm simulations to get to the protein's natural concentration 
and then, switch to Langevin dynamics (LD) simulations to simulate the dense phase of the protein chains.  From these simulations,
one can compute intermittent contact lifetime and diffusion coefficient, the codes of which are given in the GitHub respository.
The example set of files given here are for the LAF1-WT protein (150 chains), simulated using the HPS-Urry model.  See 
LAF1RGG_Variant2.dat file for the sequence of LAF1-WT.

File details:

Initial coordinate file --> 1000.0x1000.0x1000.0_box_min.gsd

Final coordinate file --> LDrestart_LAF2.gsd

HOOMD input script file --> runfile.py