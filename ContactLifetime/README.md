### This python code calculates the following:

1) Intermittent contact lifetimes between the oppositely charged residues, like charged residues, and all charged residues within an E-K sequence of a given nSCD in a file named "{SeqName}_ContactLifetimeACF.txt".

Note that the "{SeqName}_ContactLifetimeACF.txt" will have 5 columns: The first column is the time in ps, the second column is time in ns, the third column is the raw ACF between oppositely charged residues, the fourth column is the normalized ACF (goes from 0 to 1) between oppositely charged residues, the fifth column is the raw ACF between like charged residues, the sixth column is the normalized ACF (goes from 0 to 1) between like charged residues, the seventh column is the raw ACF between all charged residues, and the eigth column is the normalized ACF (goes from 0 to 1) between all charged residues.

The code was tested on version 3.9.12.

### Requirements

This code comes with 4 files: (1) two .py files, (2) one .dat sequence file, (3) one .dat amino acid diameter file ("aminoacids_vdwdiameter.dat")

### Code usage

Use the code on the GSD file (only a small file consisting of 15 trajectory stamps at specific timesteps are provided) given inside the folder GSDFile for testing.  These trajectories were generated using HOOMD-blue (version 2.9.3).

### Compile and run the python code as below

* python ninter_contactlifetimes.py GSDFile/ekv3_dumps.gsd GSDFile/ekv3_dumps.gsd 0 1 15 500 500 500 EKV3.dat EKV3

In the above line, "ninter_contactlifetimes.py" corresponds to the .py filename; "GSDFile/ekv3_dumps.gsd" corresponds to the GSD (trajectory) file location; "GSDFile/ekv3_dumps.gsd" corresponds to the GSD (topology) file location; "0" corresponds to the starting frame index; "1" corresponds to the skipping (striding) between frames; "15" corresponds to the ending frame index; "500" corresponds to the interval timesteps between trajectories used in simulations; "500" corresponds to the total number of chains in the system; "500" corresponds to the number of chains of interest; "EKV3.dat" is the sequence filename (here EKV3 with nSCD = 0.033 if given as an example); "EKV3" is the sequence name to be used in the output filename.  

Number of chains of interest means the number of chains that will be used as reference to compute rather than all chains in the system ~ typically it is the same as the total number of chains in the system.  However, to quickly run the code and check the results, use a smaller number.  If you use e.g., "3", the first 3 chains (according to GSD file indexing) will be picked and its ACF with all other chains in the system will be computed. 

Executing the above line will print a NUMPY file "{SeqName}_CLifetimeACFArray.npy" consisting of an array with the non-normalized (raw) autocorrelation function values after averaging over no. of pairs and no. of chains for each residue.
 
* python IntermittentACF.py 0 1 15 500 EKV3.dat EKV3

In the above line, "IntermittentACF.py" corresponds to the .py filename; "0" corresponds to the starting frame index; "1" corresponds to the skipping (striding) between frames; "15" corresponds to the ending frame index; "500" corresponds to the interval timesteps between trajectories used in simulations; "EKV3.dat" is the sequence filename (here EKV3 with nSCD = 0.033 if given as an example); "EKV3" is the sequence name to be used in the output filename. 

Executing the above line will read the NUMPY file "{SeqName}_CLifetimeACFArray.npy" and give out the final output in a file named "{SeqName}_ContactLifetimeACF.txt".  See above for the information contained in this file "{SeqName}_ContactLifetimeACF.txt".
