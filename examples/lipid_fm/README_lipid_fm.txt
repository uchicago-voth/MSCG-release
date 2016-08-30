# Instructions for lipid_fm example #
--------------------------------------

This example demonstrates the use of the code a biomolecule. Specifically, it is a system
of 72 lipids, where each has been mapped to 15 sites (with 8 unique types). The set-up
shows the top.in for a system with non-trivial geometry. Systems such as this demonstrate
the use of the code (and highlight the performance of sparse matrix operations) for
the coarse-graining of complex biomolecules.
As with all of the examples, the trajectory is not long enough to contain sufficient 
sampling to produce a converged result.

1) Run the range finding executable with GROMACS trajectory input:
./rangefinder.x -f cg_dopc72.trr

2) Compare the interaction range files (rmin.in and rmin_b.in) with those in the 
	"output" directory.
	
3) Run the force matching executable with a LAMMPS trajectory using either 
	dense matrix operations or sparse matrix operations.
	
3.A) For dense matrix operations, make sure that "matrix_type 0" is specified in 
	"control.in". Then, run the force matching code with GROMACS trajectory input:
./newfm.x -f cg_dopc72.trr
	Note: This may take several hours.
	Then, you can compare the results with those in the "dense_output" directory.
	
3.B) For sparse matrix operations, you must compile the Intel MKL enabled version of the 
	code (newfm_mkl.x). Make sure to select either "matrix_type 3" or "matrix_type 4" in 
	"control.in". Then, run the force matching code with GROMACS trajectory input:
./newfm_mkl.x -f cg_dopc72.trr
	Note: "matrix_type 4" runs much faster than "matrix_type 3", which takes a long time
		to solve the normal matrix after all frames have been accumulated.
	The comparisons for sparse matrix operations were generated using "matrix_type 4".
	Note: Different versions of Intel MKL may produce slightly different output 
		interaction potentials because the PARDISO solver included in Intel MKL 
		may change from version to version.