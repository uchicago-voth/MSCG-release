# Instructions for lammps_fm example #
--------------------------------------

This example demonstrates the use of the code with LAMMPS trajectory input.
The input system is a system of 1000 methanol molecules that have each been mapped
to a single, center-of-mass site. 
As with all of the examples, the trajectory is not long enough to contain sufficient 
sampling to produce a converged result.

1) The range files (rmin.in and the empty rmin_b.in) have already been prepared.
	(Optional) These ranges can be verified by running the rangefinder executable:
./rangefinder.x -l MeOH_example.dat
	
2) Run the force matching executable with LAMMPS trajectory input:
./newfm.x -l MeOH_example.dat

2) Compare the results with those in the "output" directory.