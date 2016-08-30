# Instructions for ucg_lammps_fm_dynamic_type example #
-------------------------------------------------------

This example demonstrates the use of the "dynamic_types" feature of the code, which
is only supported with LAMMPS trajectory input.
The input system is a system of 1000 1-2-dichloroethane molecules that have each been 
mapped to a single, center-of-mass site. As described in the UCG II paper, these sites
have 2 states: an anti state and a gauche state (in reference to the chloride 
orientations). Thus, the state of each site can change throughout the course of the 
mapped input trajectory as its internal state changes.
The effective state of each site is read from the "type" column of the LAMMPS trajectory.
Note: This feature is not designed to handle topology changes or unsorted inputs.
As with all of the examples, the trajectory is not long enough to contain sufficient 
sampling to produce a converged result.

1) The range files (rmin.in and the empty rmin_b.in) have already been prepared.
	(Optional) These ranges can be verified by running the rangefinder executable:
./rangefinder.x -l C2Cl2H4_Lammps.dat
	
2) Run the force matching executable with LAMMPS trajectory input:
./newfm.x -l C2Cl2H4_Lammps.dat

2) Compare the results with those in the "output" directory.

