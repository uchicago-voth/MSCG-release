# Instructions for table_fm example #
-------------------------------------

This example demonstrates the use of tabulated interactions.
The input system is a system of 250 methanol molecules that have each been mapped
to two center-of-mass sites (one for the methyl and the other for the hydroxyl). 
In this case, the bond between these two sites is tabulated (perhaps from Boltzmann
Inversion) and should be subtracted before force matching the non-bonded interactions.
As with all of the examples, the trajectory is not long enough to contain sufficient 
sampling to produce a converged result.

1) Note that the range files (rmin.in and rmin_b.in) have already been prepared.
	Since the bond between CG types 1 and 2 is tabulated, the mode for this interaction
	is "tab". The "table.in" file has already been prepared as well.
	
2) Run the force matching executable using with LAMMPS trajectory input:
./newfm.x -l cg2site.lammpstrjnew.short

3) Compare these results with the files in the "output" directory.
