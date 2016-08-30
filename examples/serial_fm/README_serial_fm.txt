# Instructions for serial_fm example #
--------------------------------------

This example demonstrates the use of the code with GROMACS trajectory input.
The input system is a system of 1000 methanol molecules that have each been mapped
to a single, center-of-mass site. 
As with all of the examples, the trajectory is not long enough to contain sufficient 
sampling to produce a converged result.

1) Run the rangefinder executable using GROMACS trajectory input to generate the range
	files (rmin.in and rmin_b.in):
./rangefinder.x -f cg.trr

2) Compare these range files with those in the "output" directory.
	
3) Run the force matching executable with GROMACS trajectory input:
./newfm.x -f cg.trr

4) Compare the results with those in the "output" directory.