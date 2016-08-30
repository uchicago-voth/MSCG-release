# Instructions for reweight_fm example #
----------------------------------------

This example demonstrates the use of the statistical reweighting feature.
The input system is a system of 1000 methanol molecules that have each been mapped
to a single, center-of-mass site. 
As with all of the examples, the trajectory is not long enough to contain sufficient 
sampling to produce a converged result.

1) To run the force matching executable using reweighting, ensure that 
	"use_statistical_reweighting 1" is in control.in. The "frame_weights.in" input
	has already been prepared. Now, you can run the force matching executable
	with GROMACS trajectory input:
./newfm.x -f example_32.trr

2) Compare these results with the files in the "output" directory that begin "yes_".

3) To run the force matching executable without reweighting, ensure that 
	"use_statistical_reweighting 0" is in control.in (or that it is not present since
	0 is the default). The "frame_weights.in" input is not needed here. Now, you can run 
	the force matching executable with GROMACS trajectory input:
./newfm.x -f example_32.trr

4) Compare these results with the files in the "output" directory that begins "no_".
