# Instructions for lammps_rem example #
--------------------------------------

This example demonstrates the use of the boltzman inversion code within rangefinder.
The input system is a system of 1000 chloromethane molecules that have been mapped
to a single, center-of-mass site.
As with all of the examples, the trajectory is not long enough to contain sufficient 
sampling to produce a converged result.

1) Prepare the range files and zeroth iteration for the REM. The zeroth iteraction is
prepared by boltmann inversion using the rangefinder code:
./rangefinder.x -l reference_trajectory.lammpstrj

2) compare the results to those in the "output" directory.
