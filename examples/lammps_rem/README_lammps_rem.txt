# Instructions for lammps_rem example #
--------------------------------------

This example demonstrates the use of the relative entropy minimization code with LAMMPS.
This example also demonstrates the use of the boltzman inversion code within rangefinder.
The input system is a system of 1000 chloromethane molecules that have been mapped
to a single, center-of-mass site.
As with all of the examples, the trajectory is not long enough to contain sufficient 
sampling to produce a converged result.

1) Prepare the range files and zeroth iteration for the REM. The zeroth iteraction is
prepared by boltmann inversion using the rangefinder code:
./rangefinder.x -l reference_trajectory.lammpstrj

1) Run a CG simulation using the resulting BI interactions given by rangefinder. (see lammps_bi example)
If you are using lammps, simulation can be run using the provided input and data files:
./lmp_mpi < CG.in
(note: the lammps input file "CG.in" will expect the name of the interactions to be in a file
named "fmec.table". This file can be created with the command:
cp 1_1.table fmec.table
)

2) Move all the files including the newrem.x excutable and BI simulation trajectory, 
but not the reference trajectory and the iterate.sh script, to a directory named "iter_0".

3) Next, the provided iterate script can be used to perform the REM routine:
source iterate.sh 10 reference_traj.lammpstrj dump.lammpstrj.cg lmp_mpi

4) compare the results to those in the "output" directory.
