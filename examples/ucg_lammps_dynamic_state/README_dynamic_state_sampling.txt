# Instructions for ucg_lammps_dynamic_state example #
-------------------------------------------------------

This example demonstrates the use of the "dynamic_state_sampling" feature of the code, 
which is only supported with LAMMPS trajectory input.
The input system is a system of 1000 methanol molecules that have each been 
mapped to a single, center-of-mass site. As described in the UCG III paper, these sites
have 2 well mixed internal states: a low-density state and high-density state. Thus, the 
each site has a probability of being in a given state (for the purposes of fitting).
The data efficiency of the input trajectory is enhanced by resampling the types of the
site a number of times (specified by "dynamic_state_samples_per_frame" in "control.in").
The probability of the lower state of each site is read from the "state" column of the 
LAMMPS trajectory. This should be a number between 0 and 1. Also, it is designed to 
sample between "types" "1" and "2". So, the top.in file and trajectory types must be
set so that there is no conflict between the switching types and other CG sites.
Note: This feature is not designed to handle topology changes or unsorted inputs.
As with all of the examples, the trajectory is not long enough to contain sufficient 
sampling to produce a converged result.

1) The range files (rmin.in and the empty rmin_b.in) have already been prepared.
   Also, "dynamic_state_sampling" is set to 1 and "dynamic_state_samples_per_frame" is
   specified in "control.in".
   Note: The top.in file lists both types 1 and 2 to allow interactions between the
   two states. Also, the range files should list the full interaction range of the 
   particles (as if they were all the same type) for all possible type combinations.
	
2) Run the force matching executable with LAMMPS trajectory input:
./newfm.x -l MeOH_state.dat

2) Compare the results with those in the "output" directory.

