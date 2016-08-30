# Instructions for ucg_lammps_fm_dynamic_type example #
-------------------------------------------------------

This example demonstrates the usage of the MS-CG static library for range finding
and force matching. Both example files make use a dummy LAMMPS reader through the
use of the "read_lammps.cpp" and "read_lammps.h" files.
This example trajectory is the same as the "lammps_fm" example and produces identical
results as that example. The input system is a system of 1000 methanol molecules that 
have each been mapped to a single, center-of-mass site. 
As with all of the examples, the trajectory is not long enough to contain sufficient 
sampling to produce a converged result.

1) The mscg.a library should be compiled following the instructions in the 
   "README_library.txt" file.
   
2) The example usage of the range finding features are demonstrated in the 
   "rangefinder_driver.cpp" file. Compile this file in the same directory as the
   "read_lammps" files and link to the MSCG library.
   e.g. icc rangefinder_driver.cpp read_lammps.cpp read_lammps.h -lmscg -o range_driver.x

3) The ranges sampled in this trajectory can be determined by running the rangefinder 
   executable:
   e.g. ./range_driver.x -l MeOH_example.dat
	
4) Compare the range files (rmin.in and rmin_b.in) to those in the output directory.

5) The example usage of the force matching features are demonstrated in the 
   "mscg_driver.cpp" file. Compile this file in the same directory as the
   "read_lammps" files and link to the MSCG library.
   e.g. icc mscg_driver.cpp read_lammps.cpp read_lammps.h -lmscg -o mscg_driver.x

6) The force matched interactions for this model can be determined by running the
   mscg executable:
   e.g. ./mscg_driver.x -l MeOH_example.dat
   
7) Compare the interactions to those in the output directory.


