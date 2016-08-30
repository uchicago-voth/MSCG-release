# Instructions for iterative_fm example #
-----------------------------------------

This example demonstrates the use of the iterative force matching feature.
The input system is the same as the one described in the "lammps_fm" example. As such,
the first run of newfm.x produces the same output as that example.
As with all of the examples, the trajectory is not long enough to contain sufficient 
sampling to produce a converged result.

1) Copy the "control0.in" file as "control.in". Then, run the rangefinder executable 
	using LAMMPS trajectory input to generate the range files (rmin.in and rmin_b.in):
./rangefinder.x -l MeOH_example.dat

2) Compare these range files with those in the "output0" directory.
	
3) For the first run, it is important that the "control.in" file has 
   "primary_output_style 2", "output_solution_flag 1", and 
   "lanyuan_iterative_method_flag 0". These  settings have already been set in 
	"control0.in". In order to use these setting, this file must be copied as 
	"control.in". Since this was already done in step 1, run the force matching executable 
   with LAMMPS trajectory input:
./newfm.x -l MeOH_example.dat

4) Compare the results with those in the "output0" directory. These are the same results
   as the "serial_fm" example.

5) For the next (and any subsequent) run, rename "result.out" to "result.in" and 
   "x.out" to "x.in". 

6) Re-analyze the trajectory using the output interactions. Note: The "MeOH_MeOH.table"
   output is a LAMMPS readable tabulated interaction. Since this code is separate from
   LAMMPS, the results re-analyzing "MeOH_example" using the interactions in
   "MeOH_MeOH.table" are provided as "MeOH_CG.dat".
   
7) It is important that the "control.in" file has "lanyuan_iterative_method_flag 1". The
   "iterative_update_rate_coeff" can also be set to a non-default value (1.0 is default).
   If additional runs are anticipated, than the other flags in mentioned in step 3 must 
   be set as in that step. Copy the "control1.in" file as "control.in". Then, run the 
   force matching executable again using LAMMPS trajectory input on the CG trajectory
   created in the previous step:
./newfm.x -f MeOH_CG.dat

8) Compare the results with those in the "output1" directory. The output interactions 
   are the sum of the previous interaction and the correction from the first iterative 
   force matching run. Note: Additional iterations require the repetition of steps 6 
   through 8.