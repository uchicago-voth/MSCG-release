# Instructions for table_fm example #
-------------------------------------

This example demonstrates the use of Bayesian MS-CG.
The input system is a system of 250 methanol molecules that have each been mapped
to two center-of-mass sites (one for the methyl and the other for the hydroxyl). 
This is the same input trajectory as the table-fm example.
In this case, all nonbonded and bonded interactions will be determined using MS-CG.
Bayesian MS-CG is used to iteratively determine the appropriate regularization to
improve the robustness of the determined interactions (especially in cases with poorly
sampled or resolved areas).
As with all of the examples, the trajectory is not long enough to contain sufficient 
sampling to produce a converged result.

1) Note that the range files (rmin.in and rmin_b.in) have already been prepared.
	
2) To generate the usual MS-CG interactions, run the force matching executable
 with LAMMPS trajectory input (after copying the control_base.in to control.in):
./newfm.x -l cg2site.lammpstrjnew.short

3) Compare these results with the files in the "base-output" directory.

4) To generate results after 10 iterations of Bayesian MS-CG, set the
bayesian_mscg_flag to 1 and bayesian_max_iterations to 10. This has been
done in control_bayesian.in. So, copy the file to control.in. Then, run the
force matching executable with LAMMPS trajectory input:
./newfm.x -l cg2site.lammpstrjnew.short

5) Compare these results with the files in the "bayesian-output" directory.
One way to determine if the iterations are convereged to look at the regularization
parameters in the alpha.out and beta.out files. Another option is look at the
basis coefficients determined for each iteration in the solution.out file.

6) To see the effect of Bayesian force matching, compare the output interactions
(the *.dat and *.table) files from these two runs. In particular, the bonded
interation 1_2_bond.dat is expected to show a difference.
