Force Matching for threebody interactions
-----------------------------------------------------------------
Reference: L. Larini, L. Lu and G. A. Voth, JCP, 132. 164107(2010)

1. Phase I: Determine the equilibrium threebody interaction angle theta0

Go to the phase1-angle directory, it should be noted that in order to scan
over the angle theta, one need to set parameter three_body_nonbonded_style=2,
and the top.in file should have a line describing the cutoff distance for
evaluating threebody interactions. Run the script and look at the output file
1_1_1.dat, find the theta that minimize the output, this would be the
equilibrium angle theta0.

2. Phase II: Parameterize the stillinger-weber interaction strength by Force
Matching

After you find the equilibrium angle theta0, calculate cos(theta0), and change
the fourth parameter of threebody interactions in top.in file to cos(theta0).
Change the three_body_nonbonded_style=3 in control.in. Run the force matching
again, you would find the output stillinger-weber interaction strength in the
file 3b.dat, and the pairwise force can be found in 1_1.dat.

