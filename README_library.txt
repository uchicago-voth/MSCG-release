MSCG, a library for performing multi-scale coarse-graining 
https://github.com/uchicago-voth/MSCG-release
By Jacob W. Wagner (UChicago in Gregory A. Voth's group)

Introduction
================================================================
MSCG is a software library for carrying out multi-scale coarse-graining
for the parameterization of coarse-grained (CG) interactions. 

MSCG comprises of several C/C++ classes/structs that can be built as a 
static library and linked to. 
Numerous examples are provided to demonstrate the library's features.

This code can also be used to generate a stand alone executable.
Please see the MSCGFM Manual for details.

Compilation - Linux / Mac OS / Windows with Cygwin
==================================================
The code is written in ANSI C++11, and compiles on many system architectures. The
package contains the C++ source code, example files, miscellaneous utilities
and documentation. On Linux, Mac OS, and Windows (using Cygwin), the
compilation and installed can be carried out using GNU Make.

To begin, the user should review the appropriate Makefile in the Make
directory inside of the src directory. The user should make sure that the compilation and 
installation settings are appropriate for their system. The GSL path and include 
variables often need to be updated. Making "lib_mscg.a" will then compile the static
library. The library will appear within the "src" directory.

Following successful compilation, the library can be used directly by
compiling a driver code (one that calls library functions) with the path for the library 
following '-L' and '-lmscg'. Alternately, the driver code can be compiled with the library 
by appending the path and name of the library to the list of code files.

MSCG will perform parallel (OpenMP) matrix calculations if the linked 
GSL or MKL libraries support it.


Compilation - Windows without Cygwin
====================================
On a Windows machine without a UNIX-like terminal environment like Cygwin, it is possible
to import and compile the library in many standard C++ development
environments. Users have reported success in building the library with
Microsoft Visual C++ Express and Code::Blocks.


Related programs
================
The required external dependency for MSCG is the GNU Scientific Library (GSL).
Optionally, sparse matrix operations require the Intel Math Kernel Library (MKL).
Also, LAPACK or MKL may be required for certain matrix operations depending on your
compilation settings. The GROMACS (GMX) variables are only used for compiling the code 
as a stand-alone executable.

Additionally, the freeware plotting program Gnuplot (available at www.gnuplot.info) 
can be used for rapid visualization of the program output interactions.

Files/Directories
=================
examples - many examples making use of the library with expected output
scripts - miscellaneous helper scripts such as a mapping script
src - source code files
MAKE - various makefiles


Developer Usage
===============
The MSCG library was developed for use with a fix in LAMMPS.
The library functions carry internal state through the void * (opaque pointer).
Thus, the order that the functions are called is important in determining their behavior.
The required order of calling functions from this library for
proper operation is the following:
mscg_startup_part1
setup_topology_and_frame
setup_bond_topology
either setup_angle_topology and setup_dihedral_topology and setup_exclusion_topology,
or setup_angle_topology and setup_dihedral_topology and generate_exclusion_topology,
or generate_angle_dihedral_and_exclusion_topology
mscg_startup_part2
mscg_process_frame for each frame
mscg_solve_and output.

Additional functions are provided for updating or modifying certain information
after the data is initially set using one of the functions before or during 
mscg_startup_part2.

The setup_*_topology functions expect 2 arrays. The first lists the number of
bonds/angles/dihedrals that begin/end at that site. The second array lists all sites
involved in the interactions mentioned in the first array. For bonds, there should be 1 
site per interaction, 2 for angles, and 3 for dihedrals.

To use the range-finding utility, change the "mscg_" prefix for the functions above to
"rangefinder_". The arguments to the "mscg_" and "rangefinder_" functions are the same.

Example usage of the force matching (mscg) and range finding capabilities are 
demonstrated in the "MSCG_library_example" sub-directory of the examples.
Eventually, the LAMMPS "mscg" fix will also be an accessible demonstration of this 
library.

At the moment the only feature in MSCGFM that is not currently supported in this
library is dynamic state sampling for use with the ultra-coarse-graining methodology.

Input Files and Usage
=====================
The required input files for running this library are the following:
* control.in 
* rmin.in/rmin_b.in
* Note: The top.in file is not used by this library. 

The following control.in parameters are largely disregarded by the library form 
of the code: n_frames and start_frame. It is only important that n_frames is divisible 
by the block_size.

Note: For dense matrices, the block_size is always 1.

Note: If the total number of frames handed to the library is not divisible by the
block_size, the last 'x' frames read are disregarded to make the total number of 
frames used to solve for interactions divisible by the block size.

Note: For three body nonbonded interactions, the information normally specified in
the "top.in" file under the "threebody" label should instead be specified in the
"top.fix" file using the same format.

Example Usage
=============

A series of files driving the range finding and force matching utilities for a simple
1-site MeOH model (the same trajectory as the usual "lammps_fm" example) can be found in
the "MSCG_library_example" sub-directory of the examples.

Copyright and Citation
======================

MSCG is released as free software through the University of Chicago - 
a detailed copyright notice is provided below, and the complete terms 
of the license can be found in the LICENSE file.

I am very interested to hear from users of the software, so if you find this
useful, please email me the group at gavoth@uchicago.edu. Also, if you plan 
to publish an academic paper using this software, please consider citing one 
of the following publications:

1) S. Izvekov and G. A. Voth, "Multiscale coarse graining of liquid-state systems", 
	J. Chem. Phys. 123, 134105 (2005). doi:10.1063/1.2038787
	
2) W. G. Noid, J-W. Chu, G. S. Ayton, V. Krishna, S. Izvekov, G. A. Voth, A. Das, and 
	H. C. Andersen, "The multiscale coarse-graining method. I. A rigorous bridge between 
	atomistic and coarse-grained models" J. Chem. Phys. 128, 244114 (2008).
	doi:10.1063/1.2938860
	
3) W. G. Noid, P. Liu, Y. Wang, J-W. Chu, G. S. Ayton, S. Izvekov, H. C. Andersen, and 
	G. A. Voth, "The multiscale coarse-graining method. II. Numerical implementation for 
	coarse-grained molecular models", J. Chem. Phys. 128, 244115 (2008). 
	doi:10.1063/1.2938857
	
4) L. Lu, S. Izvekov, A. Das, H. C. Andersen, and G. A. Voth, "Efficient, Regularized, and 
  Scalable Algorithms for Multiscale Coarse-Graining", Journal of Chemical Theory and 
  Computation, 6(3), 954-965 (2010). doi:10.1021/ct900643r

5) The collaboration between the LAMMPS developers and the Voth Group. This reference
	will be updated as more information becomes available.

The first paper describes one of the first major applications of "force matching"
atomistic data in order to determine a coarse-grained interaction. It is also 
the first use of the virial constraint.

The second paper describes the theoretical framework for the method.

The third paper describes additional implementation details.

The fourth paper describes developments in the code upto 2010. Note: Not all of those
features are supported in this release. This is largely because the increase in 
computing power and sparse matrix solvers have made the approximation used by these 
features unnecessary.

Details on the fifth reference will be provided as they become available. Primarily,
this reference will describe the use of the new MSCG fix in LAMMPS as well as 
developments in this code since the 2010 paper.

Copyright Notice
================
MSCG Copyright (c) 2016, The Voth Group at The University of Chicago. All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact the University of Chicago's Technology Transfer Office at 
info-uchicagotech@uchicago.edu.


Acknowledgments
===============
This was supported with United States Government support under and awarded
 by the Department of Defense, High Performance Computing Modernization Project
(HPCMP), under the National Defense Science and Engineering Graduate (NDSEG) 
Fellowship Program, 32 CFR 168 a. Support was also provided by the Office of 
Naval Research, under grant N00014-13-1-0058. Additional support was provided
through a Collaborative Research in Chemistry grant from the 
National Science Foundation (Grant No. CHE-0628257).