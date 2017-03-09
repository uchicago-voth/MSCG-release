//
//  mscg.cpp
//  
//
//  Created by Jacob Wagner on 04/27/16 for use with USER-MSCG package in LAMMPS.
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

// This file provides the library interface to the multi-scale coarse-graining (MS-CG)
// code, which is also referred to as the force matching (FM) code. 
// The order in which functions are intended to be called for FM is the following:
// mscg_startup_part1
// setup_topology_and_frame
// setup_bond_topology, setup_angle_topology, and setup_dihedral_topology
// either setup_exclusion_topology or generate_exclusion_topology
// mscg_startup_part2
// mscg_process_frame for each frame
// mscg_solve_and output.
// Additional functions are provided for updating or modifying certain information
// after the data is initially set using one of the functions before or during 
// mscg_startup_part2.

// The order in which functions are intended to be called for range finding is the following:
// rangefinder_startup_part1
// setup_topology_and_frame
// setup_bond_topology, setup_angle_topology, and setup_dihedral_topology
// either setup_exclusion_topology or generate_exclusion_topology
// rangefinder_startup_part2
// rangefinder_process_frame for each frame
// rangefinder_solve_and output.

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include "control_input.h"
#include "force_computation.h"
#include "fm_output.h"
#include "interaction_hashing.h"
#include "interaction_model.h"
#include "matrix.h"
#include "misc.h"
#include "trajectory_input.h"

// Solely for rangefinding capability.
#include "range_finding.h"

void* mscg_startup_part1(void* void_in);
void* rangefinder_startup_part1(void* void_in);
void* mscg_startup_part2(void* void_in);
void* rangefinder_startup_part2(void* void_in);
void* rangefinder_process_frame(void* void_in, double* const x, double* const f);
void* mscg_process_frame(void* void_in, double* const x, double* const f);
void* mscg_solve_and_output(void* void_in);
void* rangefinder_solve_and_output(void* void_in);

void* setup_frame_config(void* void_in, const int n_cg_sites, int * cg_site_types, double* box_half_lengths);
void* update_frame_config(void* void_in, const int n_cg_sites, int * cg_site_types, double* box_half_lengths);

void* setup_topology_and_frame(void* void_in, int const n_cg_sites, int const n_cg_types, char ** type_names, int* cg_site_types, double* box_half_lengths);
void* set_bond_topology(void* void_in, unsigned** bond_partners, unsigned* bond_partner_numbers); 
void* set_angle_topology(void* void_in, unsigned** angle_partners, unsigned* angle_partner_numbers);
void* set_dihedral_topology(void* void_in, unsigned** dihedral_partners, unsigned* dihedral_partner_numbers);

void* set_exclusion_topology(void* void_in, unsigned** exclusion_partners, unsigned* exclusion_partner_numbers);
void* generate_exclusion_topology(void* void_in);
void* generate_angle_dihedral_and_exclusion_topology(void* void_in);

// Prototype function definition for functions called internal to this file
void finish_fix_reading(FrameSource *const frame_source);

// Data structure holding all MSCG information.
// It is passed to the driver function (LAMMPS fix) as an opaque pointer.
struct MSCG_struct {
	int curr_frame;
	int nblocks;
	int trajectory_block_frame_index;
	int traj_frame_num;
	double start_cputime;
	FrameSource *frame_source;      // Trajectory frame data
    CG_MODEL_DATA *cg;  			// CG model parameters and data
    ControlInputs *control_input;	// Input settings read from control.in
    MATRIX_DATA *mat;				// Matrix storage structure
};