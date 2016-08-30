//
//  force_computation.h
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#ifndef _force_computation_h
#define _force_computation_h

#include <array>

#include "trajectory_input.h"
#include "interaction_model.h"

struct MATRIX_DATA;

// Initialization routines to start the FM matrix calculation

void set_up_force_computers(CG_MODEL_DATA* const cg);

// Main routine calling all other matrix element calculation routines

void calculate_frame_fm_matrix(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FrameConfig* const frame_config, PairCellList pair_cell_list, ThreeBCellList three_body_cell_list, int trajectory_block_frame_index);

// Functions for calculating the parameters that the potential basis
// functions depend on: distances, angles, and dihedrals.

double calc_squared_distance(const int k, const int l, const rvec *x, const real *simulation_box_half_lengths, double* const relative_site_position_1);
double calc_cosine_of_angle_and_intermediates(const int j, const int k, const int l, const rvec *x, const real *simulation_box_half_lengths, double* const relative_site_position_2, double* const relative_site_position_3, double* const rr1, double* const rr2);
double calc_dihedral_and_intermediates(const int i, const int j, const int k, const int l, const rvec *x, const real *simulation_box_half_lengths, double* const relative_site_position_2, double* const relative_site_position_3, double* const relative_site_position_4, double* const pb, double* const pc, double* const rpb1, double* const rpc1, double* const pbpc, double* const s);

// Free the temps used in FM matrix building, retaining only what is still needed for solution and output.

void free_fm_matrix_building_temps(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FrameSource* const frame_source);

#endif
