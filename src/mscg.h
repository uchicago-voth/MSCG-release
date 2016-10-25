#ifndef MSCG_H
#define MSCG_H

#include "batch_fm_combination.h"
#include "control_input.h"
#include "external_matrix_routines.h"
#include "fm_output.h"
#include "force_computation.h"
#include "geometry.h"
#include "interaction_hashing.h"
#include "interaction_model.h"
#include "matrix.h"
#include "misc.h"
#include "range_finding.h"
#include "splines.h"
#include "topology.h"
#include "trajectory_input.h"

void* mscg_startup_part1(void* void_in);
void* rangefinder_startup_part1(void* void_in);
void* mscg_startup_part2(void* void_in);
void* rangefinder_startup_part2(void* void_in);
void* mscg_process_frame(void* void_in, double* const x, double* const f);
void* rangefinder_process_frame(void* void_in, double* const x, double* const f);
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
int get_n_frames(void* void_in);
int get_block_size(void* void_in);

#endif
