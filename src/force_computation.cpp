//
//  force_computation.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cassert>
#include <cstdio>
#include <cmath>

#include "force_computation.h"
#include "geometry.h"
#include "interaction_model.h"
#include "matrix.h"
#include "misc.h"
#include "trajectory_input.h"
#include "splines.h"

//--------------------------------------------------------------------
// Prototypes for internal implementation-specific functions
//--------------------------------------------------------------------

// Utility functions for checking if a nonbonded interaction is excluded from the model due to bonding.

bool check_excluded_list(const TopologyData* const topo_data, const int i, const int j);

// Main routine responsible for calling single-element matrix computations,
// differing by the way that potentially interacting particles are found in 
// each frame and possibly found not to interact after.

void order_pair_nonbonded_fm_matrix_element_calculation(InteractionClassComputer* const info, calc_pair_matrix_elements calc_matrix_elements, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, const rvec* x, const real *simulation_box_half_lengths);
void order_bonded_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, const rvec* x, const real *simulation_box_half_lengths);
void order_three_body_nonbonded_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, const rvec* x, const real *simulation_box_half_lengths);

// Helper functions for the above

void process_normal_interaction_matrix_elements(InteractionClassComputer* const info, MATRIX_DATA* const mat, const int n_body, int* particle_ids, std::array<double, 3>* derivatives, const double param_value, const int virial_flag, const double param_deriv, const double distance);

// Functions for calculating individual 3-component matrix elements.

void calc_isotropic_two_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec* x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_angular_three_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec* x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_dihedral_four_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec* x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_nonbonded_1_three_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec* x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_nonbonded_2_three_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec* x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void do_nothing(InteractionClassComputer* const info, const rvec* x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);

//--------------------------------------------------------------------
// Initialization routines to start the FM matrix calculation
//--------------------------------------------------------------------

// Initialize the temps used for calculating all the matrix elements corresponding to a given class of interaction.

void set_up_force_computers(CG_MODEL_DATA* const cg)
{    
    int curr_iclass_col_index = 0;

    // Set up normal case interaction classes.
    std::list<InteractionClassSpec*>::iterator iclass_iterator;
	std::list<InteractionClassComputer*>::iterator icomp_iterator;
	for(icomp_iterator=cg->icomp_list.begin(), iclass_iterator=cg->iclass_list.begin(); icomp_iterator != cg->icomp_list.end(); icomp_iterator++, iclass_iterator++) {
        (*icomp_iterator)->set_up_computer( (*iclass_iterator), &curr_iclass_col_index);
    }

    // Set up three body nonbonded interaction classes.
    cg->three_body_nonbonded_computer.special_set_up_computer(&cg->three_body_nonbonded_interactions, &curr_iclass_col_index);
}

void InteractionClassComputer::set_up_computer(InteractionClassSpec* const ispec_pt, int *curr_iclass_col_index) 
{
    // Store the pointer to the spec.
    ispec = ispec_pt;

    // Set up spline computation for matching and tabulation
    // as needed.
    fm_s_comp = set_up_fm_spline_comp(ispec);
    if (ispec->n_to_force_match > 0) {
        fm_basis_fn_vals = std::vector<double>(fm_s_comp->get_n_coef());
    }
    table_s_comp = set_up_table_spline_comp(ispec);
    
    if (ispec->n_tabulated > 0) {
        table_basis_fn_vals = std::vector<double>(table_s_comp->get_n_coef());
    }

    // Record where this block of interaction basis functions
    // begins in the overall list.
    interaction_class_column_index = *curr_iclass_col_index;
    *curr_iclass_col_index += ispec->interaction_column_indices[ispec->n_to_force_match];

	process_interaction_matrix_elements = process_normal_interaction_matrix_elements;
    // Define the interaction class's geometric definition.
    class_set_up_computer();
}

void PairNonbondedClassComputer::class_set_up_computer(void) 
{
    cutoff2 = ispec->cutoff * ispec->cutoff;
    calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
}

void PairBondedClassComputer::class_set_up_computer(void) 
{
    cutoff2 = 1.0e6;
    calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
}

void AngularClassComputer::class_set_up_computer(void) 
{
    cutoff2 = 1.0e6;
    if (ispec->class_subtype == 1) calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
    else calculate_fm_matrix_elements = calc_angular_three_body_fm_matrix_elements;
}

void DihedralClassComputer::class_set_up_computer(void) 
{
    cutoff2 = 1.0e6;
    if (ispec->class_subtype == 1) calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
    else calculate_fm_matrix_elements = calc_dihedral_four_body_fm_matrix_elements;
}

void ThreeBodyNonbondedClassComputer::special_set_up_computer(InteractionClassSpec* const ispec_pt, int *curr_iclass_col_index)
{
    ispec = ispec_pt;
    switch(ispec->class_subtype) {
    	case 1:
    		calculate_fm_matrix_elements = calc_angular_three_body_fm_matrix_elements;
    		 break;
    		 
    	case 2:
    		if ((ispec->get_basis_type() != kBSpline) && (ispec->get_basis_type() != kBSplineAndDeriv)) {
                printf("Three body with fitted distance term can be only used with B-splines!\n");
                exit(EXIT_FAILURE);
            }
            calculate_fm_matrix_elements = calc_nonbonded_1_three_body_fm_matrix_elements;
        	break;
        
        case 3:
    		calculate_fm_matrix_elements = calc_nonbonded_2_three_body_fm_matrix_elements;
    		break;
    
    	default:
    		break;	
    }
    
    if (ispec->class_subtype > 0) {
        interaction_class_column_index = *curr_iclass_col_index;
        *curr_iclass_col_index += ispec->interaction_column_indices[ispec->n_to_force_match];
        process_interaction_matrix_elements = process_normal_interaction_matrix_elements;
    }
    fm_s_comp = new BSplineAndDerivComputer(ispec);
    fm_basis_fn_vals = std::vector<double>(fm_s_comp->get_n_coef());
}

//--------------------------------------------------------------------
// Main routine calling all other matrix element calculation routines
//--------------------------------------------------------------------

void calculate_frame_fm_matrix(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FrameConfig* const frame_config, PairCellList pair_cell_list, ThreeBCellList three_body_cell_list, int trajectory_block_frame_index)
{
    // Each frame is a set of contiguous rows in the FM matrix; get the starting row for this frame.
    int current_frame_starting_row = trajectory_block_frame_index * cg->n_cg_sites; //shift row number after each frame within one block
    
    // Wrap all coordinates to ensure they are within a single image of
    // the periodic domain and get the target forces for the calculation.
    for (unsigned l = 0; l < cg->topo_data.n_cg_sites; l++) {
        // Enforce consequences of periodic boundary conditions.
        get_minimum_image(l, frame_config->x, frame_config->simulation_box_half_lengths);
        add_target_force_from_trajectory(current_frame_starting_row, l, mat, frame_config->f);
    }
    
    // Set up a cell list and initialize the calculation temps for pair 
    // nonbonded matrix element computations.
    pair_cell_list.populateList(frame_config->current_n_sites, frame_config->x);
    if (cg->three_body_nonbonded_interactions.class_subtype > 0) {
        three_body_cell_list.populateList(frame_config->current_n_sites, frame_config->x);
    }
    
    // Calculate matrix elements by looking through interaction (cell and topology) lists to find active (and non-excluded) interactions.
    std::list<InteractionClassComputer*>::iterator icomp_iterator;
	for(icomp_iterator=cg->icomp_list.begin(); icomp_iterator != cg->icomp_list.end(); icomp_iterator++) {
        (*icomp_iterator)->calculate_interactions(mat, trajectory_block_frame_index, current_frame_starting_row, cg->n_cg_types, cg->topo_data, pair_cell_list, frame_config->x, frame_config->simulation_box_half_lengths);
    }
    cg->three_body_nonbonded_computer.calculate_3B_interactions(mat, trajectory_block_frame_index, current_frame_starting_row, cg->n_cg_types, cg->topo_data, three_body_cell_list, frame_config->x, frame_config->simulation_box_half_lengths);
}

//--------------------------------------------------------------------
// Routines for finding all active interactions to calculate FM matrix elements.
// Exclusion lists are handled in the called subroutines.
//--------------------------------------------------------------------

// Find all neighbors of all particles and call nonbonded matrix element computations for any pairs that interact. 

void PairNonbondedClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, const rvec* x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    walk_neighbor_list(mat, calculate_fm_matrix_elements, n_cg_types, topo_data, pair_cell_list, x, simulation_box_half_lengths);
}

inline void InteractionClassComputer::walk_neighbor_list(MATRIX_DATA* const mat, calc_pair_matrix_elements calc_matrix_elements, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, const rvec* x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    for (int kk = 0; kk < pair_cell_list.size; kk++) {
        k = pair_cell_list.head[kk];
        while (k >= 0) {
            l = pair_cell_list.list[k];
            while (l >= 0) {
                if (check_excluded_list(&topo_data, k, l) == false) {
                    order_pair_nonbonded_fm_matrix_element_calculation(this, calc_matrix_elements, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
                }
                l = pair_cell_list.list[l];
            }
            //do the above the 2nd time for neiboring cells
            for (int nei = 0; nei < 13; nei++) {
                int ll = pair_cell_list.stencil[13 * kk + nei];
                l = pair_cell_list.head[ll];
                while (l >= 0) {
                    if (check_excluded_list(&topo_data, k, l) == false) {
                        order_pair_nonbonded_fm_matrix_element_calculation(this, calc_matrix_elements, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
                    }
                    l = pair_cell_list.list[l];
                }
            }
            k = pair_cell_list.list[k];
        }
    }
}

// Calculate matrix elements for all bonded interactions by looping over the approriate topology lists. 

void PairBondedClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, const rvec* x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    for (k = 0; k < int(topo_data.n_cg_sites); k++) {
        for (unsigned kk = 0; kk < topo_data.bond_list->partner_numbers_[k]; kk++) {
            l = topo_data.bond_list->partners_[k][kk];
            if (k < l) order_bonded_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
        }
    }
}

void AngularClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, const rvec* x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    for (k = 0; k < int(topo_data.n_cg_sites); k++) {
        for (unsigned kk = 0; kk < topo_data.angle_list->partner_numbers_[k]; kk++) {
        	// Grab partners from angle list (organization of angle_list described in topology files).
        	// j is the "center" index while l and k are the "ends" of the angle.
        	// To avoid double counting, the interaction is only counted if the ends are
        	// ordered such that k < l.
            l = topo_data.angle_list->partners_[k][2 * kk + 1];
            j = topo_data.angle_list->partners_[k][2 * kk];
            if (k < l) order_bonded_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
        }
    }
}

void DihedralClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, const rvec* x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;

    trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    for (k = 0; k < int(topo_data.n_cg_sites); k++) {
        for (unsigned kk = 0; kk < topo_data.dihedral_list->partner_numbers_[k]; kk++) {
        	// Grab partners from dihedral list (organization of dihedral_list described in topology files).
        	// i and j are the indices for the "central bond" index while l and k are the "ends" of the dihedral.
        	// To avoid double counting, the interaction is only counted if the ends are
        	// ordered such that k < l.
            l = topo_data.dihedral_list->partners_[k][3 * kk + 2];
            i = topo_data.dihedral_list->partners_[k][3 * kk];
            j = topo_data.dihedral_list->partners_[k][3 * kk + 1];
            if (k < l) order_bonded_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
        }
    }
}

// Calculate matrix elements for three body non-bonded interactions.
// Find all pairs of neighbors of all particles and call nonbonded matrix element computations
// for any triples that interact. Exclusion lists are handled in the called subroutines.

void ThreeBodyNonbondedClassComputer::calculate_3B_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const ThreeBCellList& three_body_cell_list, const rvec* x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    if (ispec->class_subtype > 0) {                    
        trajectory_block_frame_index = traj_block_frame_index;
        current_frame_starting_row = curr_frame_starting_row;
       	walk_3B_neighbor_list(mat, n_cg_types, topo_data, three_body_cell_list, x, simulation_box_half_lengths);
	}
}

inline void InteractionClassComputer::walk_3B_neighbor_list(MATRIX_DATA* const mat, const int n_cg_types, const TopologyData& topo_data, const ThreeBCellList& three_body_cell_list, const rvec* x, const real* simulation_box_half_lengths) 
{
	for (int kk = 0; kk < three_body_cell_list.size; kk++) {
        j = three_body_cell_list.head[kk];
        while (j >= 0) {
            k = three_body_cell_list.head[kk];
            while (k >= 0) {
                if (j != k) {
                    //three body
                    l = three_body_cell_list.list[k];
                    while (l >= 0) {
                        if (l >= 0) {
	                        if (l != j) {
    	                       if ( (check_excluded_list(&topo_data, l, j) == false)  && (check_excluded_list(&topo_data, j, k) == false) ) {
        	                       order_three_body_nonbonded_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
            	               }
            	        	}
                        	l = three_body_cell_list.list[l];
						}
					}
                    for (int nei_3 = 0; nei_3 < 26; nei_3++) {
                        int ll_3 = three_body_cell_list.stencil[26 * kk + nei_3];
                        l = three_body_cell_list.head[ll_3];
                        while (l >= 0) {
                            if ( (check_excluded_list(&topo_data, l, j) == false)  && (check_excluded_list(&topo_data, j, k) == false) ) {
                                order_three_body_nonbonded_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
                            }
                            l = three_body_cell_list.list[l];
                        }
                    }
                }
                k = three_body_cell_list.list[k];
            }
            
            for (int nei = 0; nei < 26; nei++) {
                int ll = three_body_cell_list.stencil[26 * kk + nei];
                k = three_body_cell_list.head[ll];
                while (k >= 0) {
                    //three body
                    l = three_body_cell_list.list[k];
                    while (l >= 0) {
                        if ( (check_excluded_list(&topo_data, l, j) == false)  && (check_excluded_list(&topo_data, j, k) == false) ) {
                                order_three_body_nonbonded_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
                        }
                        l = three_body_cell_list.list[l];
                    }
                    for (int nei_3 = nei + 1; nei_3 < 26; nei_3++) {
                        int ll_3 = three_body_cell_list.stencil[26 * kk + nei_3];
                        l = three_body_cell_list.head[ll_3];
                        while (l >= 0) {
                            if ( (check_excluded_list(&topo_data, l, j) == false)  && (check_excluded_list(&topo_data, j, k) == false) ) {
                                order_three_body_nonbonded_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
                            }
                            l = three_body_cell_list.list[l];
                        }
                    }
                    k = three_body_cell_list.list[k];
                }
            }
            j = three_body_cell_list.list[j];
        }
    }
}

//--------------------------------------------------------------------
//  Routine for checking if nonbonded interactions should be excluded
// from the model because they are between bonded particles.
//--------------------------------------------------------------------

inline bool check_excluded_list(const TopologyData* const topo_data, const int i, const int j)
{
    // Check whether this non-bonded interaction is excluded from the model
    for (unsigned k = 0; k < topo_data->exclusion_list->partner_numbers_[i]; k++) {
        if (topo_data->exclusion_list->partners_[i][k] == unsigned(j)) return true;
    }
    return false;
}

//--------------------------------------------------------------------
// Routines responsible for calling single-interaction matrix computations,
// differing by the way that potentially interacting particles are found in 
// each frame and possibly found not to interact after.
//--------------------------------------------------------------------

void order_pair_nonbonded_fm_matrix_element_calculation(InteractionClassComputer* const info, calc_pair_matrix_elements calc_matrix_elements, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, const rvec* x, const real *simulation_box_half_lengths)
{
    // Calculate the appropriate matrix elements.
    info->index_among_defined_intrxns = info->ispec->get_index_from_hash(calc_two_body_interaction_hash(cg_site_types[info->k], cg_site_types[info->l], n_cg_types));
    info->index_among_matched_interactions = info->ispec->defined_to_matched_intrxn_index_map[info->index_among_defined_intrxns];
    info->index_among_tabulated_interactions = info->ispec->defined_to_tabulated_intrxn_index_map[info->index_among_defined_intrxns];
    
	if (info->index_among_matched_interactions == 0) return; // if the index is zero, it is not present in the model and should be ignored.
    calc_matrix_elements(info, x, simulation_box_half_lengths, mat);
}

void order_bonded_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, const rvec* x, const real *simulation_box_half_lengths)
{
     // Calculate the appropriate matrix elements.    
    info->index_among_defined_intrxns = info->ispec->get_index_from_hash(info->calculate_hash_number(cg_site_types, n_cg_types));
    info->index_among_matched_interactions = info->ispec->defined_to_matched_intrxn_index_map[info->index_among_defined_intrxns];
    info->index_among_tabulated_interactions = info->ispec->defined_to_tabulated_intrxn_index_map[info->index_among_defined_intrxns];
    
    if (info->index_among_matched_interactions == 0) return; // if the index is zero, it is not present in the model and should be ignored.

    (*info->calculate_fm_matrix_elements)(info, x, simulation_box_half_lengths, mat);
}

void order_three_body_nonbonded_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, const rvec* x, const real *simulation_box_half_lengths)
{
    ThreeBodyNonbondedClassComputer* icomp = static_cast<ThreeBodyNonbondedClassComputer*>(info);
    ThreeBodyNonbondedClassSpec* ispec = static_cast<ThreeBodyNonbondedClassSpec*>(icomp->ispec);
    
    // Calculate the appropriate matrix elements.
    icomp->index_among_defined_intrxns = info->ispec->get_index_from_hash(icomp->calculate_hash_number(cg_site_types, n_cg_types));
    if (icomp->index_among_defined_intrxns == -1) return; // if the index is -1, it is not present in the model and should be ignored.
    
    icomp->index_among_matched_interactions = ispec->defined_to_matched_intrxn_index_map[icomp->index_among_defined_intrxns];
    icomp->index_among_tabulated_interactions = ispec->defined_to_tabulated_intrxn_index_map[icomp->index_among_defined_intrxns];
    if ((icomp->index_among_matched_interactions == 0) && (icomp->index_among_tabulated_interactions == 0)) return; // if the index is zero, it is not present in the model and should be ignored.
    
    icomp->cutoff2 = ispec->three_body_nonbonded_cutoffs[icomp->index_among_defined_intrxns] * ispec->three_body_nonbonded_cutoffs[icomp->index_among_defined_intrxns];
    icomp->stillinger_weber_angle_parameter = ispec->stillinger_weber_angle_parameters_by_type[icomp->index_among_defined_intrxns];
    (*icomp->calculate_fm_matrix_elements)(icomp, x, simulation_box_half_lengths, mat); 
}

//--------------------------------------------------------------------
// Functions for calculating sets of 3-component matrix elements for each
// individual interacting set of particles
//--------------------------------------------------------------------

inline void process_normal_interaction_matrix_elements(InteractionClassComputer* const info, MATRIX_DATA* const mat, const int n_body, int* particle_ids, std::array<double, 3>* derivatives, const double param_value, const int virial_flag, const double junk = 0.0, const double junk2 = 0.0)
{
	int index_among_defined = info->index_among_defined_intrxns;
	int index_among_matched = info->index_among_matched_interactions;
    int index_among_tabulated = info->index_among_tabulated_interactions;
    int first_nonzero_basis_index;
    int temp_column_index;
    double basis_sum;
    
    if (index_among_tabulated > 0) {
		// Pull the interaction from a table. 	   
    	info->table_s_comp->calculate_basis_fn_vals(index_among_defined, param_value, first_nonzero_basis_index, info->table_basis_fn_vals);
    	basis_sum  = info->table_basis_fn_vals[0] + info->table_basis_fn_vals[1];
    	
    	// Add to force target.
		mat->accumulate_tabulated_forces(info, basis_sum, n_body, particle_ids, derivatives, mat);
    	
    	// Add to target virial if virial_flag is non-zero.
    	switch (virial_flag) {
    		case 1:
	    	    if (mat->virial_constraint_rows > 0) mat->accumulate_target_constraint_element(mat, info->trajectory_block_frame_index, -basis_sum * param_value);
        		break;
        	
        	case 0: default:
    			// These interactions do not contribute to the scalar virial.
    			// Such interactions include angles and dihedrals.
        		break;
    	}
	}

    if (index_among_matched > 0) {
	    // Compute the strength of each basis function.
	    info->fm_s_comp->calculate_basis_fn_vals(index_among_defined, param_value, first_nonzero_basis_index, info->fm_basis_fn_vals);
    	
    	// Add to the force matching.
	    mat->accumulate_matching_forces(info, first_nonzero_basis_index, info->fm_basis_fn_vals, n_body, particle_ids, derivatives, mat);
 			
    	// Add to virial matching if virial_flag is non-zero.
    	switch (virial_flag) {
    		case 1:
	    	    temp_column_index = info->interaction_class_column_index + info->ispec->interaction_column_indices[index_among_matched - 1] + first_nonzero_basis_index;
	    		for (unsigned i = 0; i < info->fm_basis_fn_vals.size(); i++) {
        			int basis_column = temp_column_index + i;
        			if (mat->virial_constraint_rows > 0)(*mat->accumulate_virial_constraint_matrix_element)(info->trajectory_block_frame_index, basis_column, info->fm_basis_fn_vals[i] * param_value, mat);
        		}
        		break;
        	
        	case 0: default:
    			// These interactions do not contribute to the scalar virial.
    			// Such interactions include angles and dihedrals.
        		break;
    	}
	}    
}

//--------------------------------------------------------------------
// Functions for calculating sets of 3-component matrix elements for each
// individual interacting set of particles
//--------------------------------------------------------------------

// Each of these functions follows the idiom of calc_isotropic_two_body_fm_matrix_elements.

void calc_isotropic_two_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec* x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    int particle_ids[2] = {info->k, info->l};
    std::array<double, 3>* derivatives = new std::array<double, 3>[1];
	double distance;
	if ( conditionally_calc_distance_and_derivatives(particle_ids, x, simulation_box_half_lengths, info->cutoff2, distance, derivatives) ) {
        int index_among_defined = info->index_among_defined_intrxns;
    	if (distance < info->ispec->lower_cutoffs[index_among_defined] ||
        	distance > info->ispec->upper_cutoffs[index_among_defined]) {
        	delete [] derivatives;
        	return;
        }
    	info->process_interaction_matrix_elements(info, mat, 2, particle_ids, derivatives, distance, 1, 0.0 , 0.0);
    }
    delete [] derivatives;
}

void calc_angular_three_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec* x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    int particle_ids[3] = {info->k, info->l, info->j}; // end indices (k, l), followed by center index (j)
    std::array<double, 3>* derivatives = new std::array<double, 3>[2];
    int index_among_defined = info->index_among_defined_intrxns;
    double angle;

    if ( conditionally_calc_angle_and_derivatives(particle_ids, x, simulation_box_half_lengths, info->cutoff2, angle, derivatives) ) {
        if (angle < info->ispec->lower_cutoffs[index_among_defined] ||
        	angle > info->ispec->upper_cutoffs[index_among_defined]) {
        	delete [] derivatives;
        	return;
        }
        info->process_interaction_matrix_elements(info, mat, 3, particle_ids, derivatives, angle, 0, 0.0, 0.0);
    }
    delete [] derivatives;
}

void calc_dihedral_four_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec* x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    int particle_ids[4] = {info->k, info->l, info->i, info->j}; // end indices (k, l) followed by central bond indices (i, j)
    std::array<double, 3>* derivatives = new std::array<double, 3>[3];
    int index_among_defined = info->index_among_defined_intrxns;
    double dihedral;

    if ( conditionally_calc_dihedral_and_derivatives(particle_ids, x, simulation_box_half_lengths, info->cutoff2, dihedral, derivatives) ) {
        if (dihedral < info->ispec->lower_cutoffs[index_among_defined] ||
        	dihedral > info->ispec->upper_cutoffs[index_among_defined]) {
        	delete [] derivatives;
        	return;
        }
		info->process_interaction_matrix_elements(info, mat, 4, particle_ids, derivatives, dihedral, 0, 0.0, 0.0);
    }
	delete [] derivatives;
}

void calc_nonbonded_1_three_body_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, 3>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    int particle_ids[3] = {info->k, info->l, info->j}; // end indices (k, l) followed by center index (j).    
    ThreeBodyNonbondedClassComputer* icomp = static_cast<ThreeBodyNonbondedClassComputer*>(info);
    ThreeBodyNonbondedClassSpec* ispec = static_cast<ThreeBodyNonbondedClassSpec*>(icomp->ispec);

	std::array<double, 3>* relative_site_position_2 = new std::array<double, 3>[1];
	std::array<double, 3>* relative_site_position_3 = new std::array<double, 3>[1];
	std::array<double, 3>* derivatives = new std::array<double, 3>[2];
	std::array<double, 3> tx1, tx2, tx;
	double theta, rr1, rr2;
    double angle_prefactor, dr1_prefactor, dr2_prefactor;
    int	this_column;
	
	bool within_cutoff = conditionally_calc_sw_angle_and_intermediates(particle_ids, x, simulation_box_half_lengths, ispec->three_body_nonbonded_cutoffs[icomp->index_among_defined_intrxns], ispec->three_body_gamma, relative_site_position_2, relative_site_position_3, derivatives, theta, rr1, rr2, angle_prefactor, dr1_prefactor, dr2_prefactor);
	if (!within_cutoff) {
		delete [] relative_site_position_2;
		delete [] relative_site_position_3;
		delete [] derivatives;
		return;
    }

    icomp->intrxn_param = theta;
  
    // Calculate the matrix elements if it's supposed to be force matched
    info->fm_s_comp->calculate_basis_fn_vals(info->index_among_defined_intrxns, info->intrxn_param, info->basis_function_column_index, info->fm_basis_fn_vals); 
    std::vector<double> basis_der_vals(info->fm_s_comp->get_n_coef());
    BSplineAndDerivComputer *fm_s_comp = static_cast<BSplineAndDerivComputer*>(icomp->fm_s_comp);
    fm_s_comp->calculate_bspline_deriv_vals(info->index_among_defined_intrxns, info->intrxn_param, info->basis_function_column_index, basis_der_vals); 
    
    int temp_row_index_1 = particle_ids[0] + icomp->current_frame_starting_row;
    int temp_row_index_2 = particle_ids[2] + icomp->current_frame_starting_row;
    int temp_row_index_3 = particle_ids[1] + icomp->current_frame_starting_row;
	int temp_column_index = icomp->interaction_class_column_index + ispec->interaction_column_indices[icomp->index_among_matched_interactions - 1] + icomp->basis_function_column_index;
        
    for (unsigned i = 0; i < info->fm_basis_fn_vals.size(); i++) {

        this_column = temp_column_index + i;
        for (int j = 0; j < DIMENSION; j++) {
        	tx1[j] = derivatives[0][j] * angle_prefactor * basis_der_vals[i] + 0.5 * dr1_prefactor * (relative_site_position_2[0][j] / rr1) * info->fm_basis_fn_vals[i]; // derivative of angle plus derivative of distance for site 0 (K)
        	tx2[j] = derivatives[1][j] * angle_prefactor * basis_der_vals[i] + 0.5 * dr2_prefactor * (relative_site_position_3[0][j] / rr2) * info->fm_basis_fn_vals[i]; // derivative of angle plust derivative of distance for site 2 (L)
        	tx[j]  = - (tx1[j] + tx2[j]); // Use Newton's third law to determine for on central site
        }
        
        (*mat->accumulate_fm_matrix_element)(temp_row_index_1, this_column, &tx1[0], mat);
        (*mat->accumulate_fm_matrix_element)(temp_row_index_2, this_column, &tx2[0], mat);
        (*mat->accumulate_fm_matrix_element)(temp_row_index_3, this_column, &tx[0], mat);
    }
    delete [] relative_site_position_2;
	delete [] relative_site_position_3;
    delete [] derivatives;
}

void calc_nonbonded_2_three_body_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, 3>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    int particle_ids[3] = {info->k, info->l, info->j}; // end indices (k, l) followed by center index (j).    
    ThreeBodyNonbondedClassComputer* icomp = static_cast<ThreeBodyNonbondedClassComputer*>(info);
    ThreeBodyNonbondedClassSpec* ispec = static_cast<ThreeBodyNonbondedClassSpec*>(icomp->ispec);
    
	std::array<double, 3>* relative_site_position_2 = new std::array<double, 3>[1];
	std::array<double, 3>* relative_site_position_3 = new std::array<double, 3>[1];
	std::array<double, 3>* derivatives = new std::array<double, 3>[2];
	std::array<double, 3> tx1, tx2, tx;
    double theta, rr1, rr2;
    double cos_theta;
    double angle_prefactor, dr1_prefactor, dr2_prefactor;
    double u, du;
    
    bool within_cutoff = conditionally_calc_sw_angle_and_intermediates(particle_ids, x, simulation_box_half_lengths, ispec->three_body_nonbonded_cutoffs[icomp->index_among_defined_intrxns], ispec->three_body_gamma, relative_site_position_2, relative_site_position_3, derivatives, theta, rr1, rr2, angle_prefactor, dr1_prefactor, dr2_prefactor);
	if (!within_cutoff) {
		delete [] relative_site_position_2;
		delete [] relative_site_position_3;
		delete [] derivatives;
		return;
    }

    icomp->intrxn_param = theta;
    theta /= DEGREES_PER_RADIAN;
    cos_theta = cos(theta);
        
    u = (cos_theta - icomp->stillinger_weber_angle_parameter) * (cos_theta - icomp->stillinger_weber_angle_parameter) * 4.184;
    du = 2.0 * (cos_theta - icomp->stillinger_weber_angle_parameter) * sin(theta) * 4.184;
    
    int temp_row_index_2 = particle_ids[2] + icomp->current_frame_starting_row;
    int temp_row_index_3 = particle_ids[1] + icomp->current_frame_starting_row;
    int temp_row_index_1 = particle_ids[0] + icomp->current_frame_starting_row;
    int temp_column_index = icomp->interaction_class_column_index + ispec->interaction_column_indices[icomp->index_among_matched_interactions - 1];
        
    for (int j = 0; j < DIMENSION; j++) {
    	tx1[j] = derivatives[0][j] * angle_prefactor * du + 0.5 * dr1_prefactor * u * (relative_site_position_2[0][j] / rr1); // derivative of angle (with harmonic cosine) plus derivative of distance for site 0 (K)
    	tx2[j] = derivatives[1][j] * angle_prefactor * du + 0.5 * dr2_prefactor * u * (relative_site_position_3[0][j] / rr1); // derivative of angle (with harmonic cosine) plus derivative of distance for site 2 (L)
    	tx[j]  = - (tx1[j] + tx2[j]); // Use Newton's third law to determine for on central site
    }
    
    (*mat->accumulate_fm_matrix_element)(temp_row_index_1, temp_column_index, &tx1[0], mat);
    (*mat->accumulate_fm_matrix_element)(temp_row_index_2, temp_column_index, &tx2[0], mat);
    (*mat->accumulate_fm_matrix_element)(temp_row_index_3, temp_column_index, &tx[0], mat); 
    
    delete [] relative_site_position_2;
	delete [] relative_site_position_3;   
	delete [] derivatives;
}

void do_nothing(InteractionClassComputer* const info, const rvec* x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat) 
{
}