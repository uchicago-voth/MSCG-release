//
//  force_computation.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include "force_computation.h"

#include <cassert>
#include <cstdio>
#include <cmath>

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

bool check_anybond_nonbonded_exclusion(const TopologyData* const topo_data, const int i, const int j);
bool check_pairbond_nonbonded_exclusion(const TopologyData* const topo_data, const int i, const int j);
bool check_excluded_list(const TopologyData* const topo_data, const int i, const int j);

// Main routine responsible for calling single-element matrix computations,
// differing by the way that potentially interacting particles are found in 
// each frame and possibly found not to interact after.

void order_pair_nonbonded_fm_matrix_element_calculation(InteractionClassComputer* const info, calc_pair_matrix_elements calc_matrix_elements, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, const rvec *x, const real *simulation_box_half_lengths);
void order_bonded_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, const rvec *x, const real *simulation_box_half_lengths);
void order_three_body_nonbonded_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, const rvec *x, const real *simulation_box_half_lengths);

// Functions for converting spline_coefficients and derivatives into matrix elements.

template<int n_body> void accumulate_table_forces(InteractionClassComputer* const info, const double &table_fn_val, const std::array<int, n_body> &particle_ids, const std::array<std::array<double, 3>, n_body - 1> &derivatives, MATRIX_DATA * const mat);
template<int n_body> void accumulate_matching_forces(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const std::array<int, n_body> &particle_ids, const std::array<std::array<double, 3>, n_body - 1> &derivatives, MATRIX_DATA * const mat);

// Functions for calculating individual 3-component matrix elements.

void calc_isotropic_two_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_angular_three_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_dihedral_four_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_nonbonded_1_three_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_nonbonded_2_three_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);

// Utility functions for calculating angles, distances, etc. in a periodic domain.

void get_minimum_image(const int l, rvec* frx, const real *simulation_box_half_lengths);
double dot_product(const double* a, const double* b);

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
    cg->three_body_nonbonded_computer.special_set_up_computer( &cg->three_body_nonbonded_interactions, &curr_iclass_col_index);
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

    // Define the interaction class's geometric definition.
    class_set_up_computer();
}

void PairNonbondedClassComputer::class_set_up_computer(void) {
    cutoff2 = ispec->cutoff * ispec->cutoff;
    calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
}

void PairBondedClassComputer::class_set_up_computer(void) {
    cutoff2 = 1.0e6;
    calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
}

void AngularClassComputer::class_set_up_computer(void) {
    cutoff2 = 1.0e6;
    if (ispec->class_subtype == 1) calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
    else calculate_fm_matrix_elements = calc_angular_three_body_fm_matrix_elements;
}

void DihedralClassComputer::class_set_up_computer(void) {
    cutoff2 = 1.0e6;
    if (ispec->class_subtype == 1) calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
    else calculate_fm_matrix_elements = calc_dihedral_four_body_fm_matrix_elements;
}

void ThreeBodyNonbondedClassComputer::special_set_up_computer(InteractionClassSpec* const ispec_pt, int *curr_iclass_col_index)
{
    ispec = ispec_pt;
    if (ispec->class_subtype > 0) {
        if (ispec->class_subtype == 1) {
            calculate_fm_matrix_elements = calc_angular_three_body_fm_matrix_elements;
        } else if (ispec->class_subtype == 2) {
            if (ispec->get_basis_type() != kBSpline) {
                printf("Three body with fitted distance term can be only used with B-splines!\n");
                exit(EXIT_FAILURE);
            }
            calculate_fm_matrix_elements = calc_nonbonded_1_three_body_fm_matrix_elements;
        } else if (ispec->class_subtype == 3) {
            calculate_fm_matrix_elements = calc_nonbonded_2_three_body_fm_matrix_elements;
        }

        interaction_class_column_index = *curr_iclass_col_index;
        *curr_iclass_col_index += ispec->interaction_column_indices[ispec->n_to_force_match];
    }
    fm_s_comp = new BSplineAndDerivComputer(ispec);
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
    for (int l = 0; l < (int)(cg->topo_data.n_cg_sites); l++) {
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
    cg->pair_nonbonded_computer.calculate_interactions(mat, trajectory_block_frame_index, current_frame_starting_row, cg->n_cg_types, cg->topo_data, pair_cell_list, frame_config->x, frame_config->simulation_box_half_lengths);
    cg->pair_bonded_computer.calculate_interactions(mat, trajectory_block_frame_index, current_frame_starting_row, cg->n_cg_types, cg->topo_data, frame_config->x, frame_config->simulation_box_half_lengths);
    cg->angular_computer.calculate_interactions(mat, trajectory_block_frame_index, current_frame_starting_row, cg->n_cg_types, cg->topo_data, frame_config->x, frame_config->simulation_box_half_lengths);
    cg->dihedral_computer.calculate_interactions(mat, trajectory_block_frame_index, current_frame_starting_row, cg->n_cg_types, cg->topo_data, frame_config->x, frame_config->simulation_box_half_lengths);
    cg->three_body_nonbonded_computer.calculate_interactions(mat, trajectory_block_frame_index, current_frame_starting_row, cg->n_cg_types, cg->topo_data, three_body_cell_list, frame_config->x, frame_config->simulation_box_half_lengths);
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

void InteractionClassComputer::walk_neighbor_list(MATRIX_DATA* const mat, calc_pair_matrix_elements calc_matrix_elements, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, const rvec* x, const real* simulation_box_half_lengths) 
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

void PairBondedClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const rvec* x, const real* simulation_box_half_lengths) 
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

void AngularClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const rvec* x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    for (k = 0; k < int(topo_data.n_cg_sites); k++) {
        for (unsigned kk = 0; kk < topo_data.angle_list->partner_numbers_[k]; kk++) {
            l = topo_data.angle_list->partners_[k][2 * kk + 1];
            j = topo_data.angle_list->partners_[k][2 * kk];
            if (k < l) order_bonded_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
        }
    }
}

void DihedralClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const rvec* x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    for (k = 0; k < int(topo_data.n_cg_sites); k++) {
        for (unsigned kk = 0; kk < topo_data.dihedral_list->partner_numbers_[k]; kk++) {
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

void ThreeBodyNonbondedClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const ThreeBCellList& three_body_cell_list, const rvec* x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    if (ispec->class_subtype > 0) {                    
        trajectory_block_frame_index = traj_block_frame_index;
        current_frame_starting_row = curr_frame_starting_row;
        for (int kk = 0; kk < three_body_cell_list.size; kk++) {
            j = three_body_cell_list.head[kk];
            while (j >= 0) {
                k = three_body_cell_list.head[kk];
                while (k >= 0) {
                    if (j != k) {
                        //three body
                        l = three_body_cell_list.list[k];
                        while (l >= 0) {
                            if (l != j) {
                                if ( (check_excluded_list(&topo_data, l, j) == false)  && (check_excluded_list(&topo_data, j, k) == false) ) {
                                    order_three_body_nonbonded_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
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
}

//--------------------------------------------------------------------
//  Routines for checking if nonbonded interactions should be excluded
// from the model because they are between bonded particles.
//--------------------------------------------------------------------

bool check_anybond_nonbonded_exclusion(const TopologyData* const topo_data, const int i, const int j)
{
    unsigned k;
    // Check whether this nonbonded interaction is excluded for the model because the pair is bonded.
    for (k = 0; k < topo_data->bond_list->partner_numbers_[i]; k++) {
        if (topo_data->bond_list->partners_[i][k] == unsigned(j)) return true;
    }
    for (k = 0; k < topo_data->angle_list->partner_numbers_[i]; k++) {
        if (topo_data->angle_list->partners_[i][2 * k + 1] == unsigned(j)) return true;
    }
    for (k = 0; k < topo_data->dihedral_list->partner_numbers_[i]; k++) {
        if (topo_data->dihedral_list->partners_[i][3 * k + 2] == unsigned(j)) return true;
    }
    return false;
}

bool check_pairbond_nonbonded_exclusion(const TopologyData* const topo_data, const int i, const int j)
{
    // Check whether this nonbonded interaction is excluded for the model because the pair is bonded.
    for (unsigned k = 0; k < topo_data->bond_list->partner_numbers_[i]; k++) {
        if (topo_data->bond_list->partners_[i][k] == unsigned(j)) return true;
    }
    return false;
}

bool check_excluded_list(const TopologyData* const topo_data, const int i, const int j)
{
    // Check whether this non-bonded interaction is excluded for the model
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

void order_pair_nonbonded_fm_matrix_element_calculation(InteractionClassComputer* const info, calc_pair_matrix_elements calc_matrix_elements, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, const rvec *x, const real *simulation_box_half_lengths)
{
    // Calculate the appropriate matrix elements.
    info->index_among_defined_intrxns = info->ispec->get_index_from_hash(calc_two_body_interaction_hash(cg_site_types[info->k], cg_site_types[info->l], n_cg_types));
    info->index_among_matched_interactions = info->ispec->defined_to_matched_intrxn_index_map[info->index_among_defined_intrxns];
    info->index_among_tabulated_interactions = info->ispec->defined_to_tabulated_intrxn_index_map[info->index_among_defined_intrxns];
    if ((info->index_among_matched_interactions == 0) && (info->index_among_tabulated_interactions == 0)) return; // if the index is zero, it is not present in the model and should be ignored.
    calc_matrix_elements(info, x, simulation_box_half_lengths, mat);
}

void order_bonded_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, const rvec *x, const real *simulation_box_half_lengths)
{
     // Calculate the appropriate matrix elements.    
    info->index_among_defined_intrxns = info->ispec->get_index_from_hash(info->calculate_hash_number(cg_site_types, n_cg_types));
    info->index_among_matched_interactions = info->ispec->defined_to_matched_intrxn_index_map[info->index_among_defined_intrxns];
    info->index_among_tabulated_interactions = info->ispec->defined_to_tabulated_intrxn_index_map[info->index_among_defined_intrxns];
    if ((info->index_among_matched_interactions == 0) && (info->index_among_tabulated_interactions == 0)) return; // if the index is zero, it is not present in the model and should be ignored.
    (*info->calculate_fm_matrix_elements)(info, x, simulation_box_half_lengths, mat);
}

void order_three_body_nonbonded_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, const rvec *x, const real *simulation_box_half_lengths)
{
    ThreeBodyNonbondedClassComputer* icomp = static_cast<ThreeBodyNonbondedClassComputer*>(info);
    ThreeBodyNonbondedClassSpec* ispec = static_cast<ThreeBodyNonbondedClassSpec*>(icomp->ispec);
    
    // Calculate the appropriate matrix elements.
    icomp->index_among_defined_intrxns = info->ispec->get_index_from_hash(icomp->calculate_hash_number(cg_site_types, n_cg_types));
    if (icomp->index_among_defined_intrxns == -1) return; // if the index is -1, it is not present in the model and should be ignored.
    
    icomp->index_among_matched_interactions = ispec->defined_to_matched_intrxn_index_map[icomp->index_among_defined_intrxns];
    icomp->index_among_tabulated_interactions = ispec->defined_to_tabulated_intrxn_index_map[icomp->index_among_defined_intrxns];
    if ((icomp->index_among_matched_interactions == 0) && (icomp->index_among_tabulated_interactions == 0)) return; // if the index is zero, it is not present in the model and should be ignored.
    
    icomp->cutoff2 = ispec->three_body_nonbonded_cutoffs[icomp->index_among_defined_intrxns];
    icomp->stillinger_weber_angle_parameter = ispec->stillinger_weber_angle_parameters_by_type[icomp->index_among_defined_intrxns];
    (*icomp->calculate_fm_matrix_elements)(icomp, x, simulation_box_half_lengths, mat); 
}

//---------------------------------------------------------------------
// Functions for adding forces on a template number of particles into
// the target vector of the matrix from their ids, derivatives, and
// a set of spline coefficients.
//---------------------------------------------------------------------

template<int n_body> void accumulate_table_forces(InteractionClassComputer* const info, const double &table_fn_val, const std::array<int, n_body> &particle_ids, const std::array<std::array<double, 3>, n_body - 1> &derivatives, MATRIX_DATA * const mat) 
{
    // Calculate the associated forces.
    // Use flat arrays for performance.
    std::array<double, 3 * n_body> forces;
    for (int j = 0; j < 3; j++) forces[3*(n_body - 1) + j] = 0.0;
    for (int i = 0; i < n_body - 1; i++) {
        for (int j = 0; j < 3; j++) {
            forces[3 * i + j] = table_fn_val * derivatives[i][j];
            forces[3 * (n_body - 1) + j] += -table_fn_val * derivatives[i][j];
        }
    }
    // Load those forces into the target vector.
    for (int i = 0; i < n_body; i++) {
        mat->accumulate_target_force_element(mat, particle_ids[i] + info->current_frame_starting_row, &forces[3 * i]);
    }
}

template<int n_body> void accumulate_matching_forces(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const std::array<int, n_body> &particle_ids, const std::array<std::array<double, 3>, n_body - 1> &derivatives, MATRIX_DATA * const mat) 
{
    // For each basis function,
    int ref_column = info->interaction_class_column_index + info->ispec->interaction_column_indices[info->index_among_matched_interactions - 1] + first_nonzero_basis_index;

    std::array<double, 3 * n_body> forces;
    for (unsigned k = 0; k < basis_fn_vals.size(); k++) {
        // Calculate the associated forces.
        // Use flat force array for performance.
        for (int j = 0; j < 3; j++) forces[3 * (n_body - 1) + j] = 0.0;
        for (int i = 0; i < n_body - 1; i++) {
            for (int j = 0; j < 3; j++) {
                forces[3 * i + j] = -basis_fn_vals[k] * derivatives[i][j];
                forces[3 * (n_body - 1) + j] += basis_fn_vals[k] * derivatives[i][j];
            }
        }
        // Load those forces into the target vector.
        for (int i = 0; i < n_body; i++) {
            (*mat->accumulate_fm_matrix_element)(particle_ids[i] + info->current_frame_starting_row, ref_column + k, &forces[3 * i], mat);
        }
    }
}

//--------------------------------------------------------------------
// Functions for calculating sets of 3-component matrix elements for each
// individual interacting set of particles
//--------------------------------------------------------------------

// Each of these functions follows the idiom of the first.

void calc_isotropic_two_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    std::array<int, 2> particle_ids = {{info->k, info->l}};
    int index_among_defined = info->index_among_defined_intrxns;
    int index_among_tabulated = info->index_among_tabulated_interactions;
    int index_among_matched = info->index_among_matched_interactions;
    double distance;
    std::array<std::array<double, 3>, 1> derivatives;
    bool within_cutoff = conditionally_calc_distance_and_derivatives(particle_ids, x, simulation_box_half_lengths, info->cutoff2, distance, derivatives);

    if (within_cutoff) {
        InteractionClassSpec *ispec = info->ispec;

        if (distance < ispec->lower_cutoffs[index_among_defined]) return;
        else if (distance > ispec->upper_cutoffs[index_among_defined]) return;

        if (index_among_tabulated > 0) {
            // Pull the interaction from a table.
            int first_nonzero_basis_index;
            info->table_s_comp->calculate_basis_fn_vals(index_among_defined, distance, first_nonzero_basis_index, info->table_basis_fn_vals);
            double basis_sum = info->table_basis_fn_vals[0] + info->table_basis_fn_vals[1];
            // Add to force target.
            accumulate_table_forces<2>(info, basis_sum, particle_ids, derivatives, mat);
            // Add to virial target.
            if (mat->virial_constraint_rows > 0) mat->accumulate_target_constraint_element(mat, info->trajectory_block_frame_index, -basis_sum * distance);
        }
        
        if (index_among_matched > 0) {
            // Compute the strength of each basis function.
            int first_nonzero_basis_index;
            info->fm_s_comp->calculate_basis_fn_vals(index_among_defined, distance, first_nonzero_basis_index, info->fm_basis_fn_vals);
            // Add to the force matching.
            accumulate_matching_forces<2>(info, first_nonzero_basis_index, info->fm_basis_fn_vals, particle_ids, derivatives, mat);
            // Add to virial matching.
            int temp_column_index = info->interaction_class_column_index + ispec->interaction_column_indices[index_among_matched - 1] + first_nonzero_basis_index;
            for (unsigned i = 0; i < info->fm_basis_fn_vals.size(); i++) {
                int basis_column = temp_column_index + i;
                if (mat->virial_constraint_rows > 0)(*mat->accumulate_virial_constraint_matrix_element)(info->trajectory_block_frame_index, basis_column, info->fm_basis_fn_vals[i] * distance, mat);
            }
        }
    } else {
        return;
    }
}

void calc_angular_three_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    std::array<int, 3> particle_ids = {{info->j, info->k, info->l}};
    int index_among_defined = info->index_among_defined_intrxns;
    int index_among_tabulated = info->index_among_tabulated_interactions;
    int index_among_matched = info->index_among_matched_interactions;
    double angle;
    std::array<std::array<double, 3>, 2> derivatives;
    bool within_cutoff = conditionally_calc_angle_and_derivatives(particle_ids, x, simulation_box_half_lengths, info->cutoff2, angle, derivatives);

    if (within_cutoff) {
        InteractionClassSpec *ispec = info->ispec;

        if (angle < ispec->lower_cutoffs[index_among_defined]) return;
        else if (angle > ispec->upper_cutoffs[index_among_defined]) return;

        if (index_among_tabulated > 0) {
            // Pull the interaction from a table.
            int first_nonzero_basis_index;
            info->table_s_comp->calculate_basis_fn_vals(index_among_defined, angle, first_nonzero_basis_index, info->table_basis_fn_vals);
            double basis_sum = info->table_basis_fn_vals[0] + info->table_basis_fn_vals[1];
            // Add to force target.
            accumulate_table_forces<3>(info, basis_sum, particle_ids, derivatives, mat);
            // Angles do not contribute to the scalar virial.
        }
        if (index_among_matched > 0) {
            // Compute the strength of each basis function.
            int first_nonzero_basis_index;
            info->fm_s_comp->calculate_basis_fn_vals(index_among_defined, angle, first_nonzero_basis_index, info->fm_basis_fn_vals);
            // Add to the force matching.
            accumulate_matching_forces<3>(info, first_nonzero_basis_index, info->fm_basis_fn_vals, particle_ids, derivatives, mat);
            // Angles do not contribute to the scalar virial.
        }
    } else {
        return;
    }
}

void calc_dihedral_four_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    std::array<int, 4> particle_ids = {{info->i, info->j, info->k, info->l}};
    int index_among_defined = info->index_among_defined_intrxns;
    int index_among_tabulated = info->index_among_tabulated_interactions;
    int index_among_matched = info->index_among_matched_interactions;
    double dihedral;
    std::array<std::array<double, 3>, 3> derivatives;
    bool within_cutoff = conditionally_calc_dihedral_and_derivatives(particle_ids, x, simulation_box_half_lengths, info->cutoff2, dihedral, derivatives);

    if (within_cutoff) {
        InteractionClassSpec *ispec = info->ispec;

        if (dihedral < ispec->lower_cutoffs[index_among_defined]) return;
        else if (dihedral > ispec->upper_cutoffs[index_among_defined]) return;

        if (index_among_tabulated > 0) {
            // Pull the interaction from a table.
            int first_nonzero_basis_index;
            info->table_s_comp->calculate_basis_fn_vals(index_among_defined, dihedral, first_nonzero_basis_index, info->table_basis_fn_vals);
            double basis_sum = info->table_basis_fn_vals[0] + info->table_basis_fn_vals[1];
            // Add to force target.
            accumulate_table_forces<4>(info, basis_sum, particle_ids, derivatives, mat);
            // Dihedrals do not contribute to the scalar virial.
        }
        if (index_among_matched > 0) {
            // Compute the strength of each basis function.
            int first_nonzero_basis_index;
            info->fm_s_comp->calculate_basis_fn_vals(index_among_defined, dihedral, first_nonzero_basis_index, info->fm_basis_fn_vals);
            // Add to the force matching.
            accumulate_matching_forces<4>(info, first_nonzero_basis_index, info->fm_basis_fn_vals, particle_ids, derivatives, mat);
            // Dihedrals do not contribute to the scalar virial.
        }
    } else {
        return;
    }
}

void calc_nonbonded_1_three_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    std::array<int, 3> particle_ids = {{info->j, info->k, info->l}};    
    ThreeBodyNonbondedClassComputer* icomp = static_cast<ThreeBodyNonbondedClassComputer*>(info);
    ThreeBodyNonbondedClassSpec* ispec = static_cast<ThreeBodyNonbondedClassSpec*>(icomp->ispec);
    unsigned i;
	int	j, tn;
    double relative_site_position_2[3], relative_site_position_3[3], rr1, rr2;
    double cos_theta;
    double theta;
    double ty1[3], ty2[3], ty[3];
    double rr_12_1, rr_11c, rr_22c;
    double sin_theta;
    double tx[3], tx1[3], tx2[3];
    double s1, s2, s1_1, s2_1, rt1, rt2;
    double c1, c2, c3;
    double tt;
    int temp_column_index;
    
    cos_theta = calc_cosine_of_angle_and_intermediates(particle_ids[0], particle_ids[1], particle_ids[2], x, simulation_box_half_lengths, relative_site_position_2, relative_site_position_3, &rr1, &rr2);
    
    if (rr1 > icomp->cutoff2 || rr2 > icomp->cutoff2) return;
    
    theta = acos(cos_theta);
    icomp->intrxn_param = theta * DEGREES_PER_RADIAN;
    
    sin_theta = sin(theta);
    rr_12_1 = 1.0 / (rr1 * rr2 * sin_theta);
    rr_11c = cos_theta / (rr1 * rr1 * sin_theta);
    rr_22c = cos_theta / (rr2 * rr2 * sin_theta);
    for (i = 0; i < 3; i++) {
        ty1[i] = -relative_site_position_3[i] * rr_12_1 + rr_11c * relative_site_position_2[i];
        ty2[i] = -relative_site_position_2[i] * rr_12_1 + rr_22c * relative_site_position_3[i];
        ty[i] = -ty1[i] - ty2[i];
    }
    
    // Calculate the matrix elements if it's supposed to be force matched
    info->fm_s_comp->calculate_basis_fn_vals(info->index_among_defined_intrxns, info->intrxn_param, info->basis_function_column_index, info->fm_basis_fn_vals); 
    std::vector<double> basis_der_vals(info->fm_s_comp->get_n_coef());
    BSplineAndDerivComputer *fm_s_comp = static_cast<BSplineAndDerivComputer*>(icomp->fm_s_comp);
    fm_s_comp->calculate_bspline_deriv_vals(info->index_among_defined_intrxns, info->intrxn_param, info->basis_function_column_index, basis_der_vals); 
    
    rt1 = rr1 - icomp->cutoff2;
    rt2 = rr2 - icomp->cutoff2;
    s1 = exp(ispec->three_body_gamma / rt1);
    s2 = exp(ispec->three_body_gamma / rt2);
    s1_1 = ispec->three_body_gamma / (rt1 * rt1) * s1;
    s2_1 = ispec->three_body_gamma / (rt2 * rt2) * s2;
    c1 = s1 * s2;
    c2 = s2 * s1_1;
    c3 = s1 * s2_1;
    
    c1 *= DEGREES_PER_RADIAN;
    
    for (i = 0; i < info->fm_basis_fn_vals.size(); i++) {
        
        temp_column_index = icomp->interaction_class_column_index + ispec->interaction_column_indices[icomp->index_among_matched_interactions - 1] + icomp->basis_function_column_index;
        tn = temp_column_index + i;
        
        tt = basis_der_vals[i] * c1;
        for (j = 0; j < 3; j++) tx1[j] = ty1[j] * tt;
        for (j = 0; j < 3; j++) tx2[j] = ty2[j] * tt;
        
        tt = info->fm_basis_fn_vals[i] * c2;
        for (j = 0; j < 3; j++) tx1[j] += relative_site_position_2[j] / rr1 * tt;
        (*mat->accumulate_fm_matrix_element)(particle_ids[1] + icomp->current_frame_starting_row, tn, tx1, mat);
        tt = info->fm_basis_fn_vals[i] * c3;
        for (j = 0; j < 3; j++) tx2[j] += relative_site_position_3[j] / rr2 * tt;
        (*mat->accumulate_fm_matrix_element)(particle_ids[2] + icomp->current_frame_starting_row, tn, tx2, mat);
        for (j = 0; j < 3; j++) tx[j] = -(tx1[j] + tx2[j]);
        (*mat->accumulate_fm_matrix_element)(particle_ids[0] + icomp->current_frame_starting_row, tn, tx, mat);
    }
}

void calc_nonbonded_2_three_body_fm_matrix_elements(InteractionClassComputer* const info, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    std::array<int, 3> particle_ids = {{info->j, info->k, info->l}};    
    ThreeBodyNonbondedClassComputer* icomp = static_cast<ThreeBodyNonbondedClassComputer*>(info);
    ThreeBodyNonbondedClassSpec* ispec = static_cast<ThreeBodyNonbondedClassSpec*>(icomp->ispec);
    
    double relative_site_position_2[3], relative_site_position_3[3], rr1, rr2;
    double cos_theta;
    double theta;
    int i, j, tn;
    double ty1[3], ty2[3], ty[3];
    double rr_12_1, rr_11c, rr_22c;
    double sin_theta;
    double tx[3], tx1[3], tx2[3];
    double s1, s2, s1_1, s2_1, rt1, rt2;
    double c1, c2, c3;
    double u, u_1;
    double tt;
    
    cos_theta = calc_cosine_of_angle_and_intermediates(particle_ids[0], particle_ids[1], particle_ids[2], x, simulation_box_half_lengths, relative_site_position_2, relative_site_position_3, &rr1, &rr2);
    if (rr1 > icomp->cutoff2 || rr2 > icomp->cutoff2) return;
    theta = acos(cos_theta);
    sin_theta = sin(theta);
    icomp->intrxn_param = theta * DEGREES_PER_RADIAN;  
    int temp_column_index;
     
    rr_12_1 = 1.0 / (rr1 * rr2 * sin_theta);
    rr_11c = cos_theta / (rr1 * rr1 * sin_theta);
    rr_22c = cos_theta / (rr2 * rr2 * sin_theta);
    for (i = 0; i < 3; i++) {
        ty1[i] = -relative_site_position_3[i] * rr_12_1 + rr_11c * relative_site_position_2[i];
        ty2[i] = -relative_site_position_2[i] * rr_12_1 + rr_22c * relative_site_position_3[i];
        ty[i] = -ty1[i] - ty2[i];
    }
    
    rt1 = rr1 - icomp->cutoff2;
    rt2 = rr2 - icomp->cutoff2;
    s1 = exp(ispec->three_body_gamma / rt1);
    s2 = exp(ispec->three_body_gamma / rt2);
    s1_1 = ispec->three_body_gamma / (rt1 * rt1) * s1;
    s2_1 = ispec->three_body_gamma / (rt2 * rt2) * s2;
    
    c1 = s1 * s2;
    c2 = s2 * s1_1;
    c3 = s1 * s2_1;
    
    u = (cos_theta - icomp->stillinger_weber_angle_parameter) * (cos_theta - icomp->stillinger_weber_angle_parameter) * 4.184;
    u_1 = 2.0 * (cos_theta - icomp->stillinger_weber_angle_parameter) * sin_theta * 4.184;
    
    temp_column_index = icomp->interaction_class_column_index + ispec->interaction_column_indices[icomp->index_among_matched_interactions - 1];
    
    tn = temp_column_index;
    
    tt = c1 * u_1;
    for (j = 0; j < 3; j++) tx1[j] = ty1[j] * tt;
    for (j = 0; j < 3; j++) tx2[j] = ty2[j] * tt;
    
    tt = c2 * u;
    for (j = 0; j < 3; j++) tx1[j] += relative_site_position_2[j] / rr1 * tt;
    (*mat->accumulate_fm_matrix_element)(particle_ids[1] + icomp->current_frame_starting_row, tn, tx1, mat);
    tt = c3 * u;
    for (j = 0; j < 3; j++) tx2[j] += relative_site_position_3[j] / rr2 * tt;
    (*mat->accumulate_fm_matrix_element)(particle_ids[2] + icomp->current_frame_starting_row, tn, tx2, mat);
    for (j = 0; j < 3; j++) tx[j] = -(tx1[j] + tx2[j]);
    (*mat->accumulate_fm_matrix_element)(particle_ids[0] + icomp->current_frame_starting_row, tn, tx, mat);    
}

//--------------------------------------------------------------------
// Utility functions for calculating the parameters that the potential
// basis functions depend on, such as distances, angles, and dihedrals
//--------------------------------------------------------------------

void get_minimum_image(const int l, rvec* frx, const real *simulation_box_half_lengths)
{
    for (int i = 0; i < 3; i++) {
        if (frx[l][i] < 0) frx[l][i] += 2.0 * simulation_box_half_lengths[i];
        else if (frx[l][i] >= 2.0 * simulation_box_half_lengths[i]) frx[l][i] -= 2.0 * simulation_box_half_lengths[i];
    }
}

double calc_squared_distance(const int k, const int l, const rvec *x, const real *simulation_box_half_lengths, double* const relative_site_position_1)
{
    double rr2 = 0.0;
    for (int i = 0; i < 3; i++) {
        relative_site_position_1[i] = x[l][i] - x[k][i];
        if (relative_site_position_1[i] > simulation_box_half_lengths[i]) relative_site_position_1[i] -= 2.0 * simulation_box_half_lengths[i];
        else if (relative_site_position_1[i] < -simulation_box_half_lengths[i]) relative_site_position_1[i] += 2.0 * simulation_box_half_lengths[i];
        
        rr2 += relative_site_position_1[i] * relative_site_position_1[i];
    }
    return rr2;
}

double dot_product(const double* a, const double* b)
{
    double t = 0.0;
    for (int i = 0; i < 3; i++) t += a[i] * b[i];
    return t;
}

double calc_cosine_of_angle_and_intermediates(const int j, const int k, const int l, const rvec *x, const real *simulation_box_half_lengths, double* const relative_site_position_2, double* const relative_site_position_3, double* const rr1, double* const rr2)
{
    *rr1 = sqrt(calc_squared_distance(j, k, x, simulation_box_half_lengths, relative_site_position_2));
    *rr2 = sqrt(calc_squared_distance(j, l, x, simulation_box_half_lengths, relative_site_position_3));
    
    double tx = (relative_site_position_2[0] * relative_site_position_3[0] + relative_site_position_2[1] * relative_site_position_3[1] + relative_site_position_2[2] * relative_site_position_3[2]) / ((*rr1) * (*rr2));
    double max = 1.0 - VERYSMALL_F;
    double min = -1.0 + VERYSMALL_F;
    if (tx > max) tx = max;
    else if (tx < min) tx = min;
    return tx;
}

//--------------------------------------------------------------------
// Free the temps used in FM matrix building, retaining only what is 
// still needed for solution and output.
//--------------------------------------------------------------------

void free_fm_matrix_building_temps(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FrameSource* const frame_source)
{    
    // Free topology information
    printf("Freeing topology information.\n");
    free_topology_data(&cg->topo_data);
    
    // Free tabulated potential data
    printf("Freeing tabulated reference potential information.\n");
	std::list<InteractionClassComputer*>::iterator icomp_iterator;
    for(icomp_iterator=cg->icomp_list.begin(); icomp_iterator != cg->icomp_list.end(); icomp_iterator++) {
        delete (*icomp_iterator)->table_s_comp;
        (*icomp_iterator)->ispec->free_force_tabulated_interaction_data();
    }    
    // Free FM matrix building temps
    printf("Freeing equation building temporaries.\n");

    if (mat->matrix_type == kDense) {
        if (mat->virial_constraint_rows > 0) delete [] mat->dense_fm_matrix;
        delete [] mat->dense_fm_rhs_vector;
    } else if (mat->matrix_type == kSparse) {
        delete [] mat->ll_sparse_matrix_row_heads;
        delete [] mat->block_fm_solution;
        delete [] mat->dense_fm_rhs_vector;
        if (mat->virial_constraint_rows > 0) delete [] mat->dense_fm_matrix;
    } else if (mat->matrix_type == kAccumulation) {
        delete [] mat->lapack_temp_workspace;
        delete [] mat->lapack_tau;
    } else if (mat->matrix_type == kSparseNormal) {
        delete [] mat->ll_sparse_matrix_row_heads;
        delete [] mat->dense_fm_rhs_vector;
        if (mat->virial_constraint_rows > 0) delete [] mat->dense_fm_matrix;
    } else if (mat->matrix_type == kSparseSparse) {
        delete [] mat->ll_sparse_matrix_row_heads;
        delete [] mat->dense_fm_rhs_vector;
        if (mat->virial_constraint_rows > 0) delete [] mat->dense_fm_matrix;
    }
}

void InteractionClassSpec::free_force_tabulated_interaction_data() 
{
    if (n_tabulated > 0) {
        for (int i = 0; i < n_tabulated; i++) {
            delete [] external_table_spline_coefficients[i];
        }
        delete [] external_table_spline_coefficients;
    }
}
