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

//#include "range_finding.cpp"

//--------------------------------------------------------------------
// Prototypes for internal implementation-specific functions
//--------------------------------------------------------------------

// Utility functions for checking if a nonbonded interaction is excluded from the model due to bonding.

bool check_excluded_list(const TopologyData* const topo_data, const int i, const int j);
bool check_density_excluded_list(const TopologyData* const topo_data, const int i, const int j);

// Main routine responsible for calling single-element matrix computations,
// differing by the way that potentially interacting particles are found in 
// each frame and possibly found not to interact after.

void order_one_body_fm_matrix_element_calculation(InteractionClassComputer* const info, calc_pair_matrix_elements calc_matrix_elements, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const & x, const real *simulation_box_half_lengths);
void order_pair_nonbonded_fm_matrix_element_calculation(InteractionClassComputer* const info, calc_pair_matrix_elements calc_matrix_elements, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths);
void order_bonded_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths);
void order_radius_of_gyration_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths);    
void order_helical_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths);
void order_three_body_nonbonded_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths);
void density_fm_matrix_element_calculation(InteractionClassComputer* const iclass, calc_pair_matrix_elements calc_matrix_elements, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths);

// Helper functions for the above

void process_completed_density(DensityClassComputer* const info, calc_pair_matrix_elements process_density, const int n_cg_types, int* const cg_site_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths);
inline void decode_density_interaction_and_calculate(DensityClassComputer* info, unsigned long interaction_flags, calc_pair_matrix_elements calc_matrix_elements, int* const cg_site_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths);
void process_normal_interaction_matrix_elements(InteractionClassComputer* const info, MATRIX_DATA* const mat, const int n_body, int* particle_ids, std::array<double, DIMENSION>* derivatives, const double param_value, const int virial_flag, const double param_deriv, const double distance);
void process_one_body_matrix_elements(InteractionClassComputer* const info, MATRIX_DATA* const mat, const int n_body, int* particle_ids, std::array<double, DIMENSION>* derivatives, const double param_value, const int virial_flag, const double junk, const double junk2);
void process_density_matrix_elements(InteractionClassComputer* const info, MATRIX_DATA* const mat, const int n_body, int* particle_ids, std::array<double, DIMENSION>* derivatives, const double density_value, const int virial_flag, const double density_derivative, const double distance);

// Functions for calculating individual 3-component matrix elements.

void calc_one_body_fm_matrix_element(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_isotropic_two_body_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_angular_three_body_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_dihedral_four_body_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_radius_of_gyration_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_helical_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_density_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_nonbonded_1_three_body_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_nonbonded_2_three_body_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
double calc_gaussian_density_derivative(DensityClassComputer* const icomp, DensityClassSpec* const ispec, const double distance);
double calc_switching_density_derivative(DensityClassComputer* const icomp, DensityClassSpec* const ispec, const double distance);
double calc_lucy_density_derivative(DensityClassComputer* const icomp, DensityClassSpec* const ispec, const double distance);
double calc_re_density_derivative(DensityClassComputer* const icomp, DensityClassSpec* const ispec, const double distance);
void do_nothing(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void accumulate_matching_order_parameter_forces(InteractionClassComputer* const info, const int first_nonzero_basis_index, double extra_derivative_value, std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);

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
    cutoff2 = ispec->cutoff * ispec->cutoff;
    class_set_up_computer();
}

void OneBodyClassComputer::class_set_up_computer(void) 
{
    cutoff2 = 1.0;
    calculate_fm_matrix_elements = calc_one_body_fm_matrix_element;
    process_interaction_matrix_elements = process_one_body_matrix_elements;
}

void PairNonbondedClassComputer::class_set_up_computer(void) 
{
	calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
}

void PairBondedClassComputer::class_set_up_computer(void) 
{
    calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
}

void AngularClassComputer::class_set_up_computer(void) 
{
    if (ispec->class_subtype == 1) calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
    else calculate_fm_matrix_elements = calc_angular_three_body_fm_matrix_elements;
}

void DihedralClassComputer::class_set_up_computer(void) 
{
    if (ispec->class_subtype == 1) calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
    else calculate_fm_matrix_elements = calc_dihedral_four_body_fm_matrix_elements;
}

void R13ClassComputer::class_set_up_computer(void) 
{
    calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
}

void R14ClassComputer::class_set_up_computer(void) 
{
    calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
}

void R15ClassComputer::class_set_up_computer(void) 
{
    calculate_fm_matrix_elements = calc_isotropic_two_body_fm_matrix_elements;
}

void HelicalClassComputer::class_set_up_computer(void) 
{
    calculate_fm_matrix_elements = calc_helical_fm_matrix_elements;
}

void RadiusofGyrationClassComputer::class_set_up_computer(void) 
{
    calculate_fm_matrix_elements = calc_radius_of_gyration_fm_matrix_elements;
}

void DensityClassComputer::class_set_up_computer(void) 
{
	DensityClassSpec* iclass = static_cast<DensityClassSpec*>(ispec);
	process_interaction_matrix_elements = process_density_matrix_elements;
	
	if(iclass->class_subtype == 0) return;
	if(iclass->class_subtype == 1) {
		printf("Will calculate density using shifted-force Gaussian weight functions.\n");
		calculate_density_values = calc_gaussian_density_values;
		calculate_density_derivative = calc_gaussian_density_derivative;
	} else if(iclass->class_subtype == 2) {
		printf("Will calculate density using shifted-force switching (tanh) weight functions.\n");
		calculate_density_values = calc_switching_density_values;
		calculate_density_derivative = calc_switching_density_derivative;
	} else if(iclass->class_subtype == 3) {
		printf("Will calculate density using Lucy-style weight functions.\n");
		calculate_density_values = calc_lucy_density_values;
		calculate_density_derivative = calc_lucy_density_derivative;
	}  else if(iclass->class_subtype == 4) {
		printf("Will calculate density using Relative Entropy-style weight functions.\n");
		calculate_density_values = calc_re_density_values;
		calculate_density_derivative = calc_re_density_derivative;
	}
	calculate_fm_matrix_elements = calc_density_fm_matrix_elements;
	process_density = do_nothing;
	
	// Allocate space to store density intermediate.
	// This approach grabs enough memory for all cg sites to have all density values for all density groups.
	// A more memory-efficient, but harder approach would be allocate each site's array based on the number of density groups needed at that site.
	// The problem with this other approach is efficiently looking up which i to use.
	density_values = new double[iclass->get_n_defined() * iclass->n_cg_sites]();
	
	// Allocate and compute constant calculation intermediates.
	denomenator = new double[iclass->get_n_defined()];
	u_cutoff = new double[iclass->get_n_defined()]();
	f_cutoff = new double[iclass->get_n_defined()]();
	
	if(iclass->class_subtype == 1) {
		for(int ii = 0; ii < iclass->get_n_defined(); ii++) {
			if (iclass->density_sigma[ii] < VERYSMALL_F) {
				printf("Density sigma parameter is too small!\n");
				exit(EXIT_FAILURE);
			}
			denomenator[ii] = 2.0 * iclass->density_sigma[ii] * iclass->density_sigma[ii];
			u_cutoff[ii] = - exp( - cutoff2 / denomenator[ii] );
			f_cutoff[ii] = - 2.0 * iclass->cutoff * u_cutoff[ii];
		}
	} else if (iclass->class_subtype == 2) {
		for(int ii = 0; ii < iclass->get_n_defined(); ii++) {
			if (iclass->density_sigma[i] < VERYSMALL_F) {
				printf("Density sigma parameter is too small!\n");
				exit(EXIT_FAILURE);
			}
			denomenator[ii] = iclass->density_sigma[ii] / 0.5;
			double arguement = (iclass->cutoff - iclass->density_switch[ii]) / iclass->density_sigma[ii];
			u_cutoff[ii] = 0.5 * tanh( arguement );
			f_cutoff[ii] =  0.5 / (iclass->density_sigma[ii] * cosh(arguement) * cosh(arguement));
			printf("%d: density_switch %lf, density_sigma %lf, cutoff %lf, u_cutoff %lf, f_cutoff %lf, denom %lf\n", ii, iclass->density_sigma[ii], iclass->density_switch[ii], iclass->cutoff, u_cutoff[ii], f_cutoff[ii], denomenator[ii]); fflush(stdout);
		}
	} else if (iclass->class_subtype == 3) {
		for(int ii = 0; ii < iclass->get_n_defined(); ii++) {
			denomenator[ii] = pow(iclass->cutoff, 4.0);
		} 
	} else if (iclass->class_subtype == 4) {
		c0 = new double[iclass->get_n_defined()]();
		c2 = new double[iclass->get_n_defined()]();
		c4 = new double[iclass->get_n_defined()]();
		c6 = new double[iclass->get_n_defined()]();
		double cutsq = iclass->cutoff * iclass->cutoff;
		for(int ii = 0; ii < iclass->get_n_defined(); ii++) {
			double x = iclass->density_sigma[ii] * iclass->density_sigma[ii] / (iclass->cutoff * iclass->cutoff);
			denomenator[ii] = (1 - x)*(1 - x)*(1 - x);
			
			c0[ii] = (1.0 - 3.0 * x)/ denomenator[ii];
			c2[ii] = 6.0 * x / (cutsq * denomenator[ii]);
			c4[ii] = 3 * (1.0 + x) / (cutsq * cutsq * denomenator[ii]);
			c6[ii] = 2.0 / (cutsq * cutsq * cutsq * denomenator[ii]);
		}
	}else {
		printf("Class_set_up_computer called for density_interactions with invalid class_subtype %d.\n", iclass->class_subtype);
		fflush(stdout);
		exit(EXIT_FAILURE);
	}
}

void DensityClassComputer::reset_density_array(void) 
{
	DensityClassSpec* iclass = static_cast<DensityClassSpec*>(ispec);
	for(int ii = 0; ii < iclass->get_n_defined() * iclass->n_cg_sites; ii++) {
		density_values[ii] = 0.0;
	}
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
   
    // Check if molecule_list should be updated.
    if (cg->molecule_flag == 1) {
    	update_molecule_list(&cg->topo_data, cg->topo_data.molecule_list);
    }
    
    // Check if molecule list should be updated.
  	if (cg->helical_interactions.class_subtype == 1) {
    	cg->helical_interactions.rebuild_helical_list(cg->topo_data.molecule_list, cg->topo_data.dihedral_list);
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

void OneBodyClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
{
	if (ispec->class_subtype == 0) return;
	trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    for (k = 0; k < (int)(topo_data.n_cg_sites); k++) {
    	order_one_body_fm_matrix_element_calculation(this, calculate_fm_matrix_elements, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
    }
}

void PairNonbondedClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    walk_neighbor_list(mat, calculate_fm_matrix_elements, n_cg_types, topo_data, pair_cell_list, x, simulation_box_half_lengths);
}

inline void InteractionClassComputer::walk_neighbor_list(MATRIX_DATA* const mat, calc_pair_matrix_elements calc_matrix_elements, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    int stencil_size = pair_cell_list.get_stencil_size();
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
            for (int nei = 0; nei < stencil_size; nei++) {
                int ll = pair_cell_list.stencil[stencil_size * kk + nei];
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

inline void DensityClassComputer::walk_density_neighbor_list(MATRIX_DATA* const mat, calc_pair_matrix_elements calc_matrix_elements, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    int stencil_size = pair_cell_list.get_stencil_size();
    for (int kk = 0; kk < pair_cell_list.size; kk++) {
        k = pair_cell_list.head[kk];
        while (k >= 0) {
            l = pair_cell_list.list[k];
            while (l >= 0) {
                if (check_density_excluded_list(&topo_data, k, l) == false) {
                    density_fm_matrix_element_calculation(this, calc_matrix_elements, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
                }
                l = pair_cell_list.list[l];
            }
            //do the above the 2nd time for neiboring cells
            for (int nei = 0; nei < stencil_size; nei++) {
                int ll = pair_cell_list.stencil[stencil_size * kk + nei];
                l = pair_cell_list.head[ll];
                while (l >= 0) {
                    if (check_density_excluded_list(&topo_data, k, l) == false) {
                        density_fm_matrix_element_calculation(this, calc_matrix_elements, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
                    }
                    l = pair_cell_list.list[l];
                }
            }
            k = pair_cell_list.list[k];
        }
    }
}

// Calculate matrix elements for all bonded interactions by looping over the approriate topology lists. 

void PairBondedClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
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

void AngularClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
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

void DihedralClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    
    if (DIMENSION != 3) {
    	printf("Dihedral calculations are currently only implemented for 3-dimensional systems.\n");
    	exit(EXIT_FAILURE);
    }

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

void R13ClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
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

void R14ClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    
    if (DIMENSION != 3) {
    	printf("Dihedral calculations are currently only implemented for 3-dimensional systems.\n");
    	exit(EXIT_FAILURE);
    }

    trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    for (k = 0; k < int(topo_data.n_cg_sites); k++) {
        for (unsigned kk = 0; kk < topo_data.dihedral_list->partner_numbers_[k]; kk++) {
        	// Grab partners from dihedral list (organization of dihedral_list described in topology files).
        	// i and j are the indices for the "central bond" index while l and k are the "ends" of the dihedral.
        	// To avoid double counting, the interaction is only counted if the ends are
        	// ordered such that k < l.
            l = topo_data.dihedral_list->partners_[k][3 * kk + 2];
            j = topo_data.dihedral_list->partners_[k][3 * kk];
            i = topo_data.dihedral_list->partners_[k][3 * kk + 1];
            if (k < l) order_bonded_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
        }
    }
}

void R15ClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    
    if (DIMENSION != 3) {
    	printf("Dihedral calculations are currently only implemented for 3-dimensional systems.\n");
    	exit(EXIT_FAILURE);
    }

    trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    for (k = 0; k < int(topo_data.n_cg_sites); k++) {
        for (unsigned kk = 0; kk < topo_data.quint_list->partner_numbers_[k]; kk++) {
        	// Grab partners from dihedral list (organization of dihedral_list described in topology files).
        	// i and j are the indices for the "central bond" index while l and k are the "ends" of the dihedral.
        	// To avoid double counting, the interaction is only counted if the ends are
        	// ordered such that k < l.
            l = topo_data.quint_list->partners_[k][4 * kk + 3];
            h = topo_data.quint_list->partners_[k][4 * kk];
            i = topo_data.quint_list->partners_[k][4 * kk + 1];
            j = topo_data.quint_list->partners_[k][4 * kk + 2];
            if (k < l) order_bonded_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
        }
    }
}

void HelicalClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0 || ispec->class_subtype == 0) return;
    
    trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    HelicalClassSpec* h_spec = static_cast<HelicalClassSpec*>(ispec);
    
    // Look through each molecule that could have active interactions
    for (k = 0; k < (int)(h_spec->topo_data_->molecule_list->n_sites_); k++) {
    	order_helical_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
    }
}

void RadiusofGyrationClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    RadiusofGyrationClassSpec* rg_spec = static_cast<RadiusofGyrationClassSpec*>(ispec);
    
    // Look through each molecule that could have active interactions
    for (k = 0; k < (int)(rg_spec->topo_data_->molecule_list->n_sites_); k++) {
    	order_radius_of_gyration_fm_matrix_element_calculation(this, topo_data.cg_site_types, n_cg_types, mat, x, simulation_box_half_lengths);
    }
}

// Calculate matrix elements for density non-bonded interactions.
// First, find the density at each site by calculating weight functions between all pairs of neighbors for all particles and call weight function calculation
// for each pair of density groups that interact. Exclusion lists are handled in the called subroutines.
// Then, calculate the matrix elements by looking through the neighbor list (for all pairs of neighbors for all particles) to calculate matrix elements.

void DensityClassComputer::calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
{
	if (ispec->get_n_defined() == 0) return;
	
	// Reset density array before accumulating weight function contributions;
	reset_density_array();
	
	// Set computer variables about matrix position
	trajectory_block_frame_index = traj_block_frame_index;
    current_frame_starting_row = curr_frame_starting_row;
    
	// First, pass through the neighbor list to compute the value of each density_group at every relavent CG site.
	walk_density_neighbor_list(mat, calculate_density_values, n_cg_types, topo_data, pair_cell_list, x, simulation_box_half_lengths);

	// Do intermediate processing (if necessary).
	process_completed_density(this, process_density, n_cg_types, topo_data.cg_site_types, mat, x, simulation_box_half_lengths);
		
	// Finally, calculate the matrix elements by combining the density, density derivative, pair distance, and pair derivative.
	walk_density_neighbor_list(mat, calculate_fm_matrix_elements, n_cg_types, topo_data, pair_cell_list, x, simulation_box_half_lengths);
}
  
// Calculate matrix elements for three body non-bonded interactions.
// Find all pairs of neighbors of all particles and call nonbonded matrix element computations
// for any triples that interact. Exclusion lists are handled in the called subroutines.

void ThreeBodyNonbondedClassComputer::calculate_3B_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const ThreeBCellList& three_body_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
{
    if (ispec->n_defined == 0) return;
    if (ispec->class_subtype > 0) {                    
        trajectory_block_frame_index = traj_block_frame_index;
        current_frame_starting_row = curr_frame_starting_row;
       	walk_3B_neighbor_list(mat, n_cg_types, topo_data, three_body_cell_list, x, simulation_box_half_lengths);
	}
}

inline void InteractionClassComputer::walk_3B_neighbor_list(MATRIX_DATA* const mat, const int n_cg_types, const TopologyData& topo_data, const ThreeBCellList& three_body_cell_list, std::array<double, DIMENSION>* const &x, const real* simulation_box_half_lengths) 
{
	int stencil_size = three_body_cell_list.get_stencil_size();
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
                    for (int nei_3 = 0; nei_3 < stencil_size; nei_3++) {
                        int ll_3 = three_body_cell_list.stencil[stencil_size * kk + nei_3];
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
            
            for (int nei = 0; nei < stencil_size; nei++) {
                int ll = three_body_cell_list.stencil[stencil_size * kk + nei];
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
                    for (int nei_3 = nei + 1; nei_3 < stencil_size; nei_3++) {
                        int ll_3 = three_body_cell_list.stencil[stencil_size * kk + nei_3];
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
//  Routines for checking if nonbonded interactions should be excluded
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

inline bool check_density_excluded_list(const TopologyData* const topo_data, const int i, const int j)
{
	// Check whetehr this non-bonded interaction is excluded from the model
	for (unsigned k = 0; k < topo_data->density_exclusion_list->partner_numbers_[i]; k++) {
        if (topo_data->density_exclusion_list->partners_[i][k] == unsigned(j)) return true;
    }
    return false;
}

//--------------------------------------------------------------------
// Routines responsible for calling single-interaction matrix computations,
// differing by the way that potentially interacting particles are found in 
// each frame and possibly found not to interact after.
//--------------------------------------------------------------------

void order_one_body_fm_matrix_element_calculation(InteractionClassComputer* const info, calc_pair_matrix_elements calc_matrix_elements, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths)
{
    // Calculate the appropriate matrix elements.
    info->index_among_defined_intrxns = cg_site_types[info->k] - 1;
    info->index_among_matched_interactions = info->ispec->defined_to_matched_intrxn_index_map[info->index_among_defined_intrxns];
    info->index_among_tabulated_interactions = info->ispec->defined_to_tabulated_intrxn_index_map[info->index_among_defined_intrxns];
    if ((info->index_among_matched_interactions == 0) && (info->index_among_tabulated_interactions == 0)) return; // if the index is zero, it is not present in the model and should be ignored.
    calc_matrix_elements(info, x, simulation_box_half_lengths, mat);
}

void order_pair_nonbonded_fm_matrix_element_calculation(InteractionClassComputer* const info, calc_pair_matrix_elements calc_matrix_elements, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths)
{
    // Calculate the appropriate matrix elements.
    info->index_among_defined_intrxns = info->ispec->get_index_from_hash(calc_two_body_interaction_hash(cg_site_types[info->k], cg_site_types[info->l], n_cg_types));
    info->set_indices();

	if (info->index_among_matched_interactions == 0) return; // if the index is zero, it is not present in the model and should be ignored.
    calc_matrix_elements(info, x, simulation_box_half_lengths, mat);
}

void order_bonded_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths)
{
     // Calculate the appropriate matrix elements.    
    info->index_among_defined_intrxns = info->ispec->get_index_from_hash(info->calculate_hash_number(cg_site_types, n_cg_types));
    info->set_indices();

    if (info->index_among_matched_interactions == 0) return; // if the index is zero, it is not present in the model and should be ignored.

    (*info->calculate_fm_matrix_elements)(info, x, simulation_box_half_lengths, mat);
}

void order_radius_of_gyration_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths)
{
	RadiusofGyrationClassSpec* rg_spec = static_cast<RadiusofGyrationClassSpec*>(info->ispec);
	int n_cg_molecules = rg_spec->topo_data_->molecule_list->n_sites_;
	// See what interaction types are active for this molecule (info->k)
	// Set the appropriate indices if it is active.
	for (int i = 0; i < rg_spec->n_molecule_groups; i++) {
		if (rg_spec->molecule_groups[i*(n_cg_molecules) + info->k] == true) {
		
			// Calculate the appropriate matrix elements.    
			info->index_among_defined_intrxns = i;
			info->set_indices();

			if (info->index_among_matched_interactions == 0) return; // if the index is zero, it is not present in the model and should be ignored.
			(*info->calculate_fm_matrix_elements)(info, x, simulation_box_half_lengths, mat);
		}
	}
}

void order_helical_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths)
{
	HelicalClassSpec* h_spec = static_cast<HelicalClassSpec*>(info->ispec);
	int n_cg_molecules = h_spec->topo_data_->molecule_list->n_sites_;
	// See what interaction types are active for this molecule (info->k)
	// Set the appropriate indices if it is active.
	for (int i = 0; i < h_spec->n_molecule_groups; i++) {
		if (h_spec->molecule_groups[i*(n_cg_molecules) + info->k] == true) {
		
			// Calculate the appropriate matrix elements.    
			info->index_among_defined_intrxns = i;
			info->set_indices();

			if (info->index_among_matched_interactions == 0) return; // if the index is zero, it is not present in the model and should be ignored.
			(*info->calculate_fm_matrix_elements)(info, x, simulation_box_half_lengths, mat);
		}
	}
}

void order_three_body_nonbonded_fm_matrix_element_calculation(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths)
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

void density_fm_matrix_element_calculation(InteractionClassComputer* const info, calc_pair_matrix_elements calc_matrix_elements, int* const cg_site_types, const int n_cg_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths)
{
	DensityClassComputer* icomp = static_cast<DensityClassComputer*>(info);
	DensityClassSpec* ispec = static_cast<DensityClassSpec*>(icomp->ispec);
    
	// Calculate the appropriate matrix elements.
	// Get the bit flags for interactions encoded for this pair of types
	unsigned long interaction_flags = ispec->site_to_density_group_intrxn_index_map[(cg_site_types[info->k] - 1) * n_cg_types + (cg_site_types[info->l] - 1)];
	decode_density_interaction_and_calculate(icomp, interaction_flags, calc_matrix_elements, cg_site_types, mat, x, simulation_box_half_lengths);
	
	// Repeat this for the reversed pair of types.
	swap_pair(info->k, info->l);
	interaction_flags = ispec->site_to_density_group_intrxn_index_map[(cg_site_types[info->k] - 1) * n_cg_types + (cg_site_types[info->l] - 1)];
	decode_density_interaction_and_calculate(icomp, interaction_flags, calc_matrix_elements, cg_site_types, mat, x, simulation_box_half_lengths);
	//restore k and l
	swap_pair(info->k, info->l);
}

//---------------------------------------------------------------------
// Helper functions of functions in the above section.
//---------------------------------------------------------------------

void process_completed_density(DensityClassComputer* const info, calc_pair_matrix_elements process_density, const int n_cg_types, int *const cg_site_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths) 
{	
	DensityClassSpec* ispec = static_cast<DensityClassSpec*>(info->ispec);
	
	// The "process_density" function pointer should be set to do_nothing for force matching and evaluate_density_sampling_range for range finding.
	// Go through all types, determine if they belong to a defined density group, then call "process_density" for each density calculated at that site.
	for(int i = 0; i < ispec->n_cg_sites; i++) {
		// Does this group belong do any defined density_group
		for(int dg1 = 0; dg1 < ispec->n_density_groups; dg1++) {
		
			if(ispec->density_groups[dg1 * ispec->n_density_groups + cg_site_types[i] - 1] == false) continue;
			
			// Go through all densities that could be calcualted at this site
			for(int dg2 = 0; dg2 < ispec->n_density_groups; dg2++) {
				info->index_among_defined_intrxns = dg1 * ispec->n_density_groups + dg2;
				info->k = i;
				process_density(info, x, simulation_box_half_lengths, mat);
			}
		}
	}
}

inline void decode_density_interaction_and_calculate(DensityClassComputer* info, unsigned long interaction_flags, calc_pair_matrix_elements calc_matrix_elements, int *const cg_site_types, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths)
{
	DensityClassSpec* ispec = static_cast<DensityClassSpec*>(info->ispec);
	// Loop through to determine which bits are non-zero.
	// For each non-zero bit, look-up the matched or tabulated interaction, and perform the requested calculation
	for(int index_counter = 0; interaction_flags != 0; index_counter++) {
		// This could easily be a while loop over interaction_flags with a manaully incremented counter.
		if(interaction_flags % 2 == 1) {
			// Look-up this index
			info->index_among_defined_intrxns = index_counter;
			std::vector<int> types = ispec->get_interaction_types(info->index_among_defined_intrxns);
			int contributing_density_group = types[1] - 1;
			info->curr_weight = ispec->density_weights[ contributing_density_group * ispec->n_density_groups + (cg_site_types[info->l] - 1)];
			(*calc_matrix_elements)(info, x, simulation_box_half_lengths, mat);
		}
		// Shift to the right and repeat the operation
		interaction_flags = interaction_flags >> 1;
	}
}

void accumulate_matching_order_parameter_forces(InteractionClassComputer* const info, const int first_nonzero_basis_index, double extra_derivative_value, std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat) 
{
	for (unsigned k = 0; k < basis_fn_vals.size(); k++) {
		basis_fn_vals[k] *= extra_derivative_value;
	}
	mat->accumulate_matching_forces(info, first_nonzero_basis_index, basis_fn_vals, n_body, particle_ids, derivatives, mat);
}

//--------------------------------------------------------------------
// Functions for calculating sets of dimension(e.g.3)-component matrix elements for each
// individual interacting set of particles
//--------------------------------------------------------------------

inline void process_normal_interaction_matrix_elements(InteractionClassComputer* const info, MATRIX_DATA* const mat, const int n_body, int* particle_ids, std::array<double, DIMENSION>* derivatives, const double param_value, const int virial_flag, const double junk = 0.0, const double junk2 = 0.0)
{
	int index_among_defined = info->index_among_defined_intrxns;
	int index_among_matched = info->index_among_matched_interactions;
    int index_among_tabulated = info->index_among_tabulated_interactions;
    int index_among_symmetric = info->index_among_symmetric_interactions;
    int index_among_symtab = info->index_among_symtab_interactions;
    int first_nonzero_basis_index;
    int temp_column_index;
    double basis_sum;
    
    if (index_among_tabulated > 0) {
		// Pull the interaction from a table. 	   
    	info->table_s_comp->calculate_basis_fn_vals(index_among_defined, param_value, first_nonzero_basis_index, info->table_basis_fn_vals);
    	if (n_body == 1) {
    		basis_sum = info->table_basis_fn_vals[0];
    	} else {
    		basis_sum  = info->table_basis_fn_vals[0] + info->table_basis_fn_vals[1];
    	}

    	// Add to force target.
		if (index_among_symtab == 0) {
			// This is an antisymmetric (i.e. force) interaction table.
    	    mat->accumulate_tabulated_forces(info, basis_sum, n_body, particle_ids, derivatives, mat);
    	} else {
		    // This is a symmetric (i.e. DOOM) interaction table.
    	    mat->accumulate_symmetric_tabulated_forces(info, basis_sum, n_body, particle_ids, derivatives, mat);
    	}
    	
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
    	if (index_among_symmetric == 0) {
    	    // This interaction is antisymmetric (i.e. force).
		    mat->accumulate_matching_forces(info, first_nonzero_basis_index, info->fm_basis_fn_vals, n_body, particle_ids, derivatives, mat);
    	} else {
        	// This interaction is symmetric (i.e. DOOM).
     		mat->accumulate_symmetric_matching_forces(info, first_nonzero_basis_index, info->fm_basis_fn_vals, n_body, particle_ids, derivatives, mat);
 		}
 			
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

inline void process_one_body_matrix_elements(InteractionClassComputer* const info, MATRIX_DATA* const mat, const int n_body, int* particle_ids, std::array<double, DIMENSION>* derivatives, const double param_value, const int virial_flag, const double junk = 0.0, const double junk2 = 0.0)
{
    int index_among_defined = info->index_among_defined_intrxns;
    int index_among_matched = info->index_among_matched_interactions;
    int index_among_tabulated = info->index_among_tabulated_interactions;
    int first_nonzero_basis_index;
    double basis_sum;

    if (index_among_tabulated > 0) {
    	// Pull the interaction from a table.
        info->table_s_comp->calculate_basis_fn_vals(index_among_defined, 1.0, first_nonzero_basis_index, info->table_basis_fn_vals);
		basis_sum = info->table_basis_fn_vals[0];
        // Add to force target.
        mat->accumulate_one_body_tabulated_force(info, basis_sum, 1, particle_ids, derivatives, mat);
        // Add to virial target.
        if (mat->virial_constraint_rows > 0) mat->accumulate_target_constraint_element(mat, info->trajectory_block_frame_index, -basis_sum);
    }
        
    if (index_among_matched > 0) {
    	// Compute the strength of each basis function.
        info->fm_s_comp->calculate_basis_fn_vals(index_among_defined, 1.0, first_nonzero_basis_index, info->fm_basis_fn_vals);
        basis_sum = info->fm_basis_fn_vals[0];
        // Add to the force matching.
        mat->accumulate_one_body_force(info, first_nonzero_basis_index, info->fm_basis_fn_vals, 1, particle_ids, derivatives, mat);
        // Add to virial matching.
        int temp_column_index = info->interaction_class_column_index + info->ispec->interaction_column_indices[index_among_matched - 1] + first_nonzero_basis_index;
        int basis_column = temp_column_index;
        if (mat->virial_constraint_rows > 0)(*mat->accumulate_virial_constraint_matrix_element)(info->trajectory_block_frame_index, basis_column, basis_sum, mat);
    }
}

inline void process_density_matrix_elements(InteractionClassComputer* const info, MATRIX_DATA* const mat, const int n_body, int* particle_ids, std::array<double, DIMENSION>* derivatives, const double density_value, const int virial_flag, const double density_derivative, const double distance)
{
    int index_among_defined = info->index_among_defined_intrxns;
    int index_among_matched = info->index_among_matched_interactions;
    int index_among_tabulated = info->index_among_tabulated_interactions;
    int first_nonzero_basis_index;
    double basis_sum;

    if (index_among_tabulated > 0) {
		// Pull the interaction from a table.
        info->table_s_comp->calculate_basis_fn_vals(index_among_defined, density_value, first_nonzero_basis_index, info->table_basis_fn_vals);
        basis_sum = info->table_basis_fn_vals[0] + info->table_basis_fn_vals[1];
        // Add to force target.
        mat->accumulate_tabulated_forces(info, basis_sum * density_derivative, 2, particle_ids, derivatives, mat);
        // Add to virial target.
        if (mat->virial_constraint_rows > 0) mat->accumulate_target_constraint_element(mat, info->trajectory_block_frame_index, -basis_sum * density_derivative * distance);
    }
    
    if (index_among_matched > 0) {
        // Compute the strength of each basis function.
       	info->fm_s_comp->calculate_basis_fn_vals(index_among_defined, density_value, first_nonzero_basis_index, info->fm_basis_fn_vals);
		// Add to the force matching.
        accumulate_matching_order_parameter_forces(info, first_nonzero_basis_index, density_derivative, info->fm_basis_fn_vals, 2, particle_ids, derivatives, mat);
        // Add to virial matching.
        int temp_column_index = info->interaction_class_column_index + info->ispec->interaction_column_indices[index_among_matched - 1] + first_nonzero_basis_index;
        for (unsigned i = 0; i < info->fm_basis_fn_vals.size(); i++) {
        	int basis_column = temp_column_index + i;
            if (mat->virial_constraint_rows > 0)(*mat->accumulate_virial_constraint_matrix_element)(info->trajectory_block_frame_index, basis_column, info->fm_basis_fn_vals[i] * distance, mat);
			// This virial expression already includes the density derivative, which was multiplied into fm_basis_fn_vals in accumulate_matching_order_parameter_forces.
        }
    }
}	

//--------------------------------------------------------------------
// Functions for calculating sets of 3-component matrix elements for each
// individual interacting set of particles
//--------------------------------------------------------------------

// Each of these functions follows the idiom of calc_isotropic_two_body_fm_matrix_elements.

void calc_one_body_fm_matrix_element(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
	int particle_ids[1] = {info->k};
    std::array<double, DIMENSION>* derivatives =  new std::array<double, DIMENSION>[1];
	for (int i = 0; i < DIMENSION; i++) derivatives[0][i] = 1.0/(double)(DIMENSION);
    info->process_interaction_matrix_elements(info, mat, 1, particle_ids, derivatives, 1.0, 1, 0.0, 0.0);
	delete [] derivatives;
}

void calc_isotropic_two_body_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    int particle_ids[2] = {info->k, info->l};
    std::array<double, DIMENSION>* derivatives = new std::array<double, DIMENSION>[1];
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

void calc_angular_three_body_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    int particle_ids[3] = {info->k, info->l, info->j}; // end indices (k, l), followed by center index (j)
    std::array<double, DIMENSION>* derivatives = new std::array<double, DIMENSION>[2];
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

void calc_dihedral_four_body_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    int particle_ids[4] = {info->k, info->l, info->i, info->j}; // end indices (k, l) followed by central bond indices (i, j)
    std::array<double, DIMENSION>* derivatives = new std::array<double, DIMENSION>[3];
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

void calc_radius_of_gyration_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
	RadiusofGyrationClassSpec* rg_spec = static_cast<RadiusofGyrationClassSpec*>(info->ispec);
	
	int index_among_defined = info->index_among_defined_intrxns;
	int n_ids = rg_spec->topo_data_->molecule_list->partner_numbers_[info->k];
	std::array<double, DIMENSION>* derivatives = new std::array<double, DIMENSION>[n_ids - 1];
	int* particle_ids = new int[n_ids];
	for (int i = 0; i < n_ids; i++) particle_ids[i] = (int)(rg_spec->topo_data_->molecule_list->partners_[info->k][i]);
	
	double radius_of_gyration;
	calc_radius_of_gyration_and_derivatives(particle_ids, x, simulation_box_half_lengths, n_ids, radius_of_gyration, derivatives);

	if (radius_of_gyration >= rg_spec->lower_cutoffs[index_among_defined] &&
	    radius_of_gyration <= rg_spec->upper_cutoffs[index_among_defined]) {
		info->process_interaction_matrix_elements(info, mat, n_ids, particle_ids, derivatives, radius_of_gyration, 0, 0.0, 0.0);
	}
	delete [] derivatives;
	delete [] particle_ids;
}

void calc_helical_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
	HelicalClassSpec* h_spec = static_cast<HelicalClassSpec*>(info->ispec);
	// k is currently the molecule id.
	int index_among_defined = info->index_among_defined_intrxns;
	
	int n_ids = h_spec->topo_data_->molecule_list->partner_numbers_[info->k];
	std::array<double, DIMENSION>* derivatives = new std::array<double, DIMENSION>[n_ids - 1];
	int* particle_ids = new int[n_ids];
	for (int i = 0; i < n_ids; i++) particle_ids[i] = (int)(h_spec->topo_data_->molecule_list->partners_[info->k][i]);
	
	int n_helical_ids = 2 * h_spec->helical_list->partner_numbers_[info->k];
	int* helical_ids = new int[n_helical_ids];
	for (int i = 0; i < n_helical_ids; i++) helical_ids[i] = h_spec->helical_list->partners_[info->k][i]; 
	
	double fraction_helical;
	calc_fraction_helical_and_derivatives(particle_ids, x, simulation_box_half_lengths, n_ids, fraction_helical, derivatives, helical_ids, n_helical_ids/2, h_spec->r0[index_among_defined], h_spec->sigma2[index_among_defined]);

	if (fraction_helical >= h_spec->lower_cutoffs[index_among_defined] &&
	    fraction_helical <= h_spec->upper_cutoffs[index_among_defined]) {
		info->process_interaction_matrix_elements(info, mat, n_ids, particle_ids, derivatives, fraction_helical, 0, 0.0, 0.0);
	}
	delete [] derivatives;
	delete [] particle_ids;
	delete [] helical_ids;
}

void calc_nonbonded_1_three_body_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    int particle_ids[3] = {info->k, info->l, info->j}; // end indices (k, l) followed by center index (j).    
    ThreeBodyNonbondedClassComputer* icomp = static_cast<ThreeBodyNonbondedClassComputer*>(info);
    ThreeBodyNonbondedClassSpec* ispec = static_cast<ThreeBodyNonbondedClassSpec*>(icomp->ispec);

	std::array<double, DIMENSION>* relative_site_position_2 = new std::array<double, DIMENSION>[1];
	std::array<double, DIMENSION>* relative_site_position_3 = new std::array<double, DIMENSION>[1];
	std::array<double, DIMENSION>* derivatives = new std::array<double, DIMENSION>[2];
	std::array<double, DIMENSION> tx1, tx2, tx;
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

void calc_nonbonded_2_three_body_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    int particle_ids[3] = {info->k, info->l, info->j}; // end indices (k, l) followed by center index (j).    
    ThreeBodyNonbondedClassComputer* icomp = static_cast<ThreeBodyNonbondedClassComputer*>(info);
    ThreeBodyNonbondedClassSpec* ispec = static_cast<ThreeBodyNonbondedClassSpec*>(icomp->ispec);
    
	std::array<double, DIMENSION>* relative_site_position_2 = new std::array<double, DIMENSION>[1];
	std::array<double, DIMENSION>* relative_site_position_3 = new std::array<double, DIMENSION>[1];
	std::array<double, DIMENSION>* derivatives = new std::array<double, DIMENSION>[2];
	std::array<double, DIMENSION> tx1, tx2, tx;
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

void calc_gaussian_density_values(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
	DensityClassComputer* icomp = static_cast<DensityClassComputer*>(info);
	DensityClassSpec* ispec = static_cast<DensityClassSpec*>(icomp->ispec);
	int particle_ids[2] = {icomp->k, icomp->l};
	int index_among_defined = icomp->index_among_defined_intrxns;
    double distance2;
    
	//Calculate the distance
	calc_squared_distance(particle_ids, x, simulation_box_half_lengths, distance2);
	if (distance2 < icomp->cutoff2) {
		// Calculate the weight function
		double distance = sqrt(distance2);
		icomp->density_values[icomp->index_among_defined_intrxns * ispec->n_cg_sites + icomp->k] +=
										icomp->curr_weight * ( exp( - distance2 / icomp->denomenator[index_among_defined]) + icomp->u_cutoff[index_among_defined]
										+ icomp->f_cutoff[index_among_defined] * (distance - ispec->cutoff) ) / icomp->denomenator[index_among_defined];
	}
}

void calc_switching_density_values(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
	DensityClassComputer* icomp = static_cast<DensityClassComputer*>(info);
	DensityClassSpec* ispec = static_cast<DensityClassSpec*>(icomp->ispec);
	int particle_ids[2] = {icomp->k, icomp->l};
    int index_among_defined = icomp->index_among_defined_intrxns;
    double distance2;
    
	//Calculate the distance
	calc_squared_distance(particle_ids, x, simulation_box_half_lengths, distance2);
	
	if (distance2 < icomp->cutoff2) {
	
		// Calculate the weight function
		double distance = sqrt(distance2);
		icomp->density_values[icomp->index_among_defined_intrxns * ispec->n_cg_sites  + icomp->k] +=
										icomp->curr_weight * -0.5 * tanh( (distance - ispec->density_switch[index_among_defined])/ispec->density_sigma[index_among_defined] )
										+ icomp->u_cutoff[index_among_defined] + icomp->f_cutoff[index_among_defined] * (distance - ispec->cutoff);
	}
}

void calc_lucy_density_values(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
	DensityClassComputer* icomp = static_cast<DensityClassComputer*>(info);
	DensityClassSpec* ispec = static_cast<DensityClassSpec*>(icomp->ispec);
	int particle_ids[2] = {icomp->k, icomp->l};
    int index_among_defined = icomp->index_among_defined_intrxns;
    double distance2;
    
	//Calculate the distance
	calc_squared_distance(particle_ids, x, simulation_box_half_lengths, distance2);
	
	if (distance2 < icomp->cutoff2) {
	
		// Calculate the weight function
		double distance = sqrt(distance2);
		double cutoff_minus_distance = ispec->cutoff - distance;
		icomp->density_values[icomp->index_among_defined_intrxns * ispec->n_cg_sites + icomp->k] +=
										icomp->curr_weight * cutoff_minus_distance * cutoff_minus_distance * cutoff_minus_distance 
										* (ispec->cutoff + 3.0*distance) / icomp->denomenator[index_among_defined];
	}
}

void calc_re_density_values(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
	DensityClassComputer* icomp = static_cast<DensityClassComputer*>(info);
	DensityClassSpec* ispec = static_cast<DensityClassSpec*>(icomp->ispec);
	int particle_ids[2] = {icomp->k, icomp->l};
    int index_among_defined = icomp->index_among_defined_intrxns;
    double distance2;
    
	//Calculate the distance
	calc_squared_distance(particle_ids, x, simulation_box_half_lengths, distance2);
	
	if (distance2 < icomp->cutoff2) {
	
		// Calculate the weight function
		if (distance2 > ispec->density_sigma[index_among_defined] * ispec->density_sigma[index_among_defined]) {
			icomp->density_values[icomp->index_among_defined_intrxns * ispec->n_cg_sites + icomp->k] +=
										icomp->curr_weight * (icomp->c0[index_among_defined] +
										distance2 * icomp->c2[index_among_defined] - 
										distance2 * distance2 * icomp->c4[index_among_defined] +
										distance2 * distance2 * distance2 * icomp->c6[index_among_defined]);
		} else {
			icomp->density_values[icomp->index_among_defined_intrxns * ispec->n_cg_sites + icomp->k] += 1.0 * icomp->curr_weight;
		}
	}
}

void calc_density_fm_matrix_elements(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
	info->index_among_matched_interactions = info->ispec->defined_to_matched_intrxn_index_map[info->index_among_defined_intrxns];
	info->index_among_tabulated_interactions = info->ispec->defined_to_tabulated_intrxn_index_map[info->index_among_defined_intrxns];
	if ((info->index_among_matched_interactions == 0) && (info->index_among_tabulated_interactions == 0)) return; // if the index is zero, it is not present in the model and should be ignored.
	
	double distance;
    int particle_ids[2] = {info->k, info->l};
    std::array<double, DIMENSION>* derivatives = new std::array<double, DIMENSION>[1];
    if ( conditionally_calc_distance_and_derivatives(particle_ids, x, simulation_box_half_lengths, info->cutoff2, distance, derivatives) ) {
            
        DensityClassComputer* icomp = static_cast<DensityClassComputer*>(info);
		DensityClassSpec* ispec = static_cast<DensityClassSpec*>(icomp->ispec);
	
		// Look-up this particular interaction's density.
		double density_value = icomp->density_values[icomp->index_among_defined_intrxns * ispec->n_cg_sites + icomp->k];
		// Calculate the weight function derivative.
		double density_derivative = (*icomp->calculate_density_derivative)(icomp, ispec, distance);
		density_derivative *= icomp->curr_weight;
		
		info->process_interaction_matrix_elements(info, mat, 2, particle_ids, derivatives, density_value, 1, density_derivative, distance);
    }	
    delete [] derivatives;
}

double calc_gaussian_density_derivative(DensityClassComputer* const icomp, DensityClassSpec* const ispec, const double distance)
{
	int index_among_defined = icomp->index_among_defined_intrxns;
	double density_derivative = - (2.0 * distance / icomp->denomenator[index_among_defined]) * exp( - distance * distance / icomp->denomenator[index_among_defined]);
	density_derivative += icomp->f_cutoff[index_among_defined];
	return density_derivative;
}

double calc_switching_density_derivative(DensityClassComputer* const icomp, DensityClassSpec* const ispec, const double distance)
{
	int index_among_defined = icomp->index_among_defined_intrxns;
	double arguement = (distance - ispec->density_switch[index_among_defined])/ ispec->density_sigma[index_among_defined];
	double density_derivative = - 0.5 / (ispec->density_sigma[index_among_defined] * cosh(arguement) * cosh(arguement));
	density_derivative += icomp->f_cutoff[index_among_defined];
	return density_derivative;
}

double calc_lucy_density_derivative(DensityClassComputer* const icomp, DensityClassSpec* const ispec, const double distance)
{
	int index_among_defined = icomp->index_among_defined_intrxns;
	double cutoff_minus_distance = ispec->cutoff - distance;
	double density_derivative = -12.0 * distance * cutoff_minus_distance * cutoff_minus_distance / icomp->denomenator[index_among_defined];
	return density_derivative;
}

double calc_re_density_derivative(DensityClassComputer* const icomp, DensityClassSpec* const ispec, const double distance)
{
	int index_among_defined = icomp->index_among_defined_intrxns;
	double distance2 = distance * distance;
	double density_derivative = 2.0 * icomp->c2[index_among_defined] - 4.0 * distance2 * icomp->c4[index_among_defined] + 6.0 * distance2 * distance2 * icomp->c6[index_among_defined];
	density_derivative *= distance;	
	return density_derivative;
}

void do_nothing(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat) 
{
}
