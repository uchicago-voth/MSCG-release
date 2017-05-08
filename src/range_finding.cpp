//
//  range_finding.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>

#include "force_computation.h"
#include "geometry.h"
#include "interaction_model.h"
#include "matrix.h"
#include "range_finding.h"
#include "misc.h"

const double VERYLARGE = 1000.0;

//----------------------------------------------------------------------------
// Prototypes for private implementation routines.
//----------------------------------------------------------------------------

// Initialization of storage for the range value arrays and their computation

void initialize_ranges(int tol, double* const lower_cutoffs, double* const upper_cutoffs, std::vector<unsigned> &num);

void initialize_single_class_range_finding_temps(InteractionClassSpec *iclass, InteractionClassComputer *icomp, TopologyData *topo_data);

// Helper functions that issues failure warnings and do special setup
void report_unrecognized_class_subtype(InteractionClassSpec *iclass);
void allocate_and_initialize_density_computer_for_range_finding(DensityClassComputer* icomp);
void setup_site_to_density_group_index_for_range(DensityClassSpec* iclass);

// Functions for computing the full range of sampling of a given class of interaction in a given trajectory.
void calc_isotropic_two_body_sampling_range(InteractionClassComputer* const icomp, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_angular_three_body_sampling_range(InteractionClassComputer* const icomp, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_dihedral_four_body_interaction_sampling_range(InteractionClassComputer* const icomp, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_radius_of_gyration_interaction_sampling_range(InteractionClassComputer* const icomp, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void evaluate_density_sampling_range(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_nothing(InteractionClassComputer* const icomp, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);

void write_interaction_range_data_to_file(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FILE* const one_body_spline_output_filep,  FILE* const nonbonded_spline_output_filep, FILE* const bonded_spline_output_filep, FILE* const density_interaction_output_filep, FILE* const radius_of_gyration_output_filep);

void write_iclass_range_specifications(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, FILE* const solution_spline_output_file);
void write_one_body_iclass_range_specifications(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, FILE* const solution_spline_output_file);
void write_single_range_specification(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, FILE* const solution_spline_output_file, const int index_among_defined);

void read_density_parameter_file(DensityClassSpec* const ispec);
void read_interaction_file_and_build_matrix(MATRIX_DATA* mat, CG_MODEL_DATA* const cg, double volume);
void read_one_param_dist_file_pair(InteractionClassComputer* const icomp, char ** const name, MATRIX_DATA* mat, const int index_among_defined_intrxns, int &counter, double num_of_pairs, double volume);
void read_one_param_dist_file_other(InteractionClassComputer* const icomp, char ** const name, MATRIX_DATA* mat, const int index_among_defined_intrxns, int &counter, double num_of_pairs);
double count_bonded_interaction(InteractionClassComputer* const icomp, char** const name, MATRIX_DATA* mat, const int index_among_defined_intrxns);
double calculate_volume(const matrix simulation_box_lengths);

// Output parameter distribution functions
void open_parameter_distribution_files_for_class(InteractionClassComputer* const icomp, char **name); 
void open_density_parameter_distribution_files_for_class(InteractionClassComputer* const icomp, char **name); 
void close_parameter_distribution_files_for_class(InteractionClassComputer* const icomp);
void generate_parameter_distribution_histogram(InteractionClassComputer* const icomp, char **name);

// Dummy implementations
void do_not_initialize_fm_matrix(MATRIX_DATA* const mat);

//------------------------------------------------------------------------
//    Implementation
//------------------------------------------------------------------------

void do_not_initialize_fm_matrix(MATRIX_DATA* const mat) {}

void initialize_ranges(int tol, double* const lower_cutoffs, double* const upper_cutoffs, std::vector<unsigned> &num)
{
    for (int i = 0; i < tol; i++) {
        lower_cutoffs[i] = VERYLARGE;
        upper_cutoffs[i] = -VERYLARGE;
        num[i] = i + 1;
    }
}

void initialize_range_finding_temps(CG_MODEL_DATA* const cg)
{   
	std::list<InteractionClassSpec*>::iterator iclass_iterator;
	std::list<InteractionClassComputer*>::iterator icomp_iterator;
	for(iclass_iterator=cg->iclass_list.begin(), icomp_iterator=cg->icomp_list.begin(); (iclass_iterator != cg->iclass_list.end()) && (icomp_iterator != cg->icomp_list.end()); iclass_iterator++, icomp_iterator++) {	
		initialize_single_class_range_finding_temps((*iclass_iterator), (*icomp_iterator), &cg->topo_data);
	}
    initialize_single_class_range_finding_temps(&cg->three_body_nonbonded_interactions, &cg->three_body_nonbonded_computer, &cg->topo_data);
	if(cg->density_interactions.class_subtype != 0) {
		read_density_parameter_file(&cg->density_interactions);
		allocate_and_initialize_density_computer_for_range_finding(&cg->density_computer);
	}
	cg->pair_nonbonded_cutoff2 = VERYLARGE * VERYLARGE;
}

void initialize_single_class_range_finding_temps(InteractionClassSpec *iclass, InteractionClassComputer *icomp, TopologyData *topo_data) 
{
	if( ((iclass->class_type == kDensity || iclass->class_type == kRadiusofGyration) && iclass->class_subtype == 0) ) {
		iclass->dummy_setup_for_defined_interactions(topo_data);
	} else {	 
		iclass->setup_for_defined_interactions(topo_data);
	}
	
    icomp->ispec = iclass;
    if (iclass->class_type == kOneBody) {
        icomp->calculate_fm_matrix_elements = calc_nothing;
    } else if (iclass->class_type == kPairNonbonded) {
        icomp->calculate_fm_matrix_elements = calc_isotropic_two_body_sampling_range;
    } else if (iclass->class_type == kPairBonded) {
        icomp->calculate_fm_matrix_elements = calc_isotropic_two_body_sampling_range;
    } else if (iclass->class_type == kAngularBonded) {
        if (iclass->class_subtype == 0) { // Angle based angular interactions
			icomp->calculate_fm_matrix_elements = calc_angular_three_body_sampling_range;
        } else if (iclass->class_subtype == 1) { // Distance based angular interactions
			icomp->calculate_fm_matrix_elements = calc_isotropic_two_body_sampling_range;
		} else {
			report_unrecognized_class_subtype(iclass);
		}
    } else if (iclass->class_type == kDihedralBonded) {
        if (iclass->class_subtype == 0) { //Angle based dihedral interactions
			icomp->calculate_fm_matrix_elements = calc_dihedral_four_body_interaction_sampling_range;
        } else if (iclass->class_subtype == 1) { // Distance based dihedral interactions
			icomp->calculate_fm_matrix_elements = calc_isotropic_two_body_sampling_range;
		} else {
			report_unrecognized_class_subtype(iclass);
		}
	} else if (iclass->class_type == kRadiusofGyration) {
		if (iclass->class_subtype == 1) { 
			icomp->calculate_fm_matrix_elements = calc_radius_of_gyration_interaction_sampling_range;
		} else if (iclass->class_subtype == 0) { // Do nothing
			icomp->calculate_fm_matrix_elements = calc_nothing;
		} else {
			report_unrecognized_class_subtype(iclass);
		}
    } else if (iclass->class_type == kDensity) {
		DensityClassComputer* dcomp = static_cast<DensityClassComputer*>(icomp);
		dcomp->cutoff2 = iclass->cutoff * iclass->cutoff;
		dcomp->process_density = evaluate_density_sampling_range;
		dcomp->calculate_fm_matrix_elements = calc_nothing;
		if (iclass->class_subtype == 1) { // Continuously varying (Gaussian) weight function
			dcomp->calculate_density_values	= calc_gaussian_density_values;
			printf("Will calculate density using shifted-force Gaussian weight functions.\n");
		} else if (iclass->class_subtype == 2) { // Switching function (tanh) weight function
			dcomp->calculate_density_values = calc_switching_density_values;
			printf("Will calculate density using shifted-force switching (tanh) weight functions.\n");
		} else if (iclass->class_subtype == 3) { // Lucy-type weight function
			dcomp->calculate_density_values = calc_lucy_density_values;
			printf("Will calculate density using Lucy-style weight functions.\n");
		}  else if (iclass->class_subtype == 4) { // Lucy-type weight function
			dcomp->calculate_density_values = calc_re_density_values;
			printf("Will calculate density using Relative-Entropy style weight functions.\n");
		} else if (iclass->class_subtype == 0) { // Do nothing
			dcomp->calculate_density_values = calc_nothing;
		} else {
			report_unrecognized_class_subtype(iclass);
		}
	}else { // For three_body_interactions
    	icomp->calculate_fm_matrix_elements = calc_nothing;
    }
    
	iclass->n_cg_types = int(topo_data->n_cg_types);
    initialize_ranges(iclass->get_n_defined(), iclass->lower_cutoffs, iclass->upper_cutoffs, iclass->defined_to_matched_intrxn_index_map);
    iclass->n_to_force_match = iclass->get_n_defined();
    iclass->interaction_column_indices = std::vector<unsigned>(iclass->n_to_force_match + 1);
	
	if(iclass->output_parameter_distribution == 1){
		if(iclass->class_type == kDensity) {
			open_density_parameter_distribution_files_for_class(icomp, topo_data->density_group_names);
 		} else if (iclass->class_type == kRadiusofGyration) {
 			RadiusofGyrationClassSpec* rg_class = static_cast<RadiusofGyrationClassSpec*>(iclass);
			open_parameter_distribution_files_for_class(icomp, rg_class->molecule_group_names);
		} else {
			open_parameter_distribution_files_for_class(icomp, topo_data->name);
		}
	}
}

void report_unrecognized_class_subtype(InteractionClassSpec *iclass)
{
	printf("Unrecognized %s class subtype!\n", iclass->get_full_name().c_str());
	fflush(stdout);
	exit(EXIT_FAILURE);
}

void allocate_and_initialize_density_computer_for_range_finding(DensityClassComputer* icomp) 
{
	DensityClassSpec* iclass = static_cast<DensityClassSpec*>(icomp->ispec);
	// Allocate space to store density intermediate.
	// This approach grabs enough memory for all cg sites to have all density values for all density groups.
	// A more memory-efficient, but harder approach would be allocate each site's array based on the number of density groups needed at that site.
	// The problem with this other approach is efficiently looking up which index to use.
	icomp->density_values = new double[iclass->get_n_defined() * iclass->n_cg_sites]();
	
	// Allocate and compute constant calculation intermediates.
	icomp->denomenator = new double[iclass->get_n_defined()]();
	icomp->u_cutoff = new double[iclass->get_n_defined()]();
	icomp->f_cutoff = new double[iclass->get_n_defined()]();
	
	if(iclass->class_subtype == 1) {
		for(int i = 0; i < iclass->get_n_defined(); i++) {
			if (iclass->density_sigma[i] < VERYSMALL) {
				printf("Density sigma parameter is too small!\n");
				exit(EXIT_FAILURE);
			}
			icomp->denomenator[i] = 2.0 * iclass->density_sigma[i] * iclass->density_sigma[i];
			icomp->u_cutoff[i] = - exp( - icomp->cutoff2 / icomp->denomenator[i] );
			icomp->f_cutoff[i] = - 2.0 * iclass->cutoff * icomp->u_cutoff[i] / icomp->denomenator[i];
			
			printf("%d: density_sigma %lf, cutoff %lf, u_cutoff %lf, f_cutoff %lf, denom %lf\n", i, iclass->density_switch[i], iclass->cutoff, icomp->u_cutoff[i], icomp->f_cutoff[i], icomp->denomenator[i]); fflush(stdout);
		}
	} else if (iclass->class_subtype == 2) {
		for(int i = 0; i < iclass->get_n_defined(); i++) {
			if (iclass->density_sigma[i] < VERYSMALL) {
				printf("Density sigma parameter is too small!\n");
				exit(EXIT_FAILURE);
			}
			icomp->denomenator[i] = iclass->density_sigma[i] / 0.5;
			double arguement = (iclass->cutoff - iclass->density_switch[i]) / iclass->density_sigma[i];
			icomp->u_cutoff[i] = 0.5 * tanh( arguement );
			icomp->f_cutoff[i] =  0.5 / (iclass->density_sigma[i] * cosh( arguement ) * cosh( arguement ));
			printf("%d: density_switch %lf, density_sigma %lf, cutoff %lf, u_cutoff %lf, f_cutoff %lf, denom %lf\n", i, iclass->density_sigma[i], iclass->density_switch[i], iclass->cutoff, icomp->u_cutoff[i], icomp->f_cutoff[i], icomp->denomenator[i]); fflush(stdout);
		}
	} else if (iclass->class_subtype == 3) {
		for(int i = 0; i < iclass->get_n_defined(); i++) {
			icomp->denomenator[i] = pow(iclass->cutoff, 4.0);
			icomp->u_cutoff[i] = 0.0;
			icomp->f_cutoff[i] =  0.0;
			printf("%d: cutoff %lf, u_cutoff %lf, f_cutoff %lf, denom %lf\n", i, iclass->cutoff, icomp->u_cutoff[i], icomp->f_cutoff[i], icomp->denomenator[i]); fflush(stdout);
		}
	} else if (iclass->class_subtype == 4) {
		icomp->c0 = new double[iclass->get_n_defined()]();
		icomp->c2 = new double[iclass->get_n_defined()]();
		icomp->c4 = new double[iclass->get_n_defined()]();
		icomp->c6 = new double[iclass->get_n_defined()]();
		double cutsq = iclass->cutoff * iclass->cutoff;
		for(int ii = 0; ii < iclass->get_n_defined(); ii++) {
			double x = iclass->density_sigma[ii] * iclass->density_sigma[ii] / (iclass->cutoff * iclass->cutoff);
			icomp->denomenator[ii] = (1 - x)*(1 - x)*(1 - x);
			icomp->u_cutoff[ii] = 0.0;
			icomp->f_cutoff[ii] =  0.0;
			
			icomp->c0[ii] = (1.0 - 3.0 * x)/ icomp->denomenator[ii];
			icomp->c2[ii] = 6.0 * x / (cutsq * icomp->denomenator[ii]);
			icomp->c4[ii] = 3 * (1.0 + x) / (cutsq * cutsq * icomp->denomenator[ii]);
			icomp->c6[ii] = 2.0 / (cutsq * cutsq * cutsq * icomp->denomenator[ii]);
		}
	} else {
		printf("Set-up called for density_interactions with invalid class_subtype %d.\n", iclass->class_subtype);
		fflush(stdout);
		exit(EXIT_FAILURE);
	}
	
	setup_site_to_density_group_index_for_range(iclass);
}

void setup_site_to_density_group_index_for_range(DensityClassSpec* iclass) 
{
	// The site_to_density_group_intrxn_index_map array indicates which pairs of sites interact for quick screening during the calculation of density.
	// The first density_group specifies where the density is calculated (at which CG sites), 
	// while the second density_group specifies which what the density is calculated of.
	// So, the density of sites belonging to the second density_group are calculated at every site belonging to the first density_group.
	
	// The value stored in site_to_density_intrxn_index_map is the sum of all "bits" (specified by the density_groups interacting) that exist between these two types.
	// Each "bit" is calculated by left-shifting the corresponding number of bits and then adding to the total.

	if(iclass->get_n_defined() <= 0) return;
	
	// Allocate the array.
	iclass->site_to_density_group_intrxn_index_map = new unsigned long [iclass->n_cg_types * iclass->n_cg_types]();
	
	// Look through types to determine which types belong to a given group.	
	// First determine which density groups each type is part of.
	
	for(int type1 = 0; type1 < iclass->n_cg_types; type1++) {
		for(int dg1 = 0; dg1 < iclass->n_density_groups; dg1++) {
			if(iclass->density_groups[dg1 * iclass->n_cg_types + type1] == false) continue;
			// This CG site type (type1) belongs to this density_group (dg1)
			// Now, look through types again to determine if the corresponding density group has an interaction with this density group.
			for(int type2 = type1; type2 < iclass->n_cg_types; type2++) {
				// Only need to look through half of the type1/type2 combinations since we can check both 1/2 and 2/1 at the same time.
				// Determine which density groups type2 belongs to.
				for(int dg2 = 0; dg2 < iclass->n_density_groups; dg2++) {
					if(iclass->density_groups[dg2 * iclass->n_cg_types + type2] == false) continue;
					// This CG site type (type2) belongs to this density_group (dg2).
					
					// Now, assume (for rangefinding) that these density groups interact with ordering type1/type2 and type2/type1.
					iclass->site_to_density_group_intrxn_index_map[type1 * iclass->n_cg_types + type2] |= 1 << (dg1 * iclass->n_density_groups + dg2);
					iclass->site_to_density_group_intrxn_index_map[type2 * iclass->n_cg_types + type1] |= 1 << (dg2 * iclass->n_density_groups + dg1);
				}
			}
		}
	}
}

//--------------------------------------------------------------------------

void calc_isotropic_two_body_sampling_range(InteractionClassComputer* const icomp, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    int particle_ids[2] = {icomp->k, icomp->l};
    double param;
    calc_distance(particle_ids, x, simulation_box_half_lengths, param);

    if (icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] > param) icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] = param;
    if (icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] < param) icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] = param;
	
	if (icomp->ispec->output_parameter_distribution == 1) {
		if (icomp->ispec->class_type == kAngularBonded || icomp->ispec->class_type == kDihedralBonded) fprintf(icomp->ispec->output_range_file_handles[icomp->index_among_defined_intrxns], "%lf\n", param);
		else if (param < icomp->ispec->cutoff) 			fprintf(icomp->ispec->output_range_file_handles[icomp->index_among_defined_intrxns], "%lf\n", param);
	}
}

void calc_angular_three_body_sampling_range(InteractionClassComputer* const icomp, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    int particle_ids[3] = {icomp->k, icomp->l, icomp->j}; // end indices (k, l) followed by center index (j)
    double param;
    calc_angle(particle_ids, x, simulation_box_half_lengths, param);

    if (icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] > param) icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] = param;
    if (icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] < param) icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] = param;
	
	if (icomp->ispec->output_parameter_distribution == 1) fprintf(icomp->ispec->output_range_file_handles[icomp->index_among_defined_intrxns], "%lf\n", param);
}

void calc_dihedral_four_body_interaction_sampling_range(InteractionClassComputer* const icomp, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
	if (mat->position_dimension != 3) {
    	printf("Dihedral calculations are currently only implemented for 3-dimensional systems.\n");
    	exit(EXIT_FAILURE);
    }

    int particle_ids[4] = {icomp->k, icomp->l, icomp->i, icomp->j}; // end indices (k, l) followed by central bond indices (i, j)
    double param;
    calc_dihedral(particle_ids, x, simulation_box_half_lengths, param);

    if (icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] > param) icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] = param;
    if (icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] < param) icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] = param;
	
	if (icomp->ispec->output_parameter_distribution == 1) fprintf(icomp->ispec->output_range_file_handles[icomp->index_among_defined_intrxns], "%lf\n", param);
}

void calc_radius_of_gyration_interaction_sampling_range(InteractionClassComputer* const icomp, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
	RadiusofGyrationClassSpec* rg_spec = static_cast<RadiusofGyrationClassSpec*>(icomp->ispec);
	int mol_id = icomp->k;
    double param;
    
    int n_ids = rg_spec->topo_data_->molecule_list->partner_numbers_[mol_id];
	int* particle_ids = new int[n_ids];
	for (int i = 0; i < n_ids; i++) particle_ids[i] = (int)(rg_spec->topo_data_->molecule_list->partners_[mol_id][i]);
	
    calc_radius_of_gyration(particle_ids, x, simulation_box_half_lengths, n_ids, param);

    if (icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] > param) icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] = param;
    if (icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] < param) icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] = param;
	
	if (icomp->ispec->output_parameter_distribution == 1) fprintf(icomp->ispec->output_range_file_handles[icomp->index_among_defined_intrxns], "%lf\n", param);

	delete [] particle_ids;
}

void evaluate_density_sampling_range(InteractionClassComputer* const info, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
	DensityClassComputer* icomp = static_cast<DensityClassComputer*>(info);
	DensityClassSpec* ispec = static_cast<DensityClassSpec*>(icomp->ispec);	
	int contributing_density_group = icomp->index_among_defined_intrxns % ispec->n_density_groups;
	double param = icomp->density_values[contributing_density_group * ispec->n_density_groups + icomp->k];
	
	if (icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] > param) icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] = param;
    if (icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] < param) icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] = param;
	
	if (icomp->ispec->output_parameter_distribution == 1) fprintf(icomp->ispec->output_range_file_handles[icomp->index_among_defined_intrxns], "%lf\n", param);
}

void calc_nothing(InteractionClassComputer* const icomp, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat) {
}

//--------------------------------------------------------------------------

void write_range_files(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat)
{
 	FILE* one_body_interaction_output_file_handle = NULL;
 	if (cg->one_body_interactions.class_subtype != 0) {
 		one_body_interaction_output_file_handle = open_file("rmin_1.in", "w");
 	}
 	
    FILE* nonbonded_interaction_output_file_handle = open_file("rmin.in", "w");
    FILE* bonded_interaction_output_file_handle = open_file("rmin_b.in", "w");
    FILE* density_interaction_output_file_handle;
    FILE* radius_of_gyration_interaction_output_file_handle;
    if (cg->density_interactions.class_subtype > 0) density_interaction_output_file_handle = open_file("rmin_den.in", "w");
	if (cg->radius_of_gyration_interactions.class_subtype > 0) radius_of_gyration_interaction_output_file_handle = open_file("rmin_rg.in", "w");
	
    write_interaction_range_data_to_file(cg, mat, one_body_interaction_output_file_handle, nonbonded_interaction_output_file_handle, bonded_interaction_output_file_handle, density_interaction_output_file_handle, radius_of_gyration_interaction_output_file_handle);
    
    if (cg->one_body_interactions.class_subtype != 0) {
 		fclose(one_body_interaction_output_file_handle);
 	}
 	fclose(nonbonded_interaction_output_file_handle);
    fclose(bonded_interaction_output_file_handle);
	if (cg->density_interactions.class_subtype > 0) fclose(density_interaction_output_file_handle);
	if (cg->radius_of_gyration_interactions.class_subtype > 0) fclose(radius_of_gyration_interaction_output_file_handle);
}

void write_interaction_range_data_to_file(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FILE* const one_body_spline_output_filep, FILE* const nonbonded_spline_output_filep, FILE* const bonded_spline_output_filep, FILE* const density_interaction_output_filep, FILE* const radius_of_gyration_interaction_output_filep)
{   
	std::list<InteractionClassSpec*>::iterator iclass_iterator;
	std::list<InteractionClassComputer*>::iterator icomp_iterator;
    for(iclass_iterator = cg->iclass_list.begin(), icomp_iterator = cg->icomp_list.begin(); (iclass_iterator != cg->iclass_list.end()) && (icomp_iterator != cg->icomp_list.end()); iclass_iterator++, icomp_iterator++) {
        if((*iclass_iterator)->class_type == kOneBody) {
        	write_one_body_iclass_range_specifications(*icomp_iterator, cg->name, mat, one_body_spline_output_filep);
        } else if ((*iclass_iterator)->class_type == kPairNonbonded) {
            write_iclass_range_specifications(*icomp_iterator, cg->name, mat, nonbonded_spline_output_filep);
        } else if ((*iclass_iterator)->class_type == kRadiusofGyration) {
			write_iclass_range_specifications(*icomp_iterator, cg->radius_of_gyration_interactions.molecule_group_names, mat, radius_of_gyration_interaction_output_filep);
		} else if ((*iclass_iterator)->class_type == kDensity) {
			write_iclass_range_specifications(*icomp_iterator, cg->density_interactions.density_group_names, mat, density_interaction_output_filep);
		} else {
            write_iclass_range_specifications(*icomp_iterator, cg->name, mat, bonded_spline_output_filep);
        }
    }
}

void write_iclass_range_specifications(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, FILE* const solution_spline_output_file) 
{
    InteractionClassSpec *iclass = icomp->ispec;
    for (int i = 0; i < iclass->get_n_defined(); i++) {
        int index_among_matched_interactions = iclass->defined_to_matched_intrxn_index_map[i];
        if (index_among_matched_interactions > 0) {
            write_single_range_specification(icomp, name, mat, solution_spline_output_file, i);
        }
    }
	
	if (iclass->output_parameter_distribution == 1) {
     	close_parameter_distribution_files_for_class(icomp);
		if(iclass->class_type == kDensity) {
			DensityClassSpec* ispec = static_cast<DensityClassSpec*>(iclass);
			generate_parameter_distribution_histogram(icomp, ispec->density_group_names);
		}
		else generate_parameter_distribution_histogram(icomp, name);
	}
}

void write_one_body_iclass_range_specifications(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, FILE* const solution_spline_output_file) 
{
    InteractionClassSpec *iclass = icomp->ispec;
    if (iclass->class_subtype == 0) return;
    for (int i = 0; i < iclass->get_n_defined(); i++) {
        fprintf(solution_spline_output_file, "%s fm\n", name[i]);
    }
}

void write_single_range_specification(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, FILE* const solution_spline_output_file, const int index_among_defined)
{
    InteractionClassSpec* ispec = icomp->ispec;
    std::string basename;
	
	if (ispec->class_type == kDensity) {
		DensityClassSpec* iclass = static_cast<DensityClassSpec*>(ispec);
		basename = iclass->get_interaction_name(iclass->density_group_names, index_among_defined, " ");
	} else if (ispec->class_type == kRadiusofGyration) {
		RadiusofGyrationClassSpec* iclass = static_cast<RadiusofGyrationClassSpec*>(ispec);
		basename = iclass->get_interaction_name(iclass->molecule_group_names, index_among_defined, " ");
	} else {
		basename = ispec->get_interaction_name(name, index_among_defined, " ");
	}
	
	fprintf(solution_spline_output_file, "%s ", basename.c_str());

    if (fabs(ispec->upper_cutoffs[index_among_defined] + VERYLARGE) < VERYSMALL_F) {
        ispec->upper_cutoffs[index_among_defined] = -1.0;
        ispec->lower_cutoffs[index_among_defined] = -1.0;
    } else if (ispec->class_type == kPairNonbonded) {
        if (ispec->lower_cutoffs[index_among_defined] > ispec->cutoff) {
            ispec->upper_cutoffs[index_among_defined] = -1.0;
            ispec->lower_cutoffs[index_among_defined] = -1.0;
        } else if (ispec->upper_cutoffs[index_among_defined] > ispec->cutoff) {
            ispec->upper_cutoffs[index_among_defined] = ispec->cutoff;
        }
    }
    fprintf(solution_spline_output_file, "%lf %lf fm", ispec->lower_cutoffs[index_among_defined], ispec->upper_cutoffs[index_among_defined]);
	
	DensityClassSpec* dspec = dynamic_cast<DensityClassSpec*>(ispec);
	if(dspec != NULL) {
		if(dspec->class_subtype == 1 || dspec->class_subtype == 4) fprintf(solution_spline_output_file, " %lf", dspec->density_sigma[index_among_defined]);
		if(dspec->class_subtype == 2) fprintf(solution_spline_output_file, " %lf %lf", dspec->density_sigma[index_among_defined], dspec->density_switch[index_among_defined]);
	}
	fprintf(solution_spline_output_file, "\n");
}

void read_density_parameter_file(DensityClassSpec* const ispec) 
{
	int num_elements;
	std::string line;
	std::string* elements = new std::string[5];
	std::ifstream prm_stream;
	prm_stream.open("den.prm", std::ifstream::in);
	if (prm_stream.fail()) {
		printf("Problem opening den.prm file!\n");
		exit(EXIT_FAILURE);
	}

	for(int i = 0; i < ispec->get_n_defined(); i++) {
		if(!std::getline(prm_stream, line)) {
			printf("More lines expected in den.prm!\n");
			exit(EXIT_FAILURE);
		}
		
		num_elements = StringSplit(line, " \t", elements);
		if(num_elements < 3) {
			printf("Each line needs to have at least 3 elements!\n");
			exit(EXIT_FAILURE);
		}
		ispec->density_sigma[i] = atof(elements[2].c_str());
		if(num_elements > 3) {
			ispec->density_switch[i] = atof(elements[3].c_str());
		}
	}
	prm_stream.close();
	delete [] elements;
}

void open_parameter_distribution_files_for_class(InteractionClassComputer* const icomp, char **name) 
{
    InteractionClassSpec* ispec = icomp->ispec;
	std::string filename;
	ispec->output_range_file_handles = new FILE*[ispec->get_n_defined()];
	
	for (int i = 0; i < ispec->get_n_defined(); i++) {
	 	filename = ispec->get_interaction_name(name, i, "_");
	 	if (!ispec->get_short_name().empty()) {
        	filename += "_" + ispec->get_short_name();
    	}
	 	filename += ".dist";
	 	ispec->output_range_file_handles[i] = open_file(filename.c_str(), "w");
	}		
}

void open_density_parameter_distribution_files_for_class(InteractionClassComputer* const icomp, char **name) 
{
	DensityClassSpec* ispec = static_cast<DensityClassSpec*>(icomp->ispec);
	std::string filename;
	ispec->output_range_file_handles = new FILE*[ispec->get_n_defined()];
	
	for (int i = 0; i < ispec->get_n_defined(); i++) {
	 	filename = ispec->get_interaction_name(ispec->density_group_names, i, "_");
	 	if (!ispec->get_short_name().empty()) {
        	filename += "_" + ispec->get_short_name();
    	}
    	filename += ".dist";
    	ispec->output_range_file_handles[i] = open_file(filename.c_str(), "w");
	}		
}

void close_parameter_distribution_files_for_class(InteractionClassComputer* const icomp) 
{
    InteractionClassSpec* ispec = icomp->ispec;	
	for (int i = 0; i < ispec->get_n_defined(); i++) {
		fclose(ispec->output_range_file_handles[i]);
	}
	delete [] ispec->output_range_file_handles;
}

void generate_parameter_distribution_histogram(InteractionClassComputer* const icomp, char **name)
{
	printf("Generating parameter distribution histogram for %s interactions.\n", icomp->ispec->get_full_name().c_str());
    InteractionClassSpec* ispec = icomp->ispec;	
	DensityClassSpec* dspec = dynamic_cast<DensityClassSpec*>(ispec);
	
	std::string filename;
	std::ifstream dist_stream;
	std::ofstream hist_stream;
	int num_bins = 0;
	int	curr_bin;
	double value; 
	double* bin_centers;
	unsigned long* bin_counts;
	for (int i = 0; i < ispec->get_n_defined(); i++) {

		// Set-up histogram based on interaction binwidth
		num_bins = (int)(ceil((ispec->upper_cutoffs[i] - ispec->lower_cutoffs[i]) / ispec->get_fm_binwidth() + 0.5));
		bin_centers = new double[num_bins]();
        bin_counts = new unsigned long[num_bins]();
        
        bin_centers[0] = ispec->lower_cutoffs[i] + 0.5 * ispec->get_fm_binwidth();
        for (int j = 1; j < num_bins; j++) {
        	bin_centers[j] = bin_centers[j - 1] + ispec->get_fm_binwidth();
        }
		
		// Open distribution file
	 	if (dspec != NULL) {
			filename = dspec->get_interaction_name(dspec->density_group_names, i, "_");
		} else {
			filename = ispec->get_interaction_name(name, i, "_");
		}
		if (!ispec->get_short_name().empty()) {
        	filename += "_" + ispec->get_short_name();
    	}
    	filename += ".dist";
		
		check_and_open_in_stream(dist_stream, filename.c_str()); 
		
		// Populate histogram by reading distribution file
		dist_stream >> value;
		while (!dist_stream.fail()) {
			curr_bin = (int)(floor((value - ispec->lower_cutoffs[i] + 0.00001) / ispec->get_fm_binwidth()));	
			if( (curr_bin < num_bins) && (curr_bin >= 0) ) {
				bin_counts[curr_bin]++;
			} else {
				printf("Warning: Bin %d is out-of-bounds. Array size: %d\n", curr_bin, num_bins);
				fflush(stdout);
			}
			dist_stream >> value;
		}

		// Write histogram to file
		if (dspec != NULL) {
			filename = dspec->get_interaction_name(dspec->density_group_names, i, "_");
		} else {
			filename = ispec->get_interaction_name(name, i, "_");
		}
		if (!ispec->get_short_name().empty()) {
        	filename += "_" + ispec->get_short_name();
    	}
    	filename += ".hist";
    
		hist_stream.open(filename, std::ofstream::out);
		
		hist_stream << "#center\tcounts\n";
		for (int j = 0; j < num_bins; j++) {
			hist_stream << bin_centers[j] << "\t" << bin_counts[j] << "\n";
		}
		
		// Close files
		hist_stream.close();
		dist_stream.close();
		delete [] bin_centers;
		delete [] bin_counts;
	}
}

void calculate_BI(CG_MODEL_DATA* const cg, MATRIX_DATA* mat, FrameSource* const fs)
{
  initialize_BI_matrix(mat, cg);

  double volume = calculate_volume(fs->simulation_box_limits);
  
  read_interaction_file_and_build_matrix(mat, cg, volume);
  mat->finish_fm(mat);
}

double calculate_volume(const matrix simulation_box_lengths)
{
  double volume = 1.0;
  int i;
  for(i = 0;i < DIMENSION; i++)
    {
      volume *= simulation_box_lengths[i][i];
    }
  return volume;
}
  

void read_interaction_file_and_build_matrix(MATRIX_DATA* mat, CG_MODEL_DATA* const cg, double volume)
{ 
  int counter = 0;
  int* sitecounter;
  sitecounter = new int[cg->n_cg_types]();
  std::list<InteractionClassComputer*>::iterator icomp_iterator;
  for(int j = 0; j < cg->n_cg_sites; j++){
    int type = cg->topo_data.cg_site_types[j];
    sitecounter[type-1]++;
  }

  for(icomp_iterator = cg->icomp_list.begin(); icomp_iterator != cg->icomp_list.end(); icomp_iterator++) {
    // For every defined interaction, 
    // if that interaction has a parameter distribution
    if ( (*icomp_iterator)->ispec->output_parameter_distribution != 1) continue;
    
    // These interactions do not generate parameter distributions
    if( (*icomp_iterator)->ispec->class_type == kOneBody|| (*icomp_iterator)->ispec->class_type == kThreeBodyNonbonded )
      {continue;}
    
    // Otherwise, process the data
    (*icomp_iterator)->table_basis_fn_vals.resize((*icomp_iterator)->ispec->get_bspline_k());
    for (unsigned i = 0; i < (*icomp_iterator)->ispec->defined_to_matched_intrxn_index_map.size(); i++) {
      if( (*icomp_iterator)->ispec->output_parameter_distribution != 1) continue;
      if( (*icomp_iterator)->ispec->class_type == kPairNonbonded ) {
	    std::vector <int> type_vector = (*icomp_iterator)->ispec->get_interaction_types(i);
	    double num_pairs = sitecounter[type_vector[0]-1] * sitecounter[type_vector[1]-1];
	    if( type_vector[0] == type_vector[1]){
	      num_pairs -= sitecounter[type_vector[0]-1];
	    }
        read_one_param_dist_file_pair((*icomp_iterator), cg->name, mat, i, counter,num_pairs, volume);
      } else if ( (*icomp_iterator)->ispec->class_type == kPairBonded ) {
	    double num_bonds = count_bonded_interaction((*icomp_iterator), cg->name, mat, i);
	    //read_one_param_dist_file_pair((*icomp_iterator), cg->name, mat, i, counter, num_bonds, volume);
        read_one_param_dist_file_pair((*icomp_iterator), cg->name, mat, i, counter, 2.0, 1.0);
      } else if( (*icomp_iterator)->ispec->class_type == kDensity ){
        DensityClassSpec* dspec = static_cast<DensityClassSpec*>((*icomp_iterator)->ispec);
	    read_one_param_dist_file_other((*icomp_iterator), dspec->density_group_names, mat, i, counter, 1);
      } else{
	    double num_angles = count_bonded_interaction((*icomp_iterator), cg->name, mat, i);
	    read_one_param_dist_file_other((*icomp_iterator), cg->name, mat, i, counter, num_angles);
      }
    }
  }  
  delete [] sitecounter;
}

void read_one_param_dist_file_pair(InteractionClassComputer* const icomp, char** const name, MATRIX_DATA* mat, const int index_among_defined_intrxns, int &counter, double num_of_pairs, double volume)
{
  std::string basename = icomp->ispec->get_interaction_name(name, index_among_defined_intrxns, "_");
  if (!icomp->ispec->get_short_name().empty())
  {
    basename += "_" + icomp->ispec->get_short_name();
  }
  std::string filename = basename + ".hist";
  FILE* curr_dist_input_file = open_file(filename.c_str(), "r");

  int i, counts;
  int *junk;
  double PI = 3.1415926;
  double r;
  double potential;
  int num_entries = (int)((icomp->ispec->upper_cutoffs[index_among_defined_intrxns] - icomp->ispec->lower_cutoffs[index_among_defined_intrxns])/icomp->ispec->get_fm_binwidth());
  std::array<double, DIMENSION>* derivatives = new std::array<double, DIMENSION>[num_entries - 1];
  char buffer[100];
  fgets(buffer,100,curr_dist_input_file); 
  for(i = 0; i < num_entries; i++)
    {
      int first_nonzero_basis_index;
      fscanf(curr_dist_input_file,"%lf %d\n",&r,&counts);
      double normalized_counts = (double)(counts) / ( 4.0*PI*r*r*(r - icomp->ispec->get_fm_binwidth()) );
      normalized_counts *= 2.0 * mat->normalization * volume / num_of_pairs;
      potential = mat->temperature*mat->boltzmann*log(normalized_counts);
      icomp->fm_s_comp->calculate_basis_fn_vals(index_among_defined_intrxns, r, first_nonzero_basis_index, icomp->table_basis_fn_vals);
      mat->accumulate_matching_forces(icomp, first_nonzero_basis_index, icomp->table_basis_fn_vals, counter, junk, derivatives, mat);
      mat->accumulate_target_force_element(mat, counter, &potential);
      counter++;
    }
  delete [] derivatives;
}

void read_one_param_dist_file_other(InteractionClassComputer* const icomp, char** const name, MATRIX_DATA* mat, const int index_among_defined_intrxns, int &counter, double num_of_pairs)
{
  std::string basename = icomp->ispec->get_interaction_name(name, index_among_defined_intrxns, "_");
  // Append the class's short name with an underscore if one exists.
  if (!icomp->ispec->get_short_name().empty())
  {
    basename += "_" + icomp->ispec->get_short_name();
  }
  std::string filename = basename + ".hist";
  FILE* curr_dist_input_file = open_file(filename.c_str(), "r");

  int i, counts;
  int *junk;
  double r;
  double potential;
  int num_entries = (int)((icomp->ispec->upper_cutoffs[index_among_defined_intrxns] - icomp->ispec->lower_cutoffs[index_among_defined_intrxns])/icomp->ispec->get_fm_binwidth());
  std::array<double, DIMENSION>* derivatives = new std::array<double, DIMENSION>[num_entries - 1];
  char buffer[100];
    fgets(buffer,100,curr_dist_input_file); 
  for(i = 0; i < num_entries; i++)
    {
      int first_nonzero_basis_index;
      fscanf(curr_dist_input_file,"%lf %d\n",&r,&counts);
      double normalized_counts = (double)(counts)*mat->normalization;
      normalized_counts *= 2.0;
      normalized_counts *= (1.0/num_of_pairs);
      potential = mat->temperature*mat->boltzmann*log(normalized_counts);
      icomp->fm_s_comp->calculate_basis_fn_vals(index_among_defined_intrxns, r, first_nonzero_basis_index, icomp->table_basis_fn_vals);
      mat->accumulate_matching_forces(icomp, first_nonzero_basis_index, icomp->table_basis_fn_vals, counter, junk, derivatives, mat);
      mat->accumulate_target_force_element(mat, counter, &normalized_counts);
      counter++;
    }
}

double count_bonded_interaction(InteractionClassComputer* const icomp, char** const name, MATRIX_DATA* mat, const int index_among_defined_intrxns)
{
  std::string basename = icomp->ispec->get_interaction_name(name, index_among_defined_intrxns, "_");
  // Append the class's short name with an underscore if one exists.
  if (!icomp->ispec->get_short_name().empty())
  {
    basename += "_" + icomp->ispec->get_short_name();
  }  
  std::string filename = basename + ".hist";
  FILE* curr_dist_input_file = open_file(filename.c_str(), "r");

  int i, counts;
  int count_summ = 0;
  double r;
  char buffer[100];
  fgets(buffer,100,curr_dist_input_file);
  int num_entries = (int)((icomp->ispec->upper_cutoffs[index_among_defined_intrxns] - icomp->ispec->lower_cutoffs[index_among_defined_intrxns])/icomp->ispec->get_fm_binwidth());
  for(i = 0; i< num_entries; i++)
    {
      fscanf(curr_dist_input_file,"%lf %d\n",&r,&counts);
      count_summ += counts;
    }
  int num_bonds = (double)( count_summ * mat->normalization );

  return num_bonds;
}

void free_name(CG_MODEL_DATA* const cg)
{
    // Free data after output.
    for (int i = 0; i < cg->n_cg_types; i++) delete [] cg->name[i];
    delete [] cg->name;
}
