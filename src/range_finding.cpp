//
//  range_finding.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cmath>
#include <cstdio>
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

// Helper function that issues failure warnings
void report_unrecognized_class_subtype(InteractionClassSpec *iclass);

// Functions for computing the full range of sampling of a given class of interaction in a given trajectory.

void calc_isotropic_two_body_sampling_range(InteractionClassComputer* const icomp, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_angular_three_body_sampling_range(InteractionClassComputer* const icomp, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_dihedral_four_body_interaction_sampling_range(InteractionClassComputer* const icomp, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
void calc_nothing(InteractionClassComputer* const icomp, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);

void write_interaction_range_data_to_file_and_free_it(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FILE* const nonbonded_spline_output_filep, FILE* const bonded_spline_output_filep);

void write_iclass_range_specifications(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, FILE* const solution_spline_output_file);
void write_single_range_specification(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, FILE* const solution_spline_output_file, const int index_among_defined);

void open_parameter_distribution_files_for_class(InteractionClassComputer* const icomp, char **name); 
void close_parameter_distribution_files_for_class(InteractionClassComputer* const icomp);
void generate_parameter_distribution_histogram(InteractionClassComputer* const icomp, char **name);

// Dummy implementations
void do_not_initialize_fm_matrix(MATRIX_DATA* const mat);

void do_not_initialize_fm_matrix(MATRIX_DATA* const mat) {}

//------------------------------------------------------------------------
//    Implementation
//------------------------------------------------------------------------

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
    cg->pair_nonbonded_cutoff2 = VERYLARGE * VERYLARGE;
}

void initialize_single_class_range_finding_temps(InteractionClassSpec *iclass, InteractionClassComputer *icomp, TopologyData *topo_data) 
{
	iclass->setup_for_defined_interactions(topo_data);
	
    icomp->ispec = iclass;
    if (iclass->class_type == kPairNonbonded) {
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
    } else { // This is implicitly for three-body interactions since no other types are defined.
    	icomp->calculate_fm_matrix_elements = calc_nothing;
    }

    initialize_ranges(iclass->get_n_defined(), iclass->lower_cutoffs, iclass->upper_cutoffs, iclass->defined_to_matched_intrxn_index_map);
    iclass->n_to_force_match = iclass->get_n_defined();
    iclass->interaction_column_indices = std::vector<unsigned>(iclass->n_to_force_match + 1);
    
    if(iclass->output_parameter_distribution == 1) open_parameter_distribution_files_for_class(icomp, topo_data->name);
}

void report_unrecognized_class_subtype(InteractionClassSpec *iclass)
{
	printf("Unrecognized %s class subtype!\n", iclass->get_full_name().c_str());
	fflush(stdout);
	exit(EXIT_FAILURE);
}

//--------------------------------------------------------------------------

void calc_isotropic_two_body_sampling_range(InteractionClassComputer* const icomp, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    std::array<int, 2> particle_ids = {{icomp->k, icomp->l}};
    double param;
    calc_distance(particle_ids, x, simulation_box_half_lengths, param);

    if (icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] > param) icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] = param;
    if (icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] < param) icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] = param;
    
    if ((icomp->ispec->output_parameter_distribution == 1) && (param < icomp->ispec->cutoff)) fprintf(icomp->ispec->output_range_file_handles[icomp->index_among_defined_intrxns], "%lf\n", param);
}

void calc_angular_three_body_sampling_range(InteractionClassComputer* const icomp, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    std::array<int, 3> particle_ids = {{icomp->j, icomp->k, icomp->l}};
    double param;
    calc_angle(particle_ids, x, simulation_box_half_lengths, param);

    if (icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] > param) icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] = param;
    if (icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] < param) icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] = param;

	if (icomp->ispec->output_parameter_distribution == 1) fprintf(icomp->ispec->output_range_file_handles[icomp->index_among_defined_intrxns], "%lf\n", param);
}

void calc_dihedral_four_body_interaction_sampling_range(InteractionClassComputer* const icomp, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat)
{
    std::array<int, 4> particle_ids = {{icomp->i, icomp->j, icomp->k, icomp->l}};
    double param;
    calc_dihedral(particle_ids, x, simulation_box_half_lengths, param);

    if (icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] > param) icomp->ispec->lower_cutoffs[icomp->index_among_defined_intrxns] = param;
    if (icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] < param) icomp->ispec->upper_cutoffs[icomp->index_among_defined_intrxns] = param;
    
    if (icomp->ispec->output_parameter_distribution == 1) fprintf(icomp->ispec->output_range_file_handles[icomp->index_among_defined_intrxns], "%lf\n", param);
}

void calc_nothing(InteractionClassComputer* const icomp, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat) {
}

//--------------------------------------------------------------------------

void write_range_files(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat)
{
    FILE* nonbonded_interaction_output_file_handle = open_file("rmin.in", "w");
    FILE* bonded_interaction_output_file_handle = open_file("rmin_b.in", "w");
    
    write_interaction_range_data_to_file_and_free_it(cg, mat, nonbonded_interaction_output_file_handle, bonded_interaction_output_file_handle);
    
    fclose(nonbonded_interaction_output_file_handle);
    fclose(bonded_interaction_output_file_handle);
}

void write_interaction_range_data_to_file_and_free_it(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FILE* const nonbonded_spline_output_filep, FILE* const bonded_spline_output_filep)
{    
	std::list<InteractionClassSpec*>::iterator iclass_iterator;
	std::list<InteractionClassComputer*>::iterator icomp_iterator;
    for(iclass_iterator = cg->iclass_list.begin(), icomp_iterator = cg->icomp_list.begin(); (iclass_iterator != cg->iclass_list.end()) && (icomp_iterator != cg->icomp_list.end()); iclass_iterator++, icomp_iterator++) {
        if ((*iclass_iterator)->class_type == kPairNonbonded) {
            write_iclass_range_specifications(*icomp_iterator, cg->name, mat, nonbonded_spline_output_filep);
        } else {
            write_iclass_range_specifications(*icomp_iterator, cg->name, mat, bonded_spline_output_filep);
        }
    }

    // Free data after output.
    printf("Done with output.\n"); fflush(stdout);
    for (int i = 0; i < cg->n_cg_types; i++) delete [] cg->name[i];
    delete [] cg->name;
}

void write_iclass_range_specifications(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, FILE* const solution_spline_output_file) {
    InteractionClassSpec *iclass = icomp->ispec;
    for (int i = 0; i < iclass->n_defined; i++) {
        int index_among_matched_interactions = iclass->defined_to_matched_intrxn_index_map[i];
        if (index_among_matched_interactions > 0) {
            write_single_range_specification(icomp, name, mat, solution_spline_output_file, i);
        }
    }
    if (iclass->output_parameter_distribution == 1) {
    	close_parameter_distribution_files_for_class(icomp);
		generate_parameter_distribution_histogram(icomp, name);
	}
}

void write_single_range_specification(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, FILE* const solution_spline_output_file, const int index_among_defined)
{
    InteractionClassSpec* ispec = icomp->ispec;
    std::vector<int> types = ispec->get_interaction_types(index_among_defined);
    for (unsigned i = 0; i < types.size(); i++) {
        fprintf(solution_spline_output_file, "%s ", name[types[i] - 1]);
    }

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
    fprintf(solution_spline_output_file, "%lf %lf fm\n", ispec->lower_cutoffs[index_among_defined], ispec->upper_cutoffs[index_among_defined]);
}

void open_parameter_distribution_files_for_class(InteractionClassComputer* const icomp, char **name) 
{
    InteractionClassSpec* ispec = icomp->ispec;	
	std::string filename;
	ispec->output_range_file_handles = new FILE*[ispec->get_n_defined()];
	
	for (int i = 0; i < ispec->get_n_defined(); i++) {
	 	filename = ispec->get_interaction_name(name, i) + ".dist";
		ispec->output_range_file_handles[i] = fopen(filename.c_str(), "w");
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
    InteractionClassSpec* ispec = icomp->ispec;	
	std::string filename;
	FILE *dist_fh, *hist_fh;
	int num_bins, curr_bin;
	double value; 
	double* bin_centers;
	unsigned long* bin_counts;
	for (int i = 0; i < ispec->get_n_defined(); i++) {

		// Set-up histogram based on interaction binwidth
		num_bins = (int)(ceil((ispec->upper_cutoffs[i] - ispec->lower_cutoffs[i]) / ispec->get_fm_binwidth()) + 0.5);
        bin_centers = new double[num_bins]();
        bin_counts = new unsigned long[num_bins]();
        
        bin_centers[0] = ispec->lower_cutoffs[i] + 0.5 * ispec->get_fm_binwidth();
        for (int j = 1; j < num_bins; j++) {
        	bin_centers[j] = bin_centers[j - 1] + ispec->get_fm_binwidth();
        }
		
		// Open distribution file
	 	filename = ispec->get_interaction_name(name, i) + ".dist";
	 	dist_fh = fopen(filename.c_str(), "r");

		// Populate histogram by reading distribution file
		while (fscanf(dist_fh, "%lf", &value) == 1) {
			curr_bin = (int)(floor((value - ispec->lower_cutoffs[i] + 0.00001) / ispec->get_fm_binwidth()));
			bin_counts[curr_bin]++;
		}

		// Write histogram to file
		filename = ispec->get_interaction_name(name, i) + ".hist";
		hist_fh = fopen(filename.c_str(), "w");
		
		fprintf(hist_fh, "center\tcounts\n");
		for (int j = 0; j < num_bins; j++) {
			fprintf(hist_fh, "%lf\t%lu\n", bin_centers[j], bin_counts[j]);
		}
		
		// Close files
		fclose(hist_fh);
		fclose(dist_fh);
		delete [] bin_centers;
		delete [] bin_counts;
	}
}

//--------------------------------------------------------------------------

void free_fm_sampling_range_calculation_temps(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat)
{
    // Free topology information.
    free_topology_data(&cg->topo_data);
    
    // Free calculation intermediate.
    delete [] mat->dense_fm_rhs_vector;
	delete [] mat->dense_fm_normal_rhs_vector;
}
