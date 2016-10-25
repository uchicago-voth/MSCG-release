//
//  splines.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cassert>
#include "splines.h"
#include "interaction_model.h"

SplineComputer* set_up_fm_spline_comp(InteractionClassSpec *ispec)
{
    if (ispec->n_to_force_match > 0) {
        if (ispec->get_basis_type() == kBSpline) {
            return new BSplineComputer(ispec);
        } else if (ispec->get_basis_type() == kLinearSpline) {
            return new LinearSplineComputer(ispec);
        } else if (ispec->get_basis_type() == kBSplineAndDeriv) {
        	return new BSplineAndDerivComputer(ispec);
        } else {
            printf("Unrecognized spline class.\n");
            exit(EXIT_FAILURE);
        }
    } else {
        return NULL;
    }
}

SplineComputer* set_up_table_spline_comp(InteractionClassSpec *ispec)
{
    if (ispec->n_tabulated > 0) {
        return new TableSplineComputer(ispec);    
    } else {
        return NULL;
    }
}


// This routine returns the input value minus the lower cutoff, rounded
// inside the allowed range adjusted by a small rounding factor.
double SplineComputer::get_param_less_lower_cutoff(const int index_among_defined, const double param_val) const 
{
    double cutoff_range = ispec_->upper_cutoffs[index_among_defined] - ispec_->lower_cutoffs[index_among_defined] - VERYSMALL;
    double param_less_lower_cutoff = param_val - ispec_->lower_cutoffs[index_among_defined];
    if (param_less_lower_cutoff < 0.0) param_less_lower_cutoff = 0.0;
    else if (param_less_lower_cutoff > cutoff_range) param_less_lower_cutoff = cutoff_range;
    return param_less_lower_cutoff;
}


BSplineComputer::BSplineComputer(InteractionClassSpec* ispec) : SplineComputer(ispec)
{
    int interaction_column_indices, n_to_print_minus_bspline_k;
    n_coef = ispec_->get_bspline_k();
    n_to_force_match = ispec_->n_to_force_match;
    n_defined = ispec_->get_n_defined();
    binwidth = ispec_->get_fm_binwidth();

    printf("Allocating b-spline temporaries for %d interactions.\n", n_to_force_match);
    bspline_workspaces = new gsl_bspline_workspace*[n_to_force_match];
    bspline_vectors = gsl_vector_alloc(n_coef);

    int counter = 0;
    for (unsigned i = 0; i < n_defined; i++) {
        if (ispec_->defined_to_matched_intrxn_index_map[i] > 0) {
            interaction_column_indices = ispec_->interaction_column_indices[counter + 1] - ispec_->interaction_column_indices[counter];
            n_to_print_minus_bspline_k = interaction_column_indices - n_coef + 2;
            bspline_workspaces[counter] = gsl_bspline_alloc(n_coef, n_to_print_minus_bspline_k);
            gsl_bspline_knots_uniform(ispec_->lower_cutoffs[i], ispec_->upper_cutoffs[i], bspline_workspaces[counter]);
            counter++;
        }
    }
}

BSplineComputer::~BSplineComputer()
{
    for (unsigned i = 0; i < n_to_force_match; i++) {
        gsl_bspline_free(bspline_workspaces[i]);
    }
    delete [] bspline_workspaces;
    gsl_vector_free(bspline_vectors);   
};

// Calculate the value of a one-parameter B-spline; direction of the corresponding
// forces is calculated in the function calling this one.

void BSplineComputer::calculate_basis_fn_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals)
{
    assert(vals.size() == n_coef);

    double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    first_nonzero_basis_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth());

    size_t junk;
    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;
    gsl_bspline_eval_nonzero(param_less_lower_cutoff + ispec_->lower_cutoffs[index_among_defined], bspline_vectors, &junk, &junk, bspline_workspaces[index_among_matched]);
    
    for (unsigned i = 0; i < n_coef; i++) {
        vals[i] = gsl_vector_get(bspline_vectors, i);
    }
}

double BSplineComputer::evaluate_spline(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
    size_t istart, iend;
    gsl_bspline_eval_nonzero(axis, bspline_vectors, &istart, &iend, bspline_workspaces[index_among_matched_interactions - 1]);
    double force = 0.0;
    for (int tn = int(istart); tn <= int(iend); tn++) {
        force += gsl_vector_get(bspline_vectors, tn - istart) * spline_coeffs[first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + tn];
    }
    return force;
}

BSplineAndDerivComputer::BSplineAndDerivComputer(InteractionClassSpec* ispec) : SplineComputer(ispec)
{
    int n_to_print_minus_bspline_k, interaction_column_indices;

    n_coef = ispec_->get_bspline_k();
    n_to_force_match = ispec_->n_to_force_match; 
    class_subtype = ispec_->class_subtype;
    n_defined = ispec_->get_n_defined();
    binwidth = ispec_->get_fm_binwidth();
    
    if ( (ispec_->class_type != kThreeBodyNonbonded) || (ispec_->class_subtype == 1) || (ispec_->class_subtype == 2) || (ispec_->class_subtype == 3) ) {
        printf("Allocating b-spline and derivative temporaries for %d interactions.\n", ispec_->get_n_defined());
        bspline_vectors = gsl_vector_alloc(n_coef);
        
       	bspline_matrices = gsl_matrix_alloc(n_coef, 2);
       	
        if (ispec_->class_type == kThreeBodyNonbonded) {
	       	bspline_workspaces = new gsl_bspline_workspace*[n_defined];
     		
     		for (unsigned counter = 0; counter < n_defined; counter++) {
            	bspline_workspaces[counter] = gsl_bspline_alloc(n_coef, n_to_print_minus_bspline_k);
            	gsl_bspline_knots_uniform(ispec_->lower_cutoffs[counter], ispec_->upper_cutoffs[counter], bspline_workspaces[counter]);
        	}
     	} else {
     		bspline_workspaces = new gsl_bspline_workspace*[n_to_force_match];
       	
      		int counter = 0;
       		for (unsigned i = 0; i < n_defined; i++) {
       			if (ispec_->defined_to_matched_intrxn_index_map[i] > 0) {
            		interaction_column_indices = ispec_->interaction_column_indices[counter + 1] - ispec_->interaction_column_indices[counter];
            		n_to_print_minus_bspline_k = interaction_column_indices - n_coef + 2;
            		bspline_workspaces[counter] = gsl_bspline_alloc(n_coef, n_to_print_minus_bspline_k);
            		gsl_bspline_knots_uniform(ispec_->lower_cutoffs[i], ispec_->upper_cutoffs[i], bspline_workspaces[counter]);
            		counter++;
        		}
        	}
        }
    }
}

BSplineAndDerivComputer::~BSplineAndDerivComputer() 
{
    if ((ispec_->class_type != kThreeBodyNonbonded) || (class_subtype == 1) || (class_subtype == 2) || (ispec_->class_subtype == 3)) {
    	if (ispec_->class_type != kThreeBodyNonbonded) {
    		for (unsigned i = 0; i < n_to_force_match; i++)	gsl_bspline_free(bspline_workspaces[i]);
    	} else {
        	for (unsigned i = 0; i < n_defined; i++)  gsl_bspline_free(bspline_workspaces[i]);
        }
        
        delete [] bspline_workspaces;
        gsl_vector_free(bspline_vectors);
        gsl_matrix_free(bspline_matrices);
    }
}

void BSplineAndDerivComputer::calculate_bspline_deriv_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals)
{
    assert(vals.size() == n_coef);

    double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    first_nonzero_basis_index = (int)(param_less_lower_cutoff / binwidth);

    size_t junk;
    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;
    gsl_bspline_deriv_eval_nonzero(param_less_lower_cutoff + ispec_->lower_cutoffs[index_among_defined], (size_t)(1), bspline_matrices, &junk, &junk, bspline_workspaces[index_among_matched]);
    
    for (unsigned i = 0; i < n_coef; i++) {
        vals[i] = -gsl_matrix_get(bspline_matrices, i, 1);
    }
}

void BSplineAndDerivComputer::calculate_basis_fn_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals) 
{
    assert(vals.size() == n_coef);
    
    double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    first_nonzero_basis_index = (int)(param_less_lower_cutoff / binwidth);
    
    size_t junk;
    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;
    gsl_bspline_eval_nonzero(param_less_lower_cutoff + ispec_->lower_cutoffs[index_among_defined], bspline_vectors, &junk, &junk, bspline_workspaces[index_among_matched]);
    
    for (unsigned i = 0; i < n_coef; i++) {
        vals[i] = gsl_vector_get(bspline_vectors, i);
    }
}

double BSplineAndDerivComputer::evaluate_spline(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
    size_t istart, iend;
    gsl_bspline_eval_nonzero(axis, bspline_vectors, &istart, &iend, bspline_workspaces[index_among_matched_interactions - 1]);
    double force = 0.0;
    
    for (int tn = int(istart); tn <= int(iend); tn++) {
    	if (spline_coeffs.size() <= first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + tn) {
	       	printf("gsl_vector_get(bspline_vectors, %d - %d) = %lf\t", tn, (int)(istart), gsl_vector_get(bspline_vectors, tn - istart)); fflush(stdout);
    		printf("spline_coeffs.size() = %u\n", (unsigned)(spline_coeffs.size())); fflush(stdout);
    		printf("index %d = ", first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + tn); fflush(stdout);
			printf("fnzbi %d + ici[%d - 1] %d + tn %d\n", first_nonzero_basis_index, index_among_matched_interactions, ispec_->interaction_column_indices[index_among_matched_interactions - 1], tn); fflush(stdout);
    	}
        force += gsl_vector_get(bspline_vectors, tn - istart) * spline_coeffs[first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + tn];
    }
    return force;
}

double BSplineAndDerivComputer::evaluate_spline_deriv(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
    size_t istart, iend;
    gsl_bspline_deriv_eval_nonzero(axis, size_t(1), bspline_matrices, &istart, &iend, bspline_workspaces[index_among_matched_interactions - 1]);
    double deriv = 0.0;
    for (int tn = int(istart); tn <= int(iend); tn++) {
    	if (spline_coeffs.size() <= first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + tn) {
	       	printf("gsl_matrix_get(bspline_vectors, %d - %d) = %lf\t", tn, (int)(istart), gsl_vector_get(bspline_vectors, tn - istart)); fflush(stdout);
    		printf("spline_coeffs.size() = %u\n", (unsigned)(spline_coeffs.size())); fflush(stdout);
    		printf("index %d = ", first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + tn); fflush(stdout);
			printf("fnzbi %d + ici[%d - 1] %d + tn %d\n", first_nonzero_basis_index, index_among_matched_interactions, ispec_->interaction_column_indices[index_among_matched_interactions - 1], tn); fflush(stdout);
    	}
        deriv += gsl_matrix_get(bspline_matrices, tn - istart, 1) * spline_coeffs[first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + tn];
    }
    return deriv;
}

LinearSplineComputer::LinearSplineComputer(InteractionClassSpec* ispec) : SplineComputer(ispec)
{
    n_coef = 2;
    n_to_force_match = ispec_->n_to_force_match;
    n_defined = ispec_->get_n_defined();
    binwidth = ispec_->get_fm_binwidth();
}

void LinearSplineComputer::calculate_basis_fn_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals)
{
    assert(vals.size() == n_coef);

    double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    
    first_nonzero_basis_index = int(param_less_lower_cutoff / binwidth);
    vals[1] = fmod(param_less_lower_cutoff / binwidth, 1.0);
    vals[0] = 1.0 - vals[1];
}

double LinearSplineComputer::evaluate_spline(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
	double param_less_lower_cutoff = axis - ispec_->lower_cutoffs[index_among_defined];
	int basis_function_column_index = (int)(param_less_lower_cutoff / binwidth);
    double remainder_after_binning = fmod(param_less_lower_cutoff / binwidth, 1.0);
    double force = 0.0;
    if (unsigned(first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + basis_function_column_index + 1) < spline_coeffs.size()) {
        force = spline_coeffs[first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + basis_function_column_index] * (1.0 - remainder_after_binning) + spline_coeffs[first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + basis_function_column_index + 1] * remainder_after_binning;
    } else if (unsigned(first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + basis_function_column_index + 1) == spline_coeffs.size()) {
        force = spline_coeffs[first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + basis_function_column_index];
    } else {
        printf("Warning: attempting to read %d columns past interaction column %d in a matrix of %lu columns, to be multiplied by a coefficient %g.", basis_function_column_index + 1, first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1], spline_coeffs.size(), remainder_after_binning);
    }
    return force;
}

TableSplineComputer::TableSplineComputer(InteractionClassSpec* ispec) : SplineComputer(ispec) 
{
    n_coef = 2;
    n_to_force_match = ispec_->n_to_force_match;
    n_defined = ispec_->get_n_defined();
    binwidth = ispec_->external_table_spline_binwidth;
}

void TableSplineComputer::calculate_basis_fn_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals)
{   
    double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    first_nonzero_basis_index = int(param_less_lower_cutoff / binwidth);
    
    vals[1] = fmod(param_less_lower_cutoff / binwidth, 1.0);
    vals[0] = 1.0 - vals[1];

    int index_among_tabulated_interactions = ispec_->defined_to_tabulated_intrxn_index_map[index_among_defined] - 1;
    vals[0] *= ispec_->external_table_spline_coefficients[index_among_tabulated_interactions][first_nonzero_basis_index];
    vals[1] *= ispec_->external_table_spline_coefficients[index_among_tabulated_interactions][first_nonzero_basis_index + 1];
}

double TableSplineComputer::evaluate_spline(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis)
{
    int dummy_first_index = 0;
    std::vector<double> vals(2);
    calculate_basis_fn_vals(index_among_defined, axis, dummy_first_index, vals);
    return vals[0] + vals[1];
}
