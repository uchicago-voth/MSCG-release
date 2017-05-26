//
//  splines.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cassert>
#include "splines.h"
#include "interaction_model.h"

// Helper function to print a sizing error message.
inline void check_bspline_size(const int control_points, const int order);
inline void check_matching_istart(const int first_nonzero_basis_index, const size_t istart);
inline double check_against_cutoffs(const double axis, const double lower_cutoff, const double upper_cutoff);
inline void check_bspline_sizing(const size_t coeffs_size, const int first_nonzero_basis_index, const int index_among_matched_interactions, const int ici_index, const int tn, const size_t istart);

SplineComputer* set_up_fm_spline_comp(InteractionClassSpec *ispec)
{
    if (ispec->n_to_force_match > 0) {
        if (ispec->get_basis_type() == kBSpline) {
            return new BSplineComputer(ispec);
        } else if (ispec->get_basis_type() == kLinearSpline) {
            return new LinearSplineComputer(ispec);
        } else if (ispec->get_basis_type() == kDelta) {
        	return new DeltaSplineComputer(ispec);
        } else if (ispec->get_basis_type() == kBSplineAndDeriv) {
        	return new BSplineAndDerivComputer(ispec);
        } else if (ispec->get_basis_type() == kNone) {
        	return NULL;
        } else {
            fprintf(stderr, "Unrecognized spline class.\n");
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
    if (param_less_lower_cutoff < 0.0) {
    	fprintf(stderr, "Value passed to spline computer (%lf) is below its lower cutoff (%lf)!\n\n", param_val, ispec_->lower_cutoffs[index_among_defined]);
    	param_less_lower_cutoff = 0.0;
    }
    else if (param_less_lower_cutoff > cutoff_range) {
    	fprintf(stderr, "Value passed to spline computer %lf is above its cutoff %lf!\n\n", param_val, ispec_->upper_cutoffs[index_among_defined]);
	    param_less_lower_cutoff = cutoff_range;
    }
    return param_less_lower_cutoff;
}


BSplineComputer::BSplineComputer(InteractionClassSpec* ispec) : SplineComputer(ispec)
{
    int interaction_column_indices, n_to_print_minus_bspline_k;
    n_coef = ispec_->get_bspline_k();
    n_to_force_match = ispec_->n_to_force_match;
    n_defined = ispec_->get_n_defined();
    binwidth = ispec_->get_fm_binwidth();

    if (n_coef < 2) {
    	printf("Spline order for BSplineAndDeriv must be at least 2!\n");
		exit(EXIT_FAILURE);
    }

    printf("Allocating b-spline temporaries for %d interactions.\n", n_to_force_match);
    bspline_workspaces = new gsl_bspline_workspace*[n_to_force_match];
    bspline_vectors = gsl_vector_alloc(n_coef);

    int counter = 0;
    for (unsigned i = 0; i < n_defined; i++) {
        if (ispec_->defined_to_matched_intrxn_index_map[i] > 0) {
            interaction_column_indices = ispec_->interaction_column_indices[counter + 1] - ispec_->interaction_column_indices[counter];
            n_to_print_minus_bspline_k = interaction_column_indices - n_coef + 2;
            check_bspline_size(n_to_print_minus_bspline_k, (int)(n_coef));
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
	size_t istart, iend;
    double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    first_nonzero_basis_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth());

    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;
    gsl_bspline_eval_nonzero(param_less_lower_cutoff + ispec_->lower_cutoffs[index_among_defined], bspline_vectors, &istart, &iend, bspline_workspaces[index_among_matched]);
    check_matching_istart(first_nonzero_basis_index, istart);
    
    for (unsigned i = 0; i < n_coef; i++) {
        vals[i] = gsl_vector_get(bspline_vectors, i);
    }
}

double BSplineComputer::evaluate_spline(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    size_t istart, iend;
    int ici_index = 0;
    double force = 0.0;
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
    double axis_val = check_against_cutoffs(axis, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    gsl_bspline_eval_nonzero(axis_val, bspline_vectors, &istart, &iend, bspline_workspaces[index_among_matched_interactions - 1]);
    if (index_among_matched_interactions > 0) {
		ici_index = ispec_->interaction_column_indices[index_among_matched_interactions - 1];
    }
    for (int tn = int(istart); tn <= int(iend); tn++) {
        check_bspline_sizing(spline_coeffs.size(), first_nonzero_basis_index, index_among_matched_interactions, ici_index, tn, istart);
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
    
    if (n_coef < 3) {
    	printf("Spline order for BSplineAndDeriv must be at least 3!\n");
		exit(EXIT_FAILURE);
    }
    
    // Exclude kThreeBody class_subtype == 0 only
    if ( (ispec_->class_type != kThreeBodyNonbonded) || (ispec_->class_subtype == 1) || (ispec_->class_subtype == 2) || (ispec_->class_subtype == 3) ) {
        printf("Allocating b-spline and derivative temporaries for %d interactions.\n", ispec_->get_n_defined());
        bspline_vectors = gsl_vector_alloc(n_coef);
        bspline_matrices = gsl_matrix_alloc(n_coef, 2);
       	bspline_workspaces = new gsl_bspline_workspace*[n_to_force_match];
     	
       	int counter = 0; // this is a stand in for index_among_matched_interxns
		for (unsigned i = 0; i < n_defined; i++) {
			if (ispec_->defined_to_matched_intrxn_index_map[i] > 0) {
				interaction_column_indices = ispec_->interaction_column_indices[counter + 1] - ispec_->interaction_column_indices[counter];
				n_to_print_minus_bspline_k = interaction_column_indices - n_coef + 2;
				check_bspline_size(n_to_print_minus_bspline_k, (int)(n_coef));
				bspline_workspaces[counter] = gsl_bspline_alloc(n_coef, n_to_print_minus_bspline_k);
				gsl_bspline_knots_uniform(ispec_->lower_cutoffs[i], ispec_->upper_cutoffs[i], bspline_workspaces[counter]);
				counter++;
			}
		}
    }
}

BSplineAndDerivComputer::~BSplineAndDerivComputer() 
{
    if ((ispec_->class_type != kThreeBodyNonbonded) || (class_subtype == 1) || (class_subtype == 2) || (ispec_->class_subtype == 3) ) {
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
	size_t istart, iend;
    double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    first_nonzero_basis_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth());

    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;
    gsl_bspline_deriv_eval_nonzero(param_less_lower_cutoff + ispec_->lower_cutoffs[index_among_defined], (size_t)(1), bspline_matrices, &istart, &iend, bspline_workspaces[index_among_matched]);
    check_matching_istart(first_nonzero_basis_index, istart);
    
    for (unsigned i = 0; i < n_coef; i++) {
    	vals[i] = -gsl_matrix_get(bspline_matrices, i, 1);
    }
}

void BSplineAndDerivComputer::calculate_basis_fn_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals) 
{
    assert(vals.size() == n_coef);
    size_t istart, iend;
    double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    first_nonzero_basis_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth());
    
    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;
    gsl_bspline_eval_nonzero(param_less_lower_cutoff + ispec_->lower_cutoffs[index_among_defined], bspline_vectors, &istart, &iend, bspline_workspaces[index_among_matched]);
    check_matching_istart(first_nonzero_basis_index, istart);

    for (unsigned i = 0; i < n_coef; i++) {
        vals[i] = gsl_vector_get(bspline_vectors, i);
    }
}

double BSplineAndDerivComputer::evaluate_spline(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    size_t istart, iend;
    int ici_index = 0;
    double force = 0.0;
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
	double axis_val = check_against_cutoffs(axis, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    gsl_bspline_eval_nonzero(axis_val, bspline_vectors, &istart, &iend, bspline_workspaces[index_among_matched_interactions - 1]);
    if (index_among_matched_interactions > 0) {
		ici_index = ispec_->interaction_column_indices[index_among_matched_interactions - 1];
    }
    for (int tn = int(istart); tn <= int(iend); tn++) {
    	check_bspline_sizing(spline_coeffs.size(), first_nonzero_basis_index, index_among_matched_interactions, ici_index, tn, istart);
    	force += gsl_vector_get(bspline_vectors, tn - istart) * spline_coeffs[first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + tn];
    }
    return force;
}

double BSplineAndDerivComputer::evaluate_spline_deriv(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    double deriv = 0.0;
    size_t istart, iend;
    int ici_index = 0;
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
    double axis_val = check_against_cutoffs(axis, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
	gsl_bspline_deriv_eval_nonzero(axis_val, size_t(1), bspline_matrices, &istart, &iend, bspline_workspaces[index_among_matched_interactions - 1]);
    if (index_among_matched_interactions > 0) {
		ici_index = ispec_->interaction_column_indices[index_among_matched_interactions - 1];
    }
    for (int tn = int(istart); tn <= int(iend); tn++) {
    	check_bspline_sizing(spline_coeffs.size(), first_nonzero_basis_index, index_among_matched_interactions, ici_index, tn, istart);
        deriv += gsl_matrix_get(bspline_matrices, tn - istart, 1) * spline_coeffs[first_nonzero_basis_index + ici_index + tn];
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

DeltaSplineComputer::DeltaSplineComputer(InteractionClassSpec* ispec) : SplineComputer(ispec)
{
    n_coef = 1;
    n_to_force_match = ispec_->n_to_force_match;
    n_defined = ispec_->get_n_defined();
    binwidth = ispec_->get_fm_binwidth();
}

// Calculate the value of the a one-parameter linear spline; direction of the 
// corresponding forces is calculated in the function calling this one.

void DeltaSplineComputer::calculate_basis_fn_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals)
{
    assert(vals.size() == n_coef);
    first_nonzero_basis_index = 0;
    vals[0] = 1.0;
}

double DeltaSplineComputer::evaluate_spline(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
	double param_less_lower_cutoff = axis - ispec_->lower_cutoffs[index_among_defined];
	int basis_function_column_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth());
    double force = 0.0;
    if ( (unsigned)(first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + basis_function_column_index) <= spline_coeffs.size()) {
	    force = spline_coeffs[first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1] + basis_function_column_index];
    } else {
        fprintf(stderr, "Warning: attempting to read %d columns past interaction column %d in a matrix of %lu columns.", basis_function_column_index + 1, first_nonzero_basis_index + ispec_->interaction_column_indices[index_among_matched_interactions - 1], spline_coeffs.size());
    }
    return force;
}

void LinearSplineComputer::calculate_basis_fn_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals)
{
    assert(vals.size() == n_coef);
    double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    first_nonzero_basis_index = int(param_less_lower_cutoff / ispec_->get_fm_binwidth());
    vals[1] = fmod(param_less_lower_cutoff / ispec_->get_fm_binwidth(), 1.0);
    vals[0] = 1.0 - vals[1];
}

double LinearSplineComputer::evaluate_spline(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    double force = 0.0;
    int ici_index = 0;
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
	double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, axis);

	int basis_function_column_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth());
    double remainder_after_binning = fmod(param_less_lower_cutoff / ispec_->get_fm_binwidth(), 1.0);
    if (index_among_matched_interactions > 0) {
		ici_index = ispec_->interaction_column_indices[index_among_matched_interactions - 1];
    }
    
    if (unsigned(first_nonzero_basis_index + ici_index + basis_function_column_index + 1) < spline_coeffs.size()) {
        force = spline_coeffs[first_nonzero_basis_index + ici_index + basis_function_column_index] * (1.0 - remainder_after_binning) + spline_coeffs[first_nonzero_basis_index + ici_index + basis_function_column_index + 1] * remainder_after_binning;
    } else if (unsigned(first_nonzero_basis_index + ici_index + basis_function_column_index + 1) == spline_coeffs.size()) {
        force = spline_coeffs[first_nonzero_basis_index + ici_index + basis_function_column_index];
    } else {
        fprintf(stderr, "Warning: attempting to read %d columns past interaction column %d in a matrix of %lu columns, to be multiplied by a coefficient %g.", basis_function_column_index + 1, first_nonzero_basis_index + ici_index, spline_coeffs.size(), remainder_after_binning);
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

inline void check_bspline_size(const int control_points, const int order)
{
	if (control_points < order) {
		fprintf(stderr, "Cannot create a BSpline with fewer control points (%d) than the spline order (%d).\n", control_points, order);
		fprintf(stderr, "Either decrease the bin width or decrease the spline order for this interaction.\n");
		exit(EXIT_FAILURE);
	}
}

inline void check_matching_istart(const int first_nonzero_basis_index, const size_t istart)
{
	if( abs(first_nonzero_basis_index - (int)(istart)) > 1) {
    	fprintf(stderr, "Internal failure! first_nonnero_basis_index %d != gsl_bspline_eval_nonzero %d!\n", first_nonzero_basis_index, (int)(istart));
    }
}

inline double check_against_cutoffs(const double axis, const double lower_cutoff, const double upper_cutoff)
{
    if (axis < lower_cutoff) {
    	if (axis + VERYSMALL_F < lower_cutoff) fprintf(stderr, "Value to evaluate (%lf) is below this spline's lower cutoff (%lf)!\n\n", axis, lower_cutoff);
    	return lower_cutoff + VERYSMALL_F;
	} else if (axis > upper_cutoff) {
    	if (axis > upper_cutoff + VERYSMALL_F) {
    		fprintf(stderr, "Difference is %d x 10E-6\n", (int)( (axis - upper_cutoff) * 10E6));
    		fprintf(stderr, "Value to evaluate (%lf) is above this spline's upper cutoff (%lf)!\n\n", axis, upper_cutoff);
    		}
    	return upper_cutoff - VERYSMALL_F;
	}
	return axis;
}

inline void check_bspline_sizing(const size_t coeffs_size, const int first_nonzero_basis_index, const int index_among_matched_interactions, const int ici_index, const int tn, const size_t istart)
{
	if ((int)(coeffs_size) <= first_nonzero_basis_index + ici_index + tn) {
		fprintf(stderr, "Internal sizing issue encountered!\n");
		fprintf(stderr, "gsl_matrix_get(bspline_vectors, %d - %d)\t", tn, (int)(istart));
		fprintf(stderr, "spline_coeffs.size() = %u\n", (unsigned)(coeffs_size));
		fprintf(stderr, "index %d = ", first_nonzero_basis_index + ici_index + tn);
		fprintf(stderr, "fnzbi %d + ici[%d - 1] %d + tn %d\n", first_nonzero_basis_index, index_among_matched_interactions, ici_index, tn);
	}
}
