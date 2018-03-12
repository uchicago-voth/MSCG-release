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
inline double check_against_cutoffs(const double axis, const double lower_cutoff, const double upper_cutoff);
inline void check_bspline_sizing(const size_t coeffs_size, const int first_nonzero_basis_index, const int index_among_matched_interactions, const int ici_index, const int tn, const size_t istart);

// Helper functions for setting up periodic splines
inline void adjust_splines_for_periodicity(const InteractionClassType class_type, const int n_coef, const std::vector<unsigned> defined_to_periodic_intrxn_index_map, std::vector<unsigned> &interaction_column_indices);
inline void shift_remaining_indices(const int start, const int bspline_k, std::vector<unsigned> &interaction_column_indices, const int size);

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
	} else if (ispec->get_basis_type() == kPower) {
	  return new PowerComputer(ispec);
	} else if (ispec->get_basis_type() == kLJ) {
	  return new LJComputer(ispec);
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

SplineComputer::SplineComputer(InteractionClassSpec* ispec) 
{
	ispec_ = ispec;
	n_coef = ispec->get_bspline_k();
	n_to_force_match = ispec->n_to_force_match; 
	n_defined = ispec->get_n_defined();
	binwidth = ispec->get_fm_binwidth();
	int size = ispec->interaction_column_indices.size();
	interaction_column_indices_ = std::vector<unsigned>(size, 0);
	for (int i = 0; i < size; i++) interaction_column_indices_[i] = ispec->interaction_column_indices[i];
}

// This routine returns the input value minus the lower cutoff, rounded
// inside the allowed range adjusted by a small rounding factor.
double SplineComputer::get_param_less_lower_cutoff(const int index_among_defined, const double param_val) const 
{
  double cutoff_range = ispec_->upper_cutoffs[index_among_defined] - ispec_->lower_cutoffs[index_among_defined]; // - VERYSMALL;
    double param_less_lower_cutoff = param_val - ispec_->lower_cutoffs[index_among_defined];
    if (param_less_lower_cutoff < 0.0) {
    	if (param_less_lower_cutoff < -VERYSMALL_F - 0.5 * ispec_->get_fm_binwidth()) fprintf(stderr, "Value passed to spline computer (%lf) is below its lower cutoff (%lf)!\n\n", param_val, ispec_->lower_cutoffs[index_among_defined]);
    	param_less_lower_cutoff = 0.0;
    }
    else if (param_less_lower_cutoff > cutoff_range) {
    	if (param_less_lower_cutoff > cutoff_range + VERYSMALL_F + 0.5 * ispec_->get_fm_binwidth()) fprintf(stderr, "Value passed to spline computer %lf is above its cutoff %lf!\n\n", param_val, ispec_->upper_cutoffs[index_among_defined]);
	    param_less_lower_cutoff = cutoff_range;
    }
    return param_less_lower_cutoff;
}


BSplineComputer::BSplineComputer(InteractionClassSpec* ispec) : SplineComputer(ispec)
{
    int ici_index, n_to_print_minus_bspline_k;

    if (n_coef < 2) {
    	printf("Spline order for BSplineAndDeriv must be at least 2!\n");
		exit(EXIT_FAILURE);
    }

    printf("Allocating b-spline temporaries for %d interactions.\n", n_to_force_match);
    bspline_workspaces = new gsl_bspline_workspace*[n_to_force_match];
    bspline_vectors = gsl_vector_alloc(n_coef);
    adjust_splines_for_periodicity(ispec->class_type, n_coef, ispec->defined_to_periodic_intrxn_index_map, interaction_column_indices_);
	
    int counter = 0;
    for (unsigned i = 0; i < n_defined; i++) {
      if (ispec_->defined_to_matched_intrxn_index_map[i] > 0) {
	ici_index = interaction_column_indices_[counter + 1] - interaction_column_indices_[counter];
	n_to_print_minus_bspline_k = ici_index - n_coef + 2;
	check_bspline_size(n_to_print_minus_bspline_k, (int)(n_coef));
	bspline_workspaces[counter] = gsl_bspline_alloc(n_coef, n_to_print_minus_bspline_k);
	if(ispec_->get_char_id() == 'a'){
	  gsl_bspline_knots_uniform((ispec_->lower_cutoffs[i]*3.14/180.0 - VERYSMALL_F), (ispec_->upper_cutoffs[i]*3.14/180.0 + VERYSMALL_F), bspline_workspaces[counter]);
	}
	else {
	  gsl_bspline_knots_uniform(ispec_->lower_cutoffs[i] - VERYSMALL_F, ispec_->upper_cutoffs[i] + VERYSMALL_F, bspline_workspaces[counter]);
	}
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
    //double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    //first_nonzero_basis_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth());
    
    double axis_val = check_against_cutoffs(param_val, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    if(ispec_->get_char_id() == 'a') axis_val *= (3.14/180.0);

    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;
    gsl_bspline_eval_nonzero(axis_val, bspline_vectors, &istart, &iend, bspline_workspaces[index_among_matched]);
    first_nonzero_basis_index = istart;
    
    for (unsigned i = 0; i < n_coef; i++) {
        vals[i] = gsl_vector_get(bspline_vectors, i);
    }
}

double BSplineComputer::evaluate_spline(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    size_t istart, iend;
    int ici_value = 0;
    double force = 0.0;
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
    double axis_val = check_against_cutoffs(axis, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    if(ispec_->get_char_id() == 'a') axis_val *= (3.14/180.0);

    gsl_bspline_eval_nonzero(axis_val, bspline_vectors, &istart, &iend, bspline_workspaces[index_among_matched_interactions - 1]);
    if (index_among_matched_interactions > 0) {
		ici_value = interaction_column_indices_[index_among_matched_interactions - 1];
    }
    for (int tn = int(istart); tn <= int(iend); tn++) {
        check_bspline_sizing(spline_coeffs.size(), first_nonzero_basis_index, index_among_matched_interactions, ici_value, tn, istart);
        force += gsl_vector_get(bspline_vectors, tn - istart) * spline_coeffs[first_nonzero_basis_index + ici_value + tn];
    }
    return force;
}

BSplineAndDerivComputer::BSplineAndDerivComputer(InteractionClassSpec* ispec) : SplineComputer(ispec)
{
    int n_to_print_minus_bspline_k, ici_index;
	//copy additional information

    
    if (n_coef < 3) {
    	printf("Spline order for BSplineAndDeriv must be at least 3!\n");
		exit(EXIT_FAILURE);
    }
 
    printf("Allocating b-spline and derivative temporaries for %d interactions.\n", ispec_->get_n_defined());
    bspline_vectors = gsl_vector_alloc(n_coef);
    bspline_matrices = gsl_matrix_alloc(n_coef, 2);

    adjust_splines_for_periodicity(ispec->class_type, n_coef, ispec->defined_to_periodic_intrxn_index_map, interaction_column_indices_);
	
    if (n_to_force_match < 1) {
      bspline_workspaces = new gsl_bspline_workspace*[1];
    } else {
      bspline_workspaces = new gsl_bspline_workspace*[n_to_force_match];
    }
	
    int counter = 0; // this is a stand in for index_among_matched_interxns
    for (unsigned i = 0; i < n_defined; i++) {
      if (ispec_->defined_to_matched_intrxn_index_map[i] > 0) {
	ici_index = interaction_column_indices_[counter + 1] - interaction_column_indices_[counter];
	n_to_print_minus_bspline_k = ici_index - n_coef + 2;
	check_bspline_size(n_to_print_minus_bspline_k, (int)(n_coef));
	bspline_workspaces[counter] = gsl_bspline_alloc(n_coef, n_to_print_minus_bspline_k);
	if(ispec_->get_char_id() == 'a'){
	  gsl_bspline_knots_uniform((ispec_->lower_cutoffs[i] * 3.14/180.0 - VERYSMALL_F), (ispec_->upper_cutoffs[i] *3.14/180.0 + VERYSMALL_F), bspline_workspaces[counter]);
	}
	else {
	  gsl_bspline_knots_uniform(ispec_->lower_cutoffs[i] - VERYSMALL_F, ispec_->upper_cutoffs[i] + VERYSMALL_F, bspline_workspaces[counter]);
	}
	counter++;
      }
    }
}

BSplineAndDerivComputer::~BSplineAndDerivComputer() 
{
	if (ispec_->class_type != kThreeBodyNonbonded) {
		for (unsigned i = 0; i < n_to_force_match; i++)	gsl_bspline_free(bspline_workspaces[i]);
	} else {
		for (unsigned i = 0; i < n_defined; i++)  gsl_bspline_free(bspline_workspaces[i]);
	}
	
	delete [] bspline_workspaces;
	gsl_vector_free(bspline_vectors);
	gsl_matrix_free(bspline_matrices);
}

void BSplineAndDerivComputer::calculate_bspline_deriv_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals)
{
    assert(vals.size() == n_coef);
    size_t istart, iend;

    //double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    //first_nonzero_basis_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth() + 0.5);

    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;

    double axis_val = check_against_cutoffs(param_val, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    gsl_bspline_deriv_eval_nonzero(axis_val, (size_t)(1), bspline_matrices, &istart, &iend, bspline_workspaces[index_among_matched]);
    first_nonzero_basis_index = istart;
    
    for (unsigned i = 0; i < n_coef; i++) {
    	vals[i] = -gsl_matrix_get(bspline_matrices, i, 1);
    }
}

void BSplineAndDerivComputer::calculate_basis_fn_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals) 
{
    assert(vals.size() == n_coef);
    size_t istart, iend;
    //double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    //first_nonzero_basis_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth() + 0.5);
    
    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;
    double axis_val = check_against_cutoffs(param_val, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    if(ispec_->get_char_id() == 'a') axis_val *= (3.14/180.0);

    gsl_bspline_eval_nonzero(axis_val, bspline_vectors, &istart, &iend, bspline_workspaces[index_among_matched]);
    first_nonzero_basis_index = istart;
    
    for (unsigned i = 0; i < n_coef; i++) {
        vals[i] = gsl_vector_get(bspline_vectors, i);
    }
}

double BSplineAndDerivComputer::evaluate_spline(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    size_t istart, iend;
    int ici_value = 0;
    double force = 0.0;
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
    double axis_val = check_against_cutoffs(axis, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    if(ispec_->get_char_id() == 'a') axis_val *= (3.14/180.0);

    gsl_bspline_eval_nonzero(axis_val, bspline_vectors, &istart, &iend, bspline_workspaces[index_among_matched_interactions - 1]);
    if (index_among_matched_interactions > 0) {
		ici_value = interaction_column_indices_[index_among_matched_interactions - 1];
    }
    for (int tn = int(istart); tn <= int(iend); tn++) {
    	check_bspline_sizing(spline_coeffs.size(), first_nonzero_basis_index, index_among_matched_interactions, ici_value, tn, istart);
    	force += gsl_vector_get(bspline_vectors, tn - istart) * spline_coeffs[first_nonzero_basis_index + ici_value + tn];
    }
    return force;
}

double BSplineAndDerivComputer::evaluate_spline_deriv(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    double deriv = 0.0;
    size_t istart, iend;
    int ici_value = 0;
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
    double axis_val = check_against_cutoffs(axis, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    if(ispec_->get_char_id() == 'a') axis_val *= (3.14/180.0);

    gsl_bspline_deriv_eval_nonzero(axis_val, size_t(1), bspline_matrices, &istart, &iend, bspline_workspaces[index_among_matched_interactions - 1]);
    if (index_among_matched_interactions > 0) {
		ici_value = interaction_column_indices_[index_among_matched_interactions - 1];
    }
    for (int tn = int(istart); tn <= int(iend); tn++) {
    	check_bspline_sizing(spline_coeffs.size(), first_nonzero_basis_index, index_among_matched_interactions, ici_value, tn, istart);
        deriv += gsl_matrix_get(bspline_matrices, tn - istart, 1) * spline_coeffs[first_nonzero_basis_index + ici_value + tn];
    }
    return deriv;
}

LinearSplineComputer::LinearSplineComputer(InteractionClassSpec* ispec) : SplineComputer(ispec)
{
	// Override generic constructor settings
    n_coef = 2;
}

DeltaSplineComputer::DeltaSplineComputer(InteractionClassSpec* ispec) : SplineComputer(ispec)
{
	// Override generic constructor settings
    n_coef = 1;
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
    if ( (unsigned)(first_nonzero_basis_index + interaction_column_indices_[index_among_matched_interactions - 1] + basis_function_column_index) <= spline_coeffs.size()) {
	    force = spline_coeffs[first_nonzero_basis_index + interaction_column_indices_[index_among_matched_interactions - 1] + basis_function_column_index];
    } else {
        fprintf(stderr, "Warning: attempting to read %d columns past interaction column %d in a matrix of %lu columns.", basis_function_column_index + 1, first_nonzero_basis_index + interaction_column_indices_[index_among_matched_interactions - 1], spline_coeffs.size());
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
    int ici_value = 0;
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
	double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, axis);

	int basis_function_column_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth());
    double remainder_after_binning = fmod(param_less_lower_cutoff / ispec_->get_fm_binwidth(), 1.0);
    if (index_among_matched_interactions > 0) {
		ici_value = interaction_column_indices_[index_among_matched_interactions - 1];
    }
    
    if (unsigned(first_nonzero_basis_index + ici_value + basis_function_column_index + 1) < spline_coeffs.size()) {
        force = spline_coeffs[first_nonzero_basis_index + ici_value + basis_function_column_index] * (1.0 - remainder_after_binning) + spline_coeffs[first_nonzero_basis_index + ici_value + basis_function_column_index + 1] * remainder_after_binning;
    } else if (unsigned(first_nonzero_basis_index + ici_value + basis_function_column_index + 1) == spline_coeffs.size()) {
        force = spline_coeffs[first_nonzero_basis_index + ici_value + basis_function_column_index];
    } else {
        fprintf(stderr, "Warning: attempting to read %d columns past interaction column %d in a matrix of %lu columns, to be multiplied by a coefficient %g.", basis_function_column_index + 1, first_nonzero_basis_index + ici_value, spline_coeffs.size(), remainder_after_binning);
    }
    return force;
}

TableSplineComputer::TableSplineComputer(InteractionClassSpec* ispec) : SplineComputer(ispec) 
{
	// Override generic constructor settings
    n_coef = 2;
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

PowerComputer::PowerComputer(InteractionClassSpec* ispec) : SplineComputer(ispec)
{
    printf("Allocating b-spline temporaries for %d interactions.\n", n_to_force_match);
}


void PowerComputer::calculate_basis_fn_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals)
{
    assert(vals.size() == n_coef);
    //size_t istart, iend;
    //double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    //first_nonzero_basis_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth() + 0.5);
    
    //    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;
    //    double axis_val = check_against_cutoffs(param_val, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    double new_param_val;
    if (ispec_->get_char_id() == 'n'){
      new_param_val = param_val * 2/ispec_->upper_cutoffs[index_among_defined];
      inverse_power_eval(param_val, vals);
    } else if (ispec_->get_char_id() == 'a'){
      new_param_val = param_val*(3.14/180);      
      power_eval(new_param_val, vals);
    } else {
      power_eval(param_val, vals);
    }
    first_nonzero_basis_index = 0;    
}

void PowerComputer::calculate_bspline_deriv_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals)
{
    assert(vals.size() == n_coef);
    //size_t istart, iend;
    //double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    //first_nonzero_basis_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth() + 0.5);
    
    //    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;
    //    double axis_val = check_against_cutoffs(param_val, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    double new_param_val;
    if (ispec_->get_char_id() == 'n'){
      new_param_val = param_val * 2/ispec_->upper_cutoffs[index_among_defined];
      inverse_deriv_eval(param_val, vals);
    } else if (ispec_->get_char_id() == 'a'){
      new_param_val = param_val*(3.14/180);      
      deriv_eval(new_param_val, vals);
    } else if (ispec_->get_char_id() == 'd'){
      new_param_val = param_val*(3.14/180);
      fourier_eval(new_param_val, vals);
    } else {
      deriv_eval(param_val, vals);
    }
    first_nonzero_basis_index = 0;
}

double PowerComputer::evaluate_spline(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis)
{
    int ici_value = 0;
    double force = 0.0;
    //   double axis_val = check_against_cutoffs(axis, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
    if (index_among_matched_interactions > 0) {
		ici_value = interaction_column_indices_[index_among_matched_interactions - 1];
    }
    double new_axis;
    if (ispec_->get_char_id() == 'n'){
      new_axis = axis * 2/ispec_->upper_cutoffs[index_among_defined];
      force = inverse_power_axis(index_among_matched_interactions, spline_coeffs,axis,ici_value,first_nonzero_basis_index);
    } else if (ispec_->get_char_id() == 'a'){
      new_axis = axis*(3.14/180);      
      force = power_axis(index_among_matched_interactions, spline_coeffs,new_axis,ici_value,first_nonzero_basis_index);
    } else if (ispec_->get_char_id() == 'd'){
      new_axis = axis*(3.14/180);
      force = fourier_axis(index_among_matched_interactions, spline_coeffs,new_axis,ici_value,first_nonzero_basis_index);
    } else {
      force = power_axis(index_among_matched_interactions, spline_coeffs,axis,ici_value,first_nonzero_basis_index);
    }
    return force;
}

double PowerComputer::evaluate_spline_deriv(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
    int ici_value = 0;
    double deriv = 0.0;

    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
    if (index_among_matched_interactions > 0) {
		ici_value = interaction_column_indices_[index_among_matched_interactions - 1];
    }
    double new_axis;
    if (ispec_->get_char_id() == 'n'){
      new_axis = axis * 2/ispec_->upper_cutoffs[index_among_defined];
      deriv = inverse_deriv_axis(index_among_matched_interactions, spline_coeffs,axis,ici_value,first_nonzero_basis_index);
    } else if (ispec_->get_char_id() == 'a'){
      new_axis = axis*(3.14/180);      
      deriv = deriv_axis(index_among_matched_interactions, spline_coeffs,new_axis,ici_value,first_nonzero_basis_index);
    } else{
      deriv = deriv_axis(index_among_matched_interactions, spline_coeffs,axis,ici_value,first_nonzero_basis_index);
    }
    return deriv;
}


void PowerComputer::power_eval(const double param_val, std::vector<double> &vals)
{
  unsigned i;
  double functemp = 1;
  for(i=0;i<n_coef;i++)
    {
      vals[i] = functemp;
      functemp *= param_val;
    }
  
}

void PowerComputer::deriv_eval(const double param_val, std::vector<double> &vals)
{
  unsigned i;
  double functemp = 1;

  vals[0] = 0.0;
  for(i=1;i<n_coef;i++)
    {
      vals[i] = i * functemp;
      functemp *= param_val;
    }
  
}

double PowerComputer::power_axis(const int index_among_defined, const std::vector<double> &spline_coeffs, const double axis_val, int ici_value, int first_nonzero_basis_index)
{
  
  unsigned i;
  double functemp = 1;
  double forcetemp = 0.0;
  
  for(i=0;i<n_coef;i++)
    {
      forcetemp += spline_coeffs[i + ici_value + first_nonzero_basis_index] * functemp;
      functemp = functemp*axis_val;
    }
  return forcetemp;
}

double PowerComputer::deriv_axis(const int index_among_defined, const std::vector<double> &spline_coeffs, const double axis_val, int ici_value, int first_nonzero_basis_index)
{
  
  unsigned i;
  double functemp = 1;
  double derivtemp = 0.0;
  
  for(i=1;i<n_coef;i++)
    {
      derivtemp += i * spline_coeffs[i + ici_value + first_nonzero_basis_index] * functemp;
      functemp = functemp*axis_val;
    }
  return derivtemp;
}

void PowerComputer::inverse_power_eval(const double param_val, std::vector<double> &vals)
{
  unsigned i;
  double functemp = 1;
  double cutofftemp = 1;
  //  printf("%lf\n",param_val);
  for(i=0;i<n_coef;i++)
    {
      functemp /= param_val;
      cutofftemp /= 2.0;
      vals[i] = functemp - cutofftemp;
    }
  
}

void PowerComputer::inverse_deriv_eval(const double param_val, std::vector<double> &vals)
{
 unsigned i;
  double functemp = 1/param_val;
  for(i=0;i<n_coef;i++)
    {
      vals[i] = -(i+1) * functemp;
      functemp /= param_val;
    }
  
}

double PowerComputer::inverse_power_axis(const int index_among_defined, const std::vector<double> &spline_coeffs, const double axis_val, int ici_value, int first_nonzero_basis_index)
{
  unsigned i;
  double functemp = 1;
  double forcetemp = 0.0;
  double cutofftemp = 1;
  for(i=0;i<n_coef;i++)
    {
      functemp = functemp/axis_val;
      cutofftemp = cutofftemp/2.0;
      forcetemp += spline_coeffs[i+ici_value + first_nonzero_basis_index] * (functemp - cutofftemp);
    }
  return forcetemp;
}

double PowerComputer::inverse_deriv_axis(const int index_among_defined, const std::vector<double> &spline_coeffs, const double axis_val, int ici_value, int first_nonzero_basis_index)
{
  unsigned i;
  double functemp = 1/axis_val;
  double derivtemp = 0.0;
  
  for(i=0;i<n_coef;i++)
    {
      derivtemp += -(i+1) * spline_coeffs[i+ici_value + first_nonzero_basis_index] * functemp;
      functemp = functemp/axis_val;
    }
  return derivtemp;
}

void PowerComputer::fourier_eval(const double param_val, std::vector<double> &vals)
{
  unsigned i;
  for(i=0;i<(n_coef/2);i++)
    {
      vals[i] = sin((i+1)*param_val);
      vals[(n_coef/2) +i] = cos((i+1)*param_val);
    }
}

void PowerComputer::fourier_deriv_eval(const double param_val, std::vector<double> &vals)
{
  unsigned i;
  for(i=0;i<(n_coef/2);i++)
    {
      vals[i] = (i+1)*cos((i+1)*param_val);
      vals[(n_coef/2) +1] = -(i+1)*sin((i+1)*param_val);
    }
}

double PowerComputer::fourier_axis(const int index_among_defined, const std::vector<double> &spline_coeffs, const double axis_val, int ici_value, int first_nonzero_basis_index)
{
  unsigned i;
  double forcetemp = 0.0;
  for(i=0;i<(n_coef/2);i++)
    {
      forcetemp += spline_coeffs[i+ici_value + first_nonzero_basis_index] * sin((i+1)*axis_val);
      forcetemp += spline_coeffs[(n_coef/2) + i+ici_value + first_nonzero_basis_index] * cos((i+1)*axis_val);
    }
  return forcetemp;
}

double PowerComputer::fourier_deriv_axis(const int index_among_defined, const std::vector<double> &spline_coeffs, const double axis_val, int ici_value, int first_nonzero_basis_index)
{
  unsigned i;
  double derivtemp = 0.0;
  for(i=0;i<(n_coef/2);i++)
    {
      derivtemp += spline_coeffs[i+ici_value + first_nonzero_basis_index] * sin((i+1)*axis_val);
      derivtemp += spline_coeffs[(n_coef/2) + i+ici_value + first_nonzero_basis_index] * cos((i+1)*axis_val);
    }
  return derivtemp;
}

 

LJComputer::LJComputer(InteractionClassSpec* ispec) : SplineComputer(ispec)
{
    printf("Allocating b-spline temporaries for %d interactions.\n", n_to_force_match);
}


void LJComputer::calculate_basis_fn_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals)
{
    assert(vals.size() == n_coef);
    //size_t istart, iend;
    //double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    //first_nonzero_basis_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth() + 0.5);
    
    //    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;
    //    double axis_val = check_against_cutoffs(param_val, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    double new_param_val;
    if (ispec_->get_char_id() == 'n'){
      new_param_val = param_val * 2/ispec_->upper_cutoffs[index_among_defined];
      inverse_power_eval(param_val, vals);
    } else if (ispec_->get_char_id() == 'a'){
      new_param_val = param_val*(3.14/180);      
      power_eval(new_param_val, vals);
    } else {
      power_eval(param_val, vals);
    }
    first_nonzero_basis_index = 0;        
}

void LJComputer::calculate_bspline_deriv_vals(const int index_among_defined, const double param_val, int &first_nonzero_basis_index, std::vector<double> &vals)
{
    assert(vals.size() == n_coef);
    //size_t istart, iend;
    //double param_less_lower_cutoff = get_param_less_lower_cutoff(index_among_defined, param_val);
    //first_nonzero_basis_index = (int)(param_less_lower_cutoff / ispec_->get_fm_binwidth() + 0.5);
    
    //    int index_among_matched = ispec_->defined_to_matched_intrxn_index_map[index_among_defined] - 1;
    //double axis_val = check_against_cutoffs(param_val, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    double new_param_val;
    if (ispec_->get_char_id() == 'n'){
      new_param_val = param_val * 2/ispec_->upper_cutoffs[index_among_defined];
      inverse_deriv_eval(param_val, vals);
    } else if (ispec_->get_char_id() == 'a'){
      new_param_val = param_val*(3.14/180);      
      deriv_eval(new_param_val, vals);
    } else {
      deriv_eval(param_val, vals);
    }
}

double LJComputer::evaluate_spline(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis)
{
    int ici_value = 0;
    double force = 0.0;
    //   double axis_val = check_against_cutoffs(axis, ispec_->lower_cutoffs[index_among_defined], ispec_->upper_cutoffs[index_among_defined]);
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
    if (index_among_matched_interactions > 0) {
		ici_value = interaction_column_indices_[index_among_matched_interactions - 1];
    }
    double new_axis;
    if (ispec_->get_char_id() == 'n'){
      new_axis = axis * 2/ispec_->upper_cutoffs[index_among_defined];
      force = inverse_power_axis(index_among_matched_interactions, spline_coeffs,axis,ici_value,first_nonzero_basis_index);
    } else if (ispec_->get_char_id() == 'a'){
      new_axis = axis*(3.14/180);      
      force = power_axis(index_among_matched_interactions, spline_coeffs,new_axis,ici_value,first_nonzero_basis_index);
    } else {
      force = power_axis(index_among_matched_interactions, spline_coeffs,axis,ici_value,first_nonzero_basis_index);
    }

    return force;
}

double LJComputer::evaluate_spline_deriv(const int index_among_defined, const int first_nonzero_basis_index, const std::vector<double> &spline_coeffs, const double axis) 
{
  int ici_value = 0;
  double deriv = 0.0;
    int index_among_matched_interactions = ispec_->defined_to_matched_intrxn_index_map[index_among_defined];
    if (index_among_matched_interactions > 0) {
		ici_value = interaction_column_indices_[index_among_matched_interactions - 1];
    }
    double new_axis;
    if (ispec_->get_char_id() == 'n'){
      new_axis = axis * 2/ispec_->upper_cutoffs[index_among_defined];
      deriv = inverse_deriv_axis(index_among_matched_interactions, spline_coeffs,axis,ici_value,first_nonzero_basis_index);
    } else if (ispec_->get_char_id() == 'a'){
      new_axis = axis*(3.14/180);      
      deriv = deriv_axis(index_among_matched_interactions, spline_coeffs,new_axis,ici_value,first_nonzero_basis_index);
    } else{
      deriv = deriv_axis(index_among_matched_interactions, spline_coeffs,axis,ici_value,first_nonzero_basis_index);
    }
    return deriv;
}
  

void LJComputer::power_eval(const double param_val, std::vector<double> &vals)
{
  unsigned i;
  double functemp = 1;
  for(i=0;i<n_coef;i++)
    {
      vals[0] = functemp;
      functemp *= param_val;
    }
}

void LJComputer::inverse_power_eval(const double param_val, std::vector<double> &vals)
{
  int i;
  double functemp = 1;
  for(i=0;i<ispec_->LJ_term1;i++)
    {
      functemp /= param_val;
    }
  vals[0] = functemp;
  for(i=ispec_->LJ_term1;i<ispec_->LJ_term2;i++)
    {
      functemp /= param_val;
    }
  vals[1] = functemp;
}



void LJComputer::deriv_eval(const double param_val, std::vector<double> &vals)
{
  unsigned i;
  double functemp = 1.0;
  vals[0] = 0;
  for(i=1;i<n_coef;i++)
    {
      vals[0]= i * functemp;
      functemp *= param_val;
    }
}

void LJComputer::inverse_deriv_eval(const double param_val, std::vector<double> &vals)
{
  int i;
  double functemp = 1/param_val;
  for(i=0;i<ispec_->LJ_term1;i++)
    {
      functemp /= param_val;
    }
  vals[0] = -ispec_->LJ_term1*functemp;
  for(i=ispec_->LJ_term1;i<ispec_->LJ_term2;i++)
    {
      functemp /= param_val;
    }
  vals[1] = -ispec_->LJ_term2*functemp;
}


double LJComputer::power_axis(const int index_among_defined, const std::vector<double> &spline_coeffs, const double axis_val, int ici_value, int first_nonzero_basis_index)
{
  unsigned i;
  double functemp = 1;
  double forcetemp = 0.0;

  for(i=0;i<n_coef;i++)
    {
      forcetemp += functemp * spline_coeffs[i + ici_value + first_nonzero_basis_index];
      functemp *= axis_val;
    }
  return forcetemp;
}

double LJComputer::inverse_power_axis(const int index_among_defined, const std::vector<double> &spline_coeffs, const double axis_val, int ici_value, int first_nonzero_basis_index)
{
  int i;
  double functemp = 1;
  double forcetemp = 0.0;
  for(i=0;i<ispec_->LJ_term1;i++)
    {
      functemp /= axis_val;
    }
  forcetemp += functemp * spline_coeffs[ici_value + first_nonzero_basis_index];
  for(i=ispec_->LJ_term1;i<ispec_->LJ_term2;i++)
    {
      functemp /= axis_val;
    }
  forcetemp += functemp * spline_coeffs[1 + ici_value + first_nonzero_basis_index];
  return forcetemp;
}

double LJComputer::deriv_axis(const int index_among_defined, const std::vector<double> &spline_coeffs, const double axis_val, int ici_value, int first_nonzero_basis_index)
{
  unsigned i;
  double derivtemp = 0.0;
  double functemp = 1.0;

  for(i=1;i<n_coef;i++)
    {
      derivtemp += i * functemp * spline_coeffs[1 + ici_value + first_nonzero_basis_index];
      functemp *= axis_val;
    }
  return derivtemp;
}

double LJComputer::inverse_deriv_axis(const int index_among_defined, const std::vector<double> &spline_coeffs, const double axis_val, int ici_value, int first_nonzero_basis_index)
{
  int i;
  double functemp = 1/axis_val;
  double derivtemp = 0.0;
  for(i=0;i<ispec_->LJ_term1;i++)
    {
      functemp /= axis_val;
    }
  derivtemp += ispec_->LJ_term1 * functemp * spline_coeffs[ici_value + first_nonzero_basis_index];
  for(i=ispec_->LJ_term1;i<ispec_->LJ_term2;i++)
    {
      functemp /= axis_val;
    }
  derivtemp += ispec_->LJ_term2 * functemp * spline_coeffs[1 + ici_value + first_nonzero_basis_index];
  
  return derivtemp;
}

void LJComputer::fourier_eval(const double param_val, std::vector<double> &vals)
{
  unsigned i;
  for(i=0;i<(n_coef/2);i++)
    {
      vals[i] = sin((i+1)*param_val);
      vals[(n_coef/2) +i] = cos((i+1)*param_val);
    }
}

void LJComputer::fourier_deriv_eval(const double param_val, std::vector<double> &vals)
{
  unsigned i;
  for(i=0;i<(n_coef/2);i++)
    {
      vals[i] = (i+1)*cos((i+1)*param_val);
      vals[(n_coef/2) +1] = -(i+1)*sin((i+1)*param_val);
    }
}

double LJComputer::fourier_axis(const int index_among_defined, const std::vector<double> &spline_coeffs, const double axis_val, int ici_value, int first_nonzero_basis_index)
{
  unsigned i;
  double forcetemp = 0.0;
  for(i=0;i<(n_coef/2);i++)
    {
      forcetemp += spline_coeffs[i+ici_value + first_nonzero_basis_index] * sin((i+1)*axis_val);
      forcetemp += spline_coeffs[(n_coef/2) + i+ici_value + first_nonzero_basis_index] * cos((i+1)*axis_val);
    }
  return forcetemp;
}

double LJComputer::fourier_deriv_axis(const int index_among_defined, const std::vector<double> &spline_coeffs, const double axis_val, int ici_value, int first_nonzero_basis_index)
{
  unsigned i;
  double derivtemp = 0.0;
  for(i=0;i<(n_coef/2);i++)
    {
      derivtemp += spline_coeffs[i+ici_value + first_nonzero_basis_index] * sin((i+1)*axis_val);
      derivtemp += spline_coeffs[(n_coef/2) + i+ici_value + first_nonzero_basis_index] * cos((i+1)*axis_val);
    }
  return derivtemp;
}

inline void check_bspline_size(const int control_points, const int order)
{
	if (control_points < order) {
		fprintf(stderr, "Cannot create a BSpline with fewer control points (%d) than the spline order (%d).\n", control_points, order);
		fprintf(stderr, "Either decrease the bin width or decrease the spline order for this interaction.\n");
		exit(EXIT_FAILURE);
	}
}

inline double check_against_cutoffs(const double axis, const double lower_cutoff, const double upper_cutoff)
{
  if (axis < lower_cutoff) {
    if (axis < lower_cutoff - VERYSMALL_F) fprintf(stderr, "Value to evaluate (%lf) is below this spline's lower cutoff (%lf)!\n\n", axis, lower_cutoff);
    return lower_cutoff - VERYSMALL_F;
  } else if (axis > upper_cutoff) {
    if (axis > upper_cutoff + VERYSMALL_F) fprintf(stderr, "Value to evaluate (%lf) is above this spline's upper cutoff (%lf)!\n\n", axis, upper_cutoff);
    return upper_cutoff + VERYSMALL_F;
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

inline void adjust_splines_for_periodicity(const InteractionClassType class_type, const int n_coef, const std::vector<unsigned> defined_to_periodic_intrxn_index_map, std::vector<unsigned> &interaction_column_indices_)
{
	if (class_type != kDihedralBonded) return;
	int n_defined = defined_to_periodic_intrxn_index_map.size();
	// Look for periodic interactions
	for (int i = 0; i < n_defined; i++) {
		// Adjust the internal copy if it is periodic (make it bigger by n_coef.
		if (defined_to_periodic_intrxn_index_map[i] == 1) {
			shift_remaining_indices(i, n_coef, interaction_column_indices_, interaction_column_indices_.size());
		}
	}
}

inline void shift_remaining_indices(const int start, const int bspline_k, std::vector<unsigned> &interaction_column_indices, const int size)
{
	for(int i = start; i < size; i++) interaction_column_indices[i] += bspline_k;
}
