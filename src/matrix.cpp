//
//  matrix.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <array>

#include "control_input.h"
#include "interaction_model.h"
#include "external_matrix_routines.h"
#include "misc.h"
#include "matrix.h"

// Matrix implementation-specific routines that are properly
// abstracted into the matrix data struct.

// Matrix initialization routines

void initialize_dense_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg);
void initialize_accumulation_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg);
void initialize_sparse_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg);
void initialize_sparse_dense_normal_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg);
void initialize_sparse_sparse_normal_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg);
void initialize_dummy_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg);
void initialize_rem_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg);

// Helper matrix initialization routines

void determine_matrix_columns_and_rows(MATRIX_DATA* const mat, CG_MODEL_DATA* const cg, int const frames_per_traj_block, int const pressure_constraint_flag);
void estimate_number_of_sparse_elements(MATRIX_DATA* const mat, CG_MODEL_DATA* const cg);
void log_n_basis_functions(InteractionClassSpec &ispec);
void determine_BI_interaction_rows_and_cols(MATRIX_DATA* mat, InteractionClassComputer* const icomp);

// Matrix reset routines

void set_dense_matrix_to_zero(MATRIX_DATA* const mat);
void set_sparse_matrix_to_zero(MATRIX_DATA* const mat);
void set_sparse_accumulation_matrix_to_zero(MATRIX_DATA* const mat);
void set_accumulation_matrix_to_zero(MATRIX_DATA* const mat);
void set_accumulation_matrix_to_zero(MATRIX_DATA* const mat, dense_matrix* const dense_fm_matrix);
void set_dummy_matrix_to_zero(MATRIX_DATA* const mat);
void set_rem_matrix_to_zero(MATRIX_DATA* const mat);

// Interface-level functions that convert force magnitude and derivatives to matrix elements.

void accumulate_vector_one_body_tabulated_force(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int *particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void accumulate_vector_tabulated_forces(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void accumulate_vector_symmetric_tabulated_forces(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void accumulate_vector_one_body_force(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void accumulate_vector_matching_forces(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void accumulate_vector_symmetric_matching_forces(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void accumulate_scalar_one_body_tabulated_force(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void accumulate_scalar_tabulated_forces(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void accumulate_scalar_one_body_force(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void accumulate_scalar_matching_forces(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void accumulate_entropy_elements(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void accumulate_tabulated_entropy(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void accumulate_BI_elements(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);

// Matrix insertion routines

void insert_sparse_matrix_element(const int i, const int j, double* const x, MATRIX_DATA* const mat);
void insert_dense_matrix_element(const int i, const int j, double* const x, MATRIX_DATA* const mat);
void insert_accumulation_matrix_element(const int i, const int j, double* const x, MATRIX_DATA* const mat);
void insert_rem_matrix_element(const int i, double const x, MATRIX_DATA* const mat);
inline void insert_BI_matrix_element(const int i, const int j, double const x, MATRIX_DATA* const mat);

void insert_scalar_sparse_matrix_element(const int i, const int j, double* const x, MATRIX_DATA* const mat);
void insert_scalar_dense_matrix_element(const int i, const int j, double* const x, MATRIX_DATA* const mat);
void insert_scalar_accumulation_matrix_element(const int i, const int j, double* const x, MATRIX_DATA* const mat);

void insert_dense_matrix_virial_element(const int m, const int n, const double x, MATRIX_DATA* const mat);
void insert_sparse_matrix_virial_element(const int m, const int n, const double x, MATRIX_DATA* const mat);
void insert_accumulation_matrix_virial_element(const int m, const int n, const double x, MATRIX_DATA* const mat);

// Vector modification routines

void accumulate_force_into_dense_target_vector(MATRIX_DATA* mat, int particle_index, double* force_element);
void accumulate_force_into_accumulation_target_vector(MATRIX_DATA* mat, int particle_index, double* force_element);
void accumulate_scalar_into_dense_target_vector(MATRIX_DATA* mat, int particle_index, double* force_element);

void accumulate_constraint_into_dense_target_vector(MATRIX_DATA *mat, int frame_index, double constraint_element);
void accumulate_constraint_into_accumulation_target_vector(MATRIX_DATA *mat, int frame_index, double constraint_element);

void calculate_target_virial_in_dense_vector(MATRIX_DATA* const mat, double *pressure_constraint_rhs_vector);
void calculate_target_virial_in_accumulation_vector(MATRIX_DATA* const mat, double *pressure_constraint_rhs_vector);

void calculate_target_force_dense_vector(int shift_i, int site_i, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &f);
void calculate_target_force_accumulation_vector(int shift_i, int site_i, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &f);

// Post-frame-block routines

void convert_dense_fm_equation_to_normal_form_and_accumulate(MATRIX_DATA* const mat);
void convert_dense_target_force_vector_to_normal_form_and_accumulate(MATRIX_DATA* const mat);
void accumulate_accumulation_matrices(MATRIX_DATA* const mat);
void solve_sparse_matrix(MATRIX_DATA* const mat);
void convert_sparse_fm_equation_to_sparse_normal_form_and_accumulate(MATRIX_DATA* const mat);
void convert_sparse_fm_equation_to_dense_normal_form_and_accumulate(MATRIX_DATA* const mat);
void calculate_frame_average_and_add_to_normal_matrix(MATRIX_DATA* const mat);
void do_nothing_to_fm_matrix(MATRIX_DATA* const mat);

// Helper solver routines

int  get_n_nonzero_matrix_elements(MATRIX_DATA* const mat);
void convert_linked_list_to_csr_matrix(MATRIX_DATA* const mat, csr_matrix& csr_fm_matrix);
void precondition_sparse_matrix(int const fm_matrix_columns, double* h, csr_matrix* csr_normal_matrix);
void sparse_matrix_addition(MATRIX_DATA* const mat, double frame_weight, int nnzmax, csr_matrix& csr_normal_matrix, csr_matrix* main_normal_matrix);
void regularize_sparse_matrix(MATRIX_DATA* const mat);
void regularize_vector_sparse_matrix(MATRIX_DATA* const mat, double* regularization_vector);
void regularize_sparse_matrix(MATRIX_DATA* const mat, csr_matrix* csr_matrix);
void regularize_vector_sparse_matrix(MATRIX_DATA* const mat, csr_matrix* csr_normal_matrix, double* regularization_vector);
void pardiso_solve(MATRIX_DATA* const mat, csr_matrix* const sparse_matrix, double* const dense_fm_normal_rhs_vector);
void solve_this_sparse_matrix(MATRIX_DATA* const mat);
inline void create_sparse_normal_form_matrix(MATRIX_DATA* const mat, const int nnzmax, csr_matrix& csr_fm_matrix, csr_matrix& csr_normal_matrix, double* const dense_fm_rhs_vector, double* const dense_rhs_normal_vector);
inline void create_dense_normal_form(MATRIX_DATA* const mat, const double frame_weight, dense_matrix* const dense_fm_matrix, dense_matrix* normal_matrix, double* const dense_fm_rhs_vector, double* dense_fm_normal_rhs_vector);
inline double calculate_dense_residual(MATRIX_DATA* const mat, dense_matrix* const dense_fm_normal_matrix, double* const dense_fm_rhs_vector, std::vector<double> &fm_solution, double normalziation);
inline double calculate_sparse_residual(MATRIX_DATA* const mat, csr_matrix* sparse_fm_normal_matrix, double* const dense_fm_rhs_vector, std::vector<double> &fm_solution, double normalization);
inline void calculate_and_apply_dense_preconditioning(MATRIX_DATA* mat, dense_matrix* dense_fm_normal_matrix, double* h);
inline void calculate_dense_svd(MATRIX_DATA* mat, int fm_matrix_columns, dense_matrix* dense_fm_normal_matrix, double* dense_fm_normal_rhs_vector, double* singular_values);
inline void calculate_dense_svd(MATRIX_DATA* mat, int fm_matrix_columns, int fm_matrix_rows, dense_matrix* dense_fm_normal_matrix, double* dense_fm_normal_rhs_vector, double* singular_values);

// After-full-trajectory routines

void average_sparse_block_fm_solutions(MATRIX_DATA* const mat);
void solve_sparse_fm_normal_equations(MATRIX_DATA* const mat);
void solve_dense_fm_normal_equations(MATRIX_DATA* const mat);
void solve_accumulation_form_fm_equations(MATRIX_DATA* const mat);

// Bootstrapping routines

void convert_dense_fm_equation_to_normal_form_and_bootstrap(MATRIX_DATA* const mat);
void solve_sparse_matrix_for_bootstrap(MATRIX_DATA* const mat);
void convert_sparse_fm_equation_to_sparse_normal_form_and_bootstrap(MATRIX_DATA* const mat);
void accumulate_accumulation_matrices_for_bootstrap(MATRIX_DATA* const mat);
void convert_sparse_fm_equation_to_dense_normal_form_and_bootstrap(MATRIX_DATA* const mat);
void average_sparse_bootstrapping_solutions(MATRIX_DATA* const mat);
void solve_sparse_fm_bootstrapping_equations(MATRIX_DATA* const mat);
void solve_dense_fm_normal_bootstrapping_equations(MATRIX_DATA* const mat);
void solve_accumulation_form_bootstrapping_equations(MATRIX_DATA* const mat);

// Matrix-implementation-dependent functions for reading 
// batches of FM matrices.

void read_binary_dense_fm_matrix(MATRIX_DATA* const mat);
void read_binary_accumulation_fm_matrix(MATRIX_DATA* const mat);
void read_binary_sparse_fm_matrix(MATRIX_DATA* const mat);
void read_regularization_vector(MATRIX_DATA* const mat);
void write_iteration(const double* alpha_vec, const double beta, std::vector<double> fm_solution, const double residual, const int iteration, FILE* alpha_fp, FILE* beta_fp, FILE* sol_fp, FILE* res_fp);

//--------------------------------------------------------------------
// Matrix initialization routines
//--------------------------------------------------------------------

// Make a matrix (constructor).

MATRIX_DATA::MATRIX_DATA(ControlInputs* const control_input, CG_MODEL_DATA *const cg)
{
    // Perform sanity checks on the new parameters 
    #if _mkl_flag == 1
	mkl_set_num_threads(control_input->num_sparse_threads);
	#else 
	if (MatrixType(control_input->matrix_type) == kSparse || MatrixType(control_input->matrix_type) == kSparseSparse) {
        printf("Cannot use sparse solving (matrix_type 1 or 4) unless compiling with MKL (use newfm_mkl.x in makefile).\n");
		exit(EXIT_FAILURE);
    }
	if (MatrixType(control_input->matrix_type) == kSparseNormal) {
        printf("Cannot use sparse accumulation (matrix_type 3) unless compiling with MKL for the time being (use newfm_mkl.x in makefile).\n");
		exit(EXIT_FAILURE);
    }
	#endif
    
    // Ignore a user's choice to output certain quantities if they will not be calculated.
    if ( (MatrixType(control_input->matrix_type) != kDense) && (MatrixType(control_input->matrix_type) != kSparseNormal) && (control_input->output_normal_equations_rhs_flag != 0) ) {
        printf("Cannot output normal equations if normal equations are not being calculated.\n");
        printf("Use a different FM matrix format.\n");
        exit(EXIT_FAILURE);
    }
    
    // Override a user's choice of block_size if it conflicts with use_statistical_reweighting flag
    if ( (control_input->use_statistical_reweighting == 1) && (control_input->frames_per_traj_block != 1) ) {
    	printf("Cannot use statistical reweighting with %d frames per trajectory block.\n", control_input->frames_per_traj_block);
    	printf("Setting block_size to 1.\n");
    	control_input->frames_per_traj_block = 1;
    }
	
	if ( (control_input->bootstrapping_flag == 1) && (control_input->frames_per_traj_block != 1) ) {
		printf("Cannot use bootstrapping with %d frames per trajectory block.\n", control_input->bootstrapping_flag);
		printf("Please change the block size to 1 and recheck your inputs before rerunning.\n");
		exit(EXIT_FAILURE);
	}
	
	if ( (control_input->volume_weighting_flag == 1) && (control_input->frames_per_traj_block != 1) ) {
		printf("Cannot use volume weighting with %d frames per trajectory block.\n", control_input->bootstrapping_flag);
		printf("Please change the block size to 1 and recheck your inputs before rerunning.\n");
		exit(EXIT_FAILURE);
	}

	if ( (MatrixType(control_input->matrix_type) == kDense) && (control_input->frames_per_traj_block != 1) ) {
		printf("Cannot use dense matrix_type (0) with %d frames per trajectory block\n", control_input->frames_per_traj_block);
		printf("Setting block_size to 1.\n");
		control_input->frames_per_traj_block = 1;
	}
	
	if( (MatrixType(control_input->matrix_type) == kREM) == 1) {
	    printf("frame block size must be 1 when doing relative entropy minimization\n");
	    printf("Setting block_size to 1.\n");
	    control_input->frames_per_traj_block = 1;
	}
	
	if (control_input->position_dimension <= 0) {
		printf("Position dimension must be a positive integer\n");
		exit(EXIT_FAILURE);
    }
    
    if (control_input->position_dimension != DIMENSION) {
    	printf("The value of position_dimension(%d) in control_input does not match the compiled dimension(%d)!\n", control_input->position_dimension, DIMENSION);
    	exit(EXIT_FAILURE);
    }
    
    // Copy over basic data members.
    output_style 					= control_input->output_style;
    output_normal_equations_rhs_flag= control_input->output_normal_equations_rhs_flag;
    output_solution_flag 			= control_input->output_solution_flag;
    rcond							= control_input->rcond;
    itnlim 							= control_input->itnlim;
    iterative_calculation_flag 		= control_input->iterative_calculation_flag;
	num_sparse_threads 				= control_input->num_sparse_threads;
	position_dimension 				= control_input->position_dimension;
	volume_weighting_flag 			= control_input->volume_weighting_flag;
	
	// Copy bootstrapping information.
	bootstrapping_flag 				= control_input->bootstrapping_flag;
	bootstrapping_full_output_flag 	= control_input->bootstrapping_full_output_flag;
	bootstrapping_num_estimates 	= control_input->bootstrapping_num_estimates;
	
	// Copy residual, regularization, and bayesian options.
	regularization_style 			= control_input->regularization_style;
    tikhonov_regularization_param 	= control_input->tikhonov_regularization_param;
	bayesian_flag					= control_input->bayesian_flag;
	bayesian_max_iter				= control_input->bayesian_max_iter;
    output_residual                 = control_input->output_residual;
    force_sq_total					= 0.0;
 
  	// Copy and  initialize output matrix options.
 	output_raw_frame_blocks			= control_input->output_raw_frame_blocks;
 	output_raw_splines				= control_input->output_raw_splines;
 	if  (output_raw_frame_blocks == 1) {
 		frame_block_fh = fopen("matrix.csr", "w");
 	}

    // Set blockwise composition weighting factors
    frames_per_traj_block 			= control_input->frames_per_traj_block;
	current_frame_weight 			= 1.0;
	use_statistical_reweighting 	= control_input->use_statistical_reweighting;
	dynamic_state_samples_per_frame = 1;
	if(control_input->dynamic_state_sampling == 1) dynamic_state_samples_per_frame = control_input->dynamic_state_samples_per_frame;
    
    // Set normalization based on the default frame weight of 1.0 now, but overwrite later if needed.
    // Will be changed in newfm.cpp if there is another value from read_frame_weights
    normalization = 1.0 /  (double) control_input->n_frames;
    
    // Set size_per_vector based on scalar_matching_flag
    size_per_vector = DIMENSION;
    if (control_input->scalar_matching_flag == 1) {
    	size_per_vector = 1;
    }
    
    // Set accumulate_*_forces function pointers for scalar versus vector matching
    if (control_input->scalar_matching_flag == 1) {
    	accumulate_matching_forces 				= accumulate_scalar_matching_forces;
		accumulate_symmetric_matching_forces 	= accumulate_scalar_matching_forces;
		accumulate_one_body_force 				= accumulate_scalar_one_body_force;
		accumulate_tabulated_forces 			= accumulate_scalar_tabulated_forces;
		accumulate_symmetric_tabulated_forces 	= accumulate_scalar_tabulated_forces;
		accumulate_one_body_tabulated_force 	= accumulate_scalar_one_body_tabulated_force;
    } else {
    	accumulate_matching_forces 				= accumulate_vector_matching_forces;
		accumulate_symmetric_matching_forces	= accumulate_vector_symmetric_matching_forces;
		accumulate_one_body_force 			 	= accumulate_vector_one_body_force;
		accumulate_tabulated_forces 			= accumulate_vector_tabulated_forces;
		accumulate_symmetric_tabulated_forces 	= accumulate_vector_symmetric_tabulated_forces;
		accumulate_one_body_tabulated_force 	= accumulate_vector_one_body_tabulated_force;
    }
   
    // Perform matrix-type-specific initializations.
    switch (control_input->matrix_type) {
    case kDense:
        matrix_type = kDense;
        initialize_dense_matrix(this, control_input, cg);
        break;
    case kSparse:
        matrix_type = kSparse;
        initialize_sparse_matrix(this, control_input, cg);
        break;
    case kAccumulation:
        matrix_type = kAccumulation;
        initialize_accumulation_matrix(this, control_input, cg);
        break;
    case kSparseNormal:
    	matrix_type = kSparseNormal;
        initialize_sparse_dense_normal_matrix(this, control_input, cg);
        break;
    case kSparseSparse:
    	matrix_type = kSparseSparse;
        initialize_sparse_sparse_normal_matrix(this, control_input, cg);
        break;
	case kDummy:
        matrix_type = kDummy;
        initialize_dummy_matrix(this, control_input, cg);
        break;
    case kREM:
        matrix_type = kREM;
        initialize_rem_matrix(this, control_input, cg);
        break;
    default:
        printf("Invalid matrix type.");
        exit(EXIT_FAILURE);
    }
    
   	if (regularization_style == 2) {
		printf("read regularization vector\n");
		read_regularization_vector(this);
	}

    printf("Finished initializing FM matrix.\n");
}

// Initialize a dense matrix.

void initialize_dense_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg)
{
    // Set the struct's pseudopolymorphic methods
    mat->set_fm_matrix_to_zero = set_dense_matrix_to_zero;
    mat->accumulate_fm_matrix_element = insert_dense_matrix_element;
    mat->accumulate_target_force_element = accumulate_force_into_dense_target_vector;
    mat->accumulate_target_constraint_element = accumulate_constraint_into_dense_target_vector;
    
    if (control_input->scalar_matching_flag == 1) {
    	mat->accumulate_fm_matrix_element = insert_scalar_dense_matrix_element;
    }
    
    if (control_input->bootstrapping_flag == 1) {
    	mat->do_end_of_frameblock_matrix_manipulations = convert_dense_fm_equation_to_normal_form_and_bootstrap;
    } else { 
	    if (control_input->iterative_calculation_flag == 0) mat->do_end_of_frameblock_matrix_manipulations = convert_dense_fm_equation_to_normal_form_and_accumulate;
	    else if (control_input->iterative_calculation_flag == 1) mat->do_end_of_frameblock_matrix_manipulations = convert_dense_target_force_vector_to_normal_form_and_accumulate;
	}
    
    if (control_input->pressure_constraint_flag == 1) mat->accumulate_virial_constraint_matrix_element = insert_dense_matrix_virial_element;

	if (control_input->bootstrapping_flag == 1) {
		mat->finish_fm = solve_dense_fm_normal_bootstrapping_equations;
	} else {
	    mat->finish_fm = solve_dense_fm_normal_equations;
	}
	
    // Determine the size of the matrix from model specifications
	determine_matrix_columns_and_rows(mat, cg, control_input->frames_per_traj_block, control_input->pressure_constraint_flag);

    // Check that the matrix dimensions are reasonable and print diagnostics.    
    if ( (unsigned)(mat->fm_matrix_rows / mat->frames_per_traj_block) * (unsigned)(control_input->n_frames) < (unsigned)(mat->fm_matrix_columns) ) {
        printf("Current number of frames in this trajectory is too low to provide a fully-determined set of FM equations. Provide more frames in the input trajectory.\n");
        exit(EXIT_FAILURE);
    }
        
    printf("Number of rows for dense matrix algorithm: %d \n", mat->fm_matrix_rows);
    printf("Number of columns for dense matrix algorithm: %d \n", mat->fm_matrix_columns);
    
    // Check that the memory usage is reasonable and print 
    // memory diagnostics if so. These are checks for integer 
    // overflow when calculating the size of the matrices.

    if ( ( (int)(INT_MAX) / mat->fm_matrix_columns) <
        (mat->fm_matrix_columns * (int)(sizeof(double)))) {
        printf("Using this number of columns will lead to integer overflow in memory allocation for the normal matrix equations. Decrease the number of basis functions.\n");
        exit(EXIT_FAILURE);
    }
    
    if ( ( (int)(INT_MAX) / mat->fm_matrix_columns) <
        (mat->fm_matrix_rows * (int)(sizeof(double)))) {
        printf("Using this number of rows and columns will lead to integer overflow in memory allocation for the framewise matrix computation. Decrease number of particles or number of basis functions.\n");
        exit(EXIT_FAILURE);
    }
    
    printf("Size of per-frame matrix: %lu bytes \n", mat->fm_matrix_rows * mat->fm_matrix_columns * sizeof(double));
    printf("Size of normal matrix: %lu bytes \n", mat->fm_matrix_columns * mat->fm_matrix_columns * sizeof(double));

    // Allocate memory for the FM matrix and target vector as well as their normal form
    mat->accumulation_matrix_columns = mat->fm_matrix_columns;
    mat->accumulation_matrix_rows = mat->fm_matrix_rows;
    mat->dense_fm_matrix = new dense_matrix(mat->fm_matrix_rows, mat->fm_matrix_columns);
    mat->dense_fm_rhs_vector = new double[mat->fm_matrix_rows]();
    
    if (control_input->bootstrapping_flag == 1) {
    	mat->bootstrapping_dense_fm_normal_rhs_vectors = new double*[control_input->bootstrapping_num_estimates];
    	mat->bootstrapping_dense_fm_normal_matrices = new dense_matrix*[control_input->bootstrapping_num_estimates];
    	for (int i = 0; i < control_input->bootstrapping_num_estimates; i++) {
    		mat->bootstrapping_dense_fm_normal_rhs_vectors[i] = new double[mat->fm_matrix_columns]();
			mat->bootstrapping_dense_fm_normal_matrices[i] = new dense_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns);
		}
    }
	mat->dense_fm_normal_matrix = new dense_matrix(mat->fm_matrix_columns , mat->fm_matrix_columns);
    mat->dense_fm_normal_rhs_vector = new double[mat->fm_matrix_columns]();
    // Initialized the matrix and vector to zero.
    printf("Initialized a dense FM matrix.\n");
}

// Initialize an accumulation matrix.

void initialize_accumulation_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg)
{
    // Set pseudopolymorphic methods
    mat->set_fm_matrix_to_zero = set_accumulation_matrix_to_zero;
    mat->accumulate_fm_matrix_element = insert_accumulation_matrix_element;
    mat->accumulate_target_force_element = accumulate_force_into_accumulation_target_vector;
    mat->accumulate_target_constraint_element = accumulate_constraint_into_accumulation_target_vector;
    
    if (control_input->bootstrapping_flag == 1) {
    	mat->do_end_of_frameblock_matrix_manipulations = accumulate_accumulation_matrices_for_bootstrap;
    } else {
		mat->do_end_of_frameblock_matrix_manipulations = accumulate_accumulation_matrices;
    }
    
    if (control_input->pressure_constraint_flag == 1) mat->accumulate_virial_constraint_matrix_element = insert_accumulation_matrix_virial_element;
    
    if (control_input->bootstrapping_flag == 1) {
		mat->finish_fm = solve_accumulation_form_bootstrapping_equations;
	} else {
	    mat->finish_fm = solve_accumulation_form_fm_equations;
	}
	
    // Determine the size of the matrix from model specifications
	determine_matrix_columns_and_rows(mat, cg, control_input->frames_per_traj_block, control_input->pressure_constraint_flag);
        
    // Check that the matrix dimensions are reasonable.
    if ( (unsigned)(mat->fm_matrix_rows / mat->frames_per_traj_block) * (unsigned)(control_input->n_frames) < (unsigned)(mat->fm_matrix_columns) ) {
        printf("Current number of frames per frame block is too low to allow use of accumulation matrices. Increase block size.\n");
        exit(EXIT_FAILURE);
    }

    mat->accumulation_matrix_columns = mat->fm_matrix_columns + 1;
    mat->accumulation_matrix_rows = mat->fm_matrix_rows + mat->accumulation_matrix_columns;
    mat->accumulation_target_forces_location = mat->fm_matrix_columns * mat->accumulation_matrix_rows;

    printf("Number of rows for accumulation matrix algorithm: %d \n", mat->accumulation_matrix_columns);
    printf("Number of columns for accumulation matrix algorithm: %d \n", mat->accumulation_matrix_rows);
    
    // Allocate memory for the FM matrix and target vector as well as temp space for the accumulation operation.
    if (control_input->bootstrapping_flag == 1) {
		printf("Bootstrapping is not currently supported for accumulation matrix type.\n");
		exit(EXIT_FAILURE);
    	mat->bootstrapping_dense_fm_normal_rhs_vectors = new double*[control_input->bootstrapping_num_estimates];
    	mat->bootstrapping_dense_fm_normal_matrices = new dense_matrix*[control_input->bootstrapping_num_estimates];
    	for (int i = 0; i < control_input->bootstrapping_num_estimates; i++) {
    		mat->bootstrapping_dense_fm_normal_rhs_vectors[i] = new double[mat->fm_matrix_columns]();
			mat->bootstrapping_dense_fm_normal_matrices[i] = new dense_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns);
		}
    }
	mat->dense_fm_matrix = new dense_matrix(mat->accumulation_matrix_rows, mat->accumulation_matrix_columns);
	mat->dense_fm_normal_rhs_vector = new double[mat->accumulation_matrix_columns]();

    mat->lapack_tau = new double[mat->accumulation_matrix_columns]();

    // Initialized the matrix to zero.
    printf("Size of per-frame matrix: %lu bytes \n", mat->accumulation_matrix_columns * mat->accumulation_matrix_rows * sizeof(double));
    printf("Initialized an accumulation algorithm FM matrix.\n");
}

// Initialize a sparse-matrix-based computation.

void initialize_sparse_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg)
{
    // Set pseudopolymorphic methods
    mat->set_fm_matrix_to_zero = set_sparse_matrix_to_zero;
    mat->accumulate_fm_matrix_element = insert_sparse_matrix_element;
    mat->accumulate_target_force_element = accumulate_force_into_dense_target_vector;
    mat->accumulate_target_constraint_element = accumulate_constraint_into_dense_target_vector;
   	mat->sparse_matrix = NULL;

    if (control_input->scalar_matching_flag == 1) {
    	mat->accumulate_fm_matrix_element = insert_scalar_sparse_matrix_element;
    }

   	if (control_input->bootstrapping_flag == 1) {
    	mat->do_end_of_frameblock_matrix_manipulations = solve_sparse_matrix_for_bootstrap;
    } else {
	    mat->do_end_of_frameblock_matrix_manipulations = solve_sparse_matrix;
    }
    
    if (control_input->pressure_constraint_flag == 1) mat->accumulate_virial_constraint_matrix_element = insert_sparse_matrix_virial_element;
    
   	if (control_input->bootstrapping_flag == 1) {
    	 mat->finish_fm = average_sparse_bootstrapping_solutions;
    } else {
	    mat->finish_fm = average_sparse_block_fm_solutions;
	}
	
    // Determine the size of the matrix blocks from model specifications
	determine_matrix_columns_and_rows(mat, cg, control_input->frames_per_traj_block, control_input->pressure_constraint_flag);
        
	// Also, determine maximum and minimium number of entries in normal form matrix
	estimate_number_of_sparse_elements(mat, cg);
	printf("Maximum number of non-zero entries is %d, Sparsity is at least %.2lf percent\n", mat->max_nonzero_normal_elements, 100.0 * (1.0 -  (double) mat->max_nonzero_normal_elements / (double) (mat->fm_matrix_columns * mat->fm_matrix_columns)) );
		
    // Check that the matrix dimensions are enough that that the equations
    // will be overdetermined (in a perfect world where all the data is 
    // linearly independent for each row).
    if (mat->fm_matrix_rows < mat->fm_matrix_columns) {
        printf("The current sparse matrix computation will lead to an underdetermined set of force matching equations. Increase the block size or use a trajectory with more particles.\n");
        exit(EXIT_FAILURE);
    }
    mat->accumulation_matrix_columns = mat->fm_matrix_columns;
    mat->accumulation_matrix_rows = mat->accumulation_matrix_rows;
 
    printf("Number of rows for dense matrix algorithm: %d \n", mat->fm_matrix_rows);
    printf("Number of columns for dense matrix algorithm: %d \n", mat->fm_matrix_columns);
 
    // Check that the memory usage is reasonable and print 
    // memory diagnostics if so. These are checks for integer 
    // overflow when calculating the size of the matrices.

    if ( (int(INT_MAX) / mat->fm_matrix_columns) <
        (mat->fm_matrix_columns * (int)(sizeof(double)))) {
        printf("Using this number of columns will lead to integer overflow in memory allocation for the normal matrix equations. Decrease the number of basis functions.\n");
        exit(EXIT_FAILURE);
    }
    
    if (mat->fm_matrix_columns > mat->fm_matrix_rows) {
    	printf("Please use a larger block size since the number of fm_matrix_columns (%d) is bigger than the fm_matrix_rows (%d).\n", mat->fm_matrix_columns, mat->fm_matrix_rows);
    	exit(EXIT_FAILURE);
    }
    
    // Allocate memory for the FM matrix in linked list format and a dense target 
    // vector as well as temp space for the solution routines and final 
    // solution averaging operation.
    mat->dense_fm_rhs_vector = new double[mat->fm_matrix_rows]();
    mat->ll_sparse_matrix_row_heads = new linked_list_sparse_matrix_row_head[mat->rows_less_virial_constraint_rows];
	for(int i = 0; i < mat->rows_less_virial_constraint_rows; i++) {
    	mat->ll_sparse_matrix_row_heads[i].n = 0;
    	mat->ll_sparse_matrix_row_heads[i].h = NULL;
    }
    if (control_input->pressure_constraint_flag == 1) mat->dense_fm_matrix = new dense_matrix(control_input->frames_per_traj_block, mat->fm_matrix_columns);
    
    // Allocate a preconditioning temp array.
    // Read "Solving lest square problems" by Lawson CL and Hanson RJ, 
    // Chapt 25 for details.
    mat->h = new double[mat->fm_matrix_columns]();
  
	// These matrices would be used for acculation of normal form before solving
	// Allocate and initialize space for the solution to a block's equations
	if (control_input->bootstrapping_flag == 1) {
		mat->bootstrap_solutions = new std::vector<double>[mat->bootstrapping_num_estimates];
		for (int i = 0; i < control_input->bootstrapping_num_estimates; i++) {
			mat->bootstrap_solutions[i] = std::vector<double>(mat->fm_matrix_columns);
		}
	}
	mat->fm_solution = std::vector<double>(mat->fm_matrix_columns);

    mat->block_fm_solution = new double[mat->fm_matrix_columns]();
    mat->fm_solution_normalization_factors = new double[mat->fm_matrix_columns]();
    printf("Initialized a sparse FM matrix.\n");
}

// Initialize a sparse-matrix-accumulation for dense normal form computation.

void initialize_sparse_dense_normal_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg)
{
    // Set pseudopolymorphic methods
    mat->set_fm_matrix_to_zero = set_sparse_accumulation_matrix_to_zero;
    mat->accumulate_fm_matrix_element = insert_sparse_matrix_element;
    mat->accumulate_target_force_element = accumulate_force_into_dense_target_vector;
    mat->accumulate_target_constraint_element = accumulate_constraint_into_dense_target_vector;
    mat->sparse_matrix = NULL;
    
    if (control_input->scalar_matching_flag == 1) {
    	mat->accumulate_fm_matrix_element = insert_scalar_sparse_matrix_element;
    }

    if (control_input->bootstrapping_flag == 1) {
    	mat->do_end_of_frameblock_matrix_manipulations = convert_sparse_fm_equation_to_dense_normal_form_and_bootstrap;
    } else {
    	mat->do_end_of_frameblock_matrix_manipulations = convert_sparse_fm_equation_to_dense_normal_form_and_accumulate;
    }
    
    if (control_input->pressure_constraint_flag == 1) mat->accumulate_virial_constraint_matrix_element = insert_sparse_matrix_virial_element;
    
    	if (control_input->bootstrapping_flag == 1) {
    	mat->finish_fm = solve_dense_fm_normal_bootstrapping_equations;
    } else {
		mat->finish_fm = solve_dense_fm_normal_equations;
	}
	
    // Determine the size of the matrix blocks from model specifications
    determine_matrix_columns_and_rows(mat, cg, control_input->frames_per_traj_block, control_input->pressure_constraint_flag);
	
	// Also, determine maximum and minimium number of entries in normal form matrix
	estimate_number_of_sparse_elements(mat, cg);
	printf("Fully dense normal matrix has %d entries.\n", mat->fm_matrix_columns * mat->fm_matrix_columns);
	printf("Maximum number of non-zero entries is %d; sparsity is at least %.2lf percent\n", mat->max_nonzero_normal_elements, 100.0 * (1.0 -  (double) mat->max_nonzero_normal_elements / (double) (mat->fm_matrix_columns * mat->fm_matrix_columns)) );
	printf("Lower bound for number of non-zero entries is %d; sparsity is at least %.2lf percent.\n", mat->min_nonzero_normal_elements, 100.0 * (1.0 -  (double) mat->min_nonzero_normal_elements / (double) (mat->fm_matrix_columns * mat->fm_matrix_columns)) );
	
    // Check that the matrix dimensions are enough that that the equations
    // will be overdetermined (in a perfect world where all the data is 
    // linearly independent for each row).
    if ( (unsigned)(mat->fm_matrix_rows / mat->frames_per_traj_block) * (unsigned)(control_input->n_frames) < (unsigned)(mat->fm_matrix_columns) ) {
        printf("Current number of frames in this trajectory is too low to provide a fully-determined set of FM equations. Provide more frames in the input trajectory.\n");
        exit(EXIT_FAILURE);
    }
    if (mat->fm_matrix_columns > mat->fm_matrix_rows) {
    	printf("Please use a larger block size since the number of fm_matrix_columns (%d) is bigger than the fm_matrix_rows (%d).\n", mat->fm_matrix_columns, mat->fm_matrix_rows);
    	exit(EXIT_FAILURE);
    }
    
    mat->accumulation_matrix_columns = mat->fm_matrix_columns;
    mat->accumulation_matrix_rows = mat->fm_matrix_rows;
 
    printf("Number of rows for dense matrix algorithm: %d \n", mat->fm_matrix_rows);
    printf("Number of columns for dense matrix algorithm: %d \n", mat->fm_matrix_columns);
 
    // Check that the memory usage is reasonable and print 
    // memory diagnostics if so. These are checks for integer 
    // overflow when calculating the size of the matrices.

    if ( (int(INT_MAX) / mat->fm_matrix_columns) <
        (mat->fm_matrix_columns * (int)(sizeof(double))) ) {
        printf("Using this number of columns will lead to integer overflow in memory allocation for the normal matrix equations. Decrease the number of basis functions.\n");
        exit(EXIT_FAILURE);
    }
    
    printf("Size of dense normal matrix: %lu bytes \n", mat->fm_matrix_columns * mat->fm_matrix_columns * sizeof(double));

    // Allocate memory for the FM matrix in linked list format and a dense target 
    // vector as well as temp space for the solution routines and final 
    // solution averaging operation.
    mat->dense_fm_rhs_vector = new double[mat->fm_matrix_rows]();
    mat->ll_sparse_matrix_row_heads = new linked_list_sparse_matrix_row_head[mat->rows_less_virial_constraint_rows];
	for(int i = 0; i < mat->rows_less_virial_constraint_rows; i++) {
    	mat->ll_sparse_matrix_row_heads[i].n = 0;
    	mat->ll_sparse_matrix_row_heads[i].h = NULL;
    }
    if (control_input->pressure_constraint_flag == 1) mat->dense_fm_matrix = new dense_matrix(control_input->frames_per_traj_block, mat->fm_matrix_columns);
	else mat->dense_fm_matrix = new dense_matrix(1, 1); // This is to line-up with memory allocation in solve_dense_matrix
	
    // These matrices are used for accumulation of normal form before solving
    if (control_input->bootstrapping_flag == 1) {
    	mat->bootstrapping_dense_fm_normal_rhs_vectors = new double*[control_input->bootstrapping_num_estimates];
    	mat->bootstrapping_dense_fm_normal_matrices = new dense_matrix*[control_input->bootstrapping_num_estimates];
    	for (int i = 0; i < control_input->bootstrapping_num_estimates; i++) {
    		mat->bootstrapping_dense_fm_normal_rhs_vectors[i] = new double[mat->fm_matrix_columns]();
			mat->bootstrapping_dense_fm_normal_matrices[i] = new dense_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns);
		}
    }
	mat->dense_fm_normal_matrix = new dense_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns);
	mat->dense_fm_normal_rhs_vector = new double[mat->fm_matrix_columns]();
	printf("Initialized a sparse normal FM matrix.\n");
}

// Initialize a sparse-matrix-accumulation for sparse normal form computation.

void initialize_sparse_sparse_normal_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg)
{
    // Set pseudopolymorphic methods
    mat->set_fm_matrix_to_zero = set_sparse_accumulation_matrix_to_zero;
    mat->accumulate_fm_matrix_element = insert_sparse_matrix_element;
    mat->accumulate_target_force_element = accumulate_force_into_dense_target_vector;
    mat->accumulate_target_constraint_element = accumulate_constraint_into_dense_target_vector;
	mat->sparse_matrix = NULL;
	
    if (control_input->scalar_matching_flag == 1) {
    	mat->accumulate_fm_matrix_element = insert_scalar_sparse_matrix_element;
    }
    

   	if (control_input->bootstrapping_flag == 1) {
    	mat->do_end_of_frameblock_matrix_manipulations = convert_sparse_fm_equation_to_sparse_normal_form_and_bootstrap;
    } else {
    	mat->do_end_of_frameblock_matrix_manipulations = convert_sparse_fm_equation_to_sparse_normal_form_and_accumulate;
	}
	
    if (control_input->pressure_constraint_flag == 1) mat->accumulate_virial_constraint_matrix_element = insert_sparse_matrix_virial_element;
    
    if (control_input->bootstrapping_flag == 1) {
    	mat->finish_fm = solve_sparse_fm_bootstrapping_equations;
    } else {
	    mat->finish_fm = solve_sparse_fm_normal_equations;
	}
	
	mat->sparse_safety_factor = control_input->sparse_safety_factor;

    // Determine the size of the matrix blocks from model specifications
    determine_matrix_columns_and_rows(mat, cg, control_input->frames_per_traj_block, control_input->pressure_constraint_flag);
	
	// Also, determine maximum and minimium number of entries in normal form matrix
	estimate_number_of_sparse_elements(mat, cg);
	printf("Fully dense normal matrix has %d entries.\n", mat->fm_matrix_columns * mat->fm_matrix_columns);
	printf("Maximum number of non-zero entries is %d; sparsity is at least %.2lf percent\n", mat->max_nonzero_normal_elements, 100.0 * (1.0 -  (double) mat->max_nonzero_normal_elements / (double) (mat->fm_matrix_columns * mat->fm_matrix_columns)) );
	printf("Lower bound for number of non-zero entries is %d; sparsity is at least %.2lf percent.\n", mat->min_nonzero_normal_elements, 100.0 * (1.0 -  (double) mat->min_nonzero_normal_elements / (double) (mat->fm_matrix_columns * mat->fm_matrix_columns)) );
    
    // Check that the matrix dimensions are enough that that the equations
    // will be overdetermined (in a perfect world where all the data is 
    // linearly independent for each row).
    if ( (unsigned)(mat->fm_matrix_rows / mat->frames_per_traj_block) * (unsigned)(control_input->n_frames) < (unsigned)(mat->fm_matrix_columns) ) {
        printf("Current number of frames in this trajectory is too low to provide a fully-determined set of FM equations. Provide more frames in the input trajectory.\n");
        exit(EXIT_FAILURE);
    }
    if (mat->fm_matrix_columns > mat->fm_matrix_rows) {
    	printf("Please use a larger block size since the number of fm_matrix_columns (%d) is bigger than the fm_matrix_rows (%d).\n", mat->fm_matrix_columns, mat->fm_matrix_rows);
    	exit(EXIT_FAILURE);
    }
    
    mat->accumulation_matrix_columns = mat->fm_matrix_columns;
    mat->accumulation_matrix_rows = mat->accumulation_matrix_rows;
 
    printf("Number of rows for dense matrix algorithm: %d \n", mat->fm_matrix_rows);
    printf("Number of columns for dense matrix algorithm: %d \n", mat->fm_matrix_columns);
 
    // Check that the memory usage is reasonable and print 
    // memory diagnostics if so. These are checks for integer 
    // overflow when calculating the size of the matrices.

    if ( (int(INT_MAX) / mat->fm_matrix_columns) <
        (mat->fm_matrix_columns * (int)(sizeof(double))) ) {
        printf("Using this number of columns will lead to integer overflow in memory allocation for the normal matrix equations. Decrease the number of basis functions.\n");
        exit(EXIT_FAILURE);
    }
    
    printf("Size of dense normal matrix: %lu bytes \n", mat->fm_matrix_columns * mat->fm_matrix_columns * sizeof(double));

    // Allocate memory for the FM matrix in linked list format and a dense target 
    // vector as well as temp space for the solution routines and final 
    // solution averaging operation.
    mat->dense_fm_rhs_vector = new double[mat->fm_matrix_rows]();
    mat->ll_sparse_matrix_row_heads = new linked_list_sparse_matrix_row_head[mat->rows_less_virial_constraint_rows];
	for(int i = 0; i < mat->rows_less_virial_constraint_rows; i++) {
    	mat->ll_sparse_matrix_row_heads[i].n = 0;
    	mat->ll_sparse_matrix_row_heads[i].h = NULL;
    }
    if (control_input->pressure_constraint_flag == 1) mat->dense_fm_matrix = new dense_matrix(control_input->frames_per_traj_block, mat->fm_matrix_columns);

    // Set up the temps and parameters for solving the sparse normal equations.
   	if (control_input->bootstrapping_flag == 1) {
		mat->bootstrap_solutions = new std::vector<double>[mat->bootstrapping_num_estimates];
		for (int i = 0; i < control_input->bootstrapping_num_estimates; i++) {
			mat->bootstrap_solutions[i] = std::vector<double>(mat->fm_matrix_columns);
		}
	}
	mat->fm_solution = std::vector<double>(mat->fm_matrix_columns);
    
    // Allocate a preconditioning temp array.
    // Read "Solving lest square problems" by Lawson CL and Hanson RJ, 
    // Chapt 25 for details.
    mat->h = new double[mat->fm_matrix_columns]();
    
	// These matrices are used for acculation of sparse normal form before solving
	if (control_input->bootstrapping_flag == 1) {
   		mat->bootstrapping_sparse_fm_normal_matrices = new csr_matrix*[mat->bootstrapping_num_estimates];
   		mat->bootstrapping_dense_fm_normal_rhs_vectors = new double*[mat->bootstrapping_num_estimates]();
   		
   		for (int i = 0; i < mat->bootstrapping_num_estimates; i++) {
   			mat->bootstrapping_sparse_fm_normal_matrices[i] = new csr_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, mat->max_nonzero_normal_elements);
			mat->bootstrapping_dense_fm_normal_rhs_vectors[i] = new double[mat->fm_matrix_columns]();
   		}
	}
	mat->dense_fm_normal_rhs_vector = new double[mat->fm_matrix_columns]();
    mat->sparse_matrix = new csr_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, mat->max_nonzero_normal_elements);
    
	printf("Initialized a sparse-sparse normal FM matrix.\n");
}

void initialize_rem_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg)
{
  //make sure to overwrite any functions pointers set in constructor
  mat->set_fm_matrix_to_zero = set_rem_matrix_to_zero;
  //  mat->accumulate_fm_matrix_element = insert_rem_matrix_element;
  mat->accumulate_matching_forces = accumulate_entropy_elements;
  mat->accumulate_tabulated_forces = accumulate_tabulated_entropy;

  determine_matrix_columns_and_rows(mat, cg, control_input->frames_per_traj_block, control_input->pressure_constraint_flag);

  mat->fm_matrix_rows = 2;

  printf("Number of rows for dense matrix algorithm: %d \n", mat->fm_matrix_rows);
  printf("Number of columns for dense matrix algorithm: %d \n", mat->fm_matrix_columns);

    mat->accumulation_matrix_columns = mat->fm_matrix_columns;
    mat->accumulation_matrix_rows = mat->fm_matrix_rows;
    mat->dense_fm_matrix = new dense_matrix(1, mat->fm_matrix_columns);
    mat->do_end_of_frameblock_matrix_manipulations = calculate_frame_average_and_add_to_normal_matrix;
    //mat->dense_fm_rhs_vector = new double[mat->fm_matrix_rows]();

    //mat->finish_fm = solve_rem_equation;

    mat->fm_solution = std::vector<double>(mat->fm_matrix_columns, 0);
    mat->previous_rem_solution = std::vector<double>(mat->fm_matrix_columns, 0);
 
    mat->temperature = control_input->temperature;
    mat->rem_chi = control_input->REM_iteration_step_size;
    mat->boltzmann = control_input->boltzmann;

    mat->dense_fm_normal_matrix = new dense_matrix(mat->fm_matrix_rows, mat->fm_matrix_columns);
  
}

// "Initialize" a dummy matrix.

void initialize_dummy_matrix(MATRIX_DATA* const mat, ControlInputs* const control_input, CG_MODEL_DATA* const cg) 
{
    mat->rows_less_virial_constraint_rows = 0;
    mat->virial_constraint_rows = 0;
    mat->dense_fm_rhs_vector = new double[DIMENSION * cg->n_cg_sites];
    mat->dense_fm_normal_rhs_vector = new double[1];
    mat->fm_solution = std::vector<double>(1);
    mat->do_end_of_frameblock_matrix_manipulations = do_nothing_to_fm_matrix;
    mat->set_fm_matrix_to_zero = set_dummy_matrix_to_zero;
    mat->temperature = control_input->temperature;
    mat->boltzmann = control_input->boltzmann;

}

void initialize_first_BI_matrix(MATRIX_DATA* const mat, CG_MODEL_DATA* const cg)
{
  determine_matrix_columns_and_rows(mat, cg, 1, 0);
  mat->fm_solution.resize(mat->fm_matrix_columns);
   
  mat->accumulate_matching_forces = accumulate_BI_elements;
  mat->accumulate_target_force_element = accumulate_scalar_into_dense_target_vector;
  
  // reset output files
  FILE* BI_matrix = fopen("BI_matrix.dat","w");
  FILE* BI_vector = fopen("BI_vector.dat","w");
  fclose(BI_matrix);
  fclose(BI_vector);
}

void initialize_next_BI_matrix(MATRIX_DATA* const mat, InteractionClassComputer* const icomp)
{
  determine_BI_interaction_rows_and_cols(mat, icomp);

  printf("%s matrix rows = %d, cols = %d\n",icomp->ispec->get_full_name().c_str(), mat->fm_matrix_rows, mat->fm_matrix_columns);fflush(stdout);
  
  // To store the dense RHS, this array only need to be of size mat->fm_matrix_rows.
  // However, the svd solver puts the solution vector into this array,
  // and  the solution vetor is of size mat->fm_matrix_columns.
  delete [] mat->dense_fm_rhs_vector;
  delete [] mat->dense_fm_normal_rhs_vector;

  // avoid allocation errors if this interactions is not actually active.
  if(mat->fm_matrix_columns == 0 || mat->fm_matrix_rows == 0) {
  	mat->dense_fm_normal_rhs_vector = new double[1]();
  	mat->dense_fm_rhs_vector = new double[1]();
  	mat->dense_fm_matrix =  new dense_matrix(1, 1);
  } else {
  	if (mat->fm_matrix_rows >= mat->fm_matrix_columns) {
		mat->dense_fm_normal_rhs_vector = new double[mat->fm_matrix_rows]();
  		mat->dense_fm_rhs_vector = new double[mat->fm_matrix_rows]();
	  } else {
  		mat->dense_fm_normal_rhs_vector = new double[mat->fm_matrix_columns]();
  		mat->dense_fm_rhs_vector = new double[mat->fm_matrix_columns]();
  	}
 	mat->dense_fm_matrix =  new dense_matrix(mat->fm_matrix_rows, mat->fm_matrix_columns);
  }
}  

//--------------------------------------------------------------------
// Bootstrapping helper routines
//--------------------------------------------------------------------

void set_bootstrapping_normalization(MATRIX_DATA* mat, double** const bootstrapping_weights, int const n_frames) 
{
	// Copy bootstrapping information from frame_source to mat.
	mat->bootstrapping_weights = bootstrapping_weights;

	// Allocate space for normalization factors.
	mat->bootstrapping_normalization = new double[mat->bootstrapping_num_estimates]();

	// Determine total frame weight for each estimate 
	double total_frame_weight;
	for(int i = 0; i < mat->bootstrapping_num_estimates; i++) {
		total_frame_weight = 0.0;
		for(int j = 0; j < n_frames; j++) {
			total_frame_weight += mat->bootstrapping_weights[i][j];
		}
		mat->bootstrapping_normalization[i] = 1.0 / total_frame_weight;
	}	
}

//--------------------------------------------------------------------
// Initialization helper routines
//--------------------------------------------------------------------

void log_n_basis_functions(InteractionClassSpec &ispec) {
    // Get the name.
    std::string name = ispec.get_full_name();
    // Capitalize it.
    name[0] = std::toupper(name[0]);
    // Print the name and number of associated basis functions.
    printf("%s: %d\n", name.c_str(), ispec.get_num_basis_func());
}

// Determine the size of the matrix blocks from model specifications
void determine_matrix_columns_and_rows( MATRIX_DATA* const mat, CG_MODEL_DATA* const cg, int const frames_per_traj_block, int const pressure_constraint_flag) 
{
	// Determine total number of columns by adding up all the columns for all classes of interaction.
	mat->fm_matrix_columns = 0;
	printf("Number of basis functions by interaction class:\n");
	std::list<InteractionClassSpec*>::iterator iclass_iterator;
	for(iclass_iterator=cg->iclass_list.begin(); iclass_iterator != cg->iclass_list.end(); iclass_iterator++) {
		mat->fm_matrix_columns += (*iclass_iterator)->get_num_basis_func();
		log_n_basis_functions(**(iclass_iterator));
	}    

	if (cg->three_body_nonbonded_interactions.class_subtype > 0) {
		mat->fm_matrix_columns += cg->three_body_nonbonded_interactions.get_num_basis_func();
		log_n_basis_functions(cg->three_body_nonbonded_interactions);
	}

    // Determine the number of rows by seeing the number of particles and the number of auxiliary scalar restraints,
	// then multiplying by the block size.
	mat->rows_less_virial_constraint_rows = cg->n_cg_sites * frames_per_traj_block;
    if (pressure_constraint_flag == 0) {
        mat->fm_matrix_rows = mat->rows_less_virial_constraint_rows * DIMENSION;
        mat->virial_constraint_rows = 0;
    } else {
        mat->fm_matrix_rows = mat->rows_less_virial_constraint_rows * DIMENSION + frames_per_traj_block;
        mat->virial_constraint_rows = frames_per_traj_block;
    }
}
    
// Estimate upper and lower bounds for the number of non-zero elements in normal matrix

void estimate_number_of_sparse_elements(MATRIX_DATA* const mat, CG_MODEL_DATA* const cg)
{
	int diff;
	std::list<InteractionClassSpec*>::iterator iclass_iterator;
	
	// Estimate maximum number of entries in normal form matrix.
	// This assumes that each interaction class is a dense block.
	mat->max_nonzero_normal_elements = 0;
	for(iclass_iterator=cg->iclass_list.begin(); iclass_iterator != cg->iclass_list.end(); iclass_iterator++) {
		mat->max_nonzero_normal_elements += (*iclass_iterator)->get_num_basis_func() * (*iclass_iterator)->get_num_basis_func();
	}
	
	if (cg->three_body_nonbonded_interactions.class_subtype > 0) mat->max_nonzero_normal_elements += cg->three_body_nonbonded_interactions.get_num_basis_func() * cg->three_body_nonbonded_interactions.get_num_basis_func();
    
	// Estimate a lower bound on the size of sparse normal form matrix.
	// This assumes that only each type-type interaction block in an interaction class is dense.
	mat->min_nonzero_normal_elements = 0;
	for(iclass_iterator=cg->iclass_list.begin(); iclass_iterator != cg->iclass_list.end(); iclass_iterator++) {
		for(int i = 1; i <= (*iclass_iterator)->n_to_force_match; i++) {
			diff = (*iclass_iterator)->interaction_column_indices[i] - (*iclass_iterator)->interaction_column_indices[i-1];
			mat->min_nonzero_normal_elements += diff * diff;
		}
	}
	
	if (cg->three_body_nonbonded_interactions.class_subtype > 0) {
		for(int i = 1; i <= cg->three_body_nonbonded_interactions.get_n_defined(); i++) {
			diff = cg->three_body_nonbonded_interactions.interaction_column_indices[i] - cg->three_body_nonbonded_interactions.interaction_column_indices[i-1];
			mat->min_nonzero_normal_elements += diff * diff;
		}
	}
}	
	
//--------------------------------------------------------------------
// Matrix reset routines
//--------------------------------------------------------------------

// Set all elements of a dense matrix to zero.

inline void set_dense_matrix_to_zero(MATRIX_DATA* const mat)
{
    mat->dense_fm_matrix->reset_matrix();
}

// Set all elements of a dense matrix to zero.
inline void set_rem_matrix_to_zero(MATRIX_DATA* const mat)
{
    mat->dense_fm_matrix->reset_matrix();
}

// Set all elements of a linked-list sparse matrix to zero.

inline void set_sparse_matrix_to_zero(MATRIX_DATA* const mat)
{
	// The row head and element information is cleared in convert_linked_list_to_csr_matrix.

    // Set the elements of the dense part of the matrix to zero.
	for (int k = 0; k < mat->virial_constraint_rows * mat->fm_matrix_columns; k++) {
        mat->dense_fm_matrix->values[k] = 0.0;
    }
}

// Set all elements of a linked-list sparse matrix to zero when accumulating normal matrix.

inline void set_sparse_accumulation_matrix_to_zero(MATRIX_DATA* const mat)
{
	// The row head and element information is cleared in convert_linked_list_to_csr_matrix.

    // Set the elements of the dense part of the matrix to zero.
   for (int k = 0; k < mat->virial_constraint_rows * mat->fm_matrix_columns; k++) {
        mat->dense_fm_matrix->values[k] = 0.0;
    }
}

// Set all elements of an accumulation matrix to zero.

void set_accumulation_matrix_to_zero(MATRIX_DATA* const mat)
{
    int k, l;
    for (k = 0; k < mat->accumulation_matrix_columns; k++) {
        for (l = 0; l < k; l++) {
            mat->dense_fm_matrix->values[l * mat->accumulation_matrix_rows + k] = 0.0;
        }
    }
    
    for (k = mat->accumulation_matrix_columns; k < mat->accumulation_matrix_rows; k++) {
        for (l = 0; l < mat->accumulation_matrix_columns; l++) {
            mat->dense_fm_matrix->values[l * mat->accumulation_matrix_rows + k] = 0.0;
        }
    }
}

void set_accumulation_matrix_to_zero(MATRIX_DATA* const mat, dense_matrix* const dense_fm_matrix)
{
    int k, l;
    for (k = 0; k < mat->accumulation_matrix_columns; k++) {
        for (l = 0; l < k; l++) {
            dense_fm_matrix->values[l * mat->accumulation_matrix_rows + k] = 0.0;
        }
    }
    
    for (k = mat->accumulation_matrix_columns; k < mat->accumulation_matrix_rows; k++) {
        for (l = 0; l < mat->accumulation_matrix_columns; l++) {
           dense_fm_matrix->values[l * mat->accumulation_matrix_rows + k] = 0.0;
        }
    }
}

void set_dummy_matrix_to_zero(MATRIX_DATA* const mat) {}

//---------------------------------------------------------------------
// Functions for adding forces on a template number of particles into
// the target vector of the matrix from their ids, derivatives, and
// a set of spline coefficients.
//---------------------------------------------------------------------

void accumulate_vector_one_body_tabulated_force(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat) 
{
    // Calculate the associated forces.
    // Use flat arrays for performance.
    std::array<double, DIMENSION> forces;
    for (int j = 0; j < DIMENSION; j++) {
        forces[j] = table_fn_val * derivatives[0][j];
    }
    // Load those forces into the target vector.
	mat->accumulate_target_force_element(mat, particle_ids[0] + info->current_frame_starting_row, &forces[0]);
}

void accumulate_vector_tabulated_forces(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat) 
{
    // Calculate the associated forces.
    // Use flat arrays for performance.
    std::vector<double> forces(DIMENSION * n_body);
    for (int j = 0; j < DIMENSION; j++) forces[DIMENSION * (n_body - 1) + j] = 0.0;
    for (int i = 0; i < n_body - 1; i++) {
        for (int j = 0; j < DIMENSION; j++) {
            forces[DIMENSION * i + j] = table_fn_val * derivatives[i][j];
            forces[DIMENSION * (n_body - 1) + j] += -table_fn_val * derivatives[i][j];
        }
    }
    // Load those forces into the target vector.
    for (int i = 0; i < n_body; i++) {
        mat->accumulate_target_force_element(mat, particle_ids[i] + info->current_frame_starting_row, &forces[DIMENSION * i]);
    }
}

void accumulate_vector_symmetric_tabulated_forces(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat) 
{
    // Calculate the associated forces.
    // Use flat arrays for performance.
    std::vector<double> forces(DIMENSION * n_body);
    for (int j = 0; j < DIMENSION; j++) forces[DIMENSION * (n_body - 1) + j] = 0.0;
    for (int i = 0; i < n_body - 1; i++) {
        for (int j = 0; j < DIMENSION; j++) {
            forces[DIMENSION * i + j] = table_fn_val * fabs(derivatives[i][j]);
            forces[DIMENSION * (n_body - 1) + j] += table_fn_val * fabs(derivatives[i][j]);
        }
    }
    // Load those forces into the target vector.
    for (int i = 0; i < n_body; i++) {
        mat->accumulate_target_force_element(mat, particle_ids[i] + info->current_frame_starting_row, &forces[DIMENSION * i]);
    }
}

void accumulate_vector_one_body_force(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat) 
{
    // For each basis function,
    int ref_column = info->interaction_class_column_index + info->ispec->interaction_column_indices[info->index_among_matched_interactions - 1] + first_nonzero_basis_index;
    std::array<double, DIMENSION> forces;
    for (unsigned k = 0; k < basis_fn_vals.size(); k++) {
        // Calculate the associated forces.
        // Use flat force array for performance.
        for (int j = 0; j < DIMENSION; j++) {
 		    forces[j] = -basis_fn_vals[k] * derivatives[0][j];
            }
        // Load those forces into the target vector.
        (*mat->accumulate_fm_matrix_element)(particle_ids[0] + info->current_frame_starting_row, ref_column + k, &forces[0], mat);
    }
}

void accumulate_vector_matching_forces(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat) 
{
    // For each basis function,
    int ref_column = info->interaction_class_column_index + info->ispec->interaction_column_indices[info->index_among_matched_interactions - 1] + first_nonzero_basis_index;
    if (ref_column >= mat->fm_matrix_columns -  1) {
		printf("ref_col %d, fm_matrix_columns %d\n", ref_column, mat->fm_matrix_columns);
		fflush(stdout);    
    }
    std::vector<double> forces(DIMENSION * n_body);
    for (unsigned k = 0; k < basis_fn_vals.size(); k++) {
        // Calculate the associated forces.
        // Use flat force array for performance.
        for (int j = 0; j < DIMENSION; j++) forces[DIMENSION * (n_body - 1) + j] = 0.0;
        for (int i = 0; i < n_body - 1; i++) {
            for (int j = 0; j < DIMENSION; j++) {
                forces[DIMENSION * i + j] = -basis_fn_vals[k] * derivatives[i][j];
                forces[DIMENSION * (n_body - 1) + j] += basis_fn_vals[k] * derivatives[i][j];
            }
        }
        // Load those forces into the target vector.
        for (int i = 0; i < n_body; i++) {
            (*mat->accumulate_fm_matrix_element)(particle_ids[i] + info->current_frame_starting_row, ref_column + k, &forces[DIMENSION * i], mat);
        }
    }
}

void accumulate_vector_symmetric_matching_forces(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat) 
{
    // For each basis function,
    int ref_column = info->interaction_class_column_index + info->ispec->interaction_column_indices[info->index_among_matched_interactions - 1] + first_nonzero_basis_index;
    std::vector<double> forces(DIMENSION * n_body);
    for (unsigned k = 0; k < basis_fn_vals.size(); k++) {
        // Calculate the associated forces.
        // Use flat force array for performance.
        for (int j = 0; j < DIMENSION; j++) forces[DIMENSION * (n_body - 1) + j] = 0.0;
        for (int i = 0; i < n_body - 1; i++) {
            for (int j = 0; j < DIMENSION; j++) {
                forces[DIMENSION * i + j] = basis_fn_vals[k] * fabs(derivatives[i][j]);
                forces[DIMENSION * (n_body - 1) + j] += basis_fn_vals[k] * fabs(derivatives[i][j]);
            }
        }
        // Load those forces into the target vector.
        for (int i = 0; i < n_body; i++) {
            (*mat->accumulate_fm_matrix_element)(particle_ids[i] + info->current_frame_starting_row, ref_column + k, &forces[DIMENSION * i], mat);
        }
    }
}

void accumulate_entropy_elements(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat)
{
  
    int ref_column = info->interaction_class_column_index + info->ispec->interaction_column_indices[info->index_among_matched_interactions - 1] + first_nonzero_basis_index; 

   for (unsigned k = 0; k < basis_fn_vals.size(); k++) {
     insert_rem_matrix_element(ref_column + k, basis_fn_vals[k], mat);    
   }
  
}

void accumulate_tabulated_entropy(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat) 
{
  printf("Tabulated interactions cannot be done through REM framework. Please remove tabulated interactions from rmin files.\n");
  exit(EXIT_FAILURE);
}

void accumulate_BI_elements(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat)
{

  int ref_column = info->interaction_class_column_index + info->ispec->interaction_column_indices[info->index_among_matched_interactions - 1] + first_nonzero_basis_index;

  for (unsigned k = 0; k < basis_fn_vals.size(); k++) {
	insert_BI_matrix_element(n_body, ref_column + k, basis_fn_vals[k], mat);
  }
}

//---------------------------------------------------------------------
// Scalar versions of the above functions (note symmetric and regular 
// interactions are the same since the derivative is ignored).
//---------------------------------------------------------------------

void accumulate_scalar_one_body_tabulated_force(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat) 
{
    // Calculate the associated forces.
    // Use flat arrays for performance.
    std::array<double, 1> forces;
    forces[0] = table_fn_val;

    // Load this force into the target vector.
	mat->accumulate_target_force_element(mat, particle_ids[0] + info->current_frame_starting_row, &forces[0]);
}

void accumulate_scalar_tabulated_forces(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat) 
{
    // Calculate the associated forces.
    // Use flat arrays for performance.
    std::vector<double> forces(n_body);
    forces[(n_body - 1)] = 0.0;
    for (int i = 0; i < n_body - 1; i++) {
        forces[i] = table_fn_val;
        forces[(n_body - 1)] += -table_fn_val;
    }
    // Load those forces into the target vector.
    for (int i = 0; i < n_body; i++) {
        mat->accumulate_target_force_element(mat, particle_ids[i] + info->current_frame_starting_row, &forces[i]);
    }
}

void accumulate_scalar_one_body_force(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat) 
{
    // For each basis function,
    int ref_column = info->interaction_class_column_index + info->ispec->interaction_column_indices[info->index_among_matched_interactions - 1] + first_nonzero_basis_index;
    std::array<double, 1> forces;
    for (unsigned k = 0; k < basis_fn_vals.size(); k++) {
        // Calculate the associated forces.
        // Use flat force array for performance.
 	    forces[0] = -basis_fn_vals[k];
        // Load those forces into the target vector.
        (*mat->accumulate_fm_matrix_element)(particle_ids[0] + info->current_frame_starting_row, ref_column + k, &forces[0], mat);
    }
}

void accumulate_scalar_matching_forces(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat) 
{
    // For each basis function,
    int ref_column = info->interaction_class_column_index + info->ispec->interaction_column_indices[info->index_among_matched_interactions - 1] + first_nonzero_basis_index;

    std::vector<double> forces(n_body);
    for (unsigned k = 0; k < basis_fn_vals.size(); k++) {
        // Calculate the associated forces.
        // Use flat force array for performance.
        forces[(n_body - 1)] = 0.0;
        for (int i = 0; i < n_body - 1; i++) {
            forces[i] = -basis_fn_vals[k];
            forces[(n_body - 1)] += basis_fn_vals[k];
        }
        // Load those forces into the target vector.
        for (int i = 0; i < n_body; i++) {
            (*mat->accumulate_fm_matrix_element)(particle_ids[i] + info->current_frame_starting_row, ref_column + k, &forces[i], mat);
        }
    }
}

//--------------------------------------------------------------------
// Matrix insertion routines
//--------------------------------------------------------------------

// Add a three-component nonzero force value to a linked list format sparse matrix.

void insert_sparse_matrix_element(const int i, const int j, double* const x, MATRIX_DATA* const mat)
{
    struct linked_list_sparse_matrix_element* curr_elem, *prev_elem, *head, *pt;
    int k;
    
    head = mat->ll_sparse_matrix_row_heads[i].h;
    
    // If the linked list is empty
    if (head == NULL) {
        pt = new linked_list_sparse_matrix_element;
        mat->ll_sparse_matrix_row_heads[i].h = pt;
        pt->col = j;
        for (k = 0; k < DIMENSION; k++) pt->valx[k] = x[k];
        pt->next = NULL;
        mat->ll_sparse_matrix_row_heads[i].n += 1;
        
    // If the linked list is not empty
    } else {
        curr_elem = head;
        prev_elem = NULL;
        while (curr_elem != NULL) {
            // If the element exists, add
            if (curr_elem->col == j) {
                for (k = 0; k < DIMENSION; k++) curr_elem->valx[k] += x[k];
                prev_elem = curr_elem;
                break;
            
            // If the new element should be ahead of the next in the list
            } else if (curr_elem->col > j) {
                
                // If the new element should be the first in the list
                if (prev_elem == NULL) {
                    pt = new linked_list_sparse_matrix_element;
                    mat->ll_sparse_matrix_row_heads[i].h = pt;
                    pt->col = j;
                    for (k = 0; k < DIMENSION; k++) pt->valx[k] = x[k];
                    pt->next = curr_elem;

                // General case
                } else {
                    pt = new linked_list_sparse_matrix_element;
                    prev_elem->next = pt;
                    pt->col = j;
                    for (k = 0; k < DIMENSION; k++) pt->valx[k] = x[k];
                    pt->next = curr_elem;
                }
                mat->ll_sparse_matrix_row_heads[i].n += 1;
                prev_elem = curr_elem;
                break;
            }
            prev_elem = curr_elem;
            curr_elem = prev_elem->next;
        }
        
        // If the new element should be the last
        if (prev_elem->col < j) {
            pt = new linked_list_sparse_matrix_element;
            prev_elem->next = pt;
            pt->col = j;
            for (k = 0; k < DIMENSION; k++) pt->valx[k] = x[k];
            pt->next = NULL;
            mat->ll_sparse_matrix_row_heads[i].n += 1;
        }
    }   
}

// Add a nonzero scalar value to a linked list format sparse matrix.

void insert_scalar_sparse_matrix_element(const int i, const int j, double* const x, MATRIX_DATA* const mat)
{
    struct linked_list_sparse_matrix_element* curr_elem, *prev_elem, *head, *pt;    
    head = mat->ll_sparse_matrix_row_heads[i].h;
    
    // If the linked list is empty
    if (head == NULL) {
        pt = new linked_list_sparse_matrix_element;
        mat->ll_sparse_matrix_row_heads[i].h = pt;
        pt->col = j;
        pt->valx[0] = x[0];
        pt->next = NULL;
        mat->ll_sparse_matrix_row_heads[i].n += 1;
        
    // If the linked list is not empty
    } else {
        curr_elem = head;
        prev_elem = NULL;
        while (curr_elem != NULL) {
            // If the element exists, add
            if (curr_elem->col == j) {
                curr_elem->valx[0] += x[0];
                prev_elem = curr_elem;
                break;
            
            // If the new element should be ahead of the next in the list
            } else if (curr_elem->col > j) {
                
                // If the new element should be the first in the list
                if (prev_elem == NULL) {
                    pt = new linked_list_sparse_matrix_element;
                    mat->ll_sparse_matrix_row_heads[i].h = pt;
                    pt->col = j;
                    pt->valx[0] = x[0];
                    pt->next = curr_elem;

                // General case
                } else {
                    pt = new linked_list_sparse_matrix_element;
                    prev_elem->next = pt;
                    pt->col = j;
                    pt->valx[0] = x[0];
                    pt->next = curr_elem;
                }
                mat->ll_sparse_matrix_row_heads[i].n += 1;
                prev_elem = curr_elem;
                break;
            }
            prev_elem = curr_elem;
            curr_elem = prev_elem->next;
        }
        
        // If the new element should be the last
        if (prev_elem->col < j) {
            pt = new linked_list_sparse_matrix_element;
            prev_elem->next = pt;
            pt->col = j;
            pt->valx[0] = x[0];
            pt->next = NULL;
            mat->ll_sparse_matrix_row_heads[i].n += 1;
        }
    }   
}
// Add a dimension-sized force element to a dense matrix.

inline void insert_dense_matrix_element(const int i, const int j, double* const x, MATRIX_DATA* const mat)
{
    mat->dense_fm_matrix->add_vector(mat->size_per_vector * i, j, x);
}

// Add a dimension-sized force element to an accumulation matrix.

inline void insert_accumulation_matrix_element(const int i, const int j, double* const x, MATRIX_DATA* const mat)
{
    mat->dense_fm_matrix->add_vector(mat->size_per_vector * i + mat->accumulation_row_shift, j, x);
}

inline void insert_rem_matrix_element(const int i, double const x, MATRIX_DATA* const mat)
{
  mat->dense_fm_matrix->add_scalar(0,i,x);
}  

inline void insert_BI_matrix_element(const int i, const int j,double const x, MATRIX_DATA* const mat)
{
  mat->dense_fm_matrix->add_scalar(i,j,x);
}

// Add a scalar force element to a dense matrix.

inline void insert_scalar_dense_matrix_element(const int i, const int j, double* const x, MATRIX_DATA* const mat)
{
    mat->dense_fm_matrix->add_scalar(mat->size_per_vector * i, j, x[0]);
}

// Add a scalar force element to an accumulation matrix.

inline void insert_scalar_accumulation_matrix_element(const int i, const int j, double* const x, MATRIX_DATA* const mat)
{
    mat->dense_fm_matrix->add_scalar(mat->size_per_vector * i + mat->accumulation_row_shift, j, x[0]);
}

// Add a one-component virial contribution element to a dense matrix held 
// alongside a sparse matrix; The virial rows will be dense, so only the
// force rows are compressed.

// Add a scalar virial contribution element to a dense matrix.

inline void insert_dense_matrix_virial_element(const int m, const int n, const double x, MATRIX_DATA* const mat)
{
    mat->dense_fm_matrix->add_scalar(mat->rows_less_virial_constraint_rows * mat->size_per_vector, n, x);

}

// Add a scalar virial contribution to a sparse matrix.

inline void insert_sparse_matrix_virial_element(const int m, const int n, const double x, MATRIX_DATA* const mat)
{
    mat->dense_fm_matrix->add_scalar(m, n, x);
}

// Add a one-component virial contribution element to an accumulation matrix.

inline void insert_accumulation_matrix_virial_element(const int m, const int n, const double x, MATRIX_DATA* const mat)
{
    mat->dense_fm_matrix->add_scalar(mat->rows_less_virial_constraint_rows * mat->size_per_vector + mat->accumulation_row_shift, n, x);
}

//--------------------------------------------------------------------
// Target calculation from trajectory data routines (move to trajectory)
//--------------------------------------------------------------------

void add_target_virials_from_trajectory(MATRIX_DATA* const mat, double *pressure_constraint_rhs_vector)
{
    if (mat->matrix_type == kDense || mat->matrix_type == kSparse) {
        calculate_target_virial_in_dense_vector(mat, pressure_constraint_rhs_vector);
    } else if (mat->matrix_type == kAccumulation) {
        calculate_target_virial_in_accumulation_vector(mat, pressure_constraint_rhs_vector);
    } else if (mat->matrix_type == kDummy) {
        return;
    }
}

// Append the target virials to the target force vector for a dense or sparse matrix
// calculation

void calculate_target_virial_in_dense_vector(MATRIX_DATA* const mat, double *pressure_constraint_rhs_vector)
{
	int frame_sample  =  mat->trajectory_block_index * mat->virial_constraint_rows;
    for (int k = 0; k < mat->virial_constraint_rows; k++) {
        mat->dense_fm_rhs_vector[mat->rows_less_virial_constraint_rows * DIMENSION + k] = pressure_constraint_rhs_vector[(int)((frame_sample + k)/ mat->dynamic_state_samples_per_frame)];
        frame_sample++;
    }
}

// Append the target virials to the target force vector for an accumulation matrix
// calculation

void calculate_target_virial_in_accumulation_vector(MATRIX_DATA* const mat, double *pressure_constraint_rhs_vector)
{
	int frame_sample  =  mat->trajectory_block_index * mat->virial_constraint_rows;
    for (int k = 0; k < mat->virial_constraint_rows; k++) {
        mat->dense_fm_matrix->values[mat->fm_matrix_columns * mat->accumulation_matrix_rows + mat->rows_less_virial_constraint_rows * DIMENSION + mat->accumulation_row_shift + k] = pressure_constraint_rhs_vector[(int)((frame_sample + k) / mat->dynamic_state_samples_per_frame)];
    	frame_sample++;
    }
}

void add_target_force_from_trajectory(int shift_i, int site_i, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &f) 
{
    if (mat->matrix_type == kDense || mat->matrix_type == kSparse || mat->matrix_type == kSparseNormal || mat->matrix_type == kSparseSparse) {
        calculate_target_force_dense_vector(shift_i, site_i, mat, f);
    } else if (mat->matrix_type == kAccumulation) {
        calculate_target_force_accumulation_vector(shift_i, site_i, mat, f);
    }
}

// Calculate the RHS vector for dense or sparse matrix calculations

void calculate_target_force_dense_vector(int shift_i, int site_i, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &f)
{
    int tn = DIMENSION * (site_i + shift_i);
    double force_sq = 0.0;
    double curr_force;
    for (int i = 0; i < DIMENSION; i++) {
    	curr_force = f[site_i][i];
    	mat->dense_fm_rhs_vector[tn + i] = curr_force;
		force_sq += curr_force * curr_force;
	}
	mat->force_sq_total += force_sq;
}

// Calculate the RHS vector for accumulation matrix calculations

inline void calculate_target_force_accumulation_vector(int shift_i, int site_i, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &f)
{
	mat->dense_fm_matrix->assign_vector(DIMENSION * (site_i + shift_i) + mat->accumulation_row_shift, mat->fm_matrix_columns, f[site_i]);
}

//--------------------------------------------------------------------
// Target accumulation routines
//--------------------------------------------------------------------

void accumulate_force_into_dense_target_vector(MATRIX_DATA* mat, int particle_index, double* force_element) 
{
    for (int i = 0; i < DIMENSION; i++) {
        mat->dense_fm_rhs_vector[DIMENSION * particle_index + i] += force_element[i];
    }
}

void accumulate_force_into_accumulation_target_vector(MATRIX_DATA* mat, int particle_index, double* force_element)
{
	mat->dense_fm_matrix->add_vector(DIMENSION * particle_index, mat->accumulation_matrix_rows, force_element);
}

void accumulate_constraint_into_dense_target_vector(MATRIX_DATA* mat, int frame_index, double constraint_element)
{
    mat->dense_fm_rhs_vector[mat->rows_less_virial_constraint_rows * DIMENSION + frame_index] += constraint_element;
}

void accumulate_constraint_into_accumulation_target_vector(MATRIX_DATA* mat, int frame_index, double constraint_element)
{
	mat->dense_fm_matrix->add_scalar(DIMENSION * mat->rows_less_virial_constraint_rows + frame_index, mat->accumulation_matrix_rows, constraint_element); 
}

void accumulate_scalar_into_dense_target_vector(MATRIX_DATA* mat, int particle_index, double* force_element)
{
  mat->dense_fm_rhs_vector[particle_index] = *force_element;
}

//--------------------------------------------------------------------
// End-of-frame-block routines
//--------------------------------------------------------------------

// The dense matrix calculation uses single-frame blocks and proceeds
// by taking the normal form of each frame's MS-CG equations and
// adding all of those up frame by frame until the trajectory is
// exhausted.

void convert_dense_fm_equation_to_normal_form_and_accumulate(MATRIX_DATA* const mat)
{
	if (mat->output_raw_frame_blocks == 1) mat->dense_fm_matrix->print_matrix_csr(mat->frame_block_fh);
    double frame_weight = mat->get_frame_weight() * mat->normalization;
 	create_dense_normal_form(mat, frame_weight, mat->dense_fm_matrix,mat->dense_fm_normal_matrix, mat->dense_fm_rhs_vector, mat->dense_fm_normal_rhs_vector);
}

void convert_dense_fm_equation_to_normal_form_and_bootstrap(MATRIX_DATA* const mat)
{
	int onei = 1.0;
	int matrix_size = mat->fm_matrix_columns * mat->fm_matrix_columns;
	if (mat->output_raw_frame_blocks == 1) mat->dense_fm_matrix->print_matrix_csr(mat->frame_block_fh);

	// Create temp normal matrix and rhs vector.
	dense_matrix* temp_normal_matrix = new dense_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns);
	double* temp_normal_rhs_vector = new double[mat->fm_matrix_columns]();

	create_dense_normal_form(mat, 1.0, mat->dense_fm_matrix, temp_normal_matrix, mat->dense_fm_rhs_vector, temp_normal_rhs_vector);
	
	// Add the matrix to the master 
	double frame_weight = mat->get_frame_weight() * mat->normalization;
	cblas_daxpy( matrix_size, frame_weight, temp_normal_matrix->values, onei, mat->dense_fm_normal_matrix->values, onei);	    
	cblas_daxpy( mat->fm_matrix_columns, frame_weight, temp_normal_rhs_vector, onei, mat->dense_fm_normal_rhs_vector, onei);
	
	// Add the matrix and vector to each of the bootstrap samples based on the weight for that frame for each bootstrap estimate.
	for (int i = 0; i < mat->bootstrapping_num_estimates; i++) {
		
		// traj_block_frame_index is the frame number processed since bootstrapping only allows a block size of 1.
		// Skip processing frame for this estimate if weight is 0.0.
		frame_weight = mat->bootstrapping_weights[i][mat->trajectory_block_index];
		if(frame_weight == 0.0) continue;
		frame_weight *= mat->bootstrapping_normalization[i];

	   cblas_daxpy( matrix_size, frame_weight, temp_normal_matrix->values, onei, mat->bootstrapping_dense_fm_normal_matrices[i]->values, onei);

	   cblas_daxpy( mat->fm_matrix_columns, frame_weight, temp_normal_rhs_vector, onei, mat->bootstrapping_dense_fm_normal_rhs_vectors[i], onei);
	}
	
	delete temp_normal_matrix;
	delete [] temp_normal_rhs_vector;
}

// As above, but ignoring the FM matrix.
// Used for Lanyuan's iterative method, in which only the FM target vector is recalculated.

void convert_dense_target_force_vector_to_normal_form_and_accumulate(MATRIX_DATA* const mat)
{
    int onei = 1;
    double oned = 1.0;
    double frame_weight = mat->get_frame_weight();

    // Take normal form of the current frame's target vector and add to the existing normal form target vector.
    cblas_dgemv(CblasColMajor, CblasTrans, mat->fm_matrix_rows, mat->fm_matrix_columns, frame_weight, mat->dense_fm_matrix->values, mat->fm_matrix_rows, mat->dense_fm_rhs_vector, onei, oned, mat->dense_fm_normal_rhs_vector, onei);
}

// Perform the accumulation operation (QR decomposition followed by composition) to combine the
// current frame's FM matrix with the growing accumulation matrix.

void accumulate_accumulation_matrices(MATRIX_DATA* const mat)
{
    int info_in;
    if (mat->output_raw_frame_blocks == 1) mat->dense_fm_matrix->print_matrix_csr(mat->frame_block_fh);

    // Initialize the operation if this is the first block.
    if (mat->trajectory_block_index == 0) {
        mat->lapack_temp_workspace = new double[1];
        mat->lapack_setup_flag = -1;
        dgeqrf_(&mat->fm_matrix_rows, &mat->accumulation_matrix_columns, mat->dense_fm_matrix->values, &mat->accumulation_matrix_rows, mat->lapack_tau, mat->lapack_temp_workspace, &mat->lapack_setup_flag, &info_in);
        mat->lapack_setup_flag = mat->lapack_temp_workspace[0];
        delete [] mat->lapack_temp_workspace;
        mat->lapack_temp_workspace = new double[mat->lapack_setup_flag];
        dgeqrf_(&mat->fm_matrix_rows, &mat->accumulation_matrix_columns, mat->dense_fm_matrix->values, &mat->accumulation_matrix_rows, mat->lapack_tau, mat->lapack_temp_workspace, &mat->lapack_setup_flag, &info_in);
        mat->accumulation_row_shift = mat->accumulation_matrix_columns;
    } else {
        dgeqrf_(&mat->accumulation_matrix_rows, &mat->accumulation_matrix_columns, mat->dense_fm_matrix->values, &mat->accumulation_matrix_rows, mat->lapack_tau, mat->lapack_temp_workspace, &mat->lapack_setup_flag, &info_in);
    }
}

void accumulate_accumulation_matrices_for_bootstrap(MATRIX_DATA* const mat)
{
	printf("Bootstrapping is not implemented for accumulation matrices.\n");
	exit(EXIT_FAILURE);
	
    int info_in;
    if (mat->output_raw_frame_blocks == 1) mat->dense_fm_matrix->print_matrix_csr(mat->frame_block_fh);

    // Initialize the operation if this is the first block.
    if (mat->trajectory_block_index == 0) {
        mat->lapack_temp_workspace = new double[1];
        mat->lapack_setup_flag = -1;
        dgeqrf_(&mat->fm_matrix_rows, &mat->accumulation_matrix_columns, mat->dense_fm_matrix->values, &mat->accumulation_matrix_rows, mat->lapack_tau, mat->lapack_temp_workspace, &mat->lapack_setup_flag, &info_in);
        mat->lapack_setup_flag = mat->lapack_temp_workspace[0];
        delete [] mat->lapack_temp_workspace;
        mat->lapack_temp_workspace = new double[mat->lapack_setup_flag];
        dgeqrf_(&mat->fm_matrix_rows, &mat->accumulation_matrix_columns, mat->dense_fm_matrix->values, &mat->accumulation_matrix_rows, mat->lapack_tau, mat->lapack_temp_workspace, &mat->lapack_setup_flag, &info_in);
        mat->accumulation_row_shift = mat->accumulation_matrix_columns;
    } else {
        dgeqrf_(&mat->accumulation_matrix_rows, &mat->accumulation_matrix_columns, mat->dense_fm_matrix->values, &mat->accumulation_matrix_rows, mat->lapack_tau, mat->lapack_temp_workspace, &mat->lapack_setup_flag, &info_in);
    }
}

// The sparse matrix is solved after every block, and the number of times each basis function coefficient is nonzero
// is recorded to allow averaging of all the solutions after the entire trajectory is read.

void solve_sparse_matrix(MATRIX_DATA* const mat)
{
   double frame_weight = mat->get_frame_weight();
   solve_this_sparse_matrix(mat);
   
   // Add the solution to the sum of all block solutions and increment 
   // the necessary normalization factors.
   
   for (int k = 0; k < mat->fm_matrix_columns; k++) {
   	  mat->block_fm_solution[k] *= mat->h[k];
      if (mat->block_fm_solution[k] < MAX_INPUT_FORCE_VALUE
                && mat->block_fm_solution[k] > -MAX_INPUT_FORCE_VALUE) {
			if ((mat->block_fm_solution[k] > VERYSMALL) 
             || (mat->block_fm_solution[k] < -VERYSMALL)) { 
				mat->fm_solution_normalization_factors[k] += frame_weight;
				mat->fm_solution[k] += mat->block_fm_solution[k] * frame_weight * mat->normalization;
			}
		}
   }
}

void solve_sparse_matrix_for_bootstrap(MATRIX_DATA* const mat)
{
   double frame_weight = mat->get_frame_weight();
   solve_this_sparse_matrix(mat);

   // Add the solution to the sum of all block solutions for master.
   for (int k = 0; k < mat->fm_matrix_columns; k++) {
   	  mat->block_fm_solution[k] *= mat->h[k];
      if (mat->block_fm_solution[k] < MAX_INPUT_FORCE_VALUE
                && mat->block_fm_solution[k] > -MAX_INPUT_FORCE_VALUE) {
			if ((mat->block_fm_solution[k] > VERYSMALL) 
             || (mat->block_fm_solution[k] < -VERYSMALL)) { 
				mat->fm_solution_normalization_factors[k] += frame_weight;
				mat->fm_solution[k] += mat->block_fm_solution[k] * frame_weight * mat->normalization;
			}
		}
   }

   // Add the solution to the sum of all block solutions and increment 
   // the necessary normalization factors for bootstrap.
   for (int k = 0; k < mat->fm_matrix_columns; k++) {
   	  mat->block_fm_solution[k] *= mat->h[k];
      if (mat->block_fm_solution[k] < MAX_INPUT_FORCE_VALUE
                && mat->block_fm_solution[k] > -MAX_INPUT_FORCE_VALUE) {
			if ((mat->block_fm_solution[k] > VERYSMALL) 
             || (mat->block_fm_solution[k] < -VERYSMALL)) { 

				mat->fm_solution_normalization_factors[k] += 1.0;             
             	for(int i = 0; i < mat->bootstrapping_num_estimates; i++) {
					mat->bootstrap_solutions[i][k] += mat->block_fm_solution[k] * mat->bootstrapping_weights[i][mat->trajectory_block_index] * mat->bootstrapping_normalization[i];
				}
			}
		}
   }
}

// The sparse matrix is converted to sparse normal form after every block, and 
// the accumulated normal form equations are solved after the entire trajectory is read.

void convert_sparse_fm_equation_to_sparse_normal_form_and_accumulate(MATRIX_DATA* const mat)
{
    // Calculate the weight of this part of the normal equations in the overall equations
	double frame_weight = mat->get_frame_weight() * mat->normalization;

    // Convert from linked list format to CSR format
    // Note: These MKL functions use a one-based index for row_sizes and column_indices
    int n_nonzero_matrix_elements = get_n_nonzero_matrix_elements(mat);
    csr_matrix csr_fm_matrix(mat->fm_matrix_rows, mat->fm_matrix_columns, n_nonzero_matrix_elements);
    convert_linked_list_to_csr_matrix(mat, csr_fm_matrix);
	if (mat->output_raw_frame_blocks == 1) csr_fm_matrix.print_matrix_csr(mat->frame_block_fh);
	
   // Convert CSR matrix and dense RHS vector to normal-form    
   // Form sparse normal-form left-hand side matrix using mkl_dcsrmultcsr
   // rows of matrix is mat->fm_matrix_rows
   // cols of matrix is mat->fm_matrix_columns
   
   // Resize estimate based on the smaller of previous frame's normal matrix size using the safety factor or the maximum number of non-zero elements
   int nnzmax = mat->max_nonzero_normal_elements;
   // Determine size of previous matrix (or 0 if first frame)
   int new_size_est = mat->sparse_matrix->row_sizes[mat->fm_matrix_columns] * (1.0 + mat->sparse_safety_factor);
   // Determine which estimate is smaller (if new_size_est is reasonable)
   if ( (new_size_est > mat->min_nonzero_normal_elements) && (new_size_est < mat->max_nonzero_normal_elements) ) {
		nnzmax = new_size_est;
		printf("Estimating number of non-zero normal matrix elements as %d < %d\n", nnzmax, mat->max_nonzero_normal_elements);
	}
	
   // Allocate temporary dense right-hand side normal form vector.
   double* dense_rhs_normal_vector = new double[mat->fm_matrix_rows]();	// the rows of the normal-form vector is number of basis functions (number of columns in input matrix)

   // Allocate space for normal_form_matrix
   csr_matrix csr_normal_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, nnzmax); // the rows of normal matrix is number of basis functions (number of columns in input matrix)
	
   // Allocate space for normal_form_matrix.
   create_sparse_normal_form_matrix(mat, nnzmax, csr_fm_matrix, csr_normal_matrix, mat->dense_fm_rhs_vector, dense_rhs_normal_vector);	

   // Accumulate normal form matrix with previous/future normal form matrices
   // Frame weight is applied to normal matrix in this step
   sparse_matrix_addition(mat, frame_weight, nnzmax, csr_normal_matrix, mat->sparse_matrix);

   #if _mkl_flag == 1
   int onei = 1;	
   // Accumulate normal form right-hand size vector with previous/future vectors
   // Frame weight is applied to normal vector in this step
   cblas_daxpy(mat->fm_matrix_columns, frame_weight,
		dense_rhs_normal_vector, onei, mat->dense_fm_normal_rhs_vector, onei);
   #endif
		
   // CSR formatted FM and normal temp matrices are freed by destructor at end of function
   // Free the intermediate normal form matrix and vector
   delete [] dense_rhs_normal_vector;
}

void convert_sparse_fm_equation_to_sparse_normal_form_and_bootstrap(MATRIX_DATA* const mat)
{
	double frame_weight = 1.0;
    // Convert from linked list format to CSR format
    // Note: These MKL functions use a one-based index for row_sizes and column_indices
    int n_nonzero_matrix_elements = get_n_nonzero_matrix_elements(mat);
    csr_matrix csr_fm_matrix(mat->fm_matrix_rows, mat->fm_matrix_columns, n_nonzero_matrix_elements);
    convert_linked_list_to_csr_matrix(mat, csr_fm_matrix);
	if (mat->output_raw_frame_blocks == 1) csr_fm_matrix.print_matrix_csr(mat->frame_block_fh);
   
   // Convert CSR matrix and dense RHS vector to normal-form    
   // Form sparse normal-form left-hand side matrix using mkl_dcsrmultcsr
   // rows of matrix is mat->fm_matrix_rows
   // cols of matrix is mat->fm_matrix_columns
   int onei=1;

   // Resize estimate based on the smaller of previous frame's normal matrix size using the safety factor or the maximum number of non-zero elements
   int nnzmax = mat->max_nonzero_normal_elements;
   int new_size_est = mat->sparse_matrix->row_sizes[mat->fm_matrix_columns] * (1.0 + mat->sparse_safety_factor);;
   
   // Determine size of previous matrix (or 0 if first frame)
   // Get largest new_size_est of all matrices
   for (int i = 0; i < mat->bootstrapping_num_estimates; i++) {
   		int temp_size_est = mat->bootstrapping_sparse_fm_normal_matrices[i]->row_sizes[mat->fm_matrix_columns] * (1.0 + mat->sparse_safety_factor);
   		if (temp_size_est > new_size_est) new_size_est = new_size_est;
   	}
   	
   // Determine which estimate is smaller (if new_size_est is reasonable)
   if ( (new_size_est > mat->min_nonzero_normal_elements) && (new_size_est < mat->max_nonzero_normal_elements) ) {
		nnzmax = new_size_est;
		printf("Estimating number of non-zero normal matrix elements as %d < %d\n", nnzmax, mat->max_nonzero_normal_elements);
	}

   // Allocate temporary dense right-hand side normal form vector.
   double* dense_rhs_normal_vector = new double[mat->fm_matrix_rows]();	// the rows of the normal-form vector is number of basis functions (number of columns in input matrix)
   // Allocate space for normal_form_matrix.
   csr_matrix csr_normal_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, nnzmax); // the rows of normal matrix is number of basis functions (number of columns in input matrix)

   // Create normal form matrix
   create_sparse_normal_form_matrix(mat, nnzmax, csr_fm_matrix, csr_normal_matrix, mat->dense_fm_rhs_vector, dense_rhs_normal_vector);	

   // Accumulate for master.
   frame_weight = mat->get_frame_weight() * mat->normalization;
   if (frame_weight != 0.0) {
	   // Accumulate normal form matrix with previous/future normal form matrices
   	   // Frame weight is applied to normal matrix in this step
	   sparse_matrix_addition(mat, frame_weight, nnzmax, csr_normal_matrix, mat->sparse_matrix);

       // Accumulate normal form right-hand size vector with previous/future vectors
       // Frame weight is applied to normal vector in this step
	   cblas_daxpy(mat->fm_matrix_columns, frame_weight, dense_rhs_normal_vector, onei, mat->dense_fm_normal_rhs_vector, onei);
	}
   
   // Accumulate for each bootstrapping estimate.
   for (int i = 0; i < mat->bootstrapping_num_estimates; i++) {
	   // Get frame weight for this estimate.
	   frame_weight = mat->bootstrapping_weights[i][mat->trajectory_block_index];
	   if(frame_weight == 0.0) continue;
	   frame_weight *= mat->bootstrapping_normalization[i];
	   
	   // Accumulate normal form matrix with previous/future normal form matrices
   	   // Frame weight is applied to normal matrix in this step
	   sparse_matrix_addition(mat, frame_weight, nnzmax, csr_normal_matrix, mat->bootstrapping_sparse_fm_normal_matrices[i]);

       // Accumulate normal form right-hand size vector with previous/future vectors
       // Frame weight is applied to normal vector in this step
	   cblas_daxpy(mat->fm_matrix_columns, frame_weight, dense_rhs_normal_vector, onei, mat->bootstrapping_dense_fm_normal_rhs_vectors[i], onei);
   }
   
   // CSR formatted FM and normal temp matrices are freed by destructor at end of function
   // Free the intermediate normal form matrix and vector
   delete [] dense_rhs_normal_vector;
}

// The sparse matrix is converted to dense normal form after every block, and 
// the accumulated normal form equations are solved after the entire trajectory is read.

void convert_sparse_fm_equation_to_dense_normal_form_and_accumulate(MATRIX_DATA* const mat) 
{
	int k, l;
    // Calculate the weight of this part of the normal equations in the overall equations
    double frame_weight = mat->get_frame_weight() * mat->normalization; 

   // Convert from linked list format to CSR format
   // Note: These MKL functions use a one-based index for row_sizes and column_indices
   int n_nonzero_matrix_elements = get_n_nonzero_matrix_elements(mat);
   csr_matrix csr_fm_matrix(mat->fm_matrix_rows, mat->fm_matrix_columns, n_nonzero_matrix_elements);
   convert_linked_list_to_csr_matrix(mat, csr_fm_matrix);
   if (mat->output_raw_frame_blocks == 1) csr_fm_matrix.print_matrix_csr(mat->frame_block_fh);
   
   // Convert CSR matrix and dense RHS vector to normal-form    
   // Form sparse normal-form left-hand side matrix using mkl_dcsrmultcsr
   // rows of matrix is mat->fm_matrix_rows
   // cols of matrix is mat->fm_matrix_columns
   int nnzmax = mat->max_nonzero_normal_elements;
   
   // Form dense right-hand side normal form vector using mkl_dcsrgemv
   double* dense_rhs_normal_vector = new double[mat->fm_matrix_rows]();	// the rows of the normal-form vector is number of basis functions (number of columns in input matrix)
   //dense_rhs_normal_vector is oversized in order for matrix-vector operation to work without overwritting (it should be cols instead of rows)
   #if _mkl_flag == 1
   char trans='t';	// normal form needs transpose of first matrix times second matrix
   int request=0;	// calculate full product matrix
   int sort = 7;	// do NOT sort any matrices
   int info=0;		// error variable
   mkl_dcsrgemv(&trans, &(mat->fm_matrix_rows), csr_fm_matrix.values, 
		csr_fm_matrix.row_sizes, csr_fm_matrix.column_indices,
   		mat->dense_fm_rhs_vector, dense_rhs_normal_vector);

   // Accumulate normal form right-hand size vector with previous/future vectors
   // Frame weight is applied to normal vector in this step  
   cblas_daxpy( mat->fm_matrix_columns, frame_weight,
		dense_rhs_normal_vector, 1, mat->dense_fm_normal_rhs_vector, 1);
   #endif
   
	// Free the intermediate normal form vector
	delete [] dense_rhs_normal_vector;
	
   // Check if it makes more sense to create intermediate normal form matrix as sparse or dense
   // Either way frame weight is applied to normal matrix
   if (nnzmax * 2 > mat->fm_matrix_columns * mat->fm_matrix_columns) {
      // Process intermediate using dense matrix
	  // Allocate space for dense normal_form_matrix
      double* normal_matrix = new double[mat->fm_matrix_columns * mat->fm_matrix_columns]();
	
	  #if _mkl_flag == 1
	  // Form dense normal form matrix using 
	  mkl_dcsrmultd(&trans, &(mat->fm_matrix_rows), &(mat->fm_matrix_columns), &(mat->fm_matrix_columns),
	    csr_fm_matrix.values, csr_fm_matrix.column_indices, csr_fm_matrix.row_sizes, 
   		csr_fm_matrix.values, csr_fm_matrix.column_indices, csr_fm_matrix.row_sizes, 
   	    normal_matrix, &(mat->fm_matrix_columns) );
	  
	  // Accumulate normal form matrix with previous/future normal form matrices
	  // This operation also applies the frame weight
	  cblas_daxpy( mat->fm_matrix_columns * mat->fm_matrix_columns, frame_weight,
	  	normal_matrix, 1, mat->dense_fm_normal_matrix->values, 1);
	  #endif
	    
	  // Free the temp normal matrix
	  delete [] normal_matrix;
	  
   } else {
	  // Process intermediate using sparse matrix
      // Allocate space for sparse normal_form_matrix
      csr_matrix csr_normal_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, nnzmax); // the rows of normal matrix is number of basis functions (number of columns in input matrix)
   
      // Form sparse normal form matrix using mkl_dcsrmultcsr
      #if _mkl_flag == 1
      mkl_dcsrmultcsr(&trans, &request, &sort, &(mat->fm_matrix_rows), &(mat->fm_matrix_columns), &(mat->fm_matrix_columns), 
   		csr_fm_matrix.values, csr_fm_matrix.column_indices, csr_fm_matrix.row_sizes, 
   		csr_fm_matrix.values, csr_fm_matrix.column_indices, csr_fm_matrix.row_sizes, 
   		csr_normal_matrix.values, csr_normal_matrix.column_indices, csr_normal_matrix.row_sizes,
   		&nnzmax, &info);
      if(info != 0) {
   		printf("Error: Value returned from mkl_dcsrmultcsr is %d!\n", info);
   		exit(EXIT_FAILURE);
      }
	  #endif
	
	  // Accumulate normal form matrix with previous/future normal form matrices
	  // It would be nice to have a function that does mixed addition with sparse and dense matrices,
	  // but for now it is being done manually.
	  for( k = 0; k < mat->fm_matrix_columns; k++) { // k is actually rows of normal matrix is this context
		for( l = csr_normal_matrix.row_sizes[k] - 1; l < csr_normal_matrix.row_sizes[k+1] - 1; l++) {
			mat->dense_fm_normal_matrix->values[ k * mat->fm_matrix_columns + csr_normal_matrix.column_indices[l] ] += csr_normal_matrix.values[l] * frame_weight;
		}
	  } 
      // CSR formatted FM and normal temp matrices are freed by destructor at end of function
  	 }
} 

void convert_sparse_fm_equation_to_dense_normal_form_and_bootstrap(MATRIX_DATA* const mat) 
{
   int k, l;
   // Calculate the weight of this part of the normal equations in the overall equations
   double frame_weight; 
   int num_elements = mat->fm_matrix_columns * mat->fm_matrix_columns;
   int onei=1;
	
   // Convert from linked list format to CSR format
   // Note: These MKL functions use a one-based index for row_sizes and column_indices
   int n_nonzero_matrix_elements = get_n_nonzero_matrix_elements(mat);
   csr_matrix csr_fm_matrix(mat->fm_matrix_rows, mat->fm_matrix_columns, n_nonzero_matrix_elements);
   convert_linked_list_to_csr_matrix(mat, csr_fm_matrix);
   if (mat->output_raw_frame_blocks == 1) csr_fm_matrix.print_matrix_csr(mat->frame_block_fh);
   
   // Convert CSR matrix and dense RHS vector to normal-form    
   // Form sparse normal-form left-hand side matrix using mkl_dcsrmultcsr
   // rows of matrix is mat->fm_matrix_rows
   // cols of matrix is mat->fm_matrix_columns
   int nnzmax = mat->max_nonzero_normal_elements;
   
   // Form dense right-hand side normal form vector using mkl_dcsrgemv.
   double* dense_rhs_normal_vector = new double[mat->fm_matrix_rows]();	// the rows of the normal-form vector is number of basis functions (number of columns in input matrix)
   //dense_rhs_normal_vector is oversized in order for matrix-vector operation to work without overwritting (it should be cols instead of rows)
   #if _mkl_flag == 1
   char trans='t';	// normal form needs transpose of first matrix times second matrix
   int request=0;	// calculate full product matrix
   int sort = 7;	// do NOT sort any matrices
   int info=0;		// error variable
   mkl_dcsrgemv(&trans, &(mat->fm_matrix_rows), csr_fm_matrix.values, 
		csr_fm_matrix.row_sizes, csr_fm_matrix.column_indices,
   		mat->dense_fm_rhs_vector, dense_rhs_normal_vector);
	#endif
   
   // Accumulate for master.
   frame_weight = mat->get_frame_weight() * mat->normalization; 
   cblas_daxpy(mat->fm_matrix_columns, frame_weight, dense_rhs_normal_vector, onei, mat->dense_fm_normal_rhs_vector, onei);
   
   // Accumulate normal form right-hand size vector with previous/future vectors.
   // Frame weight is applied to normal vector in this step.  
	for (int i = 0; i < mat->bootstrapping_num_estimates; i++) {
		frame_weight = mat->bootstrapping_weights[i][mat->trajectory_block_index];
		if(frame_weight == 0.0) continue;
		frame_weight *= mat->bootstrapping_normalization[i];
		cblas_daxpy(mat->fm_matrix_columns, frame_weight, dense_rhs_normal_vector, onei, mat->bootstrapping_dense_fm_normal_rhs_vectors[i], onei);
	}
	// Free the intermediate normal form vector
	delete [] dense_rhs_normal_vector;
	
   // Check if it makes more sense to create intermediate normal form matrix as sparse or dense.
   // Either way frame weight is applied to normal matrix.
   if (nnzmax * 2 > mat->fm_matrix_columns * mat->fm_matrix_columns) {
      // Process intermediate using dense matrix
	  // Allocate space for dense normal_form_matrix
      double* normal_matrix = new double[mat->fm_matrix_columns * mat->fm_matrix_columns]();
	
	  #if _mkl_flag == 1
	  // Form dense normal form matrix using 
	  mkl_dcsrmultd(&trans, &(mat->fm_matrix_rows), &(mat->fm_matrix_columns), &(mat->fm_matrix_columns),
	    csr_fm_matrix.values, csr_fm_matrix.column_indices, csr_fm_matrix.row_sizes, 
   		csr_fm_matrix.values, csr_fm_matrix.column_indices, csr_fm_matrix.row_sizes, 
   	    normal_matrix, &(mat->fm_matrix_columns) );
	  #endif
	  
	  // Accumulate for master.
	  frame_weight = mat->get_frame_weight() * mat->normalization; 
	  cblas_daxpy(num_elements, frame_weight, normal_matrix, onei, mat->dense_fm_normal_matrix->values, onei);
	  
	  // Accumulate normal form matrix with previous/future normal form matrices.
	  // This operation also applies the frame weight.
	  for (int i = 0; i < mat->bootstrapping_num_estimates; i++) {
		frame_weight = mat->bootstrapping_weights[i][mat->trajectory_block_index];
		if(frame_weight == 0.0) continue;
		frame_weight *= mat->bootstrapping_normalization[i];

	  	cblas_daxpy(num_elements, frame_weight, normal_matrix, onei, mat->bootstrapping_dense_fm_normal_matrices[i]->values, onei);
	  }
	    
	  // Free the temp normal matrix
	  delete [] normal_matrix;
	  
   } else {
	  // Process intermediate using sparse matrix
      // Allocate space for sparse normal_form_matrix
      csr_matrix csr_normal_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, nnzmax); // the rows of normal matrix is number of basis functions (number of columns in input matrix)
   
      // Form sparse normal form matrix using mkl_dcsrmultcsr
      #if _mkl_flag == 1
      mkl_dcsrmultcsr(&trans, &request, &sort, &(mat->fm_matrix_rows), &(mat->fm_matrix_columns), &(mat->fm_matrix_columns), 
   		csr_fm_matrix.values, csr_fm_matrix.column_indices, csr_fm_matrix.row_sizes, 
   		csr_fm_matrix.values, csr_fm_matrix.column_indices, csr_fm_matrix.row_sizes, 
   		csr_normal_matrix.values, csr_normal_matrix.column_indices, csr_normal_matrix.row_sizes,
   		&nnzmax, &info);
      if(info != 0) {
   		printf("Error: Value returned from mkl_dcsrmultcsr is %d!\n", info);
   		exit(EXIT_FAILURE);
      }
	  #endif
	
	  // Accumulate for master.
	  frame_weight = mat->get_frame_weight() * mat->normalization; 
	  for( k = 0; k < mat->fm_matrix_columns; k++) { // k is actually rows of normal matrix is this context
		for( l = csr_normal_matrix.row_sizes[k] - 1; l < csr_normal_matrix.row_sizes[k+1] - 1; l++) {
			mat->dense_fm_normal_matrix->values[ k * mat->fm_matrix_columns + csr_normal_matrix.column_indices[l] ] += csr_normal_matrix.values[l] * frame_weight;
		}
	  } 
      
	  // Accumulate normal form matrix with previous/future normal form matrices
	  // It would be nice to have a function that does mixed addition with sparse and dense matrices,
	  // but for now it is being done manually.
	  for (int i = 0; i < mat->bootstrapping_num_estimates; i++) {
		frame_weight = mat->bootstrapping_weights[i][mat->trajectory_block_index];
		if(frame_weight == 0.0) continue;
		frame_weight *= mat->bootstrapping_normalization[i];

	    for( k = 0; k < mat->fm_matrix_columns; k++) { // k is actually rows of normal matrix is this context
			for( l = csr_normal_matrix.row_sizes[k] - 1; l < csr_normal_matrix.row_sizes[k+1] - 1; l++) {
				mat->bootstrapping_dense_fm_normal_matrices[i]->values[ k * mat->fm_matrix_columns + csr_normal_matrix.column_indices[l] ] += csr_normal_matrix.values[l] * frame_weight;
			}
	  	}
	  } 
      // CSR formatted FM and normal temp matrices are freed by destructor at end of function
  	 }
} 

void calculate_frame_average_and_add_to_normal_matrix(MATRIX_DATA* const mat)
{
  int i;
  for(i = 0;i < mat->fm_matrix_columns;i++)
    {
      mat->dense_fm_normal_matrix->values[i*2] += mat->dense_fm_matrix->values[i];
      mat->dense_fm_matrix->values[i] *= mat->dense_fm_matrix->values[i];
      mat->dense_fm_normal_matrix->values[(i*2) + 1] += mat->dense_fm_matrix->values[i];
    }
  mat->dense_fm_matrix->reset_matrix();
}

void do_nothing_to_fm_matrix(MATRIX_DATA* const mat) {}

// Helper routines for sparse matrix operations.

// This function determines the number of non-zero matrix elements by walking the linked-list matrix and dense virial constraint data

int get_n_nonzero_matrix_elements(MATRIX_DATA* const mat)
{
	int n_nonzero_matrix_elements = 0;	
	
    // Begin by calculating the total number of non-zero elements in this block
    for (int k = 0; k < mat->rows_less_virial_constraint_rows; k++) {
        n_nonzero_matrix_elements += mat->ll_sparse_matrix_row_heads[k].n;
    }
    n_nonzero_matrix_elements *= DIMENSION;
    if (mat->virial_constraint_rows > 0) {
        for (int k = 0; k < mat->virial_constraint_rows * mat->fm_matrix_columns; k++) {
            if (mat->dense_fm_matrix->values[k] > VERYSMALL 
                || mat->dense_fm_matrix->values[k] < -VERYSMALL) {
                    n_nonzero_matrix_elements++;
            }
        }
    }
   printf("Rectangular FM matrix has %d non-zero elements. The sparsity is %.2lf percent.\n", n_nonzero_matrix_elements, 100.0 * (1.0 - ((double) n_nonzero_matrix_elements / (double) (mat->fm_matrix_columns * mat->fm_matrix_rows))) );    
   return n_nonzero_matrix_elements;
}
 
// Helper function to convert the sparse matrix accumulated as a linked list to CSR format
void convert_linked_list_to_csr_matrix(MATRIX_DATA* const mat, csr_matrix& csr_fm_matrix)
{   
   int row_counter, row_size, num_in_row, rowD;
   struct linked_list_sparse_matrix_element* curr_elem, *prev_elem;
   double value;
   	
   for (int k = 0; k < mat->rows_less_virial_constraint_rows; k++) {
        curr_elem = mat->ll_sparse_matrix_row_heads[k].h;
        row_counter = 0;
        while (curr_elem != NULL) {
            row_size = csr_fm_matrix.row_sizes[DIMENSION * k];
            num_in_row = mat->ll_sparse_matrix_row_heads[k].n;

            for (int i = 0; i < DIMENSION; i++) {
	            // add to element values list (adjust for built-in one-base added in row_size[0] above
	            csr_fm_matrix.values[row_size + i * num_in_row + row_counter - 1] = curr_elem->valx[i];
    	
    	        // add to column indices list
        	    csr_fm_matrix.column_indices[row_size + i * num_in_row + row_counter - 1] = curr_elem->col + 1;					// convert to one-base for columns
			}
			
            //move on and delete previous element
            prev_elem = curr_elem;
            curr_elem = prev_elem->next;
            delete prev_elem;
            row_counter++;
        }
        // add to row size list
        rowD =  DIMENSION * k;
        // Note: one-base in taken into account at element 0, so no further modification is needed for rows
        
        for (int i = 0; i < DIMENSION; i++) {
	        csr_fm_matrix.row_sizes[rowD + 1 + i] = csr_fm_matrix.row_sizes[rowD + i] + row_counter;		
		}
				
		// reset row head information
		mat->ll_sparse_matrix_row_heads[k].h = NULL;
        mat->ll_sparse_matrix_row_heads[k].n = 0;
	}

    if (mat->virial_constraint_rows > 0) {
        row_counter = csr_fm_matrix.row_sizes[mat->rows_less_virial_constraint_rows * DIMENSION] - 1; // remove one-base for processing
        for (int k = 0; k < mat->virial_constraint_rows; k++) {
            for (int l = 0; l < mat->fm_matrix_columns; l++) {
                value = mat->dense_fm_matrix->values[l * mat->virial_constraint_rows + k];
                if (value > VERYSMALL || value < -VERYSMALL) {
                    csr_fm_matrix.values[row_counter] = value;
                    csr_fm_matrix.column_indices[row_counter] = l + 1; // re-apply one-base to columns
                    row_counter++;
                }
                csr_fm_matrix.row_sizes[mat->rows_less_virial_constraint_rows * DIMENSION + k + 1] = row_counter + 1; // re-apply one-base to rows
            }
        }
    }
}

// Helper function to precondition the sparse normal equations by 
// rescaling each of the columns by its root-of-sum-of-squares-of-elements value.

void precondition_sparse_matrix(int const fm_matrix_columns, double* h, csr_matrix* csr_normal_matrix)
{
   int k, l;
   // Apply to each row of matrix, but the normal matrix rows = fm_matrix_columns.
   printf("Preconditioning sparse FM normal equations.\n");
   for (k = 0; k < fm_matrix_columns; k++) {
      h[k] = 0.0;
   }
   for (k = 0; k < fm_matrix_columns; k++) {
   	  for (l = csr_normal_matrix->row_sizes[k] - 1; l < csr_normal_matrix->row_sizes[k + 1] - 1; l++) {
         h[csr_normal_matrix->column_indices[l] - 1] += csr_normal_matrix->values[l] * csr_normal_matrix->values[l];
      }
   }
   for (k = 0; k < fm_matrix_columns; k++) {
      if (h[k] > VERYSMALL) {
         h[k] = 1.0 / sqrt(h[k]);
      } else {
         h[k] = 1.0;
      }
   }
    
   for (k = 0; k < fm_matrix_columns; k++) {
      for (l = csr_normal_matrix->row_sizes[k] - 1; l < csr_normal_matrix->row_sizes[k + 1] - 1; l++) {
         csr_normal_matrix->values[l] *= h[csr_normal_matrix->column_indices[l] - 1];
      }
   }
}

// Helper function to perform regularization

void regularize_sparse_matrix(MATRIX_DATA* const mat)
{
   double squared_regularization_parameter = mat->tikhonov_regularization_param * mat->tikhonov_regularization_param;
   csr_matrix csr_regularization_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, mat->fm_matrix_columns);
   for (int k = 0; k < mat->fm_matrix_columns; k++) {
      csr_regularization_matrix.values[k] = squared_regularization_parameter;
      csr_regularization_matrix.column_indices[k] = k + 1;
      csr_regularization_matrix.row_sizes[k] = k + 1;
   }
   csr_regularization_matrix.row_sizes[mat->fm_matrix_columns] = mat->fm_matrix_columns + 1;
        
   // Add diagonal regularization matrix to current normal matrix using extra normal matrix as intermediate
   // Resulting matrix will have at most the sum of each matrix's number of non-zero elements
   int nnzmax = mat->sparse_matrix->row_sizes[mat->fm_matrix_columns] + mat->fm_matrix_columns;
   double beta = 1.0;	// no scaling of either matrix is needed
   sparse_matrix_addition(mat, beta, nnzmax, csr_regularization_matrix, mat->sparse_matrix);
   // temp regularization matrix is automatically deleted at end of function
}

void regularize_vector_sparse_matrix(MATRIX_DATA* const mat, double* regularization_vector)
{
   csr_matrix csr_regularization_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, mat->fm_matrix_columns);
   for (int k = 0; k < mat->fm_matrix_columns; k++) {
      csr_regularization_matrix.values[k] = regularization_vector[k];
      csr_regularization_matrix.column_indices[k] = k + 1;
      csr_regularization_matrix.row_sizes[k] = k + 1;
   }
   csr_regularization_matrix.row_sizes[mat->fm_matrix_columns] = mat->fm_matrix_columns + 1;
        
   // Add diagonal regularization matrix to current normal matrix using extra normal matrix as intermediate
   // Resulting matrix will have at most the sum of each matrix's number of non-zero elements
   int nnzmax = mat->sparse_matrix->row_sizes[mat->fm_matrix_columns] + mat->fm_matrix_columns;
   double beta = 1.0;	// no scaling of either matrix is needed
   sparse_matrix_addition(mat, beta, nnzmax, csr_regularization_matrix, mat->sparse_matrix);
   // temp regularization matrix is automatically deleted at end of function
}
 
void sparse_matrix_addition(MATRIX_DATA* const mat, double frame_weight, int nnzmax, csr_matrix& csr_normal_matrix, csr_matrix* main_normal_matrix)
{
   double* extra_csr_normal_matrix_values = new double[nnzmax]();
   int* extra_csr_normal_matrix_column_indices = new int[nnzmax]();
   int* extra_csr_normal_matrix_row_sizes = new int[mat->fm_matrix_columns + 1]();

   #if _mkl_flag == 1
   char trans='n';
   int request=0;
   int sort=0;
   int info=0;
   mkl_dcsradd(&trans, &request, &sort, &(mat->fm_matrix_columns), &(mat->fm_matrix_columns), 
		main_normal_matrix->values, main_normal_matrix->column_indices, main_normal_matrix->row_sizes,
		&frame_weight, csr_normal_matrix.values, csr_normal_matrix.column_indices, csr_normal_matrix.row_sizes, 
		extra_csr_normal_matrix_values, extra_csr_normal_matrix_column_indices, extra_csr_normal_matrix_row_sizes,	
		&nnzmax, &info);
	if(info != 0) {
   		printf("Error: Value returned from mkl_dcsradd is %d!\n", info);
   		exit(EXIT_FAILURE);
   	}
	#endif
	
   	// Switch accumulated normal matrix with extra (temp array)
	main_normal_matrix->set_csr_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, nnzmax, 
						 extra_csr_normal_matrix_values, extra_csr_normal_matrix_column_indices, extra_csr_normal_matrix_row_sizes);
}

// Helper function to perform regularization

void regularize_sparse_matrix(MATRIX_DATA* const mat, csr_matrix* csr_normal_matrix)
{
   double squared_regularization_parameter = mat->tikhonov_regularization_param * mat->tikhonov_regularization_param;
   csr_matrix csr_regularization_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, mat->fm_matrix_columns);
   for (int k = 0; k < mat->fm_matrix_columns; k++) {
      csr_regularization_matrix.values[k] = squared_regularization_parameter;
      csr_regularization_matrix.column_indices[k] = k + 1;
      csr_regularization_matrix.row_sizes[k] = k + 1;
   }
   csr_regularization_matrix.row_sizes[mat->fm_matrix_columns] = mat->fm_matrix_columns + 1;
        
   // Add diagonal regularization matrix to current normal matrix using extra normal matrix as intermediate
   // Resulting matrix will have at most the sum of each matrix's number of non-zero elements
   int nnzmax = csr_normal_matrix->row_sizes[mat->fm_matrix_columns] + mat->fm_matrix_columns;
   double beta = 1.0;	// no scaling of either matrix is needed
   sparse_matrix_addition(mat, beta, nnzmax, csr_regularization_matrix, csr_normal_matrix);
   // temp regularization matrix is automatically deleted at end of function
}

void regularize_vector_sparse_matrix(MATRIX_DATA* const mat, csr_matrix* csr_normal_matrix, double* regularization_vector)
{
   csr_matrix csr_regularization_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, mat->fm_matrix_columns);
   for (int k = 0; k < mat->fm_matrix_columns; k++) {
      csr_regularization_matrix.values[k] = regularization_vector[k];
      csr_regularization_matrix.column_indices[k] = k + 1;
      csr_regularization_matrix.row_sizes[k] = k + 1;
   }
   csr_regularization_matrix.row_sizes[mat->fm_matrix_columns] = mat->fm_matrix_columns + 1;
        
   // Add diagonal regularization matrix to current normal matrix using extra normal matrix as intermediate
   // Resulting matrix will have at most the sum of each matrix's number of non-zero elements
   int nnzmax = csr_normal_matrix->row_sizes[mat->fm_matrix_columns] + mat->fm_matrix_columns;
   double beta = 1.0;	// no scaling of either matrix is needed
   sparse_matrix_addition(mat, beta, nnzmax, csr_regularization_matrix, csr_normal_matrix);
   // temp regularization matrix is automatically deleted at end of function
}
 
// Wrapper function for PARDISO sparse matrix solver

void pardiso_solve(MATRIX_DATA* const mat, csr_matrix* const sparse_matrix, double* const dense_fm_normal_rhs_vector)
{
	printf("Solving sparse normal matrix using PARDISO.\n");
	fflush(stdout);
    // Solve the normal equations using PARDISO
	// Set-up workspace and variables for PARDISO
	#if _mkl_flag == 1
	void *pt[64];					// handle for PARDISO internal memory
	int *iparm = new int[64](); 	// PARDISO parameters set to default values on first call
	int nrhs = 1;		// number of right hand side vectors to solve for
	int maxfct = 1;		// maximal number of factors with identical nonzero sparsity structure that the user would like to keep in memory at the same time
	int mnum = 1;		// select matrix to factorize (between 1 and maxfct)
	int msglvl = 0;		// 0 silent, 1 verbose output
	int error = 0;
	int mtype = 11;		// Using real and non-symmetric matrix (11)
						// Could also use real and structurally symmetric (1)
						// It looks tempting to use to real, symmetric, indefinite matrix (-2), but this requires all diagonal elements to be zero
	int* perm = new int[mat->fm_matrix_columns]();
	
	pardisoinit(pt, &mtype, iparm);
    iparm[0] = 1;
	iparm[1] = 2;							// (2)Nested dissection algorithm from METIS; (3) is OpenMP version if iparm[33]=1
	iparm[3] = 0;							// 10*L + K, where 10^-L is CGS preconditioning tolerance for ||dx_i||/dx_0
											// and K = 0 does default operation and higher values replace factorization steps with 
											// K = 1 replaces factorization with CGS iterations
	iparm[4] = 0;							// Do not use user-input perm vector for full-in reducing permuation.
	iparm[5] = 0;							// Write solution to "x" fm_solution.
	iparm[7] = mat->itnlim;					// max number of iterative refinement steps (defualt = 0, performs 2 iterations when perturbed pivots are used)
	iparm[9] = 13;							// Pivoting pertubation = -log10(eps) for pivoting pertubation (default = 13)
	iparm[10] = 0;							// Scaling so that the diagonal elements are equal to 1 and the asbolute value of the off diagonal elements is <= 1.
	iparm[11] = 0;							// Solve Ax = b with no transposition
	iparm[12] = 1;							// Use (non)-symmetric weighted matching for improved accuracy
	iparm[20] = 1;							// Allow 1x1 and 2x2 Buch and Kauffman pivoting during factorization
	iparm[23] = 1;							// Allow 2-level factorization for improved OpenMP parallelization.
	iparm[24] = 0;							// use parallel algorithm for solve.
	iparm[26] = 1;							// 1 is matrix-checker for debugging, 0 is off
	iparm[27] = 0;							// Use double precision
	iparm[30] = 0;							// Disable partial solve feature.
	iparm[33] = mat->num_sparse_threads;	// set 34th entry equal to number of processors
	iparm[34] = 0;							// Use 1-based indexing
	iparm[35] = 0;							// Do not use Schur complement method.
	iparm[36] = 0;							// CSR-format input
	iparm[55] = 0;							// Automatic pivoting control
	iparm[59] = 0;							// In-core mode.
	
	// 1) Analysis (fill-reduction analysis and symbolic factorization)
	// 2) Numerical Factorization
	// 3) Solve (forward and backward solve including iterative refinement)
	int phase = 13;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &(mat->fm_matrix_columns), sparse_matrix->values,
			sparse_matrix->row_sizes, sparse_matrix->column_indices,
			perm, &nrhs, iparm, &msglvl, dense_fm_normal_rhs_vector, mat->block_fm_solution, &error);
    if(error != 0) {
    	printf ("\nError %d during PARDISO sparse matrix solving!\n", error);
    	exit(EXIT_FAILURE);
    } 
    // 4) Clean-up (termination and memory release)
    phase = -1;
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &(mat->fm_matrix_columns), sparse_matrix->values,
			sparse_matrix->row_sizes, sparse_matrix->column_indices,
			perm, &nrhs, iparm, &msglvl, dense_fm_normal_rhs_vector, mat->block_fm_solution, &error);
    if(error != 0) {
    	printf ("\nError %d during PARDISO clean-up!\n", error);
    	exit(EXIT_FAILURE);
    } 
	
    // Free temp variables
    delete [] iparm;
    delete [] perm;
    #endif
}

void solve_this_sparse_matrix(MATRIX_DATA* const mat)
{
    // Convert from linked list format to CSR format
    // Note: These MKL functions use a one-based index for row_sizes and column_indices
    int n_nonzero_matrix_elements = get_n_nonzero_matrix_elements(mat);
	csr_matrix csr_fm_matrix(mat->fm_matrix_rows, mat->fm_matrix_columns, n_nonzero_matrix_elements);
    convert_linked_list_to_csr_matrix(mat, csr_fm_matrix);
	if (mat->output_raw_frame_blocks == 1) csr_fm_matrix.print_matrix_csr(mat->frame_block_fh);
	
   // Convert CSR matrix and dense RHS vector to normal-form    
   // Form sparse normal-form left-hand side matrix using mkl_dcsrmultcsr
   // rows of matrix is mat->fm_matrix_rows
   // cols of matrix is mat->fm_matrix_columns
   int nnzmax = mat->max_nonzero_normal_elements;
   // allocate space for normal_form_matrix
   mat->sparse_matrix = new csr_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, nnzmax);	// the rows of normal matrix is number of basis functions (number of columns in input matrix)
   
   #if _mkl_flag == 1
   char trans='t';	// normal form needs transpose of first matrix times second matrix
   int request=0;	// calculate full product matrix
   int sort = 7;	// do NOT sort any matrices
   int info=0;		// error variable
   mkl_dcsrmultcsr(&trans, &request, &sort, &(mat->fm_matrix_rows), &(mat->fm_matrix_columns), &(mat->fm_matrix_columns), 
   		csr_fm_matrix.values, csr_fm_matrix.column_indices, csr_fm_matrix.row_sizes, 
   		csr_fm_matrix.values, csr_fm_matrix.column_indices, csr_fm_matrix.row_sizes, 
   		mat->sparse_matrix->values, mat->sparse_matrix->column_indices, mat->sparse_matrix->row_sizes,
   		&nnzmax, &info);
   	if(info != 0) {
   		printf("Error: Value returned from mkl_dcsrmultcsr is %d!\n", info);
   		exit(EXIT_FAILURE);
   	}
	#endif
   	
   	printf("Actual number of non-zero normal form matrix entries is %d.\n This is a density of %.2lf percent.\n", mat->sparse_matrix->row_sizes[mat->fm_matrix_columns] - 1, 100.0 * (double) (mat->sparse_matrix->row_sizes[mat->fm_matrix_columns] - 1)/ (double) nnzmax);
	
	//Need to accumulate this matrix with previous/future normal form matrices
	
   // Form dense right-hand side normal form vector using mkl_dcsrgemv
   mat->dense_fm_normal_rhs_vector = new double[mat->fm_matrix_rows]();	// the rows of the normal-form vector is number of basis functions (number of columns in input matrix)
   #if _mkl_flag == 1
   mkl_dcsrgemv(&trans, &(mat->fm_matrix_rows), csr_fm_matrix.values, 
		csr_fm_matrix.row_sizes, csr_fm_matrix.column_indices,
   		mat->dense_fm_rhs_vector, mat->dense_fm_normal_rhs_vector);
   #endif
   	
   // Apply vector regularization if requested by user.
   if (mat->regularization_style == 2) {
	    printf("Regularizing FM normal equations.\n");
    	fflush(stdout);
    	regularize_vector_sparse_matrix(mat, mat->regularization_vector);
   }
    
   // Precondition the normal equations by rescaling each of the columns by its 
   // root-of-sum-of-squares-of-elements value.
   precondition_sparse_matrix(mat->fm_matrix_columns, mat->h, mat->sparse_matrix);
	
	// Apply Tikhonov regularization if requested by user.
    if (mat->regularization_style == 1) {
	    printf("Regularizing FM normal equations.\n");
    	fflush(stdout);
    	regularize_sparse_matrix(mat);
    }
  
    // Solve the normal equations using PARDISO
	pardiso_solve(mat, mat->sparse_matrix, mat->dense_fm_normal_rhs_vector);
	   
   // CSR formatted FM temp matrix is freed by destructor at end of function
   // Free the CSR formatting normal matrix and rhs vector
   delete mat->sparse_matrix;
   mat->sparse_matrix = NULL;
}

inline void create_sparse_normal_form_matrix(MATRIX_DATA* const mat, const int nnzmax, csr_matrix& csr_fm_matrix, csr_matrix& csr_normal_matrix, double* const dense_fm_rhs_vector, double* const dense_rhs_normal_vector)
{  
   // Convert CSR matrix and dense RHS vector to normal-form    
   // Form sparse normal-form left-hand side matrix using mkl_dcsrmultcsr
   // rows of matrix is mat->fm_matrix_rows
   // cols of matrix is mat->fm_matrix_columns

   // Form sparse normal form matrix using mkl_dcsrmultcsr
   #if _mkl_flag == 1
   char trans='t';	// normal form needs transpose of first matrix times second matrix
   int request=0;	// calculate full product matrix
   int sort = 7;	// do NOT sort any matrices
   int info=0;		// error variable
   mkl_dcsrmultcsr(&trans, &request, &sort, &(mat->fm_matrix_rows), &(mat->fm_matrix_columns), &(mat->fm_matrix_columns), 
   		csr_fm_matrix.values, csr_fm_matrix.column_indices, csr_fm_matrix.row_sizes, 
   		csr_fm_matrix.values, csr_fm_matrix.column_indices, csr_fm_matrix.row_sizes, 
   		csr_normal_matrix.values, csr_normal_matrix.column_indices, csr_normal_matrix.row_sizes,
   		&nnzmax, &info);
   	if(info != 0) {
   		printf("Error: Value returned from mkl_dcsrmultcsr is %d!\n", info);
   		exit(EXIT_FAILURE);
   	}
	#endif
    printf("Actual number of non-zero normal form matrix entries is %d.\n This is a density of %.2lf percent.\n", csr_normal_matrix.row_sizes[mat->fm_matrix_columns] - 1, 100.0 * (double) (csr_normal_matrix.row_sizes[mat->fm_matrix_columns] - 1)/ (double) nnzmax);

   // Note: dense_rhs_normal_vector is oversized in order for matrix-vector operation to work without overwritting (it should be cols instead of rows)
   #if _mkl_flag == 1
   mkl_dcsrgemv(&trans, &(mat->fm_matrix_rows), csr_fm_matrix.values, 
		csr_fm_matrix.row_sizes, csr_fm_matrix.column_indices,
   		dense_fm_rhs_vector, dense_rhs_normal_vector);
   #endif  
}

inline void create_dense_normal_form(MATRIX_DATA* const mat, const double frame_weight, dense_matrix* const dense_fm_matrix, dense_matrix* normal_matrix, double* const dense_fm_rhs_vector, double* dense_fm_normal_rhs_vector)
{	
    double oned = 1.0;    
    // Take normal form of the current frame's matrix and add to the existing normal form matrix.
    #if _mkl_flag == 1
	char upper = 'u';
	char trans = 't';
	dsyrk_(&upper, &trans, &mat->fm_matrix_columns, &mat->fm_matrix_rows, &frame_weight, dense_fm_matrix->values, &mat->fm_matrix_rows, &oned, normal_matrix->values, &mat->fm_matrix_columns);
	#else
	cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, mat->fm_matrix_columns, mat->fm_matrix_rows, frame_weight, dense_fm_matrix->values, mat->fm_matrix_rows, oned, normal_matrix->values, mat->fm_matrix_columns);
	#endif
	// Take normal form of the current frame's target vector and add to the existing normal form target vector.
	cblas_dgemv(CblasColMajor, CblasTrans, mat->fm_matrix_rows, mat->fm_matrix_columns, frame_weight, dense_fm_matrix->values, mat->fm_matrix_rows, dense_fm_rhs_vector, 1, 1.0, dense_fm_normal_rhs_vector, 1);
}

// Calculate the residual for a dense matrix.
inline double calculate_dense_residual(MATRIX_DATA* const mat, dense_matrix* const dense_fm_normal_matrix, double* const dense_fm_normal_rhs_vector, std::vector<double> &fm_solution, double normalization)
{
	double residual, normal_matrix, vector_left, vector_right;
	int i;
	// Prepare the solution for linear algebra calls.
	int onei = 1;
	double oned = 1.0;
	double* intermediate = new double[mat->fm_matrix_columns]();
	double* solution = new double[mat->fm_matrix_columns];
	for (i = 0; i < mat->fm_matrix_columns; i++) {
		solution[i] = fm_solution[i];
		intermediate[i] = 0.0;
	}
	
	// Calculate solution^T * normal_matrix * solution
	cblas_dgemv(CblasColMajor, CblasNoTrans, mat->fm_matrix_columns, mat->fm_matrix_columns, oned,
			   dense_fm_normal_matrix->values, mat->fm_matrix_columns, solution, onei,
			   oned, intermediate, onei);  
				 	
	normal_matrix = cblas_ddot(mat->fm_matrix_columns, intermediate, onei, solution, onei);
		
	// Calculate solution^T * normal_vector
	vector_left = cblas_ddot(mat->fm_matrix_columns, solution, onei, dense_fm_normal_rhs_vector, onei);
	
	// Calculate normal_vector^T * solution
	vector_right = cblas_ddot(mat->fm_matrix_columns, dense_fm_normal_rhs_vector, onei, solution, onei);
	
	// Combine all of these terms and scale by normalization (frames)
	normal_matrix /= normalization;
	vector_right  /= normalization;
	vector_left   /= normalization;
	residual = normal_matrix - vector_right - vector_left;
	
	// Add on the force_sq_total and output.
	residual += mat->force_sq_total;

	printf("Unnormalized residual: %lf = %lf - %lf - %lf + %lf\n", residual, normal_matrix, vector_left, vector_right, mat->force_sq_total);

	delete [] intermediate;
	delete [] solution;
	
	return residual;
}

// Calculate the residual for a sparse matrix.
inline double calculate_sparse_residual(MATRIX_DATA* const mat, csr_matrix* csr_normal_matrix, double* const dense_fm_normal_rhs_vector, std::vector<double> &fm_solution, double normalization)
{
	double residual, normal_matrix, vector_left, vector_right;
	int i;
	// Prepare the solution for linear algebra calls.
	int onei = 1;
	double* intermediate = new double[mat->fm_matrix_columns]();
	double* solution = new double[mat->fm_matrix_columns];
	for (i = 0; i < mat->fm_matrix_columns; i++) {
		solution[i] = fm_solution[i];
		intermediate[i] = 0.0;
	}
    	
   // Calculate solution^T * normal_matrix * solution
   #if _mkl_flag == 1
   char none='n';  // not transpose
   mkl_dcsrgemv(&none, &mat->fm_matrix_columns, csr_normal_matrix->values, 
		csr_normal_matrix->row_sizes, csr_normal_matrix->column_indices,
   		solution, intermediate);
   #endif  
	
	normal_matrix = cblas_ddot(mat->fm_matrix_columns, intermediate, onei, solution, onei);
	
	// Calculate solution^T * normal_vector
	vector_left = cblas_ddot(mat->fm_matrix_columns, solution, onei, dense_fm_normal_rhs_vector, onei);
	
	// Calculate normal_vector^T * solution
	vector_right = cblas_ddot(mat->fm_matrix_columns, dense_fm_normal_rhs_vector, onei, solution, onei);
	
	// Combine all of these terms and scale by normalization (frames)
	normal_matrix /= normalization;
	vector_right  /= normalization;
	vector_left   /= normalization;
	residual = normal_matrix - vector_right - vector_left;
	
	// Add on the force_sq_total and output.
	residual += mat->force_sq_total;
	
	printf("Unnormalized residual: %lf = %lf - %lf - %lf + %lf\n", residual, normal_matrix, vector_left, vector_right, mat->force_sq_total);
	
	delete [] intermediate;
	delete [] solution;
	
	return residual;
}

// Calculate the residual for a sparse matrix.
inline double calculate_sparse_residual(MATRIX_DATA* const mat, csr_matrix* csr_normal_matrix, double* const dense_fm_normal_rhs_vector, std::vector<double> &fm_solution, double* const h)
{
	double residual, normal_matrix, vector_left, vector_right;
	int i, k, l;
	// Prepare the solution for linear algebra calls.
	int onei = 1;
	double* intermediate = new double[mat->fm_matrix_columns]();
	double* solution = new double[mat->fm_matrix_columns];
	for (i = 0; i < mat->fm_matrix_columns; i++) {
		solution[i] = fm_solution[i];
		intermediate[i] = 0.0;
	}
	
	// Remove pre-conditioner
    for (k = 0; k < mat->fm_matrix_columns; k++) {
      	for (l = csr_normal_matrix->row_sizes[k] - 1; l < csr_normal_matrix->row_sizes[k + 1] - 1; l++) {
        	csr_normal_matrix->values[l] /= h[csr_normal_matrix->column_indices[l] - 1];
    	}
    }

	printf("Solution\n");
	for (i = 0; i < mat->fm_matrix_columns; i++) {
		printf("%lf\n", solution[i]);
	}
	printf("\n");
	fflush(stdout);
	
	printf("RHS\n");
	for (i = 0; i < mat->fm_matrix_columns; i++) {
		printf("%lf\n", dense_fm_normal_rhs_vector[i]);
	}
	printf("\n");
	fflush(stdout);
	
   // Calculate solution^T * normal_matrix * solution
   #if _mkl_flag == 1
   char none='n';  // not transpose
   mkl_dcsrgemv(&none, &mat->fm_matrix_columns, csr_normal_matrix->values, 
		csr_normal_matrix->row_sizes, csr_normal_matrix->column_indices,
   		solution, intermediate);
   #endif  
	
	normal_matrix = cblas_ddot(mat->fm_matrix_columns, intermediate, onei, solution, onei);
	
	// Calculate solution^T * normal_vector
	vector_left = cblas_ddot(mat->fm_matrix_columns, solution, onei, dense_fm_normal_rhs_vector, onei);
	
	// Calculate normal_vector^T * solution
	vector_right = cblas_ddot(mat->fm_matrix_columns, dense_fm_normal_rhs_vector, onei, solution, onei);
	
	// Combine all of these terms and scale by normalization (frames)
	normal_matrix /= mat->normalization;
	vector_right  /= mat->normalization;
	vector_left   /= mat->normalization;
	residual = normal_matrix - vector_right - vector_left;
	
	// Add on the force_sq_total and output.
	residual += mat->force_sq_total;
	
	printf("Unnormalized residual: %lf = %lf - %lf - %lf + %lf\n", residual, normal_matrix, vector_left, vector_right, mat->force_sq_total);
	
	delete [] intermediate;
	delete [] solution;
	
	return residual;
}

inline void calculate_and_apply_dense_preconditioning(MATRIX_DATA* mat, dense_matrix* dense_fm_normal_matrix, double* h)
{
	int i, j;
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        h[i] = 0.0;
    }
    
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        for (j = 0; j < mat->fm_matrix_columns; j++) {
            h[j] = h[j] + dense_fm_normal_matrix->get_scalar(i, j) * dense_fm_normal_matrix->get_scalar(i, j);
        }
    }
    
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        if (h[i] < VERYSMALL) h[i] = 1.0;
        else h[i] = 1.0 / sqrt(h[i]);
    }

   for (i = 0; i < mat->fm_matrix_columns; i++) {
      for (j = 0; j < mat->fm_matrix_columns; j++) {
         dense_fm_normal_matrix->values[j * mat->fm_matrix_columns + i] *= h[j];
      }
   }
}

inline void calculate_dense_svd(MATRIX_DATA* mat, int fm_matrix_columns, dense_matrix* dense_fm_normal_matrix, double* dense_fm_normal_rhs_vector, double* singular_values)
{
	int space_factor = 2;
	int onei = 1;
	int lapack_setup_flag = -1;
	int smlsiz = 25;
	
	double tx = fm_matrix_columns / (smlsiz + 1.0);
	int nlvl = (int)(log10(tx) / log10(2)) + 1;
	int irank_in, info_in;
	int liwork = DIMENSION * fm_matrix_columns * nlvl + 11 * fm_matrix_columns;
	liwork *= space_factor;
	
	int* iwork = new int[liwork];
	double* lapack_temp_workspace = new double[1];

	// The SVD routine works in two modes: first, the routine is run with a dummy workspace to
	// determine the size of the needed workspace, then that workspace is allocated, then the
	// routine is run again with a sufficient workspace to perform SVD.
	dgelsd_(&fm_matrix_columns, &fm_matrix_columns, &onei, dense_fm_normal_matrix->values, &fm_matrix_columns, dense_fm_normal_rhs_vector, &fm_matrix_columns, singular_values, &mat->rcond, &irank_in, lapack_temp_workspace, &lapack_setup_flag, iwork, &info_in);
	lapack_setup_flag = lapack_temp_workspace[0];
	delete [] lapack_temp_workspace;
	lapack_temp_workspace = new double[lapack_setup_flag];
	dgelsd_(&fm_matrix_columns, &fm_matrix_columns, &onei, dense_fm_normal_matrix->values, &fm_matrix_columns, dense_fm_normal_rhs_vector, &fm_matrix_columns, singular_values, &mat->rcond, &irank_in, lapack_temp_workspace, &lapack_setup_flag, iwork, &info_in);

	// Clean up the heap-allocated temps.
	delete [] lapack_temp_workspace;
	delete [] iwork;
}

inline void calculate_dense_svd(MATRIX_DATA* mat, int fm_matrix_columns, int fm_matrix_rows, dense_matrix* dense_fm_normal_matrix, double* dense_fm_rhs_vector, double* singular_values)
{
	int space_factor = 2;
	int onei = 1;
	int lapack_setup_flag = -1;
	int smlsiz = 25;
	
	double tx = fm_matrix_columns / (smlsiz + 1.0);
	int nlvl = (int)(log10(tx) / log10(2)) + 1;
	int irank_in, info_in;
	int liwork = DIMENSION * fm_matrix_columns * nlvl + 11 * fm_matrix_columns;
	liwork *= space_factor;
	
	int* iwork = new int[liwork];
	double* lapack_temp_workspace = new double[1];

	// The SVD routine works in two modes: first, the routine is run with a dummy workspace to
	// determine the size of the needed workspace, then that workspace is allocated, then the
	// routine is run again with a sufficient workspace to perform SVD.
	dgelsd_(&fm_matrix_rows, &fm_matrix_columns, &onei, dense_fm_normal_matrix->values, &fm_matrix_rows, dense_fm_rhs_vector, &fm_matrix_columns, singular_values, &mat->rcond, &irank_in, lapack_temp_workspace, &lapack_setup_flag, iwork, &info_in);
	lapack_setup_flag = lapack_temp_workspace[0];
	delete [] lapack_temp_workspace;
	lapack_temp_workspace = new double[lapack_setup_flag];
	dgelsd_(&fm_matrix_rows, &fm_matrix_columns, &onei, dense_fm_normal_matrix->values, &fm_matrix_rows, dense_fm_rhs_vector, &fm_matrix_columns, singular_values, &mat->rcond, &irank_in, lapack_temp_workspace, &lapack_setup_flag, iwork, &info_in);

	// Clean up the heap-allocated temps.
	delete [] lapack_temp_workspace;
	delete [] iwork;
}  

//--------------------------------------------------------------------
// End-of-trajectory routines
//--------------------------------------------------------------------

// The sparse matrix calculation has set up a vector which is the sum
// the blockwise solution vectors, so dividing by the number of summed
// elements gives the final block-averaged FM solution.

void average_sparse_block_fm_solutions(MATRIX_DATA* const mat)
{
    // Write a binary output of the coefficient vector if desired
    if (mat->output_style >= 2) {
        FILE* mat_out;
        mat_out = open_file("result.out", "wb");
        fwrite(&mat->fm_solution[0], sizeof(double), mat->fm_matrix_columns, mat_out);
        fwrite(&mat->fm_solution_normalization_factors[0], sizeof(double), mat->fm_matrix_columns, mat_out);
        fclose(mat_out);
        // If no other output was desired, terminate the program successfully.
        if (mat->output_style == 3) exit(EXIT_SUCCESS);
    }
    
    for (int i = 0; i < mat->fm_matrix_columns; i++) {
        if (mat->fm_solution_normalization_factors[i] > VERYSMALL) {
            mat->fm_solution[i] /=  mat->normalization  * mat->fm_solution_normalization_factors[i];
        } else mat->fm_solution[i] = 0.0;
    }
    delete [] mat->fm_solution_normalization_factors;
    delete [] mat->h;
}

void average_sparse_bootstrapping_solutions(MATRIX_DATA* const mat)
{
    // Write a binary output of the coefficient vector if desired
    if (mat->output_style >= 2) {
        FILE* mat_out;
        mat_out = open_file("result.out", "wb");
        for (int i = 0; i < mat->bootstrapping_num_estimates; i++) {
	        fwrite(&mat->bootstrap_solutions[i][0], sizeof(double), mat->fm_matrix_columns, mat_out);
    	    fwrite(&mat->fm_solution_normalization_factors[0], sizeof(double), mat->fm_matrix_columns, mat_out);
        }
        fclose(mat_out);
        // If no other output was desired, terminate the program successfully.
        if (mat->output_style == 3) exit(EXIT_SUCCESS);
    }
    
    // Apply for master
        for (int i = 0; i < mat->fm_matrix_columns; i++) {
        if (mat->fm_solution_normalization_factors[i] > VERYSMALL) {
            mat->fm_solution[i] /=  mat->normalization  * mat->fm_solution_normalization_factors[i];
        } else mat->fm_solution[i] = 0.0;
    }

    // Apply for bootstrapping estimates.
    for (int i = 0; i < mat->fm_matrix_columns; i++) {
    	for (int j = 0; j < mat->bootstrapping_num_estimates; j++) {
    		// This checks how many frames had reasonable solutions. As long as at least one frame had reasonable solutions, then that average is maintained.
	        // The overall fm_solution has already been appropriately normalized here (excluding only the contributions from extreme values).
	        if (mat->fm_solution_normalization_factors[i] < VERYSMALL) {
	        	mat->bootstrap_solutions[j][i] = 0.0;
    		}
    	}
    }
    delete [] mat->fm_solution_normalization_factors;
    delete [] mat->h;
}

// The sparse matrix equations are now in normal form and should be solved.

void solve_sparse_fm_normal_equations(MATRIX_DATA* const mat)
{
   // Back up RHS.
   double* backup_rhs = new double[mat->fm_matrix_columns];
   int matrix_size = mat->sparse_matrix->row_sizes[mat->fm_matrix_columns];
   csr_matrix* backup_normal_matrix = new csr_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns, matrix_size);
   for (int i = 0; i < mat->fm_matrix_columns; i++) {
     backup_rhs[i] = mat->dense_fm_normal_rhs_vector[i];
   }
   
   
   // Write sparse matrix as CSR and dense for parallel runs
	if (mat->output_style >= 2) {
	
		FILE* csr_out = open_file("result_csr.out", "w");
		mat->sparse_matrix->write_file(csr_out);
		for (int i = 0; i < mat->fm_matrix_columns; i++) {
			fprintf(csr_out, "%lf ", mat->dense_fm_normal_rhs_vector[i]);
		}
		fprintf(csr_out, "\n");
        fprintf(csr_out, "%lf\n", mat->force_sq_total);
        fprintf(csr_out, "%lf\n", 1.0/mat->normalization);
		fclose(csr_out);
	
		FILE* mat_out = open_file("result.out", "wb");
		int counter = 0;
		int low, high;
		double zero = 0.0;
		for (int i = 0; i < mat->fm_matrix_columns; i++) {
			low = mat->sparse_matrix->row_sizes[i] - 1;
			high = mat->sparse_matrix->row_sizes[i + 1] - 1;
			counter = low;
			for (int j = 0; j <= i; j++){
			  if( (j + 1) == mat->sparse_matrix->column_indices[counter]) {
			 	  fwrite(&mat->sparse_matrix->values[counter], sizeof(double), 1, mat_out);
				  counter++;
			   } else {
				  fwrite(&zero, sizeof(double), 1, mat_out);
			   }
			}
		}
		double inv_norm = 1.0/mat->normalization;
		fwrite(&mat->dense_fm_normal_rhs_vector[0], sizeof(double), mat->fm_matrix_columns, mat_out);
		fwrite(&mat->force_sq_total, sizeof(double), 1, mat_out);
        fwrite(&inv_norm, sizeof(double), 1, mat_out);
        fclose(mat_out);
		
		// If no other output was desired, terminate the program successfully.
		if (mat->output_style == 3) exit(EXIT_SUCCESS);
	}
   
   // Backup the sparse normal matrix
   for (int i = 0; i < mat->fm_matrix_columns; i++) {
	  backup_normal_matrix->row_sizes[i + 1] = mat->sparse_matrix->row_sizes[i + 1];
   } 
   for (int i = 0; i < matrix_size; i++) {
      backup_normal_matrix->column_indices[i] = mat->sparse_matrix->column_indices[i];
      backup_normal_matrix->values[i]         = mat->sparse_matrix->values[i];
   }
      
   // Apply vector regularization if requested by user.
   if (mat->regularization_style == 2) {
	     printf("Regularizing FM normal equations.\n");
    	 fflush(stdout);
    	 regularize_vector_sparse_matrix(mat, mat->regularization_vector);
   }
   
   // Precondition the normal equations by rescaling each of the columns by its 
   // root-of-sum-of-squares-of-elements value
   precondition_sparse_matrix(mat->fm_matrix_columns, mat->h, mat->sparse_matrix);

    // Apply Tikhonov regularization if requested by user.
    if (mat->regularization_style == 1) {
	    printf("Regularizing FM normal equations.\n");
    	fflush(stdout);
    	regularize_sparse_matrix(mat);
    }
  
    // Solve the normal equations using PARDISO
	printf("Computing solution of FM normal equations using sparse matrix operations.\n");
	mat->block_fm_solution = &(mat->fm_solution[0]);
	pardiso_solve(mat, mat->sparse_matrix, mat->dense_fm_normal_rhs_vector);
    printf("Finished PARDISO solve.\n");
	
   // Remove preconditioning effect from solution
   for (int k = 0; k < mat->fm_matrix_columns; k++) {
      mat->fm_solution[k] *= mat->h[k];
   }
   
   if (mat->output_residual == 1) {
      double residual = calculate_sparse_residual(mat, mat->sparse_matrix, mat->dense_fm_normal_rhs_vector, mat->fm_solution, mat->normalization);
      printf("residual %lf\n", residual);
   }

    // Calculate First Bayesian Estimates
    if (mat->bayesian_flag == 1) {
    	int iteration = 0;
	    int onei = 1;
    	dense_matrix* it_normal_matrix = new dense_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns);
    	dense_matrix* backup_dense_matrix = new dense_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns);
		dense_matrix* product_reg_normal_matrix_inv_normal_matrix = new dense_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns);
		double trace_product;
    	double* alpha_vec = new double[mat->fm_matrix_columns];
	    double* solution  = new double[mat->fm_matrix_columns];
		
    	for (int i = 0; i < mat->fm_matrix_columns; i++) {
			solution[i] = mat->fm_solution[i];
	    }
	    
    	double residual = calculate_sparse_residual(mat, backup_normal_matrix, backup_rhs, mat->fm_solution, mat->normalization);
		printf("iteration %d: residual %lf\n", iteration, residual);
		
	    double n_cg_sites = (double)( mat->rows_less_virial_constraint_rows/ mat->frames_per_traj_block / mat->size_per_vector);
		double n_frames = 1.0 / mat->normalization;
		
	    double alpha = (double)(mat->fm_matrix_columns) / cblas_ddot(mat->fm_matrix_columns, solution, onei, solution, onei);
		double beta  = (double)(mat->size_per_vector) * n_cg_sites * n_frames / residual;
		
		for (int i = 0; i < mat->fm_matrix_columns; i++) {
			alpha_vec[i] = alpha;
		}
				
		FILE* alpha_fp = fopen("alpha.out", "w");
		FILE* beta_fp  = fopen("beta.out",  "w");
		FILE* sol_fp   = fopen("solution.out", "w");
		FILE* res_fp   = fopen("residual.out", "w");
		FILE* ext_fp   = fopen("ext_residual.out", "w");
		write_iteration(alpha_vec, beta, mat->fm_solution, residual, iteration, alpha_fp, beta_fp, sol_fp, res_fp);
		FILE* mat_fp   = fopen("matrix.out", "w");
		FILE* inv_fp   = fopen("inverse.out", "w");
		//backup_normal_matrix->print_matrix(mat_fp);
	
		for (int i = 0; i < mat->fm_matrix_columns; i++) {
		   for (int j = backup_normal_matrix->row_sizes[i] - 1; j < backup_normal_matrix->row_sizes[i + 1] - 1; j++) {
		      backup_dense_matrix->assign_scalar(i, backup_normal_matrix->column_indices[j], backup_normal_matrix->values[j]);
		   }
		   backup_dense_matrix->add_scalar(i, i, alpha_vec[i] * mat->normalization / beta);
		}
		
		double* reg_vec = new double[mat->fm_matrix_columns];
			
		while (iteration < mat->bayesian_max_iter) {
		
			// Initialize RHS vector and normal matrix from backups.
			for(int i = 0; i < mat->fm_matrix_columns; i++) {
				mat->dense_fm_normal_rhs_vector[i] = backup_rhs[i];
			}
   			for (int i = 0; i < mat->fm_matrix_columns; i++) {
	  		   mat->sparse_matrix->row_sizes[i + 1] = backup_normal_matrix->row_sizes[i + 1];
   			} 
   			for (int i = 0; i < matrix_size; i++) {
      		   mat->sparse_matrix->column_indices[i] = backup_normal_matrix->column_indices[i];
      		   mat->sparse_matrix->values[i]         = backup_normal_matrix->values[i];
   			}
   			
			// Solve matrix with regularization
			// // First, apply regularization alpha/beta
			for(int i = 0; i < mat->fm_matrix_columns; i++) {
			   reg_vec[i] = alpha_vec[i] * mat->normalization / beta;
			}
			regularize_vector_sparse_matrix(mat, reg_vec);
			
			// // Then, apply preconditioning
			precondition_sparse_matrix(mat->fm_matrix_columns, mat->h, mat->sparse_matrix);

    		// Solve the normal equations using PARDISO
			pardiso_solve(mat, mat->sparse_matrix, mat->dense_fm_normal_rhs_vector);
    
   			// Remove preconditioning effect from solution
   			for (int k = 0; k < mat->fm_matrix_columns; k++) {
      		   mat->fm_solution[k] *= mat->h[k];
   			   solution[k] = mat->fm_solution[k];
   			}
			
			residual = calculate_sparse_residual(mat, backup_normal_matrix, backup_rhs, mat->fm_solution, mat->normalization);
			double* alpha_solution = new double[mat->fm_matrix_columns];
			for (int k = 0; k < mat->fm_matrix_columns; k++) {
				alpha_solution[k] = alpha_vec[k] * solution[k];
			}
			double alpha_product = cblas_ddot(mat->fm_matrix_columns, solution, onei, alpha_solution, onei);
			double extended_residual = beta * 0.5 * residual + 0.5 * alpha_product;
			fprintf(ext_fp, "Iteration %d: %lf\n", iteration, -extended_residual);
			printf("negative of extended residual %lf = (%lf / 2) * %lf + 1/2 * %lf\n", extended_residual, beta, residual, alpha_product);
			delete [] alpha_solution;
		
			// Calculate the values for the next round.
			iteration++;
			printf("iteration %d: residual %lf\n", iteration, residual);
			
			// Alpha needs matrix inverse and beta needs the product of the inverse with an unregularized normal matrix
			it_normal_matrix->reset_matrix();
			product_reg_normal_matrix_inv_normal_matrix->reset_matrix();

			// Actually calculate the inverse of the regularized normal matrix
			// // First, restore the dense normal matrix with regularization
			for (int i = 0; i < mat->fm_matrix_columns; i++) {
				for (int j = backup_normal_matrix->row_sizes[i] - 1; j < backup_normal_matrix->row_sizes[i + 1] - 1; j++) {
					it_normal_matrix->assign_scalar(i, backup_normal_matrix->column_indices[j], backup_normal_matrix->values[j]);
				}
				it_normal_matrix->add_scalar(i, i, alpha_vec[i] * mat->normalization / beta);
			}
			
			// // Then, invert the matrix.
			// // This routine overwrites the input matrix (it_normal_matrix).
			int* ipiv = new int[mat->fm_matrix_columns]();
			double* work = new double[mat->fm_matrix_columns * mat->fm_matrix_columns];
			int lwork = mat->fm_matrix_columns * mat->fm_matrix_columns;
			int info = 0;			
			dgetrf_(&mat->fm_matrix_columns, &mat->fm_matrix_columns, it_normal_matrix->values, &mat->fm_matrix_columns, ipiv, &info);
			dgetri_(&mat->fm_matrix_columns, it_normal_matrix->values, &mat->fm_matrix_columns, ipiv, work, &lwork, &info);
			delete [] ipiv;
			delete [] work;
			it_normal_matrix->print_matrix(inv_fp);
			
			// Calculate product using inverse of regularized matrix
			// Note: it_dense_normal_matrix is now actually the inverse of that matrix.
			double oned = 1.0;
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, mat->fm_matrix_columns, mat->fm_matrix_columns, mat->fm_matrix_columns, oned,
					it_normal_matrix->values, mat->fm_matrix_columns, backup_dense_matrix->values, mat->fm_matrix_columns, oned,
					product_reg_normal_matrix_inv_normal_matrix->values, mat->fm_matrix_columns);
		
			trace_product = 0.0;
			for (int i = 0; i < mat->fm_matrix_columns; i++) {
				trace_product += product_reg_normal_matrix_inv_normal_matrix->get_scalar(i,i);
			}
			
			// Alpha Vector
			for (int i = 0; i < mat->fm_matrix_columns; i++) {
				// Note: it_dense_normal_matrix is now actually the inverse of that matrix.
				alpha_vec[i] = 1.0 / (solution[i] * solution[i] + it_normal_matrix->get_scalar(i,i) * mat->normalization / beta);
			}
			
			// Beta Scalar
			beta =  (mat->size_per_vector * n_cg_sites * n_frames - trace_product) / residual;	

		   write_iteration(alpha_vec, beta, mat->fm_solution, residual, iteration, alpha_fp, beta_fp, sol_fp, res_fp);
		}
		// Clean-up bayesian allocated memory
		fclose(alpha_fp);
		fclose(beta_fp);
		fclose(sol_fp);
		fclose(res_fp);
		fclose(ext_fp);
		fclose(mat_fp);
		fclose(inv_fp);
		delete [] solution;
		delete [] reg_vec;
		printf("delete alpha_vec\n"); fflush(stdout);
		//delete [] alpha_vec;
		printf("delete it_normal_matrix\n"); fflush(stdout);
		//delete it_normal_matrix;
		printf("delete backup_dense_matrix\n"); fflush(stdout);
		//delete backup_dense_matrix;
		delete product_reg_normal_matrix_inv_normal_matrix;
		printf("exit bayesian\n"); fflush(stdout);
    }
    printf("after exit bayesian\n"); fflush(stdout);
    
   // Free the CSR formatted normal matrix
   delete [] backup_rhs;
   delete [] mat->dense_fm_normal_rhs_vector;
   delete backup_normal_matrix;
   delete mat->sparse_matrix;
   mat->sparse_matrix = NULL;
   delete [] mat->h;
}

void solve_sparse_fm_bootstrapping_equations(MATRIX_DATA* const mat)
{
   // Solve for master
   solve_sparse_fm_normal_equations(mat);
   mat->h = new double[mat->fm_matrix_columns]();
   
   for (int i = 0; i < mat->bootstrapping_num_estimates; i++) {

      // Apply vector regularization if requested by user.
      if (mat->regularization_style == 2) {
	     printf("Regularizing FM normal equations (estimate %d).\n", i);
    	 fflush(stdout);
    	 regularize_vector_sparse_matrix(mat, mat->bootstrapping_sparse_fm_normal_matrices[i], mat->regularization_vector);
       }
      
      // Precondition the normal equations by rescaling each of the columns by its 
      // root-of-sum-of-squares-of-elements value.
      precondition_sparse_matrix(mat->fm_matrix_columns, mat->h, mat->bootstrapping_sparse_fm_normal_matrices[i]);

      // Apply Tikhonov regularization if requested by user.
      if (mat->regularization_style == 1) {
	     printf("Regularizing FM normal equations (estimate %d).\n", i);
    	 fflush(stdout);
    	 regularize_sparse_matrix(mat, mat->bootstrapping_sparse_fm_normal_matrices[i]);
       }
  
      // Solve the normal equations using PARDISO
	   printf("Computing solution of FM normal equations using sparse matrix operations.\n");

	   mat->block_fm_solution = &(mat->bootstrap_solutions[i][0]);

	   pardiso_solve(mat, mat->bootstrapping_sparse_fm_normal_matrices[i], mat->bootstrapping_dense_fm_normal_rhs_vectors[i]);
       printf("Finished PARDISO solve.\n");
	
       // Free the CSR formatted normal matrix
       delete mat->bootstrapping_sparse_fm_normal_matrices[i];
       delete [] mat->bootstrapping_dense_fm_normal_rhs_vectors[i];
       
      // Remove preconditioning effect from solution
      for (int k = 0; k < mat->fm_matrix_columns; k++) {
         mat->bootstrap_solutions[i][k] *= mat->h[k];
      }
      
      if (mat->output_residual == 1) {
         double residual = calculate_sparse_residual(mat, mat->bootstrapping_sparse_fm_normal_matrices[i], mat->bootstrapping_dense_fm_normal_rhs_vectors[i], mat->bootstrap_solutions[i], mat->normalization);
         printf("estimate %d: residual %lf\n", i, residual);
      }
   }
   delete [] mat->h;
   //delete [] mat->dense_fm_normal_rhs_vector;
   delete [] mat->bootstrapping_sparse_fm_normal_matrices;
   delete [] mat->bootstrapping_dense_fm_normal_rhs_vectors;
}

// The dense matrix equations are now in normal form and should be solved.
// This function also handles Lanyuan's iterative method; in that case no matrix
// has been constructed and it must be read from file with the former target
// vector as well, and the target vector becomes the difference of the current and
// former target vectors.

void solve_dense_fm_normal_equations(MATRIX_DATA* const mat)
{
	int i,j,k;

    // At the time of solution, one can be sure that the raw dense matrix
    // is no longer needed.
    printf("Freeing raw FM equations.\n"); fflush(stdout);
    delete mat->dense_fm_matrix;
    
    // Store a temporary backup of the normal form target vector
    // since it is changed by the solver.
    double* backup_rhs = new double[mat->fm_matrix_columns];
    for (i = 0; i < mat->fm_matrix_columns; i++) {
    	backup_rhs[i] = mat->dense_fm_normal_rhs_vector[i];
    }
    
    if (mat->iterative_calculation_flag == 1) {
        // Read in a stored normal form matrix and normal form target vector for
        // iterative calculations
		double ttx;
        FILE* mat_in = open_file("result.in", "rb");
        double* in_rhs = new double[mat->fm_matrix_columns];
        
        for (j = 0; j < mat->fm_matrix_columns; j++) {
            for (k = 0; k <= j; k++) {
                fread(&ttx, sizeof(double), 1, mat_in);
                mat->dense_fm_normal_matrix->assign_scalar(k, j, ttx);
            }
        }
        for (j = 0; j < mat->fm_matrix_columns; j++) {
            fread(&ttx, sizeof(double), 1, mat_in);
            in_rhs[j] = ttx;
        }
        fclose(mat_in);
        
        // The target for an iterative calculation is the difference between the targets
        // for this trajectory and the previous trajectory.
        for (i = 0; i < mat->fm_matrix_columns; i++) {
            mat->dense_fm_normal_rhs_vector[i] = in_rhs[i] - mat->dense_fm_normal_rhs_vector[i];
        }
        delete [] in_rhs;
    } else {
        // Save the results in binary form for parallel runs.
        if (mat->output_style >= 2) {
            FILE* mat_out = open_file("result.out", "wb");
            for (i = 0; i < mat->fm_matrix_columns; i++) {
                fwrite(&mat->dense_fm_normal_matrix->values[i * mat->fm_matrix_columns], sizeof(double), i + 1, mat_out);
            }
       		double inv_norm = 1.0/mat->normalization;
            fwrite(&mat->dense_fm_normal_rhs_vector[0], sizeof(double), mat->fm_matrix_columns, mat_out);
        	fwrite(&mat->force_sq_total, sizeof(double), 1, mat_out);
        	fwrite(&inv_norm, sizeof(double), 1, mat_out);
            fclose(mat_out);
            // If no other output was desired, terminate the program successfully.
            if (mat->output_style == 3) exit(EXIT_SUCCESS);
        }
    }

	// Assign the upper diagonal to the lower lower diagonal (symmetric matrix)
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        for (j = 0; j < i; j++) {
            mat->dense_fm_normal_matrix->assign_scalar(i, j, mat->dense_fm_normal_matrix->values[i * mat->fm_matrix_columns + j]);
        }
    }

	// Store a temporary backup of the normal matrix since it is changed by the solver.
	dense_matrix* backup_normal_matrix = new dense_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns);
	for (i = 0; i < mat->fm_matrix_columns; i++) {
		for (int z = 0; z < mat->fm_matrix_columns; z++) {
			backup_normal_matrix->assign_scalar(z, i, mat->dense_fm_normal_matrix->get_scalar(z, i));
		}
	}
    
    // Apply vector regularization if requested.
    if (mat->regularization_style == 2) {
    	printf("Regularizing FM normal equations.\n"); fflush(stdout);
    	for (i = 0; i < mat->fm_matrix_columns; i++) {
            mat->dense_fm_normal_matrix->add_scalar(i, i, mat->regularization_vector[i]);
    	}
    }
    
    // Precondition the normal matrix using the root-of-sum-of-squares 
    // of the columns as column scaling factors.
    printf("Preconditioning FM normal equations.\n"); fflush(stdout);
    double* h = new double[mat->fm_matrix_columns];
    calculate_and_apply_dense_preconditioning(mat, mat->dense_fm_normal_matrix, h);

    // Apply Tikhonov regularization.
    if (mat->regularization_style == 1) {
    	printf("Regularizing FM normal equations.\n"); fflush(stdout);
        double squared_regularization_parameter;
        squared_regularization_parameter = mat->tikhonov_regularization_param * mat->tikhonov_regularization_param;
        for (i = 0; i < mat->fm_matrix_columns; i++) {
            mat->dense_fm_normal_matrix->add_scalar(i, i, squared_regularization_parameter);
        }
    }
    
    // Solve the normal equation by singular value decomposition using LAPACK routines.
    printf("Computing singular value decomposition of preconditioned, regularized FM normal equations.\n"); fflush(stdout);
    double* singular_values = new double[mat->fm_matrix_columns];
    calculate_dense_svd(mat, mat->fm_matrix_columns, mat->dense_fm_normal_matrix, mat->dense_fm_normal_rhs_vector, singular_values);
    
    // Print singular values.
    printf("Printing FM singular values.\n"); fflush(stdout);
    FILE* solution_file = open_file("sol_info.out", "a");
    fprintf(solution_file, "Singular vector:\n");
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        fprintf(solution_file, "%le\n", singular_values[i]);
    }
    fclose(solution_file);
    
    // Calculate the final results from the singular values.
    printf("Calculating final FM results.\n"); fflush(stdout);
    mat->fm_solution = std::vector<double>(mat->fm_matrix_columns);
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        mat->fm_solution[i] = mat->dense_fm_normal_rhs_vector[i] * h[i];
    }
   
    // Calculate and output the residual if requested.
    if (mat->output_residual == 1) {
    	double residual = calculate_dense_residual(mat, backup_normal_matrix, backup_rhs, mat->fm_solution, mat->normalization);
	    printf ("residual %lf\n", residual);
    }
    
    // Calculate First Bayesian Estimates
    if (mat->bayesian_flag == 1) {
    	int iteration = 0;
	    int onei = 1;
    	dense_matrix* it_dense_normal_matrix = new dense_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns);
		dense_matrix* product_reg_normal_matrix_inv_normal_matrix = new dense_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns);
		double trace_product;
    	double* alpha_vec = new double[mat->fm_matrix_columns];
	    double* solution  = new double[mat->fm_matrix_columns];
		
    	for (i = 0; i < mat->fm_matrix_columns; i++) {
			solution[i] = mat->fm_solution[i];
	    }
	    
    	double residual = calculate_dense_residual(mat, backup_normal_matrix, backup_rhs, mat->fm_solution, mat->normalization);
	    
	    double n_cg_sites = (double)( mat->rows_less_virial_constraint_rows/ mat->frames_per_traj_block / mat->size_per_vector);
		double n_frames = 1.0 / mat->normalization;
		
	    double alpha = (double)(mat->fm_matrix_columns) / cblas_ddot(mat->fm_matrix_columns, solution, onei, solution, onei);
		double beta  = (double)(mat->size_per_vector) * n_cg_sites * n_frames / residual;
		
		for (i = 0; i < mat->fm_matrix_columns; i++) {
			alpha_vec[i] = alpha;
		}
				
		FILE* alpha_fp = fopen("alpha.out", "w");
		FILE* beta_fp  = fopen("beta.out",  "w");
		FILE* sol_fp   = fopen("solution.out", "w");
		FILE* res_fp   = fopen("residual.out", "w");
		FILE* ext_fp   = fopen("ext_residual.out", "w");
		write_iteration(alpha_vec, beta, mat->fm_solution, residual, iteration, alpha_fp, beta_fp, sol_fp, res_fp);
		FILE* mat_fp   = fopen("matrix.out", "w");
		FILE* inv_fp   = fopen("inverse.out", "w");
		backup_normal_matrix->print_matrix(mat_fp);
		
		while (iteration < mat->bayesian_max_iter) {
		
			// Initialize RHS vector and normal matrix from backups.
			for(i = 0; i < mat->fm_matrix_columns; i++) {
				mat->dense_fm_normal_rhs_vector[i] = backup_rhs[i];
				for (j = 0; j < mat->fm_matrix_columns; j++) {
					it_dense_normal_matrix->assign_scalar(i, j, backup_normal_matrix->get_scalar(i, j));
				}
			}
			
			// Solve matrix with regularization
			// // First, apply regularization alpha/beta
			for(i = 0; i < mat->fm_matrix_columns; i++) {
				it_dense_normal_matrix->add_scalar(i, i, alpha_vec[i] * mat->normalization/ beta);
			}
			
			// // Then, apply preconditioning
			calculate_and_apply_dense_preconditioning(mat, it_dense_normal_matrix, h);

			// // Now, compute the SVD.
			for (i = 0; i < mat->fm_matrix_columns; i++) {
				singular_values[i] = 0.0;
			}
			calculate_dense_svd(mat, mat->fm_matrix_columns, it_dense_normal_matrix, mat->dense_fm_normal_rhs_vector, singular_values);
			
			for (i = 0; i < mat->fm_matrix_columns; i++) {
        		solution[i] = mat->dense_fm_normal_rhs_vector[i] * h[i];
        		mat->fm_solution[i] = solution[i];
    		}
			
			residual = calculate_dense_residual(mat, backup_normal_matrix, backup_rhs, mat->fm_solution, 1.0);
			double* alpha_solution = new double[mat->fm_matrix_columns];
			for (k = 0; k < mat->fm_matrix_columns; k++) {
				alpha_solution[k] = alpha_vec[k] * solution[k];
			}
			double alpha_product = cblas_ddot(mat->fm_matrix_columns, solution, onei, alpha_solution, onei);
			double extended_residual = beta * 0.5 * residual + 0.5 * alpha_product;
			fprintf(ext_fp, "Iteration %d: %lf\n", iteration, -extended_residual);
			printf("negative of extended residual %lf = (%lf / 2) * %lf + 1/2 * %lf\n", extended_residual, beta, residual, alpha_product);
			delete [] alpha_solution;
		
			// Calculate the values for the next round.
			iteration++;
			residual = calculate_dense_residual(mat, backup_normal_matrix, backup_rhs, mat->fm_solution, mat->normalization);
		
			// Alpha needs matrix inverse and beta needs the product of the inverse with an unregularized normal matrix
			it_dense_normal_matrix->reset_matrix();
			product_reg_normal_matrix_inv_normal_matrix->reset_matrix();
		
			// Actually calculate the inverse of the regularized normal matrix
			// // First, restore the normal matrix with regularization
			for (i = 0; i < mat->fm_matrix_columns; i++) {
				for (j = 0; j < mat->fm_matrix_columns; j++) {
					it_dense_normal_matrix->assign_scalar(i, j, backup_normal_matrix->get_scalar(i, j));
				}
				it_dense_normal_matrix->add_scalar(i, i, alpha_vec[i] * mat->normalization / beta);
			}
	
			// // Then, invert the matrix.
			// // This routine overwrites the input matrix (it_dense_normal_matrix).
			int* ipiv = new int[mat->fm_matrix_columns]();
			double* work = new double[mat->fm_matrix_columns * mat->fm_matrix_columns];
			int lwork = mat->fm_matrix_columns * mat->fm_matrix_columns;
			int info = 0;			
			dgetrf_(&mat->fm_matrix_columns, &mat->fm_matrix_columns, it_dense_normal_matrix->values, &mat->fm_matrix_columns, ipiv, &info);
			dgetri_(&mat->fm_matrix_columns, it_dense_normal_matrix->values, &mat->fm_matrix_columns, ipiv, work, &lwork, &info);
			delete [] ipiv;
			delete [] work;
			it_dense_normal_matrix->print_matrix(inv_fp);
			
			// Calculate product using inverse of regularized matrix
			// Note: it_dense_normal_matrix is now actually the inverse of that matrix.
			double oned = 1.0;
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, mat->fm_matrix_columns, mat->fm_matrix_columns, mat->fm_matrix_columns, oned,
					it_dense_normal_matrix->values, mat->fm_matrix_columns, backup_normal_matrix->values, mat->fm_matrix_columns, oned,
					product_reg_normal_matrix_inv_normal_matrix->values, mat->fm_matrix_columns);
		
			trace_product = 0.0;
			for (i = 0; i < mat->fm_matrix_columns; i++) {
				trace_product += product_reg_normal_matrix_inv_normal_matrix->get_scalar(i,i);
			}
			
			// Alpha Vector
			for (i = 0; i < mat->fm_matrix_columns; i++) {
				// Note: it_dense_normal_matrix is now actually the inverse of that matrix.
				alpha_vec[i] = 1.0 / (solution[i] * solution[i] + it_dense_normal_matrix->get_scalar(i,i) * mat->normalization / beta);
			}
			
			// Beta Scalar
			beta =  ((double)(mat->size_per_vector) * n_cg_sites * n_frames - trace_product) / residual;				

			write_iteration(alpha_vec, beta, mat->fm_solution, residual, iteration, alpha_fp, beta_fp, sol_fp, res_fp);

		}
		
		// Clean-up bayesian allocated memory
		fclose(alpha_fp);
		fclose(beta_fp);
		fclose(sol_fp);
		fclose(res_fp);
		fclose(ext_fp);
		fclose(mat_fp);
		fclose(inv_fp);
		delete [] alpha_vec;
		delete [] solution;
		delete it_dense_normal_matrix;
		delete product_reg_normal_matrix_inv_normal_matrix;
    }
    
    // For iterative calculations, the solution is a difference, so the computed quantity
    // should be added on to the previous solution value to obtain the final solution.
    if (mat->iterative_calculation_flag == 1) {
        printf("Adding iterative increment to previous solution.\n");
        fflush(stdout);
        FILE* x_in = open_file("x.in", "r");
        double* x0 = new double[mat->fm_matrix_columns];
        for (i = 0; i < mat->fm_matrix_columns; i++) fscanf(x_in, "%le", x0 + i);
        fclose(x_in);
        //test
        for (i = 0; i < mat->fm_matrix_columns; i++) mat->fm_solution[i] = mat->fm_solution[i] * mat->iterative_update_rate_coeff + x0[i];
        delete [] x0;
    }

    printf("Completed FM.\n"); fflush(stdout);
    // Restore RHS normal vector from backup.
    if (mat->output_normal_equations_rhs_flag == 1) {
        for (i = 0; i < mat->fm_matrix_columns; i++) {
        	mat->dense_fm_normal_rhs_vector[i] = backup_rhs[i];
        }
    }
    
    // Clean up allocated memory.
    delete [] h;
    delete [] singular_values;
 	delete backup_normal_matrix;
 	delete mat->dense_fm_normal_matrix;
    delete [] backup_rhs;
 	if(mat->matrix_type == 3) delete [] mat->dense_fm_normal_rhs_vector;
}
  
void solve_this_BI_equation(MATRIX_DATA* const mat, int &solution_counter)
{
  int i, j;
  FILE* BI_matrix = fopen("BI_matrix.dat","a");
  FILE* BI_vector = fopen("BI_vector.dat","a");
  for(i = 0; i < mat->fm_matrix_rows; i++){
    for(j = 0; j < mat->fm_matrix_columns; j++){
      fprintf(BI_matrix,"%lf ",mat->dense_fm_matrix->get_scalar(i,j));
    }
    fprintf(BI_matrix,"\n");
  }
  
  for(i = 0; i < mat->fm_matrix_rows;i++)
    {
      fprintf(BI_vector,"%lf\n",mat->dense_fm_rhs_vector[i]);
    }
  fclose(BI_matrix);
  fclose(BI_vector);
  
  if (mat->fm_matrix_columns > 0 && mat->fm_matrix_rows > 0) {
  	double* singular_values = new double[mat->fm_matrix_columns];
	calculate_dense_svd(mat, mat->fm_matrix_columns, mat->fm_matrix_rows, mat->dense_fm_matrix, mat->dense_fm_rhs_vector, singular_values);
  	for (i = 0; i < mat->fm_matrix_columns; i++) {
    	mat->fm_solution[solution_counter + i] = mat->dense_fm_rhs_vector[i];
  	}
  	delete [] singular_values;
  	solution_counter += mat->fm_matrix_columns;
  }
  delete mat->dense_fm_matrix;
}

void solve_dense_fm_normal_bootstrapping_equations(MATRIX_DATA* const mat)
{
    double* dd_bak;
    double ttx;
    double* dd1;
    int i,j,k;
    
    // Solve for master
    solve_dense_fm_normal_equations(mat);
    
    //Solve for bootstrapping_estimates.
    double* backup_rhs = new double[mat->fm_matrix_columns];
    double* h = new double[mat->fm_matrix_columns];
    
    // Store a temporary backup of the normal form target vector if it
    // should be output later, since it could be changed in this routine 
    // otherwise.
    if (mat->output_normal_equations_rhs_flag == 1) {
        dd_bak = new double[mat->fm_matrix_columns];
        for (i = 0; i < mat->fm_matrix_columns; i++) dd_bak[i] = mat->dense_fm_normal_rhs_vector[i];
    }
    
    // At the time of solution, one can be sure that the raw dense matrix
    // is no longer needed.
    printf("Freeing raw FM equations.\n");
    fflush(stdout);
    //delete mat->dense_fm_matrix;
    
    if (mat->iterative_calculation_flag == 1) {
        
        // Read in a stored normal form matrix and normal form target vector for
        // iterative calculations
        
        FILE* mat_in = open_file("result.in", "rb");
        dd1 = new double[mat->fm_matrix_columns];
        
        for (j = 0; j < mat->fm_matrix_columns; j++) {
            for (k = 0; k <= j; k++) {
                fread(&ttx, sizeof(double), 1, mat_in);
                mat->dense_fm_normal_matrix->assign_scalar(k, j, ttx);
            }
        }
        for (j = 0; j < mat->fm_matrix_columns; j++) {
            fread(&ttx, sizeof(double), 1, mat_in);
            dd1[j] = ttx;
        }
        fclose(mat_in);
        
        // The target for an iterative calculation is the difference between the targets
        // for this trajectory and the previous trajectory.
        for (i = 0; i < mat->fm_matrix_columns; i++)
            mat->dense_fm_normal_rhs_vector[i] = dd1[i] - mat->dense_fm_normal_rhs_vector[i];
        delete [] dd1;
        
    } else {
    
        // Save the results in binary form for parallel runs.
        if (mat->output_style >= 2) {
            FILE* mat_out = open_file("result.out", "wb");
            for (j = 0; j < mat->bootstrapping_num_estimates; j++) {
            	for (i = 0; i < mat->fm_matrix_columns; i++) {
                	fwrite(&mat->bootstrapping_dense_fm_normal_matrices[j]->values[i * mat->fm_matrix_columns], sizeof(double), i + 1, mat_out);
            	}
            	fwrite(&mat->bootstrapping_dense_fm_normal_rhs_vectors[j][0], sizeof(double), mat->fm_matrix_columns, mat_out);
            }
            double inv_norm = 1.0/mat->normalization;
            fwrite(&mat->force_sq_total, sizeof(double), 1, mat_out);
	        fwrite(&inv_norm, sizeof(double), 1, mat_out);
            fclose(mat_out);
            // If no other output was desired, terminate the program successfully.
            if (mat->output_style == 3) exit(EXIT_SUCCESS);
        }
    }

	mat->bootstrap_solutions = new std::vector<double>[mat->bootstrapping_num_estimates];
   	for (k = 0; k < mat->bootstrapping_num_estimates; k++) {
  		mat->bootstrap_solutions[k] = std::vector<double>(mat->fm_matrix_columns);
	}
	
    for (k = 0; k < mat->bootstrapping_num_estimates; k++) {
    
		// Copy over symmetric off-diagonal values in normal matrix;
	    for (i = 0; i < mat->fm_matrix_columns; i++) {
	        for (j = 0; j < i; j++) {
    	        mat->bootstrapping_dense_fm_normal_matrices[k]->assign_scalar(i, j, mat->bootstrapping_dense_fm_normal_matrices[k]->values[i * mat->fm_matrix_columns + j]);
        	}
    	}

    	// Apply vector regularization.
    	if (mat->regularization_style == 2) {
	    	printf("Regularizing FM normal equations (estimate %d).\n", k);
    		fflush(stdout);
        	for (i = 0; i < mat->fm_matrix_columns; i++) {
    	       	mat->bootstrapping_dense_fm_normal_matrices[k]->add_scalar(i, i, mat->regularization_vector[i]);
        	}
        }

	    // Precondition the normal matrix using the root-of-sum-of-squares 
    	// of the columns as column scaling factors.
    	printf("Preconditioning FM normal equations (estimate %d).\n", k);
    	fflush(stdout);
		calculate_and_apply_dense_preconditioning(mat, mat->bootstrapping_dense_fm_normal_matrices[k], h);
	    
    	// Apply Tikhonov regularization.
    	if (mat->regularization_style == 1) {
	    	printf("Regularizing FM normal equations (estimate %d).\n", k);
    		fflush(stdout);
        	double squared_regularization_parameter;
        	squared_regularization_parameter = mat->tikhonov_regularization_param * mat->tikhonov_regularization_param;
        	for (i = 0; i < mat->fm_matrix_columns; i++) {
    	       	mat->bootstrapping_dense_fm_normal_matrices[k]->add_scalar(i, i, squared_regularization_parameter);
        	}
        }
    
    	// Solve the normal equation by singular value decomposition using LAPACK routines.
    	printf("Computing singular value decomposition of preconditioned, regularized FM normal equations (estimate %d).\n", k);
    	fflush(stdout);
    	double* singular_values = new double[mat->fm_matrix_columns];
    	calculate_dense_svd(mat, mat->fm_matrix_columns, mat->bootstrapping_dense_fm_normal_matrices[k], mat->bootstrapping_dense_fm_normal_rhs_vectors[k], singular_values);
    	
    	// Print singular values.
    	printf("Printing FM singular values (estimate %d).\n", k);
    	fflush(stdout);
    	FILE* solution_file = open_file("sol_info.out", "a");
    	fprintf(solution_file, "Singular vector %d:\n", k);
    	for (i = 0; i < mat->fm_matrix_columns; i++) {
    	    fprintf(solution_file, "%le\n", singular_values[i]);
    	}
    	fclose(solution_file);
   	
   	   	// Clean up the heap-allocated temps.
    	 delete [] singular_values;

	    // Calculate the final results from the singular values.
	    printf("Calculating final FM results (estimate %d).\n", k);
    	fflush(stdout);    
    	for (i = 0; i < mat->fm_matrix_columns; i++) {
    	    mat->bootstrap_solutions[k][i] = mat->bootstrapping_dense_fm_normal_rhs_vectors[k][i] * h[i];
    	}
    	
    	// Calculate and output the residual if requested.
    	if (mat->output_residual == 1) {
	    	double residual = calculate_dense_residual(mat, mat->bootstrapping_dense_fm_normal_matrices[k], backup_rhs, mat->bootstrap_solutions[k], mat->normalization);
	    	printf ("Estimate %d: residual %lf\n", k, residual);
    	}
	}
	
    // Free the preconditioner
    delete [] h;
    delete [] backup_rhs;
    
    // For iterative calculations, the solution is a difference, so the computed quantity
    // should be added on to the previous solution value to obtain the final solution.
    if (mat->iterative_calculation_flag == 1) {
        printf("Adding iterative increment to previous solution.\n");
        fflush(stdout);
        double* x0 = new double[mat->fm_matrix_columns];
        std::ifstream x_in;
        check_and_open_in_stream(x_in, "x.in");
        for (i = 0; i < mat->fm_matrix_columns; i++) x_in >> x0[i];
        x_in.close();
        
        for (i = 0; i < mat->fm_matrix_columns; i++) mat->fm_solution[i] = mat->fm_solution[i] * mat->iterative_update_rate_coeff + x0[i];
        delete [] x0;
    }

    printf("Completed FM.\n"); fflush(stdout);

    if (mat->output_normal_equations_rhs_flag == 1) {
        for (i = 0; i < mat->fm_matrix_columns; i++) mat->dense_fm_normal_rhs_vector[i] = dd_bak[i];
        delete [] dd_bak;
    }
    // Clean up the heap-allocated matrix temps required for 
    // optional more-detailed output.
    for (k = 0; k < mat->bootstrapping_num_estimates; k++) {
	    delete mat->bootstrapping_dense_fm_normal_matrices[k];
    	delete [] mat->bootstrapping_dense_fm_normal_rhs_vectors[k];
	}
	delete [] mat->bootstrapping_dense_fm_normal_rhs_vectors;
    delete [] mat->bootstrapping_dense_fm_normal_matrices;
}

// All the individual frame-block matrices have now been accumulated into a single matrix
// with the condition number of the full-trajectory sparse matrix, as have the
// individual target vectors, so the equations are in a final, solvable form. 

void solve_accumulation_form_fm_equations(MATRIX_DATA* const mat)
{
    int i, j;
    for (i = 0; i < mat->accumulation_matrix_columns; i++) {
        mat->dense_fm_normal_rhs_vector[i] = mat->dense_fm_matrix->values[mat->fm_matrix_columns * mat->accumulation_matrix_rows + i];
    }
    set_accumulation_matrix_to_zero(mat);
    
    // Save the results in binary form and exit if no other output is desired.
    if (mat->output_style >= 2) {
        FILE* mat_out;
        mat_out = open_file("final_equations.out", "wb");
        for (i = 0; i < mat->fm_matrix_columns; i++) {
            fwrite(&mat->dense_fm_matrix->values[i * mat->accumulation_matrix_rows], sizeof(double), i + 1, mat_out);
        }
        fwrite(&mat->dense_fm_normal_rhs_vector[0], sizeof(double), mat->accumulation_matrix_columns, mat_out);
        fclose(mat_out);
        
        if (mat->output_style == 3) exit(EXIT_SUCCESS);
    }
    
    double resid = mat->dense_fm_matrix->values[(mat->accumulation_matrix_columns - 1) * mat->accumulation_matrix_rows + mat->accumulation_matrix_columns - 1];
    
    // Precondition the accumulation matrix using the root-of-sum-of-squares of the columns
    // as column scaling factors.
    double* h = new double[mat->fm_matrix_columns];
    
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        h[i] = 0.0;
    }
    
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        for (j = 0; j < mat->fm_matrix_columns; j++) {
            h[j] = h[j] + mat->dense_fm_matrix->values[j * mat->accumulation_matrix_rows + i] * mat->dense_fm_matrix->values[j * mat->accumulation_matrix_rows + i];
        }
    }
    
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        if (h[i] < VERYSMALL) h[i] = 1.0;
        else h[i] = 1.0 / sqrt(h[i]);
    }
    
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        for (j = 0; j < mat->fm_matrix_columns; j++) {
            mat->dense_fm_matrix->values[j * mat->accumulation_matrix_rows + i] *= h[j];
        }
    }
    
    // Solve the normal equation by singular value decomposition using LAPACK routines.
    // The SVD routine works in two modes: first, the routine is run with a dummy workspace to
    // determine the size of the needed workspace, then that workspace is allocated, then the
    // routine is run again with a sufficient workspace to perform SVD.

    double* singular_values = new double[mat->fm_matrix_columns];
    int onei = 1;
    int irank_in, info_in;
    int lapack_setup_flag = -1;
    double* lapack_temp_workspace = new double[1];
    dgelss_(&mat->fm_matrix_columns, &mat->fm_matrix_columns, &onei, mat->dense_fm_matrix->values, &mat->accumulation_matrix_rows, mat->dense_fm_normal_rhs_vector, &mat->fm_matrix_columns, singular_values, &mat->rcond, &irank_in, lapack_temp_workspace, &lapack_setup_flag, &info_in);
    lapack_setup_flag = lapack_temp_workspace[0];
    delete [] lapack_temp_workspace;
    lapack_temp_workspace = new double[lapack_setup_flag];
    dgelss_(&mat->fm_matrix_columns, &mat->fm_matrix_columns, &onei, mat->dense_fm_matrix->values, &mat->accumulation_matrix_rows, mat->dense_fm_normal_rhs_vector, &mat->fm_matrix_columns, singular_values, &mat->rcond, &irank_in, lapack_temp_workspace, &lapack_setup_flag, &info_in);
    delete [] lapack_temp_workspace;
    
    // Print singular values to file.
    FILE* solution_file;
    solution_file = open_file("sol_info.out", "a");
    fprintf(solution_file, "Singular vector:\n");
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        fprintf(solution_file, "%le\n", singular_values[i]);
    }
    
    // Calculate final results from the singular value decomposition.
    mat->fm_solution = std::vector<double>(mat->fm_matrix_columns);
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        mat->fm_solution[i] = mat->dense_fm_normal_rhs_vector[i] * h[i];
    }
    
    // Print the Euclidean norm of the solution and and print the residual
    double xnorm = 0.0;
    for (i = 0; i < mat->fm_matrix_columns; i++) {
        xnorm += mat->fm_solution[i] * mat->fm_solution[i];
    }
    fprintf(solution_file, "Solution 2-norm:\n%le\n", xnorm);
    
    fprintf(solution_file, "Residual 2-norm:\n%le\n", resid);
    
    // Deallocate the remaining temps.
    fclose(solution_file);
    delete [] singular_values;
    delete [] h;
    delete [] mat->dense_fm_normal_rhs_vector;
    delete mat->dense_fm_matrix;
}

void solve_accumulation_form_bootstrapping_equations(MATRIX_DATA* const mat)
{
    int i, j, k;

 	for (k = 0; k < mat->bootstrapping_num_estimates; k++) {
	    for (i = 0; i < mat->accumulation_matrix_columns; i++) {
    	    mat->bootstrapping_dense_fm_normal_rhs_vectors[k][i] = mat->bootstrapping_dense_fm_normal_matrices[k]->values[mat->fm_matrix_columns * mat->accumulation_matrix_rows + i];
    	}
    	set_accumulation_matrix_to_zero(mat, mat->bootstrapping_dense_fm_normal_matrices[k]);
	}
	    
   	if (mat->output_style >= 2) {
       	FILE* mat_out = open_file("final_equations.out", "wb");
	
		for (k = 0; k < mat->bootstrapping_num_estimates; k++) {
			// Save the results in binary form and exit if no other output is desired.
    	   	for (i = 0; i < mat->fm_matrix_columns; i++) {
        	   	fwrite(&mat->bootstrapping_dense_fm_normal_matrices[k]->values[i * mat->accumulation_matrix_rows], sizeof(double), i + 1, mat_out);
       	 	}
        	fwrite(&mat->bootstrapping_dense_fm_normal_rhs_vectors[k][0], sizeof(double), mat->accumulation_matrix_columns, mat_out);        
    	}
    
    	fclose(mat_out);
	    if (mat->output_style == 3) exit(EXIT_SUCCESS);
	}
    
    mat->bootstrap_solutions = new std::vector<double>[mat->bootstrapping_num_estimates];
   	for (k = 0; k < mat->bootstrapping_num_estimates; k++) {
   		mat->bootstrap_solutions[k] = std::vector<double>(mat->fm_matrix_columns);
	}
	
   	for (k = 0; k < mat->bootstrapping_num_estimates; k++) {

    	double resid = mat->bootstrapping_dense_fm_normal_matrices[k]->values[(mat->accumulation_matrix_columns - 1) * mat->accumulation_matrix_rows + mat->accumulation_matrix_columns - 1];
    
    	// Precondition the accumulation matrix using the root-of-sum-of-squares of the columns
    	// as column scaling factors.
    	double* h = new double[mat->fm_matrix_columns];
    
    	for (i = 0; i < mat->fm_matrix_columns; i++) {
        	h[i] = 0.0;
    	}
    
    	for (i = 0; i < mat->fm_matrix_columns; i++) {
        	for (j = 0; j < mat->fm_matrix_columns; j++) {
            	h[j] = h[j] + mat->bootstrapping_dense_fm_normal_matrices[k]->values[j * mat->accumulation_matrix_rows + i] * mat->bootstrapping_dense_fm_normal_matrices[k]->values[j * mat->accumulation_matrix_rows + i];
        	}
    	}
    
    	for (i = 0; i < mat->fm_matrix_columns; i++) {
        	if (h[i] < VERYSMALL) h[i] = 1.0;
        	else h[i] = 1.0 / sqrt(h[i]);
    	}
    
    	for (i = 0; i < mat->fm_matrix_columns; i++) {
        	for (j = 0; j < mat->fm_matrix_columns; j++) {
            	mat->bootstrapping_dense_fm_normal_matrices[k]->values[j * mat->accumulation_matrix_rows + i] *= h[j];
        	}
    	}
    
    	// Solve the normal equation by singular value decomposition using LAPACK routines.
    	// The SVD routine works in two modes: first, the routine is run with a dummy workspace to
	    // determine the size of the needed workspace, then that workspace is allocated, then the
    	// routine is run again with a sufficient workspace to perform SVD.

    	double* singular_values = new double[mat->fm_matrix_columns];
    	int onei = 1;
    	int irank_in, info_in;
    	int lapack_setup_flag = -1;
    	double* lapack_temp_workspace = new double[1];
    	dgelss_(&mat->fm_matrix_columns, &mat->fm_matrix_columns, &onei, mat->bootstrapping_dense_fm_normal_matrices[k]->values, &mat->accumulation_matrix_rows, mat->bootstrapping_dense_fm_normal_rhs_vectors[k], &mat->fm_matrix_columns, singular_values, &mat->rcond, &irank_in, lapack_temp_workspace, &lapack_setup_flag, &info_in);
    	lapack_setup_flag = lapack_temp_workspace[0];
    	delete [] lapack_temp_workspace;
    	lapack_temp_workspace = new double[lapack_setup_flag];
    	dgelss_(&mat->fm_matrix_columns, &mat->fm_matrix_columns, &onei, mat->bootstrapping_dense_fm_normal_matrices[k]->values, &mat->accumulation_matrix_rows, mat->bootstrapping_dense_fm_normal_rhs_vectors[k], &mat->fm_matrix_columns, singular_values, &mat->rcond, &irank_in, lapack_temp_workspace, &lapack_setup_flag, &info_in);
    	delete [] lapack_temp_workspace;
    
    	// Print singular values to file.
    	FILE* solution_file;
    	solution_file = open_file("sol_info.out", "a");
    	fprintf(solution_file, "Singular vector:\n");
    	for (i = 0; i < mat->fm_matrix_columns; i++) {
        	fprintf(solution_file, "%le\n", singular_values[i]);
    	}

    	// Calculate final results from the singular value decomposition.
    	for (i = 0; i < mat->fm_matrix_columns; i++) {
    	    mat->bootstrap_solutions[k][i] = mat->bootstrapping_dense_fm_normal_rhs_vectors[k][i] * h[i];
    	}
    
    	// Print the Euclidean norm of the solution and and print the residual
    	double xnorm = 0.0;
    	for (i = 0; i < mat->fm_matrix_columns; i++) {
        	xnorm += mat->fm_solution[i] * mat->fm_solution[i];
    	}
    	fprintf(solution_file, "Solution 2-norm:\n%le\n", xnorm);
    
    	fprintf(solution_file, "Residual 2-norm:\n%le\n", resid);
    
    	// Deallocate the remaining temps.
    	fclose(solution_file);
    	delete [] singular_values;
    	delete [] h;
    	delete [] mat->bootstrapping_dense_fm_normal_rhs_vectors[k];
    	delete mat->bootstrapping_dense_fm_normal_matrices[k];
	}
	delete [] mat->bootstrapping_dense_fm_normal_matrices;
	delete [] mat->bootstrapping_dense_fm_normal_rhs_vectors;
}

//--------------------------------------------------------------------
// Binary file reading routines
//--------------------------------------------------------------------

void read_binary_matrix(MATRIX_DATA* const mat)
{
    switch (mat->matrix_type) {
    case kDense: case kSparseNormal:
        read_binary_dense_fm_matrix(mat);
        break;
    case kSparse: case kSparseSparse:
        read_binary_sparse_fm_matrix(mat);
        break;
    case kAccumulation:
        read_binary_accumulation_fm_matrix(mat);
        break;
    default:
        printf("Invalid matrix type for read_binary_matrix.");
        exit(EXIT_FAILURE);
    }
}


// Read the results of a batch of dense-matrix-based FM
// calculations and add them together as if they were the
// results of blocks of an earlier trajectory.

void read_binary_dense_fm_matrix(MATRIX_DATA* const mat)
{

    int i, j, k;
    double tx;
    double inv_norm_sum = 0.0;
    double inv_norm;
    int n_batch;
    char** file_names; //name of matrix files
    FILE* file_of_input_file_names; //input file
    
    // Read the number of files to combine in this batch
    // and the file names for each.
    file_of_input_file_names = open_file("res_av.in", "r");
    fscanf(file_of_input_file_names, "%d", &n_batch);
    file_names = new char*[n_batch];
    for (i = 0; i < n_batch; i++) {
        file_names[i] = new char[100];
    }
    for (i = 0; i < n_batch; i++) {
        fscanf(file_of_input_file_names, "%s", file_names[i]);
    }
    fclose(file_of_input_file_names);

    // Read each file's dense matrix, adding them together element-by-
    // element to get a final set of normal form equations.
    // Each matrix is "un-normalized" by its number of frames before accumulating.
    FILE* single_binary_matrix_input;
	dense_matrix* read_matrix = new dense_matrix(mat->fm_matrix_columns, mat->fm_matrix_columns);
	double* read_rhs = new double[mat->fm_matrix_columns];
    for (i = 0; i < n_batch; i++) {
        read_matrix->reset_matrix();
        // Read the new normal form matrix.
        // Stored as an upper traingular matrix because it is symmetric.
        single_binary_matrix_input = open_file(file_names[i], "rb");
        for (j = 0; j < mat->fm_matrix_columns; j++) {
            for (k = 0; k <= j; k++) {
                fread(&tx, sizeof(double), 1, single_binary_matrix_input);
                read_matrix->add_scalar(k, j, tx);
            }
        }

        // Read the new normal form vector.
        for (j = 0; j < mat->fm_matrix_columns; j++) {
            fread(&tx, sizeof(double), 1, single_binary_matrix_input);
            read_rhs[j] = tx;
        }
        
        // Read the force_sq_total value
        fread(&tx, sizeof(double), 1, single_binary_matrix_input);
        mat->force_sq_total += tx;
        
        // Read the inverse normalization
        fread(&tx, sizeof(double), 1, single_binary_matrix_input);
        inv_norm = tx;
        inv_norm_sum += inv_norm;
        
        // Add the new normal form matrix to the existing one.
        // This process "unnormalizes" each element as it is added.
        // Stored as an upper traingular matrix because it is symmetric.
         for (j = 0; j < mat->fm_matrix_columns; j++) {
            for (k = 0; k <= j; k++) {
            	mat->dense_fm_normal_matrix->add_scalar(k, j, inv_norm * read_matrix->get_scalar(k, j));
            }
        }
        
        // Add the new normal form vector to the existing one.
        // This process "unnormalizes" each element as it is added.
        for (j = 0; j < mat->fm_matrix_columns; j++) {
            mat->dense_fm_normal_rhs_vector[j] += inv_norm * read_rhs[j];
        }
        fclose(single_binary_matrix_input);
        delete [] file_names[i];
    }
    delete [] file_names;
    delete [] read_rhs;
    delete read_matrix;
     
    // Normalize the normal matrix and RHS vector by the total number of frames.
 	set_normalization(mat, 1.0/inv_norm_sum);
	for (j = 0; j < mat->fm_matrix_columns; j++) {
		for (k = 0; k <= j; k++) {
			mat->dense_fm_normal_matrix->assign_scalar(k, j, mat->normalization * mat->dense_fm_normal_matrix->get_scalar(k, j));
		}
	}
	for (j = 0; j < mat->fm_matrix_columns; j++) {
		mat->dense_fm_normal_rhs_vector[j] = mat->normalization * mat->dense_fm_normal_rhs_vector[j];
	}
	    
    // The lower triangular portions of the symmetric normal form matrix
	// filled in during during the solve routine.
}

// Read the results of a batch of accumulation-matrix-based FM
// calculations and add them together as if they were the
// results of blocks of an earlier trajectory.

// Works only for trivial batches of one block each; the primary use
// is to allow for more sophisticated regularization elsewhere in this
// program.

void read_binary_accumulation_fm_matrix(MATRIX_DATA* const mat)
{
    printf("The use of combinefm with the accumulation matrix type is not supported!\n"); 
	delete mat->dense_fm_matrix;
    mat->dense_fm_matrix = new dense_matrix(mat->accumulation_matrix_columns, mat->accumulation_matrix_columns);

    //read the matrix
    int j, k;
    double tx;
    int n_batch; //number of blocks
    char** file_names; //name of matrix files
    FILE* file_of_input_file_names; //input file
    file_of_input_file_names = open_file("res_av.in", "r");
    fscanf(file_of_input_file_names, "%d", &n_batch);
    if (n_batch > 1) {
        printf("Can not read more than one block for the sequential accumulation algorithm!\n");
        exit(EXIT_FAILURE);
    }
    file_names = new char*[n_batch];

    file_names[0] = new char[100]; // Max matrix file name length is 100.

    fscanf(file_of_input_file_names, "%s", file_names[0]);
    fclose(file_of_input_file_names);

    FILE* single_binary_matrix_input;

    single_binary_matrix_input = open_file(file_names[0], "rb");
    for (j = 0; j < mat->fm_matrix_columns; j++) {
        for (k = 0; k <= j; k++) {
            fread(&tx, sizeof(double), 1, single_binary_matrix_input);
            mat->dense_fm_matrix->assign_scalar(k, j, tx);
        }
    }

    for (j = 0; j < mat->fm_matrix_columns; j++) mat->dense_fm_matrix->assign_scalar(mat->fm_matrix_columns, j, 0.0);

    for (j = 0; j < mat->accumulation_matrix_columns; j++) {
        fread(&tx, sizeof(double), 1, single_binary_matrix_input);
        mat->dense_fm_normal_rhs_vector[j] = tx;
    }
    fclose(single_binary_matrix_input);
    delete [] file_names[0];
    delete [] file_names;
}

// Read the results of a batch of sparse-matrix-based FM
// calculations and add them together as if they were the
// results of blocks of an earlier trajectory.

void read_binary_sparse_fm_matrix(MATRIX_DATA* const mat)
{
    printf("The use of combinefm with the sparse matrix type is not supported!\n"); 
	int i, j;
    int n_batch;
    FILE* file_of_input_file_names, *single_block_solution_file;
    char single_block_solution_filename[100];
    double* single_block_normalization_factors;
    
    // Allocate memory for a single block's worth of temp data.
    single_block_normalization_factors = new double[mat->fm_matrix_columns];
    
    // Open the file defining the batch size and batch file names.
    file_of_input_file_names = open_file("res_av.in", "r");
    fscanf(file_of_input_file_names, "%d", &n_batch);
    
    // For each file in the batch, read the appropriate file
    // for the solution of that batch (not normalized) and the
    // normalization factors for that solution.
    for (i = 0; i < n_batch; i++) {
        
        // Read the file
        fscanf(file_of_input_file_names, "%s", single_block_solution_filename);
        single_block_solution_file = open_file(single_block_solution_filename, "rb");
        fread(mat->block_fm_solution, sizeof(double), mat->fm_matrix_columns, single_block_solution_file);
        fread(single_block_normalization_factors, sizeof(double), mat->fm_matrix_columns, single_block_solution_file);
        fclose(single_block_solution_file);
        
        // Add that to the accumulating solution in this program
        for (j = 0; j < mat->fm_matrix_columns; j++) {
            mat->fm_solution[j] += mat->block_fm_solution[j];
            mat->fm_solution_normalization_factors[j] += single_block_normalization_factors[j];
        }
    }
    
    fclose(file_of_input_file_names);
    delete [] single_block_normalization_factors;
}

// Read the vector of regularization coefficients.

void read_regularization_vector(MATRIX_DATA* const mat)
{
    mat->regularization_vector = new double[mat->fm_matrix_columns];
    std::ifstream lambda_in;
    check_and_open_in_stream(lambda_in, "lambda.in");
    for (int i = 0; i < mat->fm_matrix_columns; i++) lambda_in >> mat->regularization_vector[i];
    lambda_in.close();
}

void write_iteration(const double* alpha_vec, const double beta, std::vector<double> fm_solution, const double residual, const int iteration, FILE* alpha_fp, FILE* beta_fp, FILE* sol_fp, FILE* res_fp)
{
	int size = fm_solution.size();
	
	fprintf(alpha_fp, "Iteration %d: ", iteration);
	for (int i = 0; i < size; i++) {
		fprintf(alpha_fp, "%lf ", alpha_vec[i]);
	}
	fprintf(alpha_fp, "\n");
	
	fprintf(beta_fp, "Iteration %d: %lf\n", iteration, beta);
	
	fprintf(sol_fp, "Iteration %d: ", iteration);
	for (int i = 0; i < size; i++) {
		fprintf(sol_fp, "%lf ", fm_solution[i]);
	}
	fprintf(sol_fp, "\n");

	fprintf(res_fp, "Iteration %d: %lf\n", iteration, residual);
}

void determine_BI_interaction_rows_and_cols(MATRIX_DATA* mat, InteractionClassComputer* const icomp)
{
  int num_entries = 0;
  // Skip if it does not have a parameter distribution to use
  if (icomp->ispec->output_parameter_distribution !=  0) {
    // For every defined interaction,
    for (unsigned i = 0; i < icomp->ispec->defined_to_matched_intrxn_index_map.size(); i++) {
      num_entries += (int)(0.5 + (icomp->ispec->upper_cutoffs[i] - icomp->ispec->lower_cutoffs[i])/ icomp->ispec->get_fm_binwidth());
    }
  }
  mat->fm_matrix_rows = num_entries;
  mat->fm_matrix_columns = icomp->ispec->get_num_basis_func();
}
