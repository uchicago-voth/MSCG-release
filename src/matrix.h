//
//  matrix.h
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//
#ifndef _matrix_h
#define _matrix_h

#include <vector>

#include "external_matrix_routines.h"

#if _mkl_flag == 1
#include "mkl.h"
#endif

struct CG_MODEL_DATA;
struct ControlInputs;

//-------------------------------------------------------------
// Matrix-equation-related type definitions
//-------------------------------------------------------------

enum MatrixType {kDense = 0, kSparse = 1, kAccumulation = 2, kSparseNormal = 3, kSparseSparse = 4, kDummy = -1};

// Linked-list-based sparse row matrix element struct. x,y,z components are stored together.

struct linked_list_sparse_matrix_element { 
    int col;                                        // Column number
    double valx[3];                                 // x,y,z components
    struct linked_list_sparse_matrix_element* next; // Pointer to the next element in the linked list
};

// Linked-list-based sparse row matrix row head struct, one linked list for each row.

struct linked_list_sparse_matrix_row_head { 
    int n;                                          // Total number of nonzero elements in this row
    struct linked_list_sparse_matrix_element* h;    // Link to first nonzero element in the row
};

// CSR sparse matrix struct w/ constructor & destructor.

struct csr_matrix {
    int n_rows;
    int n_cols;
    int max_entries;
    double *values;
    int *column_indices;
    int *row_sizes;

    inline csr_matrix(const int new_n_rows, const int new_n_cols, const int new_max_entries) : 
        n_rows(new_n_rows), n_cols(new_n_cols), max_entries(new_max_entries) {
        
        values = new double[max_entries]();
        column_indices = new int[max_entries]();
        row_sizes = new int[n_rows + 1]();
    	row_sizes[0] = 1;
	};
	
    inline csr_matrix(const csr_matrix& copy_matrix) {
    	n_rows = copy_matrix.n_rows;
    	n_cols = copy_matrix.n_cols;
    	max_entries = copy_matrix.max_entries;
    	values = copy_matrix.values;
    	column_indices = copy_matrix.column_indices;
    	row_sizes = copy_matrix.row_sizes;
    };

    inline csr_matrix(const int new_n_rows, const int new_n_cols, const int new_max_entries, double* copy_values, int* copy_column_indices, int* copy_row_sizes) :
    	n_rows(new_n_rows), n_cols(new_n_cols), max_entries(new_max_entries), values(copy_values), column_indices(copy_column_indices), row_sizes(copy_row_sizes) {
    };

    inline void set_csr_matrix(const int new_n_rows, const int new_n_cols, const int new_max_entries, double* copy_values, int* copy_column_indices, int* copy_row_sizes) {
    	delete [] values;
    	delete [] column_indices;
    	delete [] row_sizes;
    	
    	n_rows = new_n_rows;
    	n_cols = new_n_cols; 
    	max_entries = new_max_entries;
    	
    	values = copy_values;
    	column_indices = copy_column_indices;
    	row_sizes = copy_row_sizes; 
    };

    inline ~csr_matrix() {
        delete [] values;
        delete [] column_indices;
        delete [] row_sizes;
    };
};

struct dense_matrix {
    int n_rows;
    int n_cols;
    double *values;

    inline dense_matrix(const int new_n_rows, const int new_n_cols) : 
        n_rows(new_n_rows), n_cols(new_n_cols) {
        values = new double[n_rows * n_cols]();
    }

    inline dense_matrix(const dense_matrix& copy_matrix) {
    	n_rows = copy_matrix.n_rows;
    	n_cols = copy_matrix.n_cols;
    	values = copy_matrix.values;
    }

    inline dense_matrix(const int new_n_rows, const int new_n_cols, double* copy_values) :
    	n_rows(new_n_rows), n_cols(new_n_cols), values(copy_values) {
    }

	inline void print_matrix(FILE* fh) {
		fprintf(fh, "Printing Matrix:\n");
		for (int i = 0; i < n_rows; i++) {
			fprintf(fh, "Row %d: ", i);
			for (int j = 0; j < n_cols; j++) {
				fprintf(fh, "%lf\t", values[j*n_rows + i]);
			}
			fprintf(fh, "\n");
		}
	}
	
    inline void reset_matrix() {
    	for (int k = 0; k < n_rows * n_cols; k++) {
	        values[k] = 0.0;
    	}
    }
    
    inline void add_3vector(const int row, const int col, double const* x) {
    	for(int i = 0; i < 3; i++) {
    		values[ col * n_rows + row + i] += x[i];
    	}
    }

    inline void assign_3vector(const int row, const int col, double const* x) {
    	for(int i = 0; i < 3; i++) {
    		values[ col * n_rows + row + i] = x[i];
    	}
    }

	inline void add_scalar(const int row, const int col, const double x) {
		values[ col * n_rows + row] += x;
	}

	inline void assign_scalar(const int row, const int col, const double x) {
		values[ col * n_rows + row] = x;
	}
	
    inline ~dense_matrix() {
        delete [] values;
	}
};

struct MATRIX_DATA {
    // Poor-man's polymorphism.
    MatrixType matrix_type;
    void (*do_end_of_frameblock_matrix_manipulations)(MATRIX_DATA*);        // A matrix-implementation-dependent function called at the end of each frame block
    void (*accumulate_virial_constraint_matrix_element)(const int, const int, const double, MATRIX_DATA*);      // A matrix-implementation-dependent function to add virial constraint elements to the force matching matrix
    void (*accumulate_fm_matrix_element)(const int, const int, double* const, MATRIX_DATA*);                     // A matrix-implementation-dependent function to add three-vectors to the force matching matrix
    void (*accumulate_target_force_element)(MATRIX_DATA*, int, double *);                  // A matrix-implementation-dependent function to add a single force three-vector to the target force vector.
    void (*accumulate_target_constraint_element)(MATRIX_DATA*, int, double);             // A matrix-implementation-dependent function to add a single conjugate force against a constraint to the target vector.
    void (*set_fm_matrix_to_zero)(MATRIX_DATA*);            // A matrix-implementation-dependent function to reset the current force matching matrix to zero between blocks.
    void (*finish_fm)(MATRIX_DATA*);

    int frames_per_traj_block;              // Number of frames to read in a single block of FM matrix construction

    // Basic layout implementation details
    int fm_matrix_rows;                             // Number of rows for FM matrix
    int fm_matrix_columns;                          // Number of columns for FM matrix
    int rows_less_virial_constraint_rows;           // Rows less the rows reserved for virial constraints
    int virial_constraint_rows;                     // Rows specifically for virial constraints
    
    // For dense-matrix-based calculations
    dense_matrix* dense_fm_matrix;
    double* dense_fm_rhs_vector;
    dense_matrix* dense_fm_normal_matrix;                 // Normal form of the force-matching matrix. Constructed one frame at a time.
    double* dense_fm_normal_rhs_vector;             // Normal form of the target force vector. Constructed one frame at a time.
    double normalization;
    double* fm_solution_normalization_factors;      // Weighted number of times each unknown has been found nonzero in the solution vectors of all blocks
    std::vector<double> fm_solution;                // Final answers averaged over all blocks
	
    // Optional extras for any matrix_type
    int use_statistical_reweighting;        // 1 to use per-frame statistical reweighting; 0 otherwise
    int dynamic_state_samples_per_frame;	// Number of times a frame is resampled. This is 1 unless dynamic_state_sampling is 1.

    // Optional extras for dense-matrix-based calculations
    double current_frame_weight;
    int iterative_calculation_flag;         // 0 for a non-iterative calculation; 1 to use Lanyuan's iterative force matching method
    double iterative_update_rate_coeff;                     // The parameter setting the aggressiveness of the fixed-point iteration used in Lanyuan's iterative force matching method
	
	// Optional extras for bootstrapping (dense and sparse)
	int bootstrapping_flag;
	int bootstrapping_full_output_flag;
	int bootstrapping_num_estimates;
	double* bootstrapping_normalization;
	double** bootstrapping_weights;
	double** bootstrapping_dense_fm_normal_rhs_vectors;
	dense_matrix** bootstrapping_dense_fm_normal_matrices;
	csr_matrix** bootstrapping_sparse_fm_normal_matrices;
	std::vector<double>* bootstrap_solutions;
	
    // For sparse-matrix-based calculations
    int max_nonzero_normal_elements;                // Total number of nonzero values in the sparse normal matrix
	int min_nonzero_normal_elements;				// Lower bound for safe size of sparse normal matrix
	int num_sparse_threads;							// Number of threads for sparse solver
	int itnlim;										// Maximum number of iterative refinement
	double sparse_safety_factor;					// % to oversize the next frame-block's normal matrix from the current one (matrix_type = 4)
	struct linked_list_sparse_matrix_row_head* ll_sparse_matrix_row_heads;      // A linked-list-based sparse matrix
   	csr_matrix* sparse_matrix;						// CSR matrix "object" (matrix_type = 4)
	double* block_fm_solution;                      // FM solutions from one single block
    double* h;                                      // Temp for preconditioning
	
    // For accumulation-matrix-based calculations
    double fm_residual;                             // Final MS-CG residual value
    int trajectory_block_index;                     // Index of the current frame block, needed to calculate compositions of accumulation matrices
    int accumulation_matrix_rows;
    int accumulation_matrix_columns;
    int accumulation_row_shift;
    int accumulation_target_forces_location;
    int lapack_setup_flag;                          // Temp for LAPACK SVD and QR routines
    double* lapack_temp_workspace;                  // Temp for LAPACK SVD and QR routines
    double* lapack_tau;                             // Temp for LAPACK SVD and QR routines

    // Regularization parameters
    double tikhonov_regularization_param;           // Parameter for Tikhonov regularization.
    int regularization_style;                       // 0 to use no regularization; 1 to calculate results using single scalar regularization; 2 to calculate results using a set of regularization parameters in file lambda.in

    // SVD routine parameter
    double rcond;                           // SVD condition number threshold
    
    // Output specifications for matrix-based routines
    int output_style;                       // 0 to output only tables; 2 to output tables and binary block equations; 3 to output only binary block equations
    int output_normal_equations_rhs_flag;   // 1 to output the final right hand side vector of the MS-CG normal equations as well as force tables; 0 otherwise
    int output_solution_flag;               // 0 to not output the solution vector; 1 to output the solution vector in x.out
    
    // Modification functions for matrix.
    void resize_matrix(const int new_cg_sites, const int new_cols) {
    	if (new_cols != fm_matrix_columns) {
    		printf("WARNING: Resizing the number of matrix columns!\n");
    		printf("This will alter the size of the normal matrix.\n");
    		printf("Currently, this is not allowed.\n");
    		exit(EXIT_FAILURE);
    		fm_matrix_columns = new_cols;
    	}
    	
    	int new_rows_less_virial = new_cg_sites * 3 * frames_per_traj_block;
    	if (new_rows_less_virial + virial_constraint_rows == fm_matrix_rows) {
    		printf("Resize_matrix is not doing anything since matrix is the same size as before.\n");
    		return;
    	}
        // Update mat's copy of row variables.	
    	rows_less_virial_constraint_rows = new_rows_less_virial;
    	fm_matrix_rows = new_rows_less_virial + virial_constraint_rows;

		// Update dense_fm_rhs_vector.
		delete [] dense_fm_rhs_vector;
		dense_fm_rhs_vector = new double[fm_matrix_rows]();
		
		// Update the appropriate fm_matrix based on type.
		if (matrix_type == kDense) {
		    delete dense_fm_matrix;
		    dense_fm_matrix = new dense_matrix(fm_matrix_rows, fm_matrix_columns);
		} else if ( (matrix_type == kSparse) || (matrix_type == kSparseNormal) || (matrix_type == kSparseSparse) ) {
			if (sparse_matrix != NULL) {
				int max_entries = sparse_matrix->max_entries;
				delete sparse_matrix;   		
				sparse_matrix = new csr_matrix(fm_matrix_rows, fm_matrix_columns, max_entries);
	    	}
	    } else if (matrix_type == kAccumulation) {
	    	accumulation_matrix_rows = fm_matrix_rows;
	    	printf("Resizing of accumulation matricies is not supported.\n");
	    	exit(EXIT_FAILURE);
	    }
    }
};

// Matrix initialization and reset routines

MATRIX_DATA* make_matrix(ControlInputs* const control_input, CG_MODEL_DATA* const cg);

// Set FM normalization constant.

inline void set_normalization(MATRIX_DATA* mat, const double new_normalization_constant) {
    mat->normalization = new_normalization_constant;
}

void set_bootstrapping_normalization(MATRIX_DATA* mat, double** const bootstrapping_weights, int const n_frames);

// Target (RHS) vector calculation routines

void add_target_virials_from_trajectory(MATRIX_DATA* const mat, double *pressure_constraint_rhs_vector);
void add_target_force_from_trajectory(int shift_i, int site_i, MATRIX_DATA* const mat, rvec *f);

// Read serialized, partially-completed post-frameblock matrix calculation intermediates

void read_binary_matrix(MATRIX_DATA* const mat);

#endif
