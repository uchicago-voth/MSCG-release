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


#ifndef DIMENSION
#define DIMENSION 3
#endif

struct CG_MODEL_DATA;
struct ControlInputs;

typedef void (*accumulate_forces)(InteractionClassComputer* const info, const int first_nonzero_basis_index, const std::vector<double> &basis_fn_vals, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
typedef void (*accumulate_table_forces)(InteractionClassComputer* const info, const double &table_fn_val, const int n_body, const int* particle_ids, std::array<double, DIMENSION>* const &derivatives, MATRIX_DATA * const mat);
void initialize_BI_matrix(MATRIX_DATA* const mat, CG_MODEL_DATA* const cg);

//-------------------------------------------------------------
// Matrix-equation-related type definitions
//-------------------------------------------------------------

enum MatrixType {kDense = 0, kSparse = 1, kAccumulation = 2, kSparseNormal = 3, kSparseSparse = 4, kREM = 5, kDummy = -1};

// Linked-list-based sparse row matrix element struct. x,y,z components are stored together.

struct linked_list_sparse_matrix_element { 
    int col;                                        // Column number
    double valx[DIMENSION];                                 // x,y,z components
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
    }

    inline csr_matrix(const csr_matrix& copy_matrix) {
    	n_rows = copy_matrix.n_rows;
    	n_cols = copy_matrix.n_cols;
    	max_entries = copy_matrix.max_entries;
    	values = copy_matrix.values;
    	column_indices = copy_matrix.column_indices;
    	row_sizes = copy_matrix.row_sizes;
    }

    inline csr_matrix(const int new_n_rows, const int new_n_cols, const int new_max_entries, double* copy_values, int* copy_column_indices, int* copy_row_sizes) :
    	n_rows(new_n_rows), n_cols(new_n_cols), max_entries(new_max_entries), values(copy_values), column_indices(copy_column_indices), row_sizes(copy_row_sizes) {
    }
    	
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
    }

	inline void write_file(FILE* fh) {
		fprintf(fh, "%d %d %d\n", n_rows, n_cols, max_entries);
		for (int i = 0; i < row_sizes[n_rows]; i++) {
			fprintf(fh, "%lf ", values[i]);
		}
		fprintf(fh, "\n");
		for (int i = 0; i < row_sizes[n_rows]; i++) {
			fprintf(fh, "%d ", column_indices[i]);
		}
		fprintf(fh, "\n");
		for (int i = 0; i <= n_rows; i++) {
			fprintf(fh, "%d ", row_sizes[i]);
		}
		fprintf(fh, "\n");
	}
	
	inline void print_matrix(FILE* fh) {
/*	
		fprintf(fh, "Printing Matrix:\n");
		for (int i = 0; i < n_rows; i++) {
			fprintf(fh, "Row %d: ", i);
			for (int j = 0; j < n_cols; j++) {
				fprintf(fh, "%lf\t", values[j*n_rows + i]);
			}
			fprintf(fh, "\n");
		}
*/
	}
	
	inline void print_matrix_csr(FILE* fh) {
		write_file(fh);
	}
	
    inline ~csr_matrix() {
        delete [] values;
        delete [] column_indices;
        delete [] row_sizes;
    }
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
		fprintf(fh, "Rows: %d\nColumns: %d\n", n_rows,  n_cols);
		for (int i = 0; i < n_rows; i++) {
			fprintf(fh, "Row %d: ", i);
			for (int j = 0; j < n_cols; j++) {
				fprintf(fh, "%lf\t", values[j*n_rows + i]);
			}
			fprintf(fh, "\n");
		}
	}
	
	inline void print_matrix_csr(FILE* fh) {
		
		//Allocate this CSR
		int nnz = get_nnz();
		int counter = 0;
		int column_indices[nnz];
		int row_sizes[n_rows + 1];
		row_sizes[0] = 0;
		
		fprintf(fh, "%d %d %d\n", n_rows, n_cols, nnz);
		
		//Populate this CSR and print non-zero values
		for (int i = 0; i < n_rows; i++) {
			for(int j = 0; j < n_cols; j++) {
				if (fabs(values[j*n_rows + i]) > VERYSMALL)  {
					fprintf(fh,"%lf ",  values[j*n_rows +  i]);
					column_indices[counter++] = values[j*n_rows + i];
				}
			}
			row_sizes[i+1] = counter;
		}
		fprintf(fh,"\n");
		
		// print the column indices
		for (int i = 0; i < nnz; i++) {
			fprintf(fh,"%d ", column_indices[i]);
		}
		fprintf(fh,"\n");
		
		// print the row_sizes
		for (int i = 0; i <= n_rows; i++) {
			fprintf(fh,"%d ", row_sizes[i]);
		}
		fprintf(fh,"\n");
	}
	
	inline int get_nnz() {
		int nnz = 0;
		for (int i = 0; i < n_rows; i++) {
			for (int j = 0; j < n_cols; j++) {
				if (fabs(values[j*n_rows + i]) > VERYSMALL) nnz++;
			}
		}
		return nnz;
	}
	
    inline void reset_matrix() {
    	for (int k = 0; k < n_rows * n_cols; k++) {
	        values[k] = 0.0;
    	}
    }
    
    void add_vector(const int row, const int col, const double* const x) {
    	for(int i = 0; i < DIMENSION; i++) {
    		values[ col * n_rows + row + i] += x[i];
    	}
    }

    void assign_vector(const int row, const int col, const std::array<double, DIMENSION> &x) {
    	for(int i = 0; i < DIMENSION; i++) {
    		values[ col * n_rows + row + i] = x[i];
    	}
    }


	inline void add_scalar(const int row, const int col, const double x) {
		values[ col * n_rows + row] += x;
	}

	inline void assign_scalar(const int row, const int col, const double x) {
		values[ col * n_rows + row] = x;
	}
	
	inline double get_scalar(const int row, const int col) {
		return values[ col * n_rows + row];
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

	// Poor-man's polymorphism for scalar vs vector matrix operations.
	accumulate_forces accumulate_matching_forces;
	accumulate_forces accumulate_symmetric_matching_forces;
	accumulate_forces accumulate_one_body_force;
	accumulate_table_forces accumulate_tabulated_forces;
	accumulate_table_forces accumulate_symmetric_tabulated_forces;
	accumulate_table_forces accumulate_one_body_tabulated_force;

    // Basic layout implementation details
    int fm_matrix_rows;                             // Number of rows for FM matrix
    int fm_matrix_columns;                          // Number of columns for FM matrix
    int rows_less_virial_constraint_rows;           // Rows less the rows reserved for virial constraints
    int virial_constraint_rows;                     // Rows specifically for virial constraints
    int frames_per_traj_block;              // Number of frames to read in a single block of FM matrix construction
    int position_dimension;							// The number of elements needed to specify each particle's position.
    int size_per_vector;							// Either 1 or DIMENSION based on scalar_matching_flag(1 or 0, respectively);

    // For dense-matrix-based calculations
    dense_matrix* dense_fm_matrix;
    dense_matrix* dense_fm_normal_matrix;                 // Normal form of the force-matching matrix. Constructed one frame at a time.
    double* dense_fm_rhs_vector;
    double* dense_fm_normal_rhs_vector;             // Normal form of the target force vector. Constructed one frame at a time.
    double normalization;
    double* fm_solution_normalization_factors;      // Weighted number of times each unknown has been found nonzero in the solution vectors of all blocks
    std::vector<double> fm_solution;                // Final answers averaged over all blocks
	
	// REM variables
	std::vector<double> previous_rem_solution;        // Previous splie coeffs to be used in REM iteration
    double temperature;
    double rem_chi;
    double boltzmann;
    
    // Optional extras for any matrix_type
    int use_statistical_reweighting;        // 1 to use per-frame statistical reweighting; 0 otherwise
    int dynamic_state_samples_per_frame;	// Number of times a frame is resampled. This is 1 unless dynamic_state_sampling is 1.
    int volume_weighting_flag;				// 1 to use volume weighting following the approach in MS-CG V; 0 otherwise
	
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

	// Optional extras for residual, regularization, and bayesian calculations
	int output_residual;							// 1 to calculate the residual; 0 otherwise
	double force_sq_total;							
	int bayesian_flag;								// 1 to use Bayesian MS-CG to calculate regularization and interactions
	int bayesian_max_iter;
    int regularization_style;                       // 0 to use no regularization; 1 to calculate results using single scalar regularization; 2 to calculate results using a set of regularization parameters in file lambda.in
	double tikhonov_regularization_param;           // Parameter for Tikhonov regularization. (regularization_style = 1)
	double* regularization_vector;					// Vector for regularization_style 2.

    // SVD routine parameter
    double rcond;                           // SVD condition number threshold
    
    // Output specifications for matrix-based routines
    int output_style;                       // 0 to output only tables; 2 to output tables and binary block equations; 3 to output only binary block equations
    int output_normal_equations_rhs_flag;   // 1 to output the final right hand side vector of the MS-CG normal equations as well as force tables; 0 otherwise
    int output_solution_flag;               // 0 to not output the solution vector; 1 to output the solution vector in x.out
    int output_raw_splines;					// 1 to output spline contributions for each interaction processed; 0 otherwise  (default)
   	int output_raw_frame_blocks;			// 1 to output the pre-normal form fm_matrix  at the end of each frame block; 0 otherwise (default)
   	FILE* frame_block_fh;
   	

	// Constructors and destructors
	MATRIX_DATA(ControlInputs* const control_input, CG_MODEL_DATA *const cg);
	
	~MATRIX_DATA() {
		if (bootstrapping_flag == 1) {
			delete [] bootstrap_solutions;
    		delete [] bootstrapping_normalization;
    	}
    	if (regularization_style == 2) {
	   		delete [] regularization_vector;
	   	}
	   	
   	    // Free FM matrix building temps
	    printf("Freeing equation building temporaries.\n");

    	if (matrix_type == kDense) {
			if (virial_constraint_rows > 0) delete dense_fm_matrix;
			delete [] dense_fm_rhs_vector;
		} else if (matrix_type == kSparse) {
			delete [] ll_sparse_matrix_row_heads;
			delete [] block_fm_solution;
			delete [] dense_fm_rhs_vector;
			if (virial_constraint_rows > 0) delete dense_fm_matrix;
		} else if (matrix_type == kAccumulation) {
			delete [] lapack_temp_workspace;
			delete [] lapack_tau;
		} else if (matrix_type == kSparseNormal) {
			delete [] ll_sparse_matrix_row_heads;
			delete [] dense_fm_rhs_vector;
			if (virial_constraint_rows > 0) delete dense_fm_matrix;
		} else if (matrix_type == kSparseSparse) {
			delete [] ll_sparse_matrix_row_heads;
			delete [] dense_fm_rhs_vector;
			if (virial_constraint_rows > 0) delete dense_fm_matrix;
		} else if (matrix_type == kDummy) {
		    delete [] dense_fm_rhs_vector;
			delete [] dense_fm_normal_rhs_vector;
		} else if (matrix_type == kREM) {
	    	delete dense_fm_matrix;
	 	}
	 	if  (output_raw_frame_blocks == 1) {
 			fclose(frame_block_fh);
 		}
	}
   
    // Modification function for matrix.
    void resize_matrix(const int new_cg_sites, const int new_cols) {
    	if (new_cols != fm_matrix_columns) {
    		printf("WARNING: Resizing the number of matrix columns!\n");
    		printf("This will alter the size of the normal matrix.\n");
    		printf("Currently, this is not allowed.\n");
    		exit(EXIT_FAILURE);
    		fm_matrix_columns = new_cols;
    	}
    	
    	int new_rows_less_virial = new_cg_sites * DIMENSION * frames_per_traj_block;
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
    
	inline double get_frame_weight(void) {
    	if (use_statistical_reweighting) {
        	return current_frame_weight;
    	} else {
        	return 1.0;
    	}
	}
};

// Set FM normalization constant.

inline void set_normalization(MATRIX_DATA* mat, const double new_normalization_constant) {
    mat->normalization = new_normalization_constant;
}

void set_bootstrapping_normalization(MATRIX_DATA* mat, double** const bootstrapping_weights, int const n_frames);

// Target (RHS) vector calculation routines

void add_target_virials_from_trajectory(MATRIX_DATA* const mat, double *pressure_constraint_rhs_vector);
void add_target_force_from_trajectory(int shift_i, int site_i, MATRIX_DATA* const mat, std::array<double, DIMENSION>* const &f);

// Read serialized, partially-completed post-frameblock matrix calculation intermediates

void read_binary_matrix(MATRIX_DATA* const mat);

#endif
