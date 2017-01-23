//
//  mscg.cpp
//  
//
//  Created by Jacob Wagner on 04/27/16 for use with USER-MSCG package in LAMMPS.
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

// This file provides the library interface to the multi-scale coarse-graining (MS-CG)
// code, which is also referred to as the force matching (FM) code. 
// The order in which functions are intended to be called for FM is the following:
// mscg_startup_part1
// setup_topology_and_frame
// setup_bond_topology, setup_angle_topology, and setup_dihedral_topology
// either setup_exclusion_topology or generate_exclusion_topology
// mscg_startup_part2
// mscg_process_frame for each frame
// mscg_solve_and output.
// Additional functions are provided for updating or modifying certain information
// after the data is initially set using one of the functions before or during 
// mscg_startup_part2.

// The order in which functions are intended to be called for range finding is the following:
// rangefinder_startup_part1
// setup_topology_and_frame
// setup_bond_topology, setup_angle_topology, and setup_dihedral_topology
// either setup_exclusion_topology or generate_exclusion_topology
// rangefinder_startup_part2
// rangefinder_process_frame for each frame
// rangefinder_solve_and output.

#include "mscg.h"

// Prototype function definition for functions called internal to this file
void finish_fix_reading(FrameSource *const frame_source);

// Data structure holding all MSCG information.
// It is passed to the driver function (LAMMPS fix) as an opaque pointer.
struct MSCG_struct {
	int curr_frame;
	int nblocks;
	int trajectory_block_frame_index;
	int traj_frame_num;
	double start_cputime;
	FrameSource *frame_source;      // Trajectory frame data
    CG_MODEL_DATA *cg;  			// CG model parameters and data
    ControlInputs *control_input;	// Input settings read from control.in
    MATRIX_DATA *mat;				// Matrix storage structure
};

// This function starts the MSCG process by allocating memory for the mscg_struct
// and reading information from the control.in file.
// It should be called first.
void* mscg_startup_part1(void* void_in) 
{	
    // Begin to compute the total run time
    MSCG_struct* mscg_struct = new MSCG_struct;
    mscg_struct->start_cputime = clock();	
	
	mscg_struct->frame_source = new FrameSource;
    
    FrameSource *p_frame_source = mscg_struct->frame_source;
    ControlInputs *p_control_input = mscg_struct->control_input;
    
    //----------------------------------------------------------------
    // Set up the force matching procedure
    //----------------------------------------------------------------
        
    // Read the file control.in to determine the values of many of
    // the parameters in the cg struct, including the desired
    // force-matching strategy (sparse block averaged, normal 
    // equations, accumulation matrices), the desired output styles,
    // the type of basis set to use, and many others described in
    // types.h and listed in control_input.c.
    printf("Reading high level control parameters.\n");
    mscg_struct->control_input = new ControlInputs;
    mscg_struct->cg = new CG_MODEL_DATA(p_control_input);   // CG model parameters and data; put here to initialize without default constructor
    
    copy_control_inputs_to_frd(p_control_input, p_frame_source);
    
	return (void*)(mscg_struct);
}
    
    
// This function starts the range finding utility.
// Since this is almost identical to MSCG setup,
// it is a thin wrapper of mscg_setup_part1
void* rangefinder_startup_part1(void* void_in)
{
	void_in = mscg_startup_part1(void_in);
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
	if (mscg_struct->control_input->three_body_flag != 0) {
        printf("Rangefinder does not support three body nonbonded interaction ranges.\n");
        exit(EXIT_FAILURE);
    }
    mscg_struct->control_input->bootstrapping_flag = 0;
    mscg_struct->frame_source->bootstrapping_flag = 0;
    return void_in;
}

// This function continues the setup for MSCG
// after the mscg_startup_part1, setup_topology_and_frame, and setup_*_topology
// functions have already been called.
// It reads interaction range information from the rmin*.in files,
// any tabulated interactions from table.in, any frame weights from frame_weights.in,
// and virial constraint information from p_con.in.
// Also, the sets up the "ForceComputers" and initializes the matrix.
// In addition, it sets up bootstrapping samples, if needed. 
void* mscg_startup_part2(void* void_in)
{
	int total_frame_samples = 0;
	int n_blocks = 0;

	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
	FrameSource *p_frame_source = mscg_struct->frame_source;
    ControlInputs *p_control_input = mscg_struct->control_input;
	CG_MODEL_DATA *p_cg = mscg_struct->cg;

    // Read the range files rmin.in and rmax.in to determine the
    // ranges over which the FM basis functions should be defined.
    // These ranges are also used to record which interactions
    // should be fit, which should be tabulated, and which are not 
    // present in the model.
    printf("Reading interaction ranges.\n");
    read_all_interaction_ranges(p_cg);

    // The range files have specified which interactions should be
    // read from file; read the tabulated potentials from table.in
    // now if any were found.
    if (p_cg->pair_nonbonded_interactions.n_tabulated > 0 ||
        p_cg->pair_bonded_interactions.n_tabulated > 0 ||
        p_cg->angular_interactions.n_tabulated > 0 ||
        p_cg->dihedral_interactions.n_tabulated > 0) {
        printf("Reading tabulated reference potentials.\n");
        read_tabulated_interaction_file(p_cg, p_cg->topo_data.n_cg_types);
    } 
    
    // Read statistical weights for each frame if the 
    // 'use_statistical_reweighting' flag is set in control.in.
    if (p_frame_source->use_statistical_reweighting == 1) {
        printf("Reading per-frame statistical reweighting factors.\n");
        fflush(stdout);
        read_frame_weights(p_frame_source, p_control_input->starting_frame, p_control_input->n_frames); 
    }
    
    // Generate bootstrapping weights if the
    // 'bootstrapping_flag' is set in control.in.
    if (p_frame_source->bootstrapping_flag == 1) {
    	printf("Generating bootstrapping frame weights.\n");
    	fflush(stdout);
    	generate_bootstrapping_weights(p_frame_source, p_control_input->n_frames);
    }

    // Read the input virials if the correct flag was set in control.in.
    if (p_frame_source->pressure_constraint_flag == 1) {
        printf("Reading virial constraint target.\n");
        read_virial_constraint_vector(p_frame_source, p_control_input->starting_frame, p_control_input->n_frames);
    }
    
    // Use the trajectory type inferred from trajectory file 
    // extensions to specify how the trajectory files should be 
    // read.
    if (p_frame_source->bootstrapping_flag == 1) {
		p_frame_source->mt_rand_gen = std::mt19937(p_frame_source->random_num_seed);
	}
	if (p_frame_source->dynamic_state_sampling == 1) p_frame_source->sampleTypesFromProbs();
	
    // Assign a host of function pointers in 'cg' new definitions
    // based on matrix implementation, basis set type, etc.
    set_up_force_computers(p_cg);

    // Initialize the force-matching matrix.
    printf("Initializing FM matrix.\n");
    mscg_struct->mat = new MATRIX_DATA(p_control_input, p_cg);
    if (p_frame_source->use_statistical_reweighting == 1) {
        set_normalization(mscg_struct->mat, 1.0 / p_frame_source->total_frame_weights);
    }
    if (p_frame_source->bootstrapping_flag == 1) {
    	// Multiply the reweighting frame weights by the bootstrapping weights to determine the appropriate
    	// net frame weights and normalizations.
    	if(p_frame_source->use_statistical_reweighting == 1) {
    		combine_reweighting_and_boostrapping_weights(p_frame_source);
    	}
    	set_bootstrapping_normalization(mscg_struct->mat, p_frame_source->bootstrapping_weights, p_frame_source->n_frames);
    }
    
    // Set up the loop index limits for the inner and outer loops.
	if (p_frame_source->dynamic_state_sampling == 1) {
		total_frame_samples = p_frame_source->n_frames * p_frame_source->dynamic_state_samples_per_frame;
	}
	
	printf("Check matrix type, block size, and samples per frame\n"); fflush(stdout);
    if (mscg_struct->mat->matrix_type == kDense) {
        n_blocks = total_frame_samples;
        mscg_struct->mat->frames_per_traj_block = 1;
    } else {
		// Check if number of frames is divisible by frames per trajectory block.
		if (total_frame_samples % mscg_struct->mat->frames_per_traj_block != 0) {
			printf("Total number of frame samples %d is not divisible by block size %d.\n", total_frame_samples, mscg_struct->mat->frames_per_traj_block);
			exit(EXIT_FAILURE);
		}
		n_blocks = total_frame_samples / mscg_struct->mat->frames_per_traj_block;
	}
	(*mscg_struct->mat->set_fm_matrix_to_zero)(mscg_struct->mat);

	//Initialize other data types
	p_frame_source->cleanup = finish_fix_reading;
	p_frame_source->current_frame_n = 1;
    
	mscg_struct->nblocks = n_blocks;
	mscg_struct->curr_frame = 0;
	mscg_struct->mat->accumulation_row_shift = 0;
	mscg_struct->mat->trajectory_block_index = 0;
    mscg_struct->trajectory_block_frame_index = 0;
	mscg_struct->traj_frame_num = 0;

    // Record the dimensions of the matrix after initialization in a
    // solution file.
    FILE* solution_file = fopen("sol_info.out", "w");
    fprintf(solution_file, "fm_matrix_rows:%d; fm_matrix_columns:%d;\n",
            mscg_struct->mat->fm_matrix_rows, mscg_struct->mat->fm_matrix_columns);
    fclose(solution_file);
    
	return(void*)(mscg_struct);
}

// This function continues the setup for range finding
// after the rangefinder_startup_part1, setup_topology_and_frame, and setup_*_topology
// functions have already been called.

// It reads any frame weights from frame_weights.in.
// Also, the sets up the "ForceComputers" and initializes the matrix.
// In addition, it sets up bootstrapping samples, if needed. 

void* rangefinder_startup_part2(void* void_in)
{
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
	int total_frame_samples = mscg_struct->control_input->n_frames;
	int n_blocks = 0;

    // Read the range files rmin.in and rmax.in to determine which
    // interactions should be searched for interaction ranges.
    printf("Reading interaction ranges.\n");
    initialize_range_finding_temps(mscg_struct->cg);

    printf("Allocating dummy force matching matrix temps.\n");
    mscg_struct->control_input->matrix_type = kDummy;
    mscg_struct->mat = new MATRIX_DATA(mscg_struct->control_input, mscg_struct->cg);

    // Read statistical weights for each frame if the 
    // 'use_statistical_reweighting' flag is set in control.in.
    if (mscg_struct->frame_source->use_statistical_reweighting == 1) {
        printf("Reading per-frame statistical reweighting factors.\n");
        fflush(stdout);
        read_frame_weights(mscg_struct->frame_source, mscg_struct->control_input->starting_frame, mscg_struct->control_input->n_frames); 
    }
    
    // Initialize the force-matching matrix.
    printf("Initializing FM matrix.\n");
    mscg_struct->mat = new MATRIX_DATA(mscg_struct->control_input, mscg_struct->cg);
	
	//Initialize other data types
	mscg_struct->frame_source->cleanup = finish_fix_reading;
	mscg_struct->frame_source->current_frame_n = 1;
    
	mscg_struct->nblocks = n_blocks;
	mscg_struct->mat->accumulation_row_shift = 0;
	mscg_struct->mat->trajectory_block_index = 0;
    mscg_struct->trajectory_block_frame_index = 0;
	mscg_struct->traj_frame_num = 0;

	// Check and modify matrix settings.
    if (mscg_struct->mat->matrix_type == kDense) {
        mscg_struct->nblocks = total_frame_samples;
        mscg_struct->mat->frames_per_traj_block = 1;
    } else {
		// Check if number of frames is divisible by frames per trajectory block.
		if (total_frame_samples % mscg_struct->mat->frames_per_traj_block != 0) {
			printf("Total number of frame samples %d is not divisible by block size %d.\n", total_frame_samples, mscg_struct->mat->frames_per_traj_block);
			exit(EXIT_FAILURE);
		}
		mscg_struct->nblocks = total_frame_samples / mscg_struct->mat->frames_per_traj_block;
	}
    mscg_struct->mat->accumulation_row_shift = 0;
	(*mscg_struct->mat->set_fm_matrix_to_zero)(mscg_struct->mat);

	return(void*)(mscg_struct);
}

//----------------------------------------------------------------
// Do the force matching
//----------------------------------------------------------------

// Process each frame of the trajectory to build the MSCG matrix.
// This function should be called after mscg_setup_part2, but before
// mscg_solve_and_output.
void* mscg_process_frame(void* void_in, double* const x, double* const f)
{       		
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
	FrameSource *p_frame_source = mscg_struct->frame_source;
	CG_MODEL_DATA *p_cg = mscg_struct->cg;
    
	// Convert 1D x and f arrays into rvec array
	FrameConfig* p_frame_config = p_frame_source->frame_config;
	for(int i = 0; i < p_frame_config->current_n_sites; i++) {
		for(int j = 0; j < 3; j++) {
			p_frame_config->x[i][j] = x[i*3 + j];
			p_frame_config->f[i][j] = f[i*3 + j];
		}
	}	
	
	// Initialize the cell linked lists for finding neighbors in the provided frames;
    // NVT trajectories are assumed, so this only needs to be done once.
    PairCellList pair_cell_list = PairCellList();
    ThreeBCellList three_body_cell_list = ThreeBCellList();
    pair_cell_list.init(p_cg->pair_nonbonded_interactions.cutoff, p_frame_source);
    if (p_cg->three_body_nonbonded_interactions.class_subtype > 0) {
        double max_cutoff = 0.0;
        for (int i = 0; i < p_cg->three_body_nonbonded_interactions.get_n_defined(); i++) {
            max_cutoff = fmax(max_cutoff, p_cg->three_body_nonbonded_interactions.three_body_nonbonded_cutoffs[i]);
        }
        three_body_cell_list.init(max_cutoff, p_frame_source);
    }
    
    // The trajectory_block_frame_index is incremented for each frame-sample processed.
    // When this index reaches the block size (frames_per_traj_block),
    // The end-of-frame-block routines are called.
    // Then, the trajectory_block_index is incremented.
    mscg_struct->curr_frame++;
 	int traj_frame_num = mscg_struct->traj_frame_num;
	int trajectory_block_frame_index = mscg_struct->trajectory_block_frame_index;
	int times_sampled = 1;
    
    // If reweighting is being used, scale the block of the FM matrix for this frame
    // by the appropriate weighting factor
    if (p_frame_source->use_statistical_reweighting) {
       printf("Reweighting entries for trajectory frame %d. ", traj_frame_num);
       mscg_struct->mat->current_frame_weight = p_frame_source->frame_weights[traj_frame_num];
	}
            
    //Skip processing frame if frame weight is 0.
    if (p_frame_source->use_statistical_reweighting && mscg_struct->mat->current_frame_weight == 0.0) {
    	traj_frame_num++;
    	if (p_frame_source->dynamic_state_sampling != 0) times_sampled = p_frame_source->dynamic_state_samples_per_frame;
    	mscg_struct->trajectory_block_frame_index += times_sampled;
    	
    	if(trajectory_block_frame_index >= mscg_struct->mat->frames_per_traj_block) {
    		// Print status and do end-of-block computations before wiping the blockwise matrix and beginning anew.
        	printf("\r%d (%d) frames have been sampled. ", p_frame_source->current_frame_n, (mscg_struct->mat->trajectory_block_index + 1) * mscg_struct->mat->frames_per_traj_block);
        	fflush(stdout);
        	(*mscg_struct->mat->do_end_of_frameblock_matrix_manipulations)(mscg_struct->mat);
        	(*mscg_struct->mat->set_fm_matrix_to_zero)(mscg_struct->mat);
        	trajectory_block_frame_index=0;
        	mscg_struct->mat->trajectory_block_index++;
        }

    } else {
    
    	// Apply virial constraint, if appropriate.
        add_target_virials_from_trajectory(mscg_struct->mat, mscg_struct->frame_source->pressure_constraint_rhs_vector);

    	// Otherwise, process this frame once unless dynamic_state_sampling is used, in which case it is resampled in this do-while loop.
    	do {    
			if (p_frame_source->dynamic_state_sampling != 0) p_frame_source->sampleTypesFromProbs();
	    	calculate_frame_fm_matrix(p_cg, mscg_struct->mat, p_frame_config, pair_cell_list, three_body_cell_list, trajectory_block_frame_index);
    		times_sampled++;
    		traj_frame_num++;
    		trajectory_block_frame_index++;
    		
     		if(trajectory_block_frame_index >= mscg_struct->mat->frames_per_traj_block) {
    			// Print status and do end-of-block computations before wiping the blockwise matrix and beginning anew.
        		printf("\r%d (%d) frames have been sampled. ", p_frame_source->current_frame_n, (mscg_struct->mat->trajectory_block_index + 1) * mscg_struct->mat->frames_per_traj_block);
        		fflush(stdout);
        		(*mscg_struct->mat->do_end_of_frameblock_matrix_manipulations)(mscg_struct->mat);
        		(*mscg_struct->mat->set_fm_matrix_to_zero)(mscg_struct->mat);
        		trajectory_block_frame_index=0;
        		mscg_struct->mat->trajectory_block_index++;
        	}
		} while ( (p_frame_source->dynamic_state_sampling != 0) && (times_sampled < p_frame_source->dynamic_state_samples_per_frame) );
    }
	p_frame_source->current_frame_n++;
	mscg_struct->traj_frame_num = traj_frame_num;
	mscg_struct->trajectory_block_frame_index = trajectory_block_frame_index;
	
	return (void*)(mscg_struct);
}

// Process each frame of the trajectory to search interaction ranges.
// This function should be called after rangefinder_setup_part2, \
// but before rangefinder_solve_and_output.

void* rangefinder_process_frame(void* void_in, double* const x, double* const f)
{
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
	FrameSource *p_frame_source = mscg_struct->frame_source;
	CG_MODEL_DATA *p_cg = mscg_struct->cg;
    
	// Convert 1D x and f arrays into rvec array
	FrameConfig* p_frame_config = p_frame_source->frame_config;
	for(int i = 0; i < p_frame_config->current_n_sites; i++) {
		for(int j = 0; j < 3; j++) {
			p_frame_config->x[i][j] = x[i*3 + j];
			p_frame_config->f[i][j] = f[i*3 + j];
		}
	}	
	
	// Initialize the cell linked lists for finding neighbors in the provided frames;
    // NVT trajectories are assumed, so this only needs to be done once.
    PairCellList pair_cell_list = PairCellList();
    ThreeBCellList three_body_cell_list = ThreeBCellList();
    pair_cell_list.init(p_cg->pair_nonbonded_interactions.cutoff, p_frame_source);
    if (p_cg->three_body_nonbonded_interactions.class_subtype > 0) {
        double max_cutoff = 0.0;
        for (int i = 0; i < p_cg->three_body_nonbonded_interactions.get_n_defined(); i++) {
            max_cutoff = fmax(max_cutoff, p_cg->three_body_nonbonded_interactions.three_body_nonbonded_cutoffs[i]);
        }
        three_body_cell_list.init(max_cutoff, p_frame_source);
    }
        
    // The trajectory_block_frame_index is incremented for each frame-sample processed.
    // When this index reaches the block size (frames_per_traj_block),
    // The end-of-frame-block routines are called.
    // Then, the trajectory_block_index is incremented.
	int traj_frame_num = mscg_struct->traj_frame_num;
	int trajectory_block_frame_index = mscg_struct->trajectory_block_frame_index;
	int times_sampled = 1;
      
  	// If reweighting is being used, scale the block of the FM matrix for this frame
    // by the appropriate weighting factor
    if (p_frame_source->use_statistical_reweighting) {
       printf("Reweighting entries for trajectory frame %d. ", traj_frame_num);
       mscg_struct->mat->current_frame_weight = p_frame_source->frame_weights[traj_frame_num];
	}
            
    //Skip processing frame if frame weight is 0.
    if (p_frame_source->use_statistical_reweighting && mscg_struct->mat->current_frame_weight == 0.0) {
    	traj_frame_num++;
    	if (p_frame_source->dynamic_state_sampling != 0) times_sampled = p_frame_source->dynamic_state_samples_per_frame;
    	mscg_struct->trajectory_block_frame_index += times_sampled;
    	
    	if(trajectory_block_frame_index >= mscg_struct->mat->frames_per_traj_block) {
    		// Print status and do end-of-block computations before wiping the blockwise matrix and beginning anew.
        	printf("\r%d (%d) frames have been sampled. ", p_frame_source->current_frame_n, (mscg_struct->mat->trajectory_block_index + 1) * mscg_struct->mat->frames_per_traj_block);
        	fflush(stdout);
        	(*mscg_struct->mat->do_end_of_frameblock_matrix_manipulations)(mscg_struct->mat);
        	(*mscg_struct->mat->set_fm_matrix_to_zero)(mscg_struct->mat);
        	trajectory_block_frame_index=0;
        	mscg_struct->mat->trajectory_block_index++;
        }

    } else {
    
    	// Apply virial constraint, if appropriate.
        add_target_virials_from_trajectory(mscg_struct->mat, mscg_struct->frame_source->pressure_constraint_rhs_vector);

    	// Otherwise, process this frame once unless dynamic_state_sampling is used, in which case it is resampled in this do-while loop.
    	do {    
			if (p_frame_source->dynamic_state_sampling != 0) p_frame_source->sampleTypesFromProbs();
	    	calculate_frame_fm_matrix(p_cg, mscg_struct->mat, p_frame_config, pair_cell_list, three_body_cell_list, trajectory_block_frame_index);
    		times_sampled++;
    		traj_frame_num++;
    		trajectory_block_frame_index++;
    		
     		if(trajectory_block_frame_index >= mscg_struct->mat->frames_per_traj_block) {
    			// Print status and do end-of-block computations before wiping the blockwise matrix and beginning anew.
        		printf("\r%d (%d) frames have been sampled. ", p_frame_source->current_frame_n, (mscg_struct->mat->trajectory_block_index + 1) * mscg_struct->mat->frames_per_traj_block);
        		fflush(stdout);
        		(*mscg_struct->mat->do_end_of_frameblock_matrix_manipulations)(mscg_struct->mat);
        		(*mscg_struct->mat->set_fm_matrix_to_zero)(mscg_struct->mat);
        		trajectory_block_frame_index=0;
        		mscg_struct->mat->trajectory_block_index++;
        	}
		} while ( (p_frame_source->dynamic_state_sampling != 0) && (times_sampled < p_frame_source->dynamic_state_samples_per_frame) );
    }
	p_frame_source->current_frame_n++;
	mscg_struct->traj_frame_num = traj_frame_num;
	mscg_struct->trajectory_block_frame_index = trajectory_block_frame_index;
	
	return (void*)(mscg_struct);
}

// Solve the MSCG matrix to generate interactions
// and output those interactions.
// This function should be called last, after all frames
// have been processed using mscg_process_frame.
void* mscg_solve_and_output(void* void_in) 
{
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
	FrameSource *p_frame_source = mscg_struct->frame_source;
    ControlInputs *p_control_input = mscg_struct->control_input;
    MATRIX_DATA *mat = mscg_struct->mat;
    CG_MODEL_DATA *p_cg = mscg_struct->cg;
    
    // Free the space used to build the force-matching matrix that is
    // not necessary for finding a solution to the final matrix
    // equations.

    // Close the trajectory and free the relevant temp variables.
    p_frame_source->cleanup(p_frame_source);

    printf("Finished constructing FM equations.\n");
    if (p_frame_source->bootstrapping_flag == 1) {
		free_bootstrapping_weights(p_frame_source);
	}
	
	// Check that the actual number of frames read matches
	// the number specified in control.in
	if (mscg_struct->curr_frame != p_control_input->n_frames) {
		printf("Warning: The number of frames processed does not match the number of frames specified in the control.in file!\n");
		printf("Please set the number of frames in the control.in file (%d) to match the actual number of frames read (%d).\n", mscg_struct->curr_frame, p_control_input->n_frames);
		fflush(stdout);
		
		// See if this caused some frames not to be converted to normal form.
		int total_frame_samples = mscg_struct->curr_frame;
		if (p_frame_source->dynamic_state_sampling == 1) {
			total_frame_samples *= p_frame_source->dynamic_state_samples_per_frame;
		}
		if (total_frame_samples % mscg_struct->mat->frames_per_traj_block != 0) {
			printf("Warning: Total number of actual frame samples %d is not divisible by block size %d.\n", total_frame_samples, mscg_struct->mat->frames_per_traj_block);
			printf("This can cause some frames to be excluded from the calculation of interactions.\n"); 
			fflush(stdout);
		}

	}
	
    // Find the solution to the force-matching equations set up in
    // previous steps. The solution routines may also print out
    // singular values, residuals, raw matrix equations, etc. as
    // necessary.
    printf("Finishing FM.\n");
    mat->finish_fm(mat);

    // Write tabulated interaction files resulting from the basis set
    // coefficients found in the solution step.
    printf("Writing final output.\n"); fflush(stdout);
    write_fm_interaction_output_files(p_cg, mat);
	
	if (p_frame_source->bootstrapping_flag == 1) {
		delete [] mat->bootstrap_solutions;
    	delete [] mat->bootstrapping_normalization;
    }

    delete mat;
    
    // Record the time and print total elapsed time for profiling purposes.
    double end_cputime = clock();
    double elapsed_cputime = ((double)(end_cputime - mscg_struct->start_cputime)) / CLOCKS_PER_SEC;
    printf("%lf seconds used.\n", elapsed_cputime);
	return (void*)(mscg_struct);
}

// Finalize interaction ranges and generate output.
// This function should be called last, after all frames
// have been processed using rangefinder_process_frame.
void* rangefinder_solve_and_output(void* void_in) 
{
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
    
    // Close the trajectory and free the relevant temp variables.
    mscg_struct->frame_source->cleanup(mscg_struct->frame_source);

     // Free the space used to build the force-matching matrix that is
    // not necessary for the interaction range outputs.
  	printf("Ending range finding.\n");
    printf("Writing final output.\n");
    write_range_files(mscg_struct->cg, mscg_struct->mat);
  	if (mscg_struct->frame_source->bootstrapping_flag == 1) {
		free_bootstrapping_weights(mscg_struct->frame_source);
	}
	
    delete mscg_struct->mat;
    
    // Record the time and print total elapsed time for profiling purposes.
    double end_cputime = clock();
    double elapsed_cputime = ((double)(end_cputime - mscg_struct->start_cputime)) / CLOCKS_PER_SEC;
    printf("%lf seconds used.\n", elapsed_cputime);
	return (void*)(mscg_struct);
}

//-------------------------------------------------------------
// Supporting functions for FrameConfig and FrameSource structs 
// (adapted from trajectory_input.cpp).
//-------------------------------------------------------------
// Initializes the FrameSource struct.
// This function is not intended to be called externally.
// It should only be called by setup_topology_and_frame.
// Otherwise the alias between the two cg_site_types arrays could be disrupted.
void* setup_frame_config(void* void_in, const int n_cg_sites, int * cg_site_types, double* box_half_lengths)
{
	assert(n_cg_sites > 0);	
	
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
	FrameSource *p_frame_source = mscg_struct->frame_source;

	FrameConfig* p_frame_config = new FrameConfig(n_cg_sites);
	p_frame_source->frame_config = p_frame_config;
	
	// Set number of sites types.
	p_frame_source->frame_config->cg_site_types = cg_site_types;
	
	// Set box size.
	for (int i = 0; i < 3; i++) mscg_struct->frame_source->frame_config->simulation_box_half_lengths[i] = (real)(box_half_lengths[i]);
    
    // Check misc. control settings.
    if (mscg_struct->frame_source->dynamic_state_sampling == 1) {
		printf("Error: Dynamic State Sampling is not currently support when using the LAMMPS fix!\n");
		exit(EXIT_FAILURE);
	}
    if ( (mscg_struct->frame_source->dynamic_state_sampling == 1) || (mscg_struct->frame_source->bootstrapping_flag == 1) ) {
        mscg_struct->frame_source->mt_rand_gen = std::mt19937(mscg_struct->frame_source->random_num_seed);
    }
	return (void*)(mscg_struct);
}

// This function updates site types, particle numbers, and box size changes. 
// It can be used if the number of particles changes.
// Use of this function with the rest of the library is unsupported since the code was 
// not designed to handle a change in the number of particles.
// Use at your own risk.
void* update_frame_config(void* void_in, const int n_cg_sites, int * cg_site_types, double* box_half_lengths)
{
	assert(n_cg_sites > 0);	
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
	FrameConfig* p_frame_config = mscg_struct->frame_source->frame_config;
	
	// Resize frame config and the arrays only if the number of sites has changed.
	if(p_frame_config->current_n_sites != n_cg_sites) {
		delete p_frame_config;
		p_frame_config = new FrameConfig(n_cg_sites);
		
		// Also, update other copies of n_cg_sites.
		mscg_struct->cg->topo_data.n_cg_sites = n_cg_sites;
		mscg_struct->cg->n_cg_sites = n_cg_sites;
		
		// Resize the pre-normal matrix for the correct number of rows.
		// Note: This does not affect the normal matrix, which is still
		// the same size (# basis sets = columns by # of basis sets).
		mscg_struct->mat->resize_matrix(n_cg_sites, mscg_struct->mat->fm_matrix_columns);
		// Note: The bond/angle/dihedral/exclusion topology still needs to be updated.
		// Likely, one would call set_*_topology again to update the array pointer.
	}
	
	// Update number of particle types.
	p_frame_config->cg_site_types = cg_site_types;
	
	// Update box size.
	for (int i = 0; i < 3; i++) p_frame_config->simulation_box_half_lengths[i] = box_half_lengths[i];

	return (void*)(mscg_struct);
}

// Clean-up function after all frames are read.
void finish_fix_reading(FrameSource *const frame_source)
{
    // Cleanup allocated memory.
    if (frame_source->use_statistical_reweighting == 1) delete [] frame_source->frame_weights;
    if (frame_source->pressure_constraint_flag == 1) delete [] frame_source->pressure_constraint_rhs_vector;
    delete frame_source->frame_config;
}

//-------------------------------------------------------------
// Supporting functions for TopologyData structs 
// (using topology.h).
//-------------------------------------------------------------

// This function initializes the TopologyData struct.
// This function is intended to be called between 
// part1 and part2 of the mscg_startup functions.
void* setup_topology_and_frame(void* void_in, int const n_cg_sites, int const n_cg_types, char ** type_names, int* cg_site_types, double* box_half_lengths)
{
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
	CG_MODEL_DATA *p_cg = mscg_struct->cg;
    TopologyData* p_topo_data = &(p_cg->topo_data);

	// Copy basic information to CGMD and TopologyData.
	mscg_struct->cg->topo_data.n_cg_sites = n_cg_sites;
	mscg_struct->cg->n_cg_sites = mscg_struct->cg->topo_data.n_cg_sites;
	
	mscg_struct->cg->topo_data.n_cg_types = n_cg_types;	
	mscg_struct->cg->n_cg_types = mscg_struct->cg->topo_data.n_cg_types;

	// Run initializer method for TopologyData.
	// This allocates bond, angle, and dihedral lists and activation_flags
	// as well as the topo_list copy of type names.
	initialize_topology_data(&(mscg_struct->cg->topo_data));
	
	// Allocate space for names and copy names
    mscg_struct->cg->name = new char*[mscg_struct->cg->n_cg_types];
    for (unsigned i = 0; i < unsigned(mscg_struct->cg->n_cg_types); i++) {
        mscg_struct->cg->name[i] = new char[MAX_CG_TYPE_NAME_LENGTH + 1];
        sscanf(type_names[i], "%s", mscg_struct->cg->name[i]);
    }
    
    // Clear the topo_data struct's copy of cg_site_types.
    // Point it to frame_config copy instead.
    delete [] p_topo_data->cg_site_types;
    
    // Process three-body nonbonded interaction information, if appropriate.
    if (p_cg->three_body_nonbonded_interactions.class_subtype > 0) {
    	FILE* top_in = fopen("top.fix", "r");
		char buff[100], parameter_name[50];
		int line = 0;
		unsigned tbtype, i;
		int* tb_i;
		int* tb_j;
		int* tb_k;

		// Check that label matches three body flag.
		fgets(buff, 100, top_in);
		line++;
		sscanf(buff, "%s%d", parameter_name, &tbtype);
		if (strcmp(parameter_name, "threebody") != 0) {
			printf("Error: Expected threebody label instead of %s on line %d of top.fix file.\n", parameter_name, line);
			exit(EXIT_FAILURE);
		}

		// Allocate necessary arrays based on tbtype.
		p_cg->three_body_nonbonded_interactions.stillinger_weber_angle_parameters_by_type = new double[tbtype];
		p_cg->three_body_nonbonded_interactions.three_body_nonbonded_cutoffs = new double[tbtype];

		tb_i = new int[tbtype];
		tb_j = new int[tbtype];
		tb_k = new int[tbtype];
	
		// Read interaction specifications.
		for (i = 0; i < tbtype; i++) {
			fgets(buff, 100, top_in);
			line++;
			sscanf(buff, "%d%d%d%lf%lf", tb_i + i, tb_j + i, tb_k + i, p_cg->three_body_nonbonded_interactions.stillinger_weber_angle_parameters_by_type + i, p_cg->three_body_nonbonded_interactions.three_body_nonbonded_cutoffs + i);
			tb_i[i]--;
			tb_j[i]--;
			tb_k[i]--;
		}
	
		p_cg->three_body_nonbonded_interactions.tb_n = new int[p_topo_data->n_cg_types]();
		p_cg->three_body_nonbonded_interactions.tb_list = new int*[p_topo_data->n_cg_types];
	
		for (i = 0; i < tbtype; i++) {
			p_cg->three_body_nonbonded_interactions.tb_n[tb_i[i]]++;
		}	
	
		for (i = 0; i < p_topo_data->n_cg_types; i++) {
			p_cg->three_body_nonbonded_interactions.tb_list[i] = new int[p_cg->three_body_nonbonded_interactions.tb_n[i] * 2];
			p_cg->three_body_nonbonded_interactions.tb_n[i] = 0;
		}
	
		for (i = 0; i < tbtype; i++) {
			p_cg->three_body_nonbonded_interactions.tb_list[tb_i[i]][p_cg->three_body_nonbonded_interactions.tb_n[tb_i[i]] * 2] = tb_j[i] + 1;
			p_cg->three_body_nonbonded_interactions.tb_list[tb_i[i]][p_cg->three_body_nonbonded_interactions.tb_n[tb_i[i]] * 2 + 1] = tb_k[i] + 1;
			p_cg->three_body_nonbonded_interactions.tb_n[tb_i[i]]++;
		}

		p_cg->three_body_nonbonded_interactions.set_n_defined(tbtype);
		delete [] tb_i;
		delete [] tb_j;
		delete [] tb_k;
	}
	
    void_in = setup_frame_config( (void*)(mscg_struct), n_cg_sites, cg_site_types, box_half_lengths);
	mscg_struct = (MSCG_struct*)(void_in);
	mscg_struct->cg->topo_data.cg_site_types = mscg_struct->frame_source->frame_config->cg_site_types;

	return (void*)(mscg_struct);
}

// This function sets topology information for bonds.
// This function is intended to be called between 
// the setup_topology and the mscg_startup_part2 functions.
void* set_bond_topology(void* void_in, unsigned** bond_partners, unsigned* bond_partner_numbers) 
{
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
    CG_MODEL_DATA *p_cg = mscg_struct->cg;
    TopologyData* p_topo_data = &(p_cg->topo_data);

	// Delete pre-allocated arrays before assigning inputs to these pointers.
	// This should only be done if it is the first time that this function is called.
	if (p_topo_data->bond_list->modified == 0) {
		for (int i = 0; i < p_cg->n_cg_sites; i++) delete [] p_topo_data->bond_list->partners_[i];
		delete [] p_topo_data->bond_list->partners_;
		delete [] p_topo_data->bond_list->partner_numbers_;
	}
	
	// Indicate that these arrays should no longer be freed by the object delete.
	p_topo_data->bond_list->modified = 1;	
	// The number of partners each CG site has for this topological attribute.
	p_topo_data->bond_list->partners_ = bond_partners;
	// For a given CG site, it lists the CG site indices of all partnered particles (for this topological attribute).
	p_topo_data->bond_list->partner_numbers_ = bond_partner_numbers;
	
	// Now, look through bond data to determine activation_flags
	int cg_site1, cg_site2;
	int cg_type1, cg_type2;
	int hash_val;
	for (int i = 0; i < p_cg->n_cg_sites; i++) {
        // Only process this particle if it is part of a bond.
        // Look at each set of angle IDs stored for this particle.
 		for(unsigned j = 0; j < p_topo_data->bond_list->partner_numbers_[i]; j++) {
			// Set the site indices for the sites in this bond.
        	cg_site1 = i;
	        cg_site2 = p_topo_data->bond_list->partners_[i][j];
        
        	// Determine the bond type and activate it (set bond_type_activation_flags).
        	cg_type1 = p_topo_data->cg_site_types[cg_site1];
        	cg_type2 = p_topo_data->cg_site_types[cg_site2];
        	hash_val = calc_two_body_interaction_hash(cg_type1, cg_type2, p_topo_data->n_cg_types);
        	p_topo_data->bond_type_activation_flags[hash_val] = 1;
    	}
    }

	return (void*)(mscg_struct);
}	

// This function sets topology information for angles.
// This function is intended to be called between 
// the setup_topology and the mscg_startup_part2 functions.
void* set_angle_topology(void* void_in, unsigned** angle_partners, unsigned* angle_partner_numbers) 
{
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
    CG_MODEL_DATA *p_cg = mscg_struct->cg;
    TopologyData* p_topo_data = &(p_cg->topo_data);

	// Delete pre-allocated arrays before assigning inputs to these pointers.
	// This should only be done if it is the first time that this function is called.
	if (p_topo_data->angle_list->modified == 0) {
		for (int i = 0; i < p_cg->n_cg_sites; i++) delete [] p_topo_data->angle_list->partners_[i];
		delete [] p_topo_data->angle_list->partners_;
		delete [] p_topo_data->angle_list->partner_numbers_;
	}
	
	// Indicate that these arrays should no longer be freed by the object delete.
	p_topo_data->angle_list->modified = 1;
	// The number of partners each CG site has for this topological attribute.
	p_topo_data->angle_list->partners_ = angle_partners;
	// For a given CG site, it lists the CG site indices of all partnered particles (for this topological attribute).
	p_topo_data->angle_list->partner_numbers_ = angle_partner_numbers;
	
	// Determine the angle type and set angle_type_activation_flags
	int cg_site1, cg_site2, cg_site3;
	int cg_type1, cg_type2, cg_type3;
	int hash_val;
	for (int i = 0; i < p_cg->n_cg_sites; i++) {
        // Only process this particle if it is stores information about an angle
        // Note: Angles only appear for a given angle if the particle i is at the end of the angle. 
 		// Look at each set of angle IDs stored for this particle.
 		for(unsigned j = 0; j < p_topo_data->angle_list->partner_numbers_[i]; j++) {
 			// Grab the particle IDs.
 			// Note: The particle i and the 2nd index in parter are the ends of the angle while the first index is the middle.
 			// However, this code takes the angle as Midle - End - End.
 			// So the transposition of 1 and 2 is intentional
 			cg_site2 = i;
 			cg_site1 = p_topo_data->angle_list->partners_[i][j*2];
 			cg_site3 = p_topo_data->angle_list->partners_[i][j*2 + 1];
 		
 			// Look up the interaction has based on the types in the angle
			cg_type1 = p_topo_data->cg_site_types[cg_site1];
        	cg_type2 = p_topo_data->cg_site_types[cg_site2];
        	cg_type3 = p_topo_data->cg_site_types[cg_site3];
        	hash_val = calc_three_body_interaction_hash(cg_type1, cg_type2, cg_type3, p_topo_data->n_cg_types);
        
        	// Activate this particular type combination angle.
        	p_topo_data->angle_type_activation_flags[hash_val] = 1;
        }
    }
	
	return (void*)(mscg_struct);
}	

// This function sets topology information for dihedrals.
// This function is intended to be called between 
// the setup_topology and the mscg_startup_part2 functions.
void* set_dihedral_topology(void* void_in, unsigned** dihedral_partners, unsigned* dihedral_partner_numbers) 
{
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
    CG_MODEL_DATA *p_cg = mscg_struct->cg;
    TopologyData* p_topo_data = &(p_cg->topo_data);

	// Delete pre-allocated arrays before assigning inputs to these pointers.
	// This should only be done if it is the first time that this function is called.
	if (p_topo_data->dihedral_list->modified == 0) {
		for (int i = 0; i < p_cg->n_cg_sites; i++) delete [] p_topo_data->dihedral_list->partners_[i];
		delete [] p_topo_data->dihedral_list->partners_;
		delete [] p_topo_data->dihedral_list->partner_numbers_;
	}
	
	// Indicate that these arrays should no longer be freed by the object delete.
	p_topo_data->dihedral_list->modified = 1;
	// The number of partners each CG site has for this topological attribute.
	p_topo_data->dihedral_list->partners_ = dihedral_partners;
	p_topo_data->dihedral_list->partner_numbers_ = dihedral_partner_numbers;
	
	// Determine the dihedral type and set dihedral_type_activation_flags
	int cg_site1, cg_site2, cg_site3, cg_site4;
	int cg_type1, cg_type2, cg_type3, cg_type4;
	int hash_val;
	for (int i = 0; i < p_cg->n_cg_sites; i++) {
        // Only process this particle if it is stores information about an angle
        // Note: Dihedrals only appear for a given angle if the particle i is at the end of the dihedral. 
 		// Look at each set of dihedral IDs stored for this particle.
 		for(unsigned j = 0; j < p_topo_data->dihedral_list->partner_numbers_[i]; j++) {
 			// Grab the particle IDs.
 			// Note: The particle i and the 3nd index in parter are the ends of the angle while the first and second indices are in the middle.
 			// However, this code takes the angle as Midle - Middle - End - End.
 			// So the transposition of 1 and 3 is intentional
 			cg_site3 = i;
 			cg_site2 = p_topo_data->angle_list->partners_[i][j*3];
 			cg_site1 = p_topo_data->angle_list->partners_[i][j*3 + 1];
 			cg_site4 = p_topo_data->angle_list->partners_[i][j*3 + 2];
 		
 			// Look up the interaction has based on the types in the angle
			cg_type1 = p_topo_data->cg_site_types[cg_site1];
        	cg_type2 = p_topo_data->cg_site_types[cg_site2];
        	cg_type3 = p_topo_data->cg_site_types[cg_site3];
        	cg_type4 = p_topo_data->cg_site_types[cg_site4];
        	hash_val = calc_four_body_interaction_hash(cg_type1, cg_type2, cg_type3, cg_type4, p_topo_data->n_cg_types);
           
        	// Activate this particular type combination angle.
        	p_topo_data->dihedral_type_activation_flags[hash_val] = 1;
        }
    }
	
	return (void*)(mscg_struct);
}	

// This function sets exclusion information for interactions.
// If the exclusion type should be generated automatically
// or is different in the control_input.in file than the LAMMPS
// setting, generate_exclusion_topology_should be called instead.
// This function is intended to be called between 
// the setup_topology and the mscg_startup_part2 functions.
void* set_exclusion_topology(void* void_in, unsigned** exclusion_partners, unsigned* exclusion_partner_numbers) 
{
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
    CG_MODEL_DATA *p_cg = mscg_struct->cg;
    TopologyData* p_topo_data = &(p_cg->topo_data);

	// Delete pre-allocated arrays before assigning inputs to these pointers.
	// This should only be done if it is the first time that this function is called.
	if (p_topo_data->exclusion_list->modified == 0) {
		for (int i = 0; i < p_cg->n_cg_sites; i++) delete [] p_topo_data->exclusion_list->partners_[i];
		delete [] p_topo_data->exclusion_list->partners_;
		delete [] p_topo_data->exclusion_list->partner_numbers_;
	}
	
	// Indicate that these arrays should no longer be freed by the object delete.
	p_topo_data->exclusion_list->modified = 1;
	// The number of partners each CG site has for this topological attribute.
	p_topo_data->exclusion_list->partners_ = exclusion_partners;
	// For a given CG site, it lists the CG site indices of all partnered particles (for this topological attribute).
	p_topo_data->exclusion_list->partner_numbers_ = exclusion_partner_numbers;
	
	return (void*)(mscg_struct);
}

// This function generates exclusion information for interactions.
// If the exclusion type should be generated automatically
// or is different in the control_input.in file than the LAMMPS
// setting, generate_exclusion_topology_should be called instead.
// This function is intended to be called between 
// the setup_topology (after bond, angle, and/or dihedrals are set)
// and the mscg_startup_part2 functions.
void* generate_exclusion_topology(void* void_in)
{
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
    CG_MODEL_DATA *p_cg = mscg_struct->cg;
    TopologyData* p_topo_data = &(p_cg->topo_data);
	setup_excluded_list(p_topo_data, p_topo_data->exclusion_list, p_topo_data->excluded_style);
	
	return (void*)(mscg_struct);
}

// This function generates angle, dihedral, and exclusion information for interactions.
// The bond topology must already be set before calling this function.
// This function is intended to be called between 
// the setup_topology (setup_bond_topology) and the mscg_startup_part2 functions.
void* generate_angle_dihedral_and_exclusion_topology(void* void_in)
{
	MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
    CG_MODEL_DATA *p_cg = mscg_struct->cg;
    TopologyData* p_topo_data = &(p_cg->topo_data);

	unsigned i, j, k, l;
	int cg_site1, cg_site2;
	int cg_type1, cg_type2, cg_type3;
	int hash_val, n_angles, n_dihedrals;
	int total_angles = 0;
	int total_dihedrals = 0;
	
	// Generate angle topology.
	// Loop over all CG sites.
    for (i = 0; i < p_topo_data->n_cg_sites; i++) {
        n_angles = 0;
        // Loop over sites bonded to that site.
        for (j = 0; j < p_topo_data->bond_list->partner_numbers_[i]; j++) {
                
            cg_site1 = p_topo_data->bond_list->partners_[i][j];

            // Check the sites bonded to that one
            // and consider adding the angle to the angle table.
                for (l = 0; l < p_topo_data->bond_list->partner_numbers_[cg_site1]; l++) {
                    
                    cg_site2 = p_topo_data->bond_list->partners_[cg_site1][l];
                    if (unsigned(cg_site2) == i) continue;
                    
                    p_topo_data->angle_list->partners_[i][2 * n_angles] = cg_site1;
                    p_topo_data->angle_list->partners_[i][2 * n_angles + 1] = cg_site2;
                    n_angles++;
                    
                    cg_type1 = p_topo_data->cg_site_types[cg_site1];
                    cg_type2 = p_topo_data->cg_site_types[i];
                    cg_type3 = p_topo_data->cg_site_types[cg_site2];
                    
                    hash_val = calc_three_body_interaction_hash(cg_type1, cg_type2, cg_type3, p_topo_data->n_cg_types);
                    p_topo_data->angle_type_activation_flags[hash_val] = 1;
            }
            p_topo_data->angle_list->partner_numbers_[i] = n_angles;
        }
    }
	for (i = 0; i < p_topo_data->n_cg_sites; i++) {
	    total_angles += p_topo_data->angle_list->partner_numbers_[i];
	}
    printf("Automatically generated angle topology; %d angles of %d angle types.\n", total_angles/2, calc_n_active_interactions(p_topo_data->angle_type_activation_flags, calc_n_distinct_triples(p_topo_data->n_cg_types)));

	// Generate dihedral topology.
	    // Loop over CG sites in the molecule 
    for (i = 0; i < p_topo_data->n_cg_sites; i++) {
            
        // Infer dihedrals for this site using a combination of the 
        // angle topology and the pair bond topology.
        n_dihedrals = 0;
        for (j = 0; j < p_topo_data->angle_list->partner_numbers_[i]; j++) {
        	for (k = 0; k < p_topo_data->bond_list->partner_numbers_[p_topo_data->angle_list->partners_[i][2 * j + 1]]; k++) {
                if (p_topo_data->bond_list->partners_[p_topo_data->angle_list->partners_[i][2 * j + 1]][k] == p_topo_data->angle_list->partners_[i][2 * j]) continue;
                    if (p_topo_data->bond_list->partners_[p_topo_data->angle_list->partners_[i][2 * j + 1]][k] == i) continue;
                    p_topo_data->dihedral_list->partners_[i][3 * n_dihedrals] = p_topo_data->angle_list->partners_[i][2 * j];
                    p_topo_data->dihedral_list->partners_[i][3 * n_dihedrals + 1] = p_topo_data->angle_list->partners_[i][2 * j + 1];
                    p_topo_data->dihedral_list->partners_[i][3 * n_dihedrals + 2] = p_topo_data->bond_list->partners_[p_topo_data->angle_list->partners_[i][2 * j + 1]][k];
                    hash_val = calc_four_body_interaction_hash(p_topo_data->cg_site_types[p_topo_data->dihedral_list->partners_[i][3 * n_dihedrals]], p_topo_data->cg_site_types[p_topo_data->dihedral_list->partners_[i][3 * n_dihedrals + 1]], p_topo_data->cg_site_types[i], p_topo_data->cg_site_types[p_topo_data->dihedral_list->partners_[i][3 * n_dihedrals + 2]], p_topo_data->n_cg_types);
                    p_topo_data->dihedral_type_activation_flags[hash_val] = 1;
                    n_dihedrals++;
            }
        }
        p_topo_data->dihedral_list->partner_numbers_[i] = n_dihedrals;
    }
	for (i = 0; i < p_topo_data->n_cg_sites; i++) {
	    total_dihedrals += p_topo_data->dihedral_list->partner_numbers_[i];
	}
    printf("Automatically generated dihedral topology; %d dihedrals of %d dihedral types.\n", total_dihedrals/2, calc_n_active_interactions(p_topo_data->dihedral_type_activation_flags, calc_n_distinct_quadruples(p_topo_data->n_cg_types)));

	// Generate exclusion topology.
	setup_excluded_list(p_topo_data, p_topo_data->exclusion_list, p_topo_data->excluded_style);
	
	return (void*)(mscg_struct);
}

// Returns the n_frames parameter
int get_n_frames(void* void_in)
{
  MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
  return mscg_struct->control_input->n_frames;
}

// Returns the block_size parameter
int get_block_size(void* void_in)
{
  MSCG_struct* mscg_struct = (MSCG_struct*)(void_in);
  return mscg_struct->control_input->frames_per_traj_block;
}
