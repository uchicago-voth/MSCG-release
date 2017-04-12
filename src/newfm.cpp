//
//  newfm.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "control_input.h"
#include "force_computation.h"
#include "fm_output.h"
#include "interaction_hashing.h"
#include "interaction_model.h"
#include "matrix.h"
#include "misc.h"
#include "trajectory_input.h"

void construct_full_fm_matrix(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FrameSource* const frame_source);

int main(int argc, char* argv[])
{
    // Begin to compute the total run time
    double start_cputime = clock();
    FrameSource frame_source;      // Trajectory frame data; see types.h
    
    //----------------------------------------------------------------
    // Set up the force matching procedure
    //----------------------------------------------------------------
    
    // Parse the command-line arguments of the program; these are 
    // only used to set the trajectory input files and do nothing 
    // else.
    printf("Parsing command line arguments.\n");
    parse_command_line_arguments(argc, argv, &frame_source); 
    
    // Read the file control.in to determine the values of many of
    // the parameters in the cg struct, including the desired
    // force-matching strategy (sparse block averaged, normal 
    // equations, accumulation matrices), the desired output styles,
    // the type of basis set to use, and many others described in
    // types.h and listed in control_input.c.
    printf("Reading high level control parameters.\n");
    ControlInputs control_input; 		// Control parameters read from control.in
	CG_MODEL_DATA cg(&control_input);   // CG model parameters and data (InteractionClasses and Computers)
    copy_control_inputs_to_frd(&control_input, &frame_source);

    // Read the topology file top.in to determine the definitions of
    // all molecules in the system and their topologies, then to 
    // construct all bond lists, angle lists, and dihedral lists. 
    // The format is based on Gromacs conventions.
    printf("Reading topology file.\n");
    read_topology_file(&cg.topo_data, &cg);
    
    // Read the range files rmin.in and rmax.in to determine the
    // ranges over which the FM basis functions should be defined.
    // These ranges are also used to record which interactions
    // should be fit, which should be tabulated, and which are not 
    // present in the model.
    printf("Reading interaction ranges.\n");
    read_all_interaction_ranges(&cg);

    // The range files have specified which interactions should be
    // read from file; read the tabulated potentials from table.in
    // now if any were found.
    if (cg.one_body_interactions.n_tabulated > 0 ||
    	cg.pair_nonbonded_interactions.n_tabulated > 0 ||
        cg.pair_bonded_interactions.n_tabulated > 0 ||
        cg.angular_interactions.n_tabulated > 0 ||
        cg.dihedral_interactions.n_tabulated > 0 ||
        cg.r13_interactions.n_tabulated > 0 ||
        cg.r14_interactions.n_tabulated > 0 ||
        cg.r15_interactions.n_tabulated > 0 ||
        cg.radius_of_gyration_interactions.n_tabulated > 0 ||
		cg.density_interactions.n_tabulated > 0) {
        printf("Reading tabulated reference potentials.\n");
        read_tabulated_interaction_file(&cg, cg.topo_data.n_cg_types);
    } 
    
    // Read statistical weights for each frame if the 
    // 'use_statistical_reweighting' flag is set in control.in.
    if (frame_source.use_statistical_reweighting == 1) {
        printf("Reading per-frame statistical reweighting factors.\n");
        fflush(stdout);
        read_frame_weights(&frame_source, control_input.starting_frame, control_input.n_frames); 
    }
        
    // Generate bootstrapping weights if the
    // 'bootstrapping_flag' is set in control.in.
    if (frame_source.bootstrapping_flag == 1) {
    	printf("Generating bootstrapping frame weights.\n");
    	fflush(stdout);
    	generate_bootstrapping_weights(&frame_source, control_input.n_frames);
    }

    // Read the input virials if the correct flag was set in control.in.
    if (frame_source.pressure_constraint_flag == 1) {
        printf("Reading virial constraint target.\n");
        read_virial_constraint_vector(&frame_source, control_input.starting_frame, control_input.n_frames);
    }
    
    // Use the trajectory type inferred from trajectory file 
    // extensions to specify how the trajectory files should be 
    // read.
    printf("Beginning to read frames.\n");
    frame_source.get_first_frame(&frame_source, cg.topo_data.n_cg_sites, cg.topo_data.cg_site_types, cg.topo_data.molecule_ids);
	if (frame_source.dynamic_state_sampling == 1) frame_source.sampleTypesFromProbs();
	
    // Assign a host of function pointers in 'cg' new definitions
    // based on matrix implementation, basis set type, etc.
    set_up_force_computers(&cg);

    // Initialize the force-matching matrix.
    printf("Initializing FM matrix.\n");
    MATRIX_DATA mat(&control_input, &cg);
    if (frame_source.use_statistical_reweighting == 1) {
        set_normalization(&mat, 1.0 / frame_source.total_frame_weights);
    }
    if (frame_source.bootstrapping_flag == 1) {
    	// Multiply the reweighting frame weights by the bootstrapping weights to determine the appropriate
    	// net frame weights and normalizations.
    	if(frame_source.use_statistical_reweighting == 1) {
    		combine_reweighting_and_boostrapping_weights(&frame_source);
    	}
    	set_bootstrapping_normalization(&mat, frame_source.bootstrapping_weights, frame_source.n_frames);
    }
        
    // Record the dimensions of the matrix after initialization in a
    // solution file.
    FILE* solution_file = open_file("sol_info.out", "w");
    fprintf(solution_file, "fm_matrix_rows:%d; fm_matrix_columns:%d;\n",
            mat.fm_matrix_rows, mat.fm_matrix_columns);
    fclose(solution_file);

    //----------------------------------------------------------------
    // Do the force matching
    //----------------------------------------------------------------

    // Process the whole trajectory to build the force-matching matrix
    // of the appropriate type.
    printf("Constructing FM equations.\n");
    construct_full_fm_matrix(&cg, &mat, &frame_source);

    // Free the space used to build the force-matching matrix that is
    // not necessary for finding a solution to the final matrix
    // equations.
    printf("Finished constructing FM equations.\n");
    if (frame_source.bootstrapping_flag == 1) {
		free_bootstrapping_weights(&frame_source);
	}
	
    // Find the solution to the force-matching equations set up in
    // previous steps. The solution routines may also print out
    // singular values, residuals, raw matrix equations, etc. as
    // necessary.
    printf("Finishing FM.\n");
    mat.finish_fm(&mat);

    // Write tabulated interaction files resulting from the basis set
    // coefficients found in the solution step.
    printf("Writing final output.\n"); fflush(stdout);
    write_fm_interaction_output_files(&cg, &mat);
	
    // Record the time and print total elapsed time for profiling purposes.
    double end_cputime = clock();
    double elapsed_cputime = ((double)(end_cputime - start_cputime)) / CLOCKS_PER_SEC;
    printf("%f seconds used.\n", elapsed_cputime);
    return 0;
}

void construct_full_fm_matrix(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FrameSource* const frame_source)
{
    int n_blocks;
    int read_stat = 1;
    int total_frame_samples = frame_source->n_frames;
	int traj_frame_num = 0;
	int times_sampled = 1;
	double* ref_box_half_lengths = new double[frame_source->position_dimension];
    
    // Skip the desired number of frames before starting the matrix building loops.
    frame_source->move_to_start_frame(frame_source);
    
    // Perform initial generation of cell lists user for generating neighbor lists.
    // This list will only be rebuilt if the box dimensions change.
    
    // Initialize the cell linked lists for finding neighbors in the provided frames;
    PairCellList pair_cell_list = PairCellList();
    ThreeBCellList three_body_cell_list = ThreeBCellList();
    pair_cell_list.init(cg->pair_nonbonded_interactions.cutoff, frame_source);
    if (cg->three_body_nonbonded_interactions.class_subtype > 0) {
    	double max_cutoff = 0.0;
        for (int i = 0; i < cg->three_body_nonbonded_interactions.get_n_defined(); i++) {
        	max_cutoff = fmax(max_cutoff, cg->three_body_nonbonded_interactions.three_body_nonbonded_cutoffs[i]);
        }
    three_body_cell_list.init(max_cutoff, frame_source);
    }
    
	// Record this box's dimensions.
	for (int i = 0; i < frame_source->position_dimension; i++) {
		ref_box_half_lengths[i] = frame_source->frame_config->simulation_box_half_lengths[i];
	}
    
    // Begin the main building loops. This routine operates as a for loop
    // over frame blocks wrapped around a loop over frames within each block.
    // In the inner loop, frames are read every iteration and new matrix elements are computed.
    // In the outer loop, the blockwise matrix is incorporated into the total equations, 
    // then wiped for the process to start again with the next iteration.

    // Set up the loop index limits for the inner and outer loops.
    if (frame_source->dynamic_state_sampling == 1) {
		total_frame_samples = frame_source->n_frames * frame_source->dynamic_state_samples_per_frame;
	}	
    if (mat->matrix_type == kDense) {
        n_blocks = total_frame_samples;
	    mat->frames_per_traj_block = 1;
    } else {
		// Check if number of frames is divisible by frames per trajectory block.
		if (total_frame_samples % mat->frames_per_traj_block != 0) {
			printf("Total number of frame samples %d is not divisible by block size %d.\n", total_frame_samples, mat->frames_per_traj_block);
			exit(EXIT_FAILURE);
		}
		n_blocks = total_frame_samples / mat->frames_per_traj_block;
	}

    mat->accumulation_row_shift = 0;

    // For each block of frame samples.
    printf("Entering primary matrix-building loop.\n"); fflush(stdout);
    for (mat->trajectory_block_index = 0; mat->trajectory_block_index < n_blocks; mat->trajectory_block_index++) {
        
        // Wipe the matrix, then calculate the target virial for all frames in this block.
        (*mat->set_fm_matrix_to_zero)(mat);
        add_target_virials_from_trajectory(mat, frame_source->pressure_constraint_rhs_vector);

        // For each frame sample in this block
        for (int trajectory_block_frame_index = 0; trajectory_block_frame_index < mat->frames_per_traj_block; trajectory_block_frame_index++) {
	
		    // Check that the last frame was read successfully (read at end of each iteration)
    		if (read_stat == 0) {
        		printf("Failure reading frame %d (%d). Check trajectory for errors.\n", frame_source->current_frame_n, mat->trajectory_block_index * mat->frames_per_traj_block + trajectory_block_frame_index);
        		exit(EXIT_FAILURE);
    		}

            // If reweighting is being used, scale the block of the FM matrix for this frame
            // by the appropriate weighting factor
            if (frame_source->use_statistical_reweighting) {
                int frame_index = mat->trajectory_block_index * mat->frames_per_traj_block + trajectory_block_frame_index;
                printf("Reweighting entries for frame %d. ", frame_index);
                mat->current_frame_weight = frame_source->frame_weights[frame_index];
            }
            
            //Skip processing frame if frame weight is 0.
            if (frame_source->use_statistical_reweighting && mat->current_frame_weight == 0.0) {
            } else {
            
            	// Check if the simulation box has changed.
            	int box_change = 0;
            	for (int i = 0; i < frame_source->position_dimension; i++) {
					if ( fabs(ref_box_half_lengths[i] - frame_source->frame_config->simulation_box_half_lengths[i]) > VERYSMALL_F ) {
						box_change = 1;
						break;
					}
				}
				
				// Redo cell list set-up and update reference box size if box has changed.
				if (box_change == 1) {
	            	// Re-initialize the cell linked lists for finding neighbors in the provided frames;
  					pair_cell_list = PairCellList();
    				three_body_cell_list = ThreeBCellList();
    				pair_cell_list.init(cg->pair_nonbonded_interactions.cutoff, frame_source);
    				if (cg->three_body_nonbonded_interactions.class_subtype > 0) {
        				double max_cutoff = 0.0;
        				for (int i = 0; i < cg->three_body_nonbonded_interactions.get_n_defined(); i++) {
            				max_cutoff = fmax(max_cutoff, cg->three_body_nonbonded_interactions.three_body_nonbonded_cutoffs[i]);
        				}
        				three_body_cell_list.init(max_cutoff, frame_source);
    				}
    			
    				// Update the reference_box_half_lengths for this new box size.
    				for (int i = 0; i < frame_source->position_dimension; i++) {
    					ref_box_half_lengths[i] = frame_source->frame_config->simulation_box_half_lengths[i];
    				}
    			}
    			
    			// Modify frame weight if using volume weighting.
    			if (mat->volume_weighting_flag == 1) {
    				real* half_lengths = frame_source->frame_config->simulation_box_half_lengths;
    				double volume = 1.0;
    				for (int i = 0; i < mat->position_dimension; i++) volume *= 2.0 * half_lengths[i];
    				mat->current_frame_weight *= volume * volume;
    			}
				
				// Process frame information.
                FrameConfig* frame_config = frame_source->getFrameConfig();
    			calculate_frame_fm_matrix(cg, mat, frame_config, pair_cell_list, three_body_cell_list, trajectory_block_frame_index);
            }
			
            // Read the next frame; the success of this read will be
            // checked at the start of the next iteration of the loop.
            if (frame_source->dynamic_state_sampling == 0) {
				// Read next frame.
				// Only do this if we are not currently process the last frame.
				if ( ((trajectory_block_frame_index + 1) < mat->frames_per_traj_block) ||
			         ((mat->trajectory_block_index + 1) < n_blocks) ) {
					read_stat = (*frame_source->get_next_frame)(frame_source);  
				}
				traj_frame_num++;
				
			} else if (times_sampled < frame_source->dynamic_state_samples_per_frame) {
				// Resample this frame.
				frame_source->sampleTypesFromProbs();
				times_sampled++;
				
			} else {
				// Read next frame, sample frame, and reset sampling counter.
				// Only do this if we are not currently process the last frame.
				if ( ((trajectory_block_frame_index + 1) < mat->frames_per_traj_block) ||
			         ((mat->trajectory_block_index + 1) < n_blocks) ) {
					read_stat = (*frame_source->get_next_frame)(frame_source);  
				}
				frame_source->sampleTypesFromProbs();
				times_sampled = 1;
				traj_frame_num++;
			}
		}
		
        // Print status and do end-of-block computations before wiping the blockwise matrix and beginning anew
        printf("\r%d (%d) frames have been sampled. ", frame_source->current_frame_n, (mat->trajectory_block_index + 1) * mat->frames_per_traj_block);
        fflush(stdout);
        (*mat->do_end_of_frameblock_matrix_manipulations)(mat);
	}

    printf("\nFinishing frame parsing.\n");
    
    // Close the trajectory and free the relevant temp variables.
    frame_source->cleanup(frame_source);
    delete [] ref_box_half_lengths;
}
