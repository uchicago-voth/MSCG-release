//
//  rangefinder.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "control_input.h"
#include "force_computation.h"
#include "interaction_hashing.h"
#include "interaction_model.h"
#include "matrix.h"
#include "range_finding.h"
#include "misc.h"
#include "trajectory_input.h"
#include "fm_output.h"

void construct_full_fm_matrix(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FrameSource* const frame_source);
inline bool any_active_parameter_distributions(CG_MODEL_DATA* const cg);

int main(int argc, char* argv[])
{
    double start_cputime = clock();

    FrameSource fs;

    printf("Parsing command line arguments.\n");
    parse_command_line_arguments(argc, argv, &fs);
    printf("Reading high level control parameters.\n");
    ControlInputs control_input;
    CG_MODEL_DATA cg(&control_input);   // CG model parameters and data; put here to initialize without default constructor
    copy_control_inputs_to_frd(&control_input, &fs);
    if (control_input.three_body_flag != 0) {
        printf("Rangefinder does not support three body nonbonded interaction ranges.\n");
        exit(EXIT_FAILURE);
    }
    printf("Reading topology file.\n");
    read_topology_file(&cg.topo_data, &cg);

    // Read statistical weights for each frame if the 
    // 'use_statistical_reweighting' flag is set in control.in.
	// For rangefinder, this is only used to exclude frames that 0.0 weight
    if (fs.use_statistical_reweighting == 1) {
        printf("Reading per-frame statistical reweighting factors.\n");
        fflush(stdout);
        read_frame_weights(&fs, control_input.starting_frame, control_input.n_frames); 
    }

    printf("Reading first frame.\n");
    fs.get_first_frame(&fs, cg.n_cg_sites, cg.topo_data.cg_site_types, cg.topo_data.molecule_ids);

    printf("Reading interaction ranges.\n");
    initialize_range_finding_temps(&cg);

    printf("Allocating dummy force matching matrix temps.\n");
    control_input.matrix_type = kDummy;
    MATRIX_DATA mat(&control_input, &cg);

    printf("Beginning range finding.\n");
    construct_full_fm_matrix(&cg, &mat, &fs);
    printf("Ending range finding.\n");
    
    printf("Writing final output.\n");
    write_range_files(&cg, &mat);

    //This is part of the BI routine
    //Only calculate if at least one parameter distribution exists
	if (any_active_parameter_distributions(&cg) == true) {

		read_all_interaction_ranges(&cg);

		set_up_force_computers(&cg);

		printf("hello after force computers\n");fflush(stdout);
	
		calculate_BI(&cg,&mat,&fs);

		printf("hello after calculate BI\n");fflush(stdout);
	
		write_fm_interaction_output_files(&cg,&mat);
	} else {
		// Clean-up allocated memory
		free_name(&cg);
	}
	
    //print cpu time used
    double end_cputime = clock();
    double elapsed_cputime = ((double)(end_cputime - start_cputime)) / CLOCKS_PER_SEC;
    printf("%f seconds used.\n", elapsed_cputime);

    return 0;
}

void construct_full_fm_matrix(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FrameSource* const frame_source)
{
    int n_blocks;
	int total_frame_samples = frame_source->n_frames;
	int traj_frame_num = 0;
	int times_sampled = 1;
    int read_stat = 1;
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
                printf("Reweighting entries for trajectory frame %d. ", traj_frame_num);
                mat->current_frame_weight = frame_source->frame_weights[traj_frame_num];
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

inline bool any_active_parameter_distributions(CG_MODEL_DATA* const cg) {
	std::list<InteractionClassSpec*>::iterator iclass_iterator;
	for(iclass_iterator = cg->iclass_list.begin(); iclass_iterator != cg->iclass_list.end(); iclass_iterator++) {
        if((*iclass_iterator)->output_parameter_distribution == 1) {
        	return true;
        }
    }
    return false;
}
