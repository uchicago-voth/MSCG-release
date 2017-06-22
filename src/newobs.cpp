//
//  newobs.cpp
//
//  The driver implements relative entropy for (scalar) frame-wise observables.
//  It uses least squares fitting.
//
//  Copyright (c) 2017 The Voth Group at The University of Chicago. All rights reserved.
//

#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "control_input.h"
#include "misc.h"
#include "trajectory_input.h"
#include "interaction_model.h"
#include "matrix.h"
#include "force_computation.h"
#include "fm_output.h"
#include "misc.h"

void construct_full_fm_matrix(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FrameSource* const fs);

int main(int argc, char* argv[])
{
    // Begin to compute the total run time
    double start_cputime = clock();

    // Trajectory frame data; see types.h
    FrameSource fs_cg;      // CG trajectory
    ControlInputs control_input;
    
    //----------------------------------------------------------------
    // Set up the entropy minimizing procedure
    //----------------------------------------------------------------
    
    // Parse the command-line arguments of the program; these are 
    // only used to set the trajectory input files and do nothing 
    // else.
    printf("Parsing command line arguments.\n");
    parse_command_line_arguments(argc,argv, &fs_cg);
    
    // Read the file control.in to determine the values of many of
    // the parameters in the cg struct, including the desired
    // force-matching strategy (sparse block averaged, normal 
    // equations, accumulation matrices), the desired output styles,
    // the type of basis set to use, and many others described in
    // types.h and listed in control_input.c.
    printf("Reading high level control parameters.\n");
    CG_MODEL_DATA cg(&control_input);   // CG model parameters and data; put here to initialize without default constructor
    copy_control_inputs_to_frd(&control_input, &fs_cg);
    
    // Read the topology file top.in to determine the definitions of
    // all molecules in the system and their topologies, then to 
    // construct all bond lists, angle lists, and dihedral lists. 
    // The format is based on Gromacs conventions.
    printf("Reading topology file.\n");
    read_topology_file(&cg.topo_data, &cg);

    // Read the range files rmin.in to determine the
    // ranges over which the FM basis functions should be defined.
    // These ranges are also used to record which interactions
    // should be fit, which should be tabulated, and which are not 
    // present in the model.
    printf("Reading interaction ranges.\n");
	screen_interaction_basis(&cg);
    read_all_interaction_ranges(&cg);

    // Read statistical weights for each frame if the 
    // 'use_statistical_reweighting' flag is set in control.in.
    if (fs_cg.use_statistical_reweighting == 1) {
    	printf("Reading per-frame statistical reweighting factors for trajectory.\n");
    	fflush(stdout);
    	read_frame_weights(&fs_cg, control_input.starting_frame, control_input.n_frames, "in"); 
    }

    // Generate bootstrapping weights if the
    // 'bootstrapping_flag' is set in control.in.
	if (fs_cg.bootstrapping_flag == 1) { // REM bootstrapping currently is only implemented to act on CG NOT REF.
    	printf("Generating bootstrapping frame weights.\n");
    	fflush(stdout);
    	generate_bootstrapping_weights(&fs_cg, control_input.n_frames);
    }

    // Use the trajectory type inferred from trajectory file 
    // extensions to specify how the trajectory files should be 
    // read.
    printf("Beginning to read frames.\n");
    fs_cg.get_first_frame(&fs_cg, cg.topo_data.n_cg_sites, cg.topo_data.cg_site_types, cg.topo_data.molecule_ids);
    
	// Perhaps we should have a control.in option that enables reading/use of the
	// observables.cg file (where only observables.ref would be read by default).

    // Read in reference and CG observable values (1 value per frame).
    read_frame_values("observables.ref", control_input.starting_frame, control_input.n_frames, fs_cg.ref_observables);
    if (control_input->cg_observable_flag == 1) {
    	read_frame_values("observables.cg", control_input.starting_frame, control_input.n_frames, fs_cg.cg_observables);
    } else {
    	fs_cg.cg_observables = new double[control_input.n_frames]();
    }

    //  dynamic state sampling for REM systems
	if (fs_cg.dynamic_state_sampling == 1) {
		fs_cg.sampleTypesFromProbs();
	}

    // Assign a host of function pointers in 'cg' new definitions
    // based on matrix implementation, basis set type, etc.
    set_up_force_computers(&cg);

    // Initialize the entropy minimizing matrix.
    printf("setting up RE observable matrix\n");
    control_input.matrix_type = kObs;  
    MATRIX_DATA mat_cg(&control_input, &cg); /*CHECK*/
    
    if (fs_cg.use_statistical_reweighting == 1) {
        set_normalization(&mat_cg, 1.0 / fs_cg.total_frame_weights);
	}

    if (fs_cg.bootstrapping_flag == 1) {
    	// Allocate for bootstrapping only for cg and only if appropriate
  		allocate_bootstrapping(&mat_cg, &control_input, mat_cg.fm_matrix_columns, mat_cg.fm_matrix_columns);

    	// Multiply the reweighting frame weights by the bootstrapping weights to determine the appropriate
    	// net frame weights and normalizations.
    	if(fs_cg.use_statistical_reweighting == 1) {
    		combine_reweighting_and_boostrapping_weights(&fs_cg);
    	}
    	set_bootstrapping_normalization(&mat_cg, fs_cg.bootstrapping_weights, fs_cg.n_frames);
    }
    
    printf("Constructing RE observable equations.\n");
    construct_full_fm_matrix(&cg, &mat_cg, &fs_cg);
        
    // Find the solution to the entropy observable equations set up in the previous steps. 
    // This uses the usual dense FM normal matrix solving since it is actually least squares.
    printf("Finishing RE Observable matching.\n");
    mat_cg.finish_fm(&mat_cg);
    if (fs_cg.bootstrapping_flag == 1) {
		free_bootstrapping_weights(&fs_cg);
	}


    // Write tabulated interaction files resulting from the basis set
    // coefficients found in the solution step.
    printf("Writing final output.\n");
    write_fm_interaction_output_files(&cg, &mat_cg);
    
    // Clean-up special data
    delete [] fs_cg.ref_observables;
    delete [] fs_cg.cg_observables;
    
    // Record the time and print total elapsed time for profiling purposes.
    double end_cputime = clock();
    double elapsed_cputime = ((double)(end_cputime - start_cputime)) / CLOCKS_PER_SEC;
    printf("%f seconds used.\n", elapsed_cputime);
    return 0;
}

void construct_full_fm_matrix(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FrameSource* const fs)
{
  int n_blocks;
  int total_frame_samples = fs->n_frames;
  int traj_frame_num = 0;
  int times_sampled=1;
  int read_stat = 1;

  //skip the desired number of frames before starting the matrix building loop
  fs->move_to_start_frame(fs);

  //Begin the main loop.
  if (fs->dynamic_state_sampling == 1) {
		total_frame_samples = fs->n_frames * fs->dynamic_state_samples_per_frame;
  }

  if (total_frame_samples % mat->frames_per_traj_block != 0) {
    printf("Total number of frames samples %d is not divisible by block size %d.\n",total_frame_samples, mat->frames_per_traj_block);
    exit(EXIT_FAILURE);
  }
  n_blocks = total_frame_samples / mat->frames_per_traj_block;
  
  mat->accumulation_row_shift = 0;

  //initialize the cell linked list
  //NVT assumed so that this only needs to be done once
  PairCellList pair_cell_list = PairCellList();
  ThreeBCellList three_body_cell_list = ThreeBCellList();
  pair_cell_list.init(cg->pair_nonbonded_interactions.cutoff, fs);
  
  if( cg->three_body_nonbonded_interactions.class_subtype > 0){
    printf("Three body non-bonded interactions are not yet supported by newrem\n");
    exit(EXIT_FAILURE);
  }
  
  if (fs->pressure_constraint_flag != 0){
      printf("Pressure constraints are not supported by newobs!\n");
      exit(EXIT_FAILURE);
  }

  printf("Entering expectation calculator.\n");fflush(stdout);

  for (mat->trajectory_block_index = 0; mat->trajectory_block_index < n_blocks; mat->trajectory_block_index++) {
  
	// reset the matrix to 0.
	(*mat->set_fm_matrix_to_zero)(mat);
	
    //for each frame in this block
    for (int trajectory_block_frame_index = 0; trajectory_block_frame_index < mat->frames_per_traj_block; trajectory_block_frame_index++){
	  
	  // check that the last frame read was successful
      if (read_stat == 0){
	    printf("Failure reading frame %d (%d). Check trajectory for errors.\n", fs->current_frame_n, mat->trajectory_block_index * mat->frames_per_traj_block + trajectory_block_frame_index);
	    exit(EXIT_FAILURE);
      }
      
	  // Accumulate the difference in observable values to the RHS vector
	  double observable_difference = fs->ref_observables[traj_frame_num] - fs->cg_observables[traj_frame_num];
	  mat->accumulate_target_force_element(mat, trajectory_block_frame_index, &observable_difference);
	  
	  // If reweighting is being used, scale the block of the FM matrix for this frame
	  // by the appropriate weighting factor
	  if (fs->use_statistical_reweighting == 1) {
		  int frame_index = mat->trajectory_block_index * mat->frames_per_traj_block + trajectory_block_frame_index;
		  printf("Reweighting entries for frame %d. ", frame_index);
		  mat->current_frame_weight = fs->frame_weights[frame_index];
	  }

	  //Skip processing frame if frame weight is 0.
      if (fs->use_statistical_reweighting == 1 && mat->current_frame_weight == 0.0) {
      } else {
      	FrameConfig* frame_config = fs->getFrameConfig();
        calculate_frame_fm_matrix(cg, mat, frame_config, pair_cell_list, three_body_cell_list, trajectory_block_frame_index);
      }
      
      if(fs->dynamic_state_sampling == 0){
	    //Read next frame
	    // Only do this if we are not currently process the last frame.
	    if( ((trajectory_block_frame_index + 1) < mat->frames_per_traj_block) || ((mat->trajectory_block_index + 1) < n_blocks) ){
	      read_stat = (*fs->get_next_frame)(fs);
	    }
	    traj_frame_num++;
	  } else if (times_sampled < fs->dynamic_state_samples_per_frame) {
		// Resample this frame.
		fs->sampleTypesFromProbs();
		times_sampled++;
		
	  } else {
		// Read next frame, sample frame, and reset sampling counter.
		// Only do this if we are not currently process the last frame.
		if ( ((trajectory_block_frame_index + 1) < mat->frames_per_traj_block) ||
			 ((mat->trajectory_block_index + 1) < n_blocks) ) {
			read_stat = (*fs->get_next_frame)(fs);  
		}
		fs->sampleTypesFromProbs();
		times_sampled = 1;
		traj_frame_num++; 
      }
    }
    printf("\r%d (%d) frames have been sampled. ", fs->current_frame_n, (mat->trajectory_block_index + 1) * mat->frames_per_traj_block);
    fflush(stdout);
    (*mat->do_end_of_frameblock_matrix_manipulations)(mat);
  }

  printf("Finishing frame parsing.\n");
  fs->cleanup(fs);
}