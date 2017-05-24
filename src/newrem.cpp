//
//  newrem.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

// Note: This currently only does reweighting of REF trajectory.

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

void read_previous_rem_solution(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat);
void calculate_new_rem_parameters(MATRIX_DATA* const mat_cg, MATRIX_DATA* const mat_ref);
void construct_full_fm_matrix(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat, FrameSource* const fs);
inline void screen_interaction_basis(CG_MODEL_DATA* const cg);

int main(int argc, char* argv[])
{
    // Begin to compute the total run time
    double start_cputime = clock();

    // Trajectory frame data; see types.h
    FrameSource fs_cg;      //CG trajectory
    FrameSource fs_ref;      //reference trajectory
    ControlInputs control_input;
    
    //----------------------------------------------------------------
    // Set up the entropy minimizing procedure
    //----------------------------------------------------------------
    
    // Parse the command-line arguments of the program; these are 
    // only used to set the trajectory input files and do nothing 
    // else.
    printf("Parsing command line arguments.\n");
    parse_entropy_command_line_arguments(argc,argv, &fs_cg, &fs_ref);
    
    // Read the file control.in to determine the values of many of
    // the parameters in the cg struct, including the desired
    // force-matching strategy (sparse block averaged, normal 
    // equations, accumulation matrices), the desired output styles,
    // the type of basis set to use, and many others described in
    // types.h and listed in control_input.c.
    printf("Reading high level control parameters.\n");
    CG_MODEL_DATA cg(&control_input);   // CG model parameters and data; put here to initialize without default constructor
    copy_control_inputs_to_frd(&control_input, &fs_cg);
    copy_control_inputs_to_frd(&control_input, &fs_ref);

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
	screen_interaction_basis(&cg);
    read_all_interaction_ranges(&cg);

    // Read statistical weights for each frame if the 
    // 'use_statistical_reweighting' flag is set in control.in.
    if (fs_ref.use_statistical_reweighting == 1) {
    	fs_cg.use_statistical_reweighting = 0; // REM reweighting currently is only implemented to act on REF NOT CG.
    	printf("Reading per-frame statistical reweighting factors for reference trajectory.\n");
    	fflush(stdout);
    	read_frame_weights(&fs_ref, control_input.starting_frame, control_input.n_frames); 
    }

    //bootstraping??? Read in of weights would happen here.

    // Use the trajectory type inferred from trajectory file 
    // extensions to specify how the trajectory files should be 
    // read.
    printf("Beginning to read frames.\n");
    printf("In CG trajectory:");
    fs_cg.get_first_frame(&fs_cg, cg.topo_data.n_cg_sites, cg.topo_data.cg_site_types, cg.topo_data.molecule_ids);
    printf("In reference trajectory:");
    fs_ref.get_first_frame(&fs_ref, cg.topo_data.n_cg_sites, cg.topo_data.cg_site_types, cg.topo_data.molecule_ids);
	
	//  dynamic state sampling for REM systems
	if (fs_ref.dynamic_state_sampling == 1) {
		// Sample the types independently
		fs_cg.sampleTypesFromProbs();
		fs_ref.sampleTypesFromProbs();
	}

    // Assign a host of function pointers in 'cg' new definitions
    // based on matrix implementation, basis set type, etc.
    set_up_force_computers(&cg);

    // Initialize the entropy minimizing matrix.
    printf("setting up REM\n");
    control_input.matrix_type = kREM;
    
    MATRIX_DATA mat_cg(&control_input, &cg);
    MATRIX_DATA mat_ref(&control_input, &cg);
    
    if (fs_ref.use_statistical_reweighting == 1) {
        set_normalization(&mat_ref, 1.0 / fs_ref.total_frame_weights);
	}

	//Bootstraping weights would be set here
    
    printf("starting read-in of CG data\n");
   construct_full_fm_matrix(&cg,&mat_cg,&fs_cg);
    
    printf("starting read-in of reference data\n");
    construct_full_fm_matrix(&cg,&mat_ref,&fs_ref);    
    
    //Read in spline coefficents used in the previous iteration.

    printf("Reading in previous iteration's solution\n");
    read_previous_rem_solution(&cg, &mat_cg);

    //Find the solution to the entropy minimization equations set up in
    //previous steps. The solution routine will be moved to matrix after
    //testing has been completed.
    printf("Calculating new REM parameters\n");
    calculate_new_rem_parameters(&mat_cg, &mat_ref);

    // Write tabulated interaction files resulting from the basis set
    // coefficients found in the solution step.
    printf("Writing final output.\n");
    write_fm_interaction_output_files(&cg,&mat_cg);
    
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
    printf("three body non-bonded interactions are not yet supported by newrem\n");
    exit(EXIT_FAILURE);
  }

  printf("entering expectation calculator\n");fflush(stdout);

  (*mat->set_fm_matrix_to_zero)(mat);
  
  for (mat->trajectory_block_index = 0; mat->trajectory_block_index < n_blocks; mat->trajectory_block_index++) {    

    if(fs->pressure_constraint_flag != 0){
      printf("pressure constraints are not supported by newrem\n");
      exit(EXIT_FAILURE);
    }
    
    //for each frame in this block
    for(int trajectory_block_frame_index = 0; trajectory_block_frame_index < mat->frames_per_traj_block; trajectory_block_frame_index++){

      //chck that the last frame read was successful
      if (read_stat == 0){
	    printf("Failure reading frame %d (%d). Check trajectory for errors.\n", fs->current_frame_n, mat->trajectory_block_index * mat->frames_per_traj_block + trajectory_block_frame_index);
	    exit(EXIT_FAILURE);
      }

	  // If reweighting is being used, scale the block of the FM matrix for this frame
	  // by the appropriate weighting factor
	  if (fs->use_statistical_reweighting) {
		  int frame_index = mat->trajectory_block_index * mat->frames_per_traj_block + trajectory_block_frame_index;
		  printf("Reweighting entries for frame %d. ", frame_index);
		  mat->current_frame_weight = fs->frame_weights[frame_index];
	  }

	  //Skip processing frame if frame weight is 0.
      if (fs->use_statistical_reweighting && mat->current_frame_weight == 0.0) {
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

void read_previous_rem_solution(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat)
{
  FILE* spline_input_file = open_file("b-spline-previous.out","r");

  std::string type_name;
  char stringjunk[100];
  int counter = 0;

  std::list<InteractionClassComputer*>::iterator icomp_iterator;
  for(icomp_iterator = cg->icomp_list.begin(); icomp_iterator != cg->icomp_list.end(); icomp_iterator++) {
    // For every defined interaction,
    for (unsigned i = 0; i < (*icomp_iterator)->ispec->defined_to_matched_intrxn_index_map.size(); i++) {
        // If that interaction is being matched,
        if ((*icomp_iterator)->ispec->defined_to_matched_intrxn_index_map[i] != 0) {
	      int index_among_matched = (*icomp_iterator)->ispec->defined_to_matched_intrxn_index_map[i];
          int n_basis_funcs = (*icomp_iterator)->ispec->interaction_column_indices[index_among_matched] - (*icomp_iterator)->ispec->interaction_column_indices[index_among_matched - 1];
	      fgets(stringjunk, 100, spline_input_file);
	      for (int k = 0; k < n_basis_funcs; k++) {
		    fscanf(spline_input_file, "%lf ", &mat->previous_rem_solution[counter]);
		    counter += 1;
	      }
	    fscanf(spline_input_file, "\n");
      }
	}
  }
  fclose(spline_input_file);
}		    	 

void calculate_new_rem_parameters(MATRIX_DATA* const mat_cg, MATRIX_DATA* const mat_ref)
{
  double beta = 1.0 / (mat_cg->temperature * mat_cg->boltzmann);
  double chi = mat_cg->rem_chi;
  const double SMALL = 0.001;

  // Apply normalization to matrices
  for(int i = 0; i < (mat_cg->fm_matrix_columns * mat_cg->fm_matrix_rows); i++)
    {
      mat_ref->dense_fm_normal_matrix->values[i] *= mat_ref->normalization;
      mat_cg->dense_fm_normal_matrix->values[i] *= mat_cg->normalization;
    }
  
  for(int j = 0; j < mat_cg->fm_matrix_columns; j++) 
    {
      mat_ref->previous_rem_solution[j] += mat_ref->dense_fm_normal_matrix->values[j * 2] * beta;
      mat_ref->previous_rem_solution[j] -= mat_cg->dense_fm_normal_matrix->values[j * 2] * beta;
      mat_ref->fm_solution[j] += mat_cg->dense_fm_normal_matrix->values[(j*2) + 1] * beta * beta;
      mat_ref->fm_solution[j] -= mat_cg->dense_fm_normal_matrix->values[j * 2] * mat_cg->dense_fm_normal_matrix->values[j * 2] * beta * beta;
    }

  for(int k = 0; k < mat_cg->fm_matrix_columns; k++) {
      if(mat_ref->fm_solution[k] == 0) {
	      mat_ref->fm_solution[k] = SMALL;	  	
      }
      //This is the gradient decent equation
      //lamda_new = lamda_old - chi * dS/dlamda / Hessian(i,i)
      //lamda_new = mat_cg->fm_solution
      //lamda_old = mat_cg->previous_rem_solution
      //dS/dlamda = mat_ref->previous_rem_solution
      //Hessian   = mat-ref->fm_solution
      mat_cg->fm_solution[k] = mat_cg->previous_rem_solution[k] - mat_ref->previous_rem_solution[k] / mat_ref->fm_solution[k] * chi;

      if((mat_cg->fm_solution[k] - mat_cg->previous_rem_solution[k]) > 100.0) {
	    mat_cg->fm_solution[k] = mat_cg->previous_rem_solution[k] + 100.0;
	  }
      if((mat_cg->fm_solution[k] - mat_cg->previous_rem_solution[k]) < -100.0) {
	    mat_cg->fm_solution[k] = mat_cg->previous_rem_solution[k] - 100.0;
	  }
  }
}

inline void screen_interaction_basis(CG_MODEL_DATA* const cg) {
	std::list<InteractionClassSpec*>::iterator iclass_iterator;
	for(iclass_iterator = cg->iclass_list.begin(); iclass_iterator != cg->iclass_list.end(); iclass_iterator++) {
		if ((*iclass_iterator)->class_type != kOneBody) {
	        (*iclass_iterator)->set_basis_type(kBSplineAndDeriv);
    	}
    }
}