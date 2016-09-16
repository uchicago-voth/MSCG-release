//
//  trajectory_input.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <stdint.h>

#include "control_input.h"
#include "misc.h"
#include "trajectory_input.h"

extern "C" {
#if _exclude_gromacs == 1
#else
#include "xdrfile/xdrfile_xtc.h"
#include "xdrfile/xdrfile_trr.h"
#endif
}

// Local structs only used in trajectory_input.cpp
// These are PIMPLs of FrameSource

//-------------------------------------------------------------
// struct for keeping track of LAMMPS frame data
//-------------------------------------------------------------

struct LammpsData {
	std::ifstream trajectory_stream;
	int type_pos;			// Index for type element in frame body
	int x_pos;				// Starting index for position elements in frame body
	int f_pos;				// Starting index for force elements in frame body
	int state_pos;			// Starting index for state probabilities in frame_body
	int header_size;		// Number of columns for header/body of frame
	std::string* elements; 	// Array to store tokenized header elements 
	double* cg_site_state_probabilities;   // A list of the probabilities for all states of all CG particles (used if dynamic_state_sampling = 1) (currently only for 2 states)
};

//-------------------------------------------------------------
// struct for keeping track of GROMACS frame data
//-------------------------------------------------------------

struct XRDData {
#if _exclude_gromacs == 1
#else
    char extra_trajectory_filename[1000];           // Second trajectory file name (forces for .xtc, not present for .trr)
    XDRFILE* trajectory_filepointer;
    XDRFILE* extra_trajectory_filepointer;
    int read_fr;                                    // Return value holder for the gromacs xtc library functions; see gromacs documentation.
#endif
};

// Prototypes for exclusively internal functions.

// Misc. small helpers.
void report_traj_input_suffix_error(const char *suffix);
void report_usage_error(const char *exe_name);
void check_file_extension(const char* name, const char* suffix);

// Read the initial frame of a trajectory.
void read_initial_trr_frame(FrameSource* const frame_source, const int n_cg_sites, int* cg_site_types);
void read_initial_xtc_frame(FrameSource* const frame_source, const int n_cg_sites,  int* cg_site_types);
void read_initial_lammps_frame(FrameSource* const frame_source, const int n_cg_sites, int* cg_site_types);

// Read a frame of a trajectory after the first has been read.
int read_next_trr_frame(FrameSource* const frame_source);
int read_next_xtc_frame(FrameSource* const frame_source);
int read_next_lammps_frame(FrameSource* const frame_source);

// Read all frames up until a starting frame.
void default_move_to_starting_frame(FrameSource* const frame_source);

// Finish reading a trajectory by closing relevant files and cleaning up temps.
void finish_trr_reading(FrameSource* const frame_source);
void finish_xtc_reading(FrameSource* const frame_source);
void finish_lammps_reading(FrameSource* const frame_source);

// Additional helper functions.
void read_lammps_header(LammpsData* const lammps_data, int* const current_n_sites, int* const timestep, real* const time, matrix box, const int dynamic_types, const int dynamic_state_sampling);
int read_lammps_body(LammpsData* const lammps_data, FrameConfig* const frame_config, const int dynamic_types, const int dynamic_state_sampling);
inline void set_random_number_seed(const uint_fast32_t random_num_seed);

//-------------------------------------------------------------
// Misc. small file-reading helper functions.
//-------------------------------------------------------------

// Report if the command line arguments are not of the correct format.

void report_usage_error(const char *exe_name)
{
    printf("Usage: %s -f file.trr OR %s -f file.xtc -f1 file1.xtc\n", exe_name, exe_name);
    exit(EXIT_SUCCESS);
}

// Report if the trajectories do not have the expected file extension.

void report_traj_input_suffix_error(const char *suffix)
{
    printf("Failed to find a valid file extension (.%s) in the supplied trajectory filename.\n", suffix);
    exit(EXIT_SUCCESS);
}

// Check if the trajectories have the expected file extension.

void check_file_extension(const char* name, const char* suffix)
{
    int len, pos;
    len = strlen(name);
    pos = -1;
    for (int i = 0; i < len; i++) {
        if (name[i] == '.') pos = i;
    }
    if (pos < 0) report_traj_input_suffix_error(suffix);
    if (strcmp(&name[pos + 1], suffix) != 0) report_traj_input_suffix_error(suffix);
}

//-------------------------------------------------------------
// Generic trajectory reading and control functions
//-------------------------------------------------------------

// Parse the command line arguments to determine the trajectory to
// read from and the format of that trajectory.

void parse_command_line_arguments(const int num_arg, char** arg, FrameSource* const frame_source)
{
    if (num_arg != 3 && num_arg != 5) report_usage_error(arg[0]);
    else if (num_arg == 3) {
        if (strcmp(arg[1], "-f") == 0) {
            sscanf(arg[2], "%s", frame_source->trajectory_filename);
            check_file_extension(arg[2], "trr");
            frame_source->trajectory_type = kGromacsTRR;
            frame_source->get_first_frame = read_initial_trr_frame;
            frame_source->get_next_frame = read_next_trr_frame;
            frame_source->cleanup = finish_trr_reading;
            #if _exclude_gromacs == 1
            printf("Cannot read TRR files when _exclude_gromacs is 1. Please recompile without this option and try again.\n");
			fflush(stdout);
			exit(EXIT_FAILURE);
			#endif
        } else if (strcmp(arg[1], "-l") == 0) {
            sscanf(arg[2], "%s", frame_source->trajectory_filename);
            frame_source->trajectory_type = kLAMMPSDump;
            frame_source->get_first_frame = read_initial_lammps_frame;
            frame_source->get_next_frame = read_next_lammps_frame;
            frame_source->cleanup = finish_lammps_reading;
        } else {
            report_usage_error(arg[0]);
        }
    } else if (num_arg == 5) {
        if (strcmp(arg[1], "-f") != 0 || strcmp(arg[3], "-f1") != 0) report_usage_error(arg[0]);
        sscanf(arg[2], "%s", frame_source->trajectory_filename);
        #if _exclude_gromacs == 1
       	printf("Cannot read XTC files when _exclude_gromacs is 1. Please recompile without this option and try again.\n");
		fflush(stdout);
		exit(EXIT_FAILURE);
        #else
        sscanf(arg[4], "%s", frame_source->gromacs_data->extra_trajectory_filename);
        #endif
        check_file_extension(arg[2], "xtc");
        check_file_extension(arg[4], "xtc");

        frame_source->trajectory_type = kGromacsXTC;
        frame_source->get_first_frame = read_initial_xtc_frame;
        frame_source->get_next_frame = read_next_xtc_frame;
        frame_source->cleanup = finish_xtc_reading;
    }
    frame_source->move_to_start_frame = default_move_to_starting_frame;
}

void copy_control_inputs_to_frd(ControlInputs* const control_input, FrameSource* const frame_source)
{
    frame_source->use_statistical_reweighting = control_input->use_statistical_reweighting;
    frame_source->pressure_constraint_flag = control_input->pressure_constraint_flag;
    frame_source->dynamic_types = control_input->dynamic_types;
    frame_source->dynamic_state_sampling = control_input->dynamic_state_sampling;
	frame_source->dynamic_state_samples_per_frame = control_input->dynamic_state_samples_per_frame;
	frame_source->bootstrapping_flag = control_input->bootstrapping_flag;
	frame_source->bootstrapping_num_subsamples = control_input->bootstrapping_num_subsamples;
	frame_source->bootstrapping_num_estimates = control_input->bootstrapping_num_estimates;
    frame_source->random_num_seed = control_input->random_num_seed;
    frame_source->starting_frame = control_input->starting_frame;
    frame_source->n_frames = control_input->n_frames;
}

void finish_trr_reading(FrameSource *const frame_source)
{
 	#if _exclude_gromacs == 1
	#else
    xdrfile_close(frame_source->gromacs_data->trajectory_filepointer);
    delete frame_source->frame_config;
    delete frame_source->gromacs_data;
    if (frame_source->use_statistical_reweighting == 1) delete [] frame_source->frame_weights;
    if (frame_source->pressure_constraint_flag == 1) delete [] frame_source->pressure_constraint_rhs_vector;
    #endif
}

void finish_xtc_reading(FrameSource *const frame_source)
{
 	#if _exclude_gromacs == 1
	#else
	xdrfile_close(frame_source->gromacs_data->trajectory_filepointer);
    xdrfile_close(frame_source->gromacs_data->extra_trajectory_filepointer);
    delete frame_source->frame_config;
    delete frame_source->gromacs_data;
    if (frame_source->use_statistical_reweighting == 1) delete [] frame_source->frame_weights;
    if (frame_source->pressure_constraint_flag == 1) delete [] frame_source->pressure_constraint_rhs_vector;
    #endif
}

void finish_lammps_reading(FrameSource *const frame_source)
{
    //close trajectory file
    frame_source->lammps_data->trajectory_stream.close();
    
    //cleanup allocated memory
    if (frame_source->use_statistical_reweighting == 1) delete [] frame_source->frame_weights;
    if (frame_source->pressure_constraint_flag == 1) delete [] frame_source->pressure_constraint_rhs_vector;
    if ( (frame_source->dynamic_types == 1) || (frame_source->dynamic_state_sampling == 1) ) frame_source->frame_config->cg_site_types = NULL; //undo alias of cg.topo_data.cg_site_types
    if (frame_source->dynamic_state_sampling == 1) delete [] frame_source->lammps_data->cg_site_state_probabilities;
    delete [] frame_source->lammps_data->elements;
	delete frame_source->frame_config;
	delete frame_source->lammps_data;
}


//-------------------------------------------------------------
// Frame-by-frame trajectory reading functions
//-------------------------------------------------------------

// Read the initial frame of a .trr-format trajectory.

void read_initial_trr_frame(FrameSource* const frame_source, const int n_cg_sites, int* cg_site_types)
{
	#if _exclude_gromacs == 1
	#else
    real junk_floating_point;
    int n_sites;
	
 	if (frame_source->dynamic_types == 1) {
    	printf("Use of dynamic_types is only supported for LAMMPS frames (not TRR frames)!\n");
    	exit(EXIT_FAILURE);
    }
    if (frame_source->dynamic_state_sampling == 1) {
    	printf("Use of dynamic_state_sampling is only supported for LAMMPS frames (not TRR frames)!\n");
    	exit(EXIT_FAILURE);
    }
    
    // Get the number of sites in this initial frame and allocate memory to store their forces and positions.
    read_trr_natoms(frame_source->trajectory_filename, &n_sites);
    frame_source->frame_config = new FrameConfig(n_sites);
    frame_source->gromacs_data = new XRDData;

    // Check that the trajectory is consistent with the desired CG model.
    if (n_cg_sites != frame_source->frame_config->current_n_sites) {
        printf("Warning: Number of CG sites defined in top.in is not consistent with trajectory!\n");
    } else {
        printf("Number of CG sites defined in top.in (%d) is consistent with trajectory.\n", n_cg_sites);
    }
    
    // Use Gromacs xtc library routines to read all the data stored in the .trr file for a single frame.
    frame_source->gromacs_data->trajectory_filepointer = xdrfile_open(frame_source->trajectory_filename, "r");
    if (read_trr(frame_source->gromacs_data->trajectory_filepointer, frame_source->frame_config->current_n_sites, &frame_source->current_timestep, &frame_source->time, &junk_floating_point, frame_source->simulation_box_limits, frame_source->frame_config->x, NULL, frame_source->frame_config->f) != exdrOK) {
        printf("Can not read the first frame!\n");
        exit(EXIT_FAILURE);
    }
    frame_source->current_frame_n = 1;

	if (frame_source->bootstrapping_flag == 1) {
		frame_source->mt_rand_gen = std::mt19937(frame_source->random_num_seed);
	}
	
    // Finish up by changing information simply determined by the data just read.
    for (int i = 0; i < 3; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
    return;
    #endif
}

// Read the initial frame of a .xtc-format trajectory.

void read_initial_xtc_frame(FrameSource* const frame_source, const int n_cg_sites, int* cg_site_types)
{
	#if _exclude_gromacs == 1
	#else
    int junk_integer;
    real junk_floating_point;
    int n_atoms, n_sites;
    
    if (frame_source->dynamic_types == 1) {
    	printf("Use of dynamic_types is only supported for LAMMPS frames (not XTC frames)!\n");
    	exit(EXIT_FAILURE);
    }
	if (frame_source->dynamic_state_sampling == 1) {
    	printf("Use of dynamic_state_sampling is only supported for LAMMPS frames (not XTC frames)!\n");
    	exit(EXIT_FAILURE);
    }
    
    // Get the number of sites in this initial frame and allocate memory to store their forces and positions.
    read_xtc_natoms(frame_source->trajectory_filename, &n_sites);
    read_xtc_natoms(frame_source->gromacs_data->extra_trajectory_filename, &n_atoms);
    if (n_sites != n_atoms) {
        printf("Atom numbers are not consistent between two xtc files!\n");
        exit(EXIT_FAILURE);
    }
    frame_source->frame_config = new FrameConfig(n_sites);
    frame_source->gromacs_data = new XRDData;
    
    // Check that the trajectory is consistent with the desired CG model.
    if (n_cg_sites != frame_source->frame_config->current_n_sites) {
        printf("Warning: Number of CG sites defined in top.in is not consistent with trajectory!\n");
    }  else {
        printf("Number of CG sites defined in top.in (%d) is consistent with trajectory.\n", n_cg_sites);

    }
    
    // Use Gromacs xtc library routines to read all the data stored in the .xtc files for a single frame.
    frame_source->gromacs_data->trajectory_filepointer = xdrfile_open(frame_source->trajectory_filename, "r");
    frame_source->gromacs_data->extra_trajectory_filepointer = xdrfile_open(frame_source->gromacs_data->extra_trajectory_filename, "r");
    if ((read_xtc(frame_source->gromacs_data->extra_trajectory_filepointer, frame_source->frame_config->current_n_sites, &frame_source->current_timestep, &frame_source->time, frame_source->simulation_box_limits, frame_source->frame_config->f, &junk_floating_point)
         != exdrOK) ||
        (read_xtc(frame_source->gromacs_data->trajectory_filepointer, frame_source->frame_config->current_n_sites, &junk_integer, &junk_floating_point, frame_source->simulation_box_limits, frame_source->frame_config->x, &junk_floating_point)
         != exdrOK)) {
            printf("Can not read the first frame!\n");
            exit(EXIT_FAILURE);
    }
    frame_source->current_frame_n = 1;

	if (frame_source->bootstrapping_flag == 1) {
		frame_source->mt_rand_gen = std::mt19937(frame_source->random_num_seed);
	}
    
    // Finish up by changing information simply determined by the data just read.
    for (int i = 0; i < 3; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
    return;
    #endif
}

// Read the initial frame of a lammps-dump-format trajectory.

void read_initial_lammps_frame(FrameSource* const frame_source, const int n_cg_sites, int* cg_site_types)
{
	assert(n_cg_sites > 0);	
	frame_source->lammps_data = new LammpsData;
	frame_source->lammps_data->header_size = 0;
	int n_sites = 0;
	
    // Get the number of sites in this initial frame and allocate memory to store their forces and positions.
	frame_source->lammps_data->trajectory_stream.open(frame_source->trajectory_filename, std::ifstream::in);
	if (frame_source->lammps_data->trajectory_stream.fail()) {
		printf("Problem opening lammps trajcetory %s\n", frame_source->trajectory_filename);
		exit(EXIT_FAILURE);
	}
	
	//read header for first frame 
	read_lammps_header(frame_source->lammps_data, &n_sites, &frame_source->current_timestep, &frame_source->time, frame_source->simulation_box_limits, frame_source->dynamic_types, frame_source->dynamic_state_sampling);
	if(n_sites <= 0) {
		exit(EXIT_FAILURE);
    }
    
    //allocate position and force vectors
    frame_source->frame_config = new FrameConfig(n_sites);
    frame_source->lammps_data->elements = new std::string[frame_source->lammps_data->header_size];
    if (frame_source->dynamic_state_sampling == 1) frame_source->lammps_data->cg_site_state_probabilities = new double[n_sites];
	else frame_source->lammps_data->state_pos = -1;
    if ( (frame_source->dynamic_types == 1) || (frame_source->dynamic_state_sampling == 1) ) frame_source->frame_config->cg_site_types = cg_site_types;
    else frame_source->lammps_data->type_pos = -1;
        
    // Check that the trajectory is consistent with the desired CG model.
    if (n_cg_sites != frame_source->frame_config->current_n_sites) {
        printf("Warning: Number of CG sites defined in top.in is not consistent with trajectory!\n");
    } else {
        printf("Number of CG sites defined in top.in (%d) is consistent with trajectory.\n", n_cg_sites);
    }
    
    if ( (frame_source->dynamic_types == 1) && (frame_source->dynamic_state_sampling == 1) ) {
		printf("Warning: Dynamic_state_sampling will override dynamic_types!\n");
	}
    
    //read the body of the frame into memory
    if ( read_lammps_body(frame_source->lammps_data, frame_source->frame_config, frame_source->dynamic_types, frame_source->dynamic_state_sampling) != 1 ) {
    	printf("Cannot read the first frame!\n");			
    	if ( (frame_source->dynamic_types == 1) || (frame_source->dynamic_state_sampling == 1) ) frame_source->frame_config->cg_site_types = NULL; //undo aliasing to cg.topo_data.cg_site_types
		if (frame_source->dynamic_state_sampling == 1) delete [] frame_source->lammps_data->cg_site_state_probabilities;			
		delete [] frame_source->lammps_data->elements;			
		delete frame_source->frame_config;
    	exit(EXIT_FAILURE);
    }
    frame_source->current_frame_n = 1;

   	// Setup random number generator, if appropriate.
 	if ( (frame_source->dynamic_state_sampling == 1) || (frame_source->bootstrapping_flag == 1) ) {
        frame_source->mt_rand_gen = std::mt19937(frame_source->random_num_seed);
    }
    
    // Finish up by changing information simply determined by the data just read.
    for (int i = 0; i < 3; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
    return;
}

// Read a frame of a .trr-format trajectory after the first has been read.

int read_next_trr_frame(FrameSource* const frame_source)
{
	real junk_floating_point;
    int return_val;
    
    #if _exclude_gromacs == 1
	#else
    // Use Gromacs xtc library routines to read all the data stored in the .trr file for a single frame.
    if (read_trr(frame_source->gromacs_data->trajectory_filepointer, frame_source->frame_config->current_n_sites, &frame_source->current_timestep, &frame_source->time, &junk_floating_point, frame_source->simulation_box_limits, frame_source->frame_config->x, NULL, frame_source->frame_config->f) == exdrOK) return_val = 1;
    else return_val = 0;
    
    // Finish up by changing information simply determined by the data just read.
    for (int i = 0; i < 3; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
    frame_source->current_frame_n += 1;
	#endif
    
    // Return any relevant error codes.
    return return_val;
}

// Read a frame of a .xtc-format trajectory after the first has been read

int read_next_xtc_frame(FrameSource* const frame_source)
{
    int junk_integer;
    real junk_floating_point;
    int return_val;
    
    #if _exclude_gromacs == 1
    #else
    // Use Gromacs xtc library routines to read all the data stored in the .xtc files for a single frame.
    if ((read_xtc(frame_source->gromacs_data->extra_trajectory_filepointer, frame_source->frame_config->current_n_sites, &frame_source->current_timestep, &frame_source->time, frame_source->simulation_box_limits, frame_source->frame_config->f, &junk_floating_point)
         == exdrOK) &&
        (read_xtc(frame_source->gromacs_data->trajectory_filepointer, frame_source->frame_config->current_n_sites, &junk_integer, &junk_floating_point, frame_source->simulation_box_limits, frame_source->frame_config->x, &junk_floating_point)
         == exdrOK)) return_val = 1;
    else return_val = 0;
    
    // Finish up by changing information simply determined by the data just read.
    for (int i = 0; i < 3; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
    frame_source->current_frame_n += 1;
	#endif
	
    // Return any relevant error codes.
    return return_val;
}

// Read a frame of a lammps-dump-format trajectory after the first has been read.

int read_next_lammps_frame(FrameSource* const frame_source)
{
	int return_value = 1;  
	int reference_atoms  = frame_source->frame_config->current_n_sites;
	
	read_lammps_header(frame_source->lammps_data, &frame_source->frame_config->current_n_sites, &frame_source->current_timestep, &frame_source->time, frame_source->simulation_box_limits, frame_source->dynamic_types, frame_source->dynamic_state_sampling);    

 	if (reference_atoms != frame_source->frame_config->current_n_sites) {
 		printf("Warning: Number of CG sites defined in top.in is not consistent with trajectory!\n");
 		return_value = 0;
 	} else if ( read_lammps_body(frame_source->lammps_data, frame_source->frame_config, frame_source->dynamic_types, frame_source->dynamic_state_sampling) != 1) {
    	printf("Cannot read the frame at time %lf!\n", frame_source->time);
    	return_value = 0;
    }
 
    // Finish up by changing information simply determined by the data just read.
	for (int i = 0; i < 3; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
    frame_source->current_frame_n += 1;

 	// Return 1 if successful, 0 otherwise.
 	return return_value;
}

void default_move_to_starting_frame(FrameSource* const frame_source) {
    for (int i = 0; i < frame_source->starting_frame - 1; i++) {
        if ((*frame_source->get_next_frame)(frame_source) == 0) {
            printf("Failure attempting to skip frame %d. Check the trajectory file for errors.\n", i);
            exit(EXIT_FAILURE);
        }
    }
}

//-------------------------------------------------------------
// Helper functions for reading LAMMPS header and body
//-------------------------------------------------------------

void read_lammps_header(LammpsData *const lammps_data, int* const current_n_sites, int *const timestep, real *const time, matrix box, const int dynamic_types, const int dynamic_state_sampling)
{
	double low = 0.0;
	double high = 0.0;
	std::string line;
	int flag = 1; 
	
	while(flag == 1) {
		//read next line of header (and wrap-up if end-of-file)
		if(!std::getline(lammps_data->trajectory_stream, line)) {
			printf("\nIt appears that the file is no longer open.\n");
			printf("Please check that you are not attempting to read past the end of the file and try again.\n");
			exit(EXIT_FAILURE);
		}	
	
		//test if it is a labeled line (all LAMMPS labels start with "ITEM:")
		if( line.compare(0, 5, "ITEM:") == 0 ) {
			//find out which label matched (skip space after ITEM:)
			if( line.compare(6, 15, "NUMBER OF ATOMS") == 0) {
				
				//read number of atoms
				lammps_data->trajectory_stream >> *current_n_sites;
				
			} else if( line.compare(6, 10, "BOX BOUNDS") == 0) {
					
				//read in bounds (low high) for 3 dimensions
				for(int pos=0; pos<3; pos++) {
					lammps_data->trajectory_stream >> low >> high;
					box[pos][pos] = high - low;
				}	
				
			} else if( line.compare(6, 8, "TIMESTEP") == 0) {
				
				//read in timestep value
				lammps_data->trajectory_stream >> *time;
				(*timestep)++;
			
			} else if( line.compare(6, 5, "ATOMS") == 0) {
				
				//read labels for body of frame
				flag = 0; 
				size_t prev = 11;
				size_t next = 0;
				int set_x = 0;
				int set_f = 0;
				int set_type = 0;
				int set_state = 0;
				lammps_data->header_size = 0;
				
				//check for xpos and fpos as we tokenize string to determine number of columns in body
				while ((next = line.find_first_of(" ", prev)) != std::string::npos) {
					if( (next - prev) == 0 ) { //check if empty
						prev++;
						continue;
					} else if( (line.compare(prev, 1, "x") == 0 )||
					           (line.compare(prev, 1, "xu") == 0) ) {
						lammps_data->x_pos = lammps_data->header_size;
						set_x = 1;
					} else if( line.compare(prev, 2, "fx") == 0 ) {
						lammps_data->f_pos = lammps_data->header_size;
						set_f = 1;
					} else if( line.compare(prev, 4, "type") == 0 ) {		
						lammps_data->type_pos = lammps_data->header_size;
						set_type = 1;
					} else if( line.compare(prev, 5, "state") == 0 ) {
						lammps_data->state_pos = lammps_data->header_size;
						set_state = 1;
					}
					lammps_data->header_size++;
					prev = next;
				}
		
				if (prev < line.size()) {
                	if( line.compare(prev, 1, "x") == 0 ) {
                        lammps_data->x_pos = lammps_data->header_size;
                    	set_x = 1;
                    } else if( line.compare(prev, 2, "fx") == 0 ) {
                        lammps_data->f_pos = lammps_data->header_size;
                    	set_f = 1;
                    } else if( line.compare(prev, 4, "type") == 0 ) {
                        lammps_data->type_pos = lammps_data->header_size;
                    	set_type = 1;
					} else if( line.compare(prev, 5, "state") == 0 ) {
						lammps_data->state_pos = lammps_data->header_size;
						set_state = 1;
					}
                	lammps_data->header_size++;
                }
				
				//verify that necessary information was extracted to input
				if( (set_x == 0) || (set_f == 0) ) {
					printf("Warning: Was not able to find either x position or f position when parsing LAMMPS frame header!\n");
					exit(EXIT_FAILURE);
				}
				if ( (dynamic_types == 1) && (set_type == 0) ) {
					printf("Warning: Type information not detected in header when parsing LAMMPS frame header!\n");
					exit(EXIT_FAILURE);
				}
				if ( (dynamic_state_sampling == 1) && (set_state == 0) ) {
					printf("Warning: State probability information not detected in header when parsing LAMMPS frame header!\n");
					exit(EXIT_FAILURE);
				}
				
			} else {
				printf("Unrecognized line in frame header: %s", line.c_str() );
				
			}	//close inner if/else if structure
		}		//close ITEM match
	}			//close while	
	return;
}

int read_lammps_body(LammpsData *const lammps_data, FrameConfig *const frame_config, const int dynamic_types, const int dynamic_state_sampling)
{
	//read in current_n_sites lines to extract position and force information
	int j = 0;
	int return_value = 1;
	std::string line;
	
	//allocate space for array of strings based on lammps_data->header_size
	for(int i=0; i < frame_config->current_n_sites; i++)
	{
		//read in next line and tokenize 
		std::getline(lammps_data->trajectory_stream, line, '\n');	
		if(  (j = StringSplit(line, " \t", lammps_data->elements)) != lammps_data->header_size ) {	//allow for trailing white space
			printf("Warning: Number of fields detected in frame body");
			printf(" (%d) does not agree with number expected from frame header (%d)!\n", j, lammps_data->header_size);
			return_value = -1;
			}
			
		//extract position information
		for(j = 0; j < 3; j++) {
			frame_config->x[i][j] = atof( lammps_data->elements[j + lammps_data->x_pos].c_str() );
		}
		
		//extract force information
		for(j = 0; j < 3; j++) {
			frame_config->f[i][j] = atof( lammps_data->elements[j + lammps_data->f_pos].c_str() );
		}

		//extract type information
		if(dynamic_types == 1) { //check if dynamic_type is set
			frame_config->cg_site_types[i] = atoi( lammps_data->elements[lammps_data->type_pos].c_str() );
		}
		if(dynamic_state_sampling == 1) { // check if dynamic_state_sampling is set
			lammps_data->cg_site_state_probabilities[i] = atof( lammps_data->elements[lammps_data->state_pos].c_str() );
		}
	}
	return return_value;
}

void FrameSource::sampleTypesFromProbs()
{
	double rand;
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);
	// Determine each site's type/state by comparing the probability against a random number
	for(int i = 0; i < frame_config->current_n_sites; i++) {
		// Generate random number [0,1] using Mersenne Twister.
		rand = uniform_dist(mt_rand_gen);
		// Make state assignment based on comparison.
		if (rand > lammps_data->cg_site_state_probabilities[i]) frame_config->cg_site_types[i] = 2;
		else frame_config->cg_site_types[i] = 1;
	}
}


//-------------------------------------------------------------
// Whole-trajectory reading functions
//-------------------------------------------------------------

// Read information to reweight all frames from an auxiliary file 'frame_weights.in'.

void read_frame_weights(FrameSource* const frame_source, const int start_frame, const int n_frames)
{
    double junk;
    double total = 0.0;
    int i;
    
    frame_source->frame_weights = new double[n_frames];
    FILE* pp = open_file("frame_weights.in", "r");
    
    for (i = 0; i < start_frame - 1; i++) {
        fscanf(pp, "%lf", &junk);
    }
    for (i = 0; i < n_frames; i++) {
        fscanf(pp, "%lf", frame_source->frame_weights + i);
        total += frame_source->frame_weights[i];
    }
    fclose(pp);
	frame_source->total_frame_weights = 1.0 / total;
}

// Read information relating to the virial constraint for all frames from an auxiliary file 'p_con.in'.
// pressure_constraint_rhs_vector is the RHS of eq. (12) in JCP,123,134105,2005

void read_virial_constraint_vector(FrameSource* const frame_source, const int start_frame, const int n_frames)
{
    double junk;
    int i;
    
    frame_source->pressure_constraint_rhs_vector = new double [n_frames];
    FILE* pp = open_file("p_con.in", "r");
    
    for (i = 0; i < start_frame - 1; i++) {
        fscanf(pp, "%lf", &junk);
    }
    for (i = 0; i < n_frames; i++) {
        fscanf(pp, "%lf", frame_source->pressure_constraint_rhs_vector + i);
    }
    fclose(pp);
}

//-------------------------------------------------------------
// Whole-trajectory generating functions
//-------------------------------------------------------------

void generate_bootstrapping_weights(FrameSource* const frame_source, const int num_frames)
{
	std::uniform_int_distribution<int> uniform_dist(0, num_frames - 1);
	int rand_frame;
	
	if (frame_source->bootstrapping_num_estimates < 1) {
		printf("Cannot request 0 or negative bootstrapping estimates (%d).\n", frame_source->bootstrapping_num_estimates);
		fflush(stdout);
		exit(EXIT_FAILURE);
	}
	if (frame_source->bootstrapping_num_subsamples < 1) {
		printf("Cannot request 0 or negative bootstrapping subsamples per estimate (%d).\n", frame_source->bootstrapping_num_subsamples);
		fflush(stdout);
		exit(EXIT_FAILURE);
	}
	
	// Allocate and initialize bootstrapping_weights
	frame_source->bootstrapping_weights = new double*[frame_source->bootstrapping_num_estimates];
	for (int i = 0; i < frame_source->bootstrapping_num_estimates; i++) {
		frame_source->bootstrapping_weights[i] = new double[num_frames]();
	}
	
	// For each estimate, generate num_subsamples discrete weights.
	// Weights are assigned using random number generator.
	for (int estimate = 0; estimate < frame_source->bootstrapping_num_estimates; estimate++) {
		for (int sample = 0; sample < frame_source->bootstrapping_num_subsamples; sample++) {
			rand_frame = uniform_dist(frame_source->mt_rand_gen);
			frame_source->bootstrapping_weights[estimate][rand_frame] += 1.0;
		}
	}
}

void combine_reweighting_and_boostrapping_weights(FrameSource* const frame_source) {
	for(int i = 0; i < frame_source->bootstrapping_num_estimates; i++) {
		for(int j = 0; j < frame_source->n_frames; j++) {
			frame_source->bootstrapping_weights[i][j] *= frame_source->frame_weights[j];
		}
	}
}

void free_bootstrapping_weights(FrameSource* const frame_source)
{
	for (int i = 0; i < frame_source->bootstrapping_num_estimates; i++) {
		delete [] frame_source->bootstrapping_weights[i];
	}
	delete [] frame_source->bootstrapping_weights;
}
//--------------------------------------------------------------------
// Cell list routines for two- or three-body nonbonded interactions
//--------------------------------------------------------------------

// Initializer for cell lists, using derived class's stencil set up routine.

void BaseCellList::init(const double cutoff, const FrameSource* const fr)
{
    setUpCellListCells(cutoff, fr->frame_config->simulation_box_half_lengths, fr->frame_config->current_n_sites);
    setUpCellListStencil();
}

// Set up the lists and the spatial decomposition.

void BaseCellList::setUpCellListCells(const double cutoff, const real*  simulation_box_half_lengths, const int current_n_sites)
{
    assert(cutoff > 0);
    if (cutoff > simulation_box_half_lengths[0] || cutoff > simulation_box_half_lengths[1] || cutoff > simulation_box_half_lengths[2]) {
        printf("Cutoff is larger than half of the simulation box size!\n");
        exit(EXIT_FAILURE);
    }

	// Determine the number of cells in the box first by calculateng the number of cells needed to span each dimension.
    x_cell_number = (int)(2.0 * simulation_box_half_lengths[0] / cutoff);
    y_cell_number = (int)(2.0 * simulation_box_half_lengths[1] / cutoff);
    z_cell_number = (int)(2.0 * simulation_box_half_lengths[2] / cutoff);
    if (x_cell_number < 3 || y_cell_number < 3 || z_cell_number < 3) {
		// This is a special case where the number of cells is smaller than the stencil in at least one dimension.
		x_cell_number = y_cell_number = z_cell_number = 3;
        x_cell_size = y_cell_size = z_cell_size = -1.0;
    } else {
        x_cell_size = 2.0 * simulation_box_half_lengths[0] / x_cell_number;
        y_cell_size = 2.0 * simulation_box_half_lengths[1] / y_cell_number;
        z_cell_size = 2.0 * simulation_box_half_lengths[2] / z_cell_number;
    }
	// Allocate arrays based on the total number of cells needed to cover the entire simulation box.
    size = x_cell_number * y_cell_number * z_cell_number;
    head = std::vector<int>(size);
    list = std::vector<int>(current_n_sites);
}

// Populate the cell lists.

void BaseCellList::populateList(const int n_particles, const rvec* particle_positions)
{
    assert(n_particles > 0);
    int i, icell;
    int xy_cell_number = x_cell_number * y_cell_number;
	// Calculate the inverse of the size of a cell in each dimension.
    double cellx1 = 1.0 / x_cell_size;
    double celly1 = 1.0 / y_cell_size;
    double cellz1 = 1.0 / z_cell_size;
    
    if (x_cell_size > 0.0) {
		// Assign the particles to cells and build the neighbor list for each cell (backwards).
        // Initialize the cell heads.
		for (i = 0; i < size; i++) {
            head[i] = -1;
        }
        for (i = 0; i < n_particles; i++) {
			// Determine each particles cell. 
			// This "hash" for each cell refers to the cell's index since that data is stored in a flat array (x + y * x_offset + z * x_offset * y_offset).
            icell = (int)(particle_positions[i][0] * cellx1) &
                  + (int)(particle_positions[i][1] * celly1) * x_cell_number &
                  + (int)(particle_positions[i][2] * cellz1) * xy_cell_number;
			// Make this particle the head of its cell's list.
			// This particle now "points" to the index that was previously the head of this cell's list (-1 if it is the first particle).
            list[i] = head[icell];
            head[icell] = i;
        }
    } else {
		// In this special case, it does not make sense to use actual cells.
		// So, this is simply a list of all particles with each particle connected to the particles adjacent to it (by index) in the list.
        for (i = 1; i < size; i++) {
            head[i] = -1;
        }
        head[0] = 0;
        for (i = 0; i < n_particles - 1; i++) {
            list[i] = i + 1;
        }
        list[n_particles - 1] = -1;
    }
}

// Set up a pair list stencil.

void PairCellList::setUpCellListStencil()
{
    int ix, iy, iz, imap;
    int m01 = x_cell_number * y_cell_number;
    stencil = std::vector<int>(size * 13);
    for (iz = 0; iz < z_cell_number; iz++) {
        for (iy = 0; iy < y_cell_number; iy++) {
            for (ix = 0; ix < x_cell_number; ix++) {
                imap = (ix + iy * x_cell_number + iz * m01) * 13;
                stencil[imap] = (ix + 1) % x_cell_number + iy * x_cell_number + iz * m01;
                stencil[imap + 1] = (ix + 1) % x_cell_number + (iy + 1) % y_cell_number * x_cell_number + iz * m01;
                stencil[imap + 2] = (ix + 1) % x_cell_number + iy * x_cell_number + (iz + 1) % z_cell_number * m01;
                stencil[imap + 3] = (ix + 1) % x_cell_number + (iy - 1 + y_cell_number) % y_cell_number * x_cell_number + iz * m01;
                stencil[imap + 4] = (ix + 1) % x_cell_number + iy * x_cell_number + (iz - 1 + z_cell_number) % z_cell_number * m01;
                stencil[imap + 5] = (ix + 1) % x_cell_number + (iy + 1) % y_cell_number * x_cell_number + (iz - 1 + z_cell_number) % z_cell_number * m01;
                stencil[imap + 6] = (ix + 1) % x_cell_number + (iy - 1 + y_cell_number) % y_cell_number * x_cell_number + (iz + 1) % z_cell_number * m01;
                stencil[imap + 7] = (ix + 1) % x_cell_number + (iy - 1 + y_cell_number) % y_cell_number * x_cell_number + (iz - 1 + z_cell_number) % z_cell_number * m01;
                stencil[imap + 8] = (ix + 1) % x_cell_number + (iy + 1) % y_cell_number * x_cell_number + (iz + 1) % z_cell_number * m01;
                stencil[imap + 9] = ix + (iy + 1) % y_cell_number * x_cell_number + (iz + 1) % z_cell_number * m01;
                stencil[imap + 10] = ix + (iy + 1) % y_cell_number * x_cell_number + iz * m01;
                stencil[imap + 11] = ix + (iy + 1) % y_cell_number * x_cell_number + (iz - 1 + z_cell_number) % z_cell_number * m01;
                stencil[imap + 12] = ix + iy * x_cell_number + (iz + 1) % z_cell_number * m01;
            }
        }
    }
}

// Set up a three body list stencil.

void ThreeBCellList::setUpCellListStencil()
{
    int ix, iy, iz, i, j, k, imap;
    int m01 = x_cell_number * y_cell_number;
    stencil = std::vector<int>(size * 26);
    for (iz = 0; iz < z_cell_number; iz++) {
        for (iy = 0; iy < y_cell_number; iy++) {
            for (ix = 0; ix < x_cell_number; ix++) {
                imap = (ix + iy * x_cell_number + iz * m01) * 26;
                for (i = -1; i <= 1; i++) {
                    for (j = -1; j <= 1; j++) {
                        for (k = -1; k <= 1; k++) {
                            if (i != 0 || j != 0 || k != 0) {
                                stencil[imap] = 0;
                                stencil[imap] += (ix + i + x_cell_number) % x_cell_number;
                                stencil[imap] += (iy + j + y_cell_number) % y_cell_number * x_cell_number;
                                stencil[imap] += (iz + k + z_cell_number) % z_cell_number * m01;
                                imap++;
                            }
                        }
                    }
                }   
            }
        }
    }
}

