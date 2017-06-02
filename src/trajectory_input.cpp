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
#include <random>
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
	int mol_pos;			// Index for molecule id in frame body (if molecule_flag == 1)
	int header_size;		// Number of columns for header/body of frame
	std::string* elements; 	// Array to store tokenized header elements 
	double* cg_site_state_probabilities;   // A list of the probabilities for all states of all CG particles (used if dynamic_state_sampling = 1) (currently only for 2 states)
	int (*read_lammps_body)(LammpsData *const lammps_data, FrameConfig *const frame_config, const int dynamic_types, const int molecule_flag, const int dynamic_state_sampling);
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
    rvec* x;
    rvec* f;
    int read_fr;                                    // Return value holder for the gromacs xtc library functions; see gromacs documentation.
	
	inline XRDData(int n_sites) {
		x = new rvec[n_sites + 1];
		f = new rvec[n_sites + 1];
	};
	
	inline ~XRDData() {
		delete [] x;
		delete [] f;
	};
	
	inline void convert_rvec_to_vector(std::array<double, DIMENSION>* &position, std::array<double, DIMENSION>* &force, const int n_cg_sites) {
		for (int i = 0; i < n_cg_sites; i++) {
			for(int j = 0; j < DIMENSION; j++) {
				position[i][j] = x[i][j];
				force[i][j] = f[i][j];
			}
		}
	}
	
#endif
};

// Prototypes for exclusively internal functions.

// Helper for command line to file type setup
void trr_setup(FrameSource* const frame_source, const char* filename);
void lammps_setup(FrameSource* const frame_source, const char* filename);
void xtc_setup(FrameSource* const frame_source, const char* filename1, const char* filename2);

// Misc. small helpers.
inline void report_traj_input_suffix_error(const char *suffix);
inline void report_usage_error(const char *exe_name);
inline void report_rem_usage_error(const char *exe_name);
inline void report_invalid_setting(const char *flag, const char* suffix);
inline void check_molecule_sites(const int n_expected, const int n_read);
inline void check_file_extension(const char* name, const char* suffix);

// Read the initial frame of a trajectory.
void read_initial_trr_frame(FrameSource* const frame_source, const int n_cg_sites, int* cg_site_types, int* mol_ids);
void read_initial_xtc_frame(FrameSource* const frame_source, const int n_cg_sites,  int* cg_site_types, int* mol_ids);
void read_initial_lammps_frame(FrameSource* const frame_source, const int n_cg_sites, int* cg_site_types, int* mol_ids);
void initial_nothing(FrameSource* const frame_source, const int n_cg_sites, int* cg_site_types, int* mol_ids);

// Read a frame of a trajectory after the first has been read.
int read_next_trr_frame(FrameSource* const frame_source);
int read_next_xtc_frame(FrameSource* const frame_source);
int read_next_lammps_frame(FrameSource* const frame_source);
int read_junk_lammps_frame(FrameSource* const frame_source);
int next_nothing(FrameSource* const frame_source);

// Read all frames up until a starting frame.
void default_move_to_starting_frame(FrameSource* const frame_source);

// Read frame-wise entries into an array.
inline void read_stream_into_array(std::ifstream &in_file, const int start_frame, const int n_frames, double* &values);

// Finish reading a trajectory by closing relevant files and cleaning up temps.
inline void finish_general_reading(FrameSource *const frame_source);
void finish_trr_reading(FrameSource* const frame_source);
void finish_xtc_reading(FrameSource* const frame_source);
void finish_lammps_reading(FrameSource* const frame_source);

// Additional helper functions.
void read_lammps_header(LammpsData* const lammps_data, int* const current_n_sites, int* const timestep, real* const time, matrix box, const int dynamic_types, const int molecule_flag, const int dynamic_state_sampling);
int read_dimension_lammps_body(LammpsData* const lammps_data, FrameConfig* const frame_config, const int dynamic_types, const int molecule_flag, const int dynamic_state_sampling);
int read_scalar_lammps_body(LammpsData* const lammps_data, FrameConfig* const frame_config, const int dynamic_types, const int molecule_flag, const int dynamic_state_sampling);
inline void set_random_number_seed(const uint_fast32_t random_num_seed);

//-------------------------------------------------------------
// Misc. small file-reading helper functions.
//-------------------------------------------------------------

// Report if the command line arguments are not of the correct format.

inline void report_usage_error(const char *exe_name)
{
    printf("Usage: %s -f file.trr OR %s -f file.xtc -f1 file1.xtc OR %s -l file.lammpstrj\n", exe_name, exe_name, exe_name);
    exit(EXIT_SUCCESS);
}

inline void report_rem_usage_error(const char *exe_name)
{
    printf("Usage: %s -f_ref file1.trr -f_cg file2.trr OR %s -l_ref file1.lammpstrj -l_cg file2.lammpstrj\n", exe_name, exe_name);
    exit(EXIT_SUCCESS);
}

// Report if the trajectories do not have the expected file extension.

inline void report_traj_input_suffix_error(const char *suffix)
{
    printf("Failed to find a valid file extension (.%s) in the supplied trajectory filename.\n", suffix);
    exit(EXIT_SUCCESS);
}

// Report if a flag setting is compatible with the trajectory type
inline void report_invalid_setting(const char *flag, const char* suffix)
{
	printf("Use of %s is only supported for LAMMPS frames (not %s frames)!\n", flag, suffix);
    exit(EXIT_FAILURE);
}

// Determine output by comparing the expected and read number of cg sites.
inline void check_molecule_sites(const int n_expected, const int n_read) 
{
	if (n_expected != n_read) {
        printf("Warning: Number of CG sites defined in top.in is not consistent with trajectory!\n");
    }  else {
        printf("Number of CG sites defined in top.in (%d) is consistent with trajectory.\n", n_expected);
    }
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
        	trr_setup(frame_source, arg[2]); 
        } else if (strcmp(arg[1], "-l") == 0) {
            lammps_setup(frame_source, arg[2]);
        } else {
            report_usage_error(arg[0]);
        }
    } else if (num_arg == 5) {
        if (strcmp(arg[1], "-f") != 0 || strcmp(arg[3], "-f1") != 0) report_usage_error(arg[0]);
        xtc_setup(frame_source, arg[2], arg[4]);
    }
    frame_source->move_to_start_frame = default_move_to_starting_frame;
}

void parse_entropy_command_line_arguments(const int num_arg, char** arg, FrameSource* const frame_source_cg, FrameSource* const frame_source_ref)
{
  int checker_cg = 0;
  int checker_ref = 0;
  if (num_arg != 5) report_rem_usage_error(arg[0]);
  else {
    parse_command_line_set(arg[1], arg[2], frame_source_cg, frame_source_ref, checker_cg, checker_ref);
    parse_command_line_set(arg[3], arg[4], frame_source_cg, frame_source_ref, checker_cg, checker_ref);
  } 
  if (frame_source_cg->trajectory_filename == frame_source_ref->trajectory_filename) {
	printf("CG and reference trajectory are the same file: %s!", frame_source_cg->trajectory_filename); fflush(stdout);
	exit(EXIT_FAILURE);
      }
  if (checker_cg == 0) {
    printf("CG trajectory is not set. Please set this using -l_cg or -f_cg\n"); fflush(stdout);
    exit(EXIT_FAILURE);
  }
  if (checker_ref == 0) {
    printf("Reference trajectory is not set. Please set this using -l_ref or -f_ref\n"); fflush(stdout);
    exit(EXIT_FAILURE);
  }
    frame_source_cg->move_to_start_frame = default_move_to_starting_frame;
    frame_source_ref->move_to_start_frame = default_move_to_starting_frame;
}

void parse_command_line_set(const char* arg1, const char* arg2, FrameSource* const frame_source_cg, FrameSource* const frame_source_ref, int& checker_cg, int& checker_ref)
{
    if (strcmp(arg1, "-f_cg") == 0) {
      checker_cg = 1;
      trr_setup(frame_source_cg, arg2);
    } else if (strcmp(arg1, "-f_ref") ==0){
      checker_ref = 1;
	  trr_setup(frame_source_ref, arg2);
    } else if (strcmp(arg1, "-l_cg") == 0) {
      checker_cg = 1;
	  lammps_setup(frame_source_cg, arg2);
    } else if (strcmp(arg1, "-l_ref") == 0) {
      checker_ref = 1;
      lammps_setup(frame_source_ref, arg2);
    } else {
      exit(EXIT_FAILURE);
    }
}
 
void trr_setup(FrameSource* const frame_source, const char* filename)
{
	sscanf(filename, "%s", frame_source->trajectory_filename);
	check_file_extension(filename, "trr");
	frame_source->trajectory_type = kGromacsTRR;
	frame_source->get_first_frame = read_initial_trr_frame;
	frame_source->get_next_frame = read_next_trr_frame;
	frame_source->get_junk_frame = read_next_trr_frame;
	frame_source->cleanup = finish_trr_reading;
	#if _exclude_gromacs == 1
	printf("Cannot read TRR files when _exclude_gromacs is 1. Please recompile without this option and try again.\n");
	fflush(stdout);
	exit(EXIT_FAILURE);
	#endif 
}

void lammps_setup(FrameSource* const frame_source, const char* filename)
{
	sscanf(filename, "%s", frame_source->trajectory_filename);
	frame_source->trajectory_type = kLAMMPSDump;
	frame_source->get_first_frame = read_initial_lammps_frame;
	frame_source->get_next_frame = read_next_lammps_frame;
	frame_source->get_junk_frame = read_junk_lammps_frame;
	frame_source->cleanup = finish_lammps_reading;
}

void xtc_setup(FrameSource* const frame_source, const char* filename1, const char* filename2)
{
	sscanf(filename1, "%s", frame_source->trajectory_filename);
	#if _exclude_gromacs == 1
	printf("Cannot read XTC files when _exclude_gromacs is 1. Please recompile without this option and try again.\n");
	fflush(stdout);
	exit(EXIT_FAILURE);
	#else
	sscanf(filename2, "%s", frame_source->gromacs_data->extra_trajectory_filename);
	#endif
	check_file_extension(filename1, "xtc");
	check_file_extension(filename2, "xtc");

	frame_source->trajectory_type = kGromacsXTC;
	frame_source->get_first_frame = read_initial_xtc_frame;
	frame_source->get_next_frame = read_next_xtc_frame;
	frame_source->get_junk_frame = read_next_xtc_frame;
	frame_source->cleanup = finish_xtc_reading;
}

void copy_control_inputs_to_frd(ControlInputs* const control_input, FrameSource* const frame_source)
{
    frame_source->use_statistical_reweighting = control_input->use_statistical_reweighting;
    frame_source->pressure_constraint_flag = control_input->pressure_constraint_flag;
    frame_source->dynamic_types = control_input->dynamic_types;
    frame_source->molecule_flag = control_input->molecule_flag;
    frame_source->dynamic_state_sampling = control_input->dynamic_state_sampling;
	frame_source->dynamic_state_samples_per_frame = control_input->dynamic_state_samples_per_frame;
	frame_source->bootstrapping_flag = control_input->bootstrapping_flag;
	frame_source->bootstrapping_num_subsamples = control_input->bootstrapping_num_subsamples;
	frame_source->bootstrapping_num_estimates = control_input->bootstrapping_num_estimates;
    frame_source->random_num_seed = control_input->random_num_seed;
    frame_source->position_dimension = control_input->position_dimension;
    frame_source->starting_frame = control_input->starting_frame;
    frame_source->n_frames = control_input->n_frames;
    frame_source->scalar_matching_flag = control_input->scalar_matching_flag;
    
    if(frame_source->position_dimension != DIMENSION) {
    	printf("The value of position_dimension(%d) in control_input does not match the compiled dimension(%d)!\n", control_input->position_dimension, DIMENSION);
    	exit(EXIT_FAILURE);
    }
}

void copy_control_inputs_to_frd(ControlInputs * const control_input, FrameSource* fs_ref, FrameSource* fs_cg)
{
	// Do basic calls for each frame source copy
    copy_control_inputs_to_frd(control_input, fs_cg);
    copy_control_inputs_to_frd(control_input, fs_ref);
    
    // The use_statistical_reweighting option acts on the CG frame source and is already copied
    // The reference_statistical_reweighting option acts on the reference frame source 
    fs_ref->use_statistical_reweighting = control_input->reference_statistical_reweighting;
	
	// Bootstrapping is only for CG trajectory (which is already copied)
	// so set the ref bootstrapping flag to 0
	fs_ref->bootstrapping_flag = 0;
	
	// Handle non-trajectory reference input.
	if (control_input->REM_reference_style == 1) {
		fs_ref;
		fs_ref->trajectory_type = kRef;
		fs_ref->get_first_frame = initial_nothing;
		fs_ref->get_next_frame = next_nothing;
		fs_ref->get_junk_frame = next_nothing;
		fs_ref->cleanup = finish_general_reading;
	}
}

inline void finish_general_reading(FrameSource *const frame_source)
{
    delete frame_source->frame_config;
    if (frame_source->use_statistical_reweighting == 1) delete [] frame_source->frame_weights;
    if (frame_source->pressure_constraint_flag == 1) delete [] frame_source->pressure_constraint_rhs_vector;
}

void finish_trr_reading(FrameSource *const frame_source)
{
 	#if _exclude_gromacs == 1
	#else
    xdrfile_close(frame_source->gromacs_data->trajectory_filepointer);
    delete frame_source->gromacs_data;
    finish_general_reading(frame_source);
    #endif
}

void finish_xtc_reading(FrameSource *const frame_source)
{
 	#if _exclude_gromacs == 1
	#else
	xdrfile_close(frame_source->gromacs_data->trajectory_filepointer);
    xdrfile_close(frame_source->gromacs_data->extra_trajectory_filepointer);
    delete frame_source->gromacs_data;
    finish_general_reading(frame_source);
    #endif
}

void finish_lammps_reading(FrameSource *const frame_source)
{
    //close trajectory file
    frame_source->lammps_data->trajectory_stream.close();
    
    //cleanup allocated memory
    if ( (frame_source->dynamic_types == 1) || (frame_source->dynamic_state_sampling == 1) ) frame_source->frame_config->cg_site_types = NULL; //undo alias of cg.topo_data.cg_site_types
    if (frame_source->molecule_flag == 1) {
    	frame_source->frame_config->molecule_ids = NULL; //undo alias of cg.topo_data.molecule_ids
    }
    if (frame_source->dynamic_state_sampling == 1) delete [] frame_source->lammps_data->cg_site_state_probabilities;
    delete [] frame_source->lammps_data->elements;
	delete frame_source->lammps_data;
	
	finish_general_reading(frame_source);
}

//-------------------------------------------------------------
// Frame-by-frame trajectory reading functions
//-------------------------------------------------------------

// Read the initial frame of a .trr-format trajectory.

void read_initial_trr_frame(FrameSource* const frame_source, const int n_cg_sites, int* cg_site_types, int* mol_ids)
{
	#if _exclude_gromacs == 1
	#else
    real junk_floating_point;
    int n_sites;
	
 	if (frame_source->dynamic_types == 1) report_invalid_setting("dynamic_types", "TRR");
    if (frame_source->molecule_flag == 1) report_invalid_setting("molecule_flag", "TRR");
    if (frame_source->dynamic_state_sampling == 1) report_invalid_setting("dynamic_state_sampling", "TRR");
    
    if (frame_source->position_dimension != 3) {
    	printf("TRR frames only support 3 dimensional positions!\n");
    	exit(EXIT_FAILURE);
    };
    
    // Get the number of sites in this initial frame and allocate memory to store their forces and positions.
    read_trr_natoms(frame_source->trajectory_filename, &n_sites);
    frame_source->frame_config = new FrameConfig(n_sites);
    frame_source->gromacs_data = new XRDData(n_sites);

    // Check that the trajectory is consistent with the desired CG model.
    check_molecule_sites(n_cg_sites, frame_source->frame_config->current_n_sites);
    
    // Use Gromacs xtc library routines to read all the data stored in the .trr file for a single frame.
    frame_source->gromacs_data->trajectory_filepointer = xdrfile_open(frame_source->trajectory_filename, "r");
    if (read_trr(frame_source->gromacs_data->trajectory_filepointer, frame_source->frame_config->current_n_sites, &frame_source->current_timestep, &frame_source->time, &junk_floating_point, frame_source->simulation_box_limits, frame_source->gromacs_data->x, NULL, frame_source->gromacs_data->f) != exdrOK) {
        printf("Can not read the first frame!\n");
        exit(EXIT_FAILURE);
    }
    frame_source->current_frame_n = 1;

	if (frame_source->bootstrapping_flag == 1) {
		frame_source->mt_rand_gen = std::mt19937(frame_source->random_num_seed);
	}
	
    // Finish up by changing information simply determined by the data just read.
    frame_source->gromacs_data->convert_rvec_to_vector(frame_source->frame_config->x, frame_source->frame_config->f, frame_source->frame_config->current_n_sites);
    for (int i = 0; i < DIMENSION; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
    return;
    #endif
}

// Read the initial frame of a .xtc-format trajectory.

void read_initial_xtc_frame(FrameSource* const frame_source, const int n_cg_sites, int* cg_site_types, int* mol_ids)
{
	#if _exclude_gromacs == 1
	#else
    int junk_integer;
    real junk_floating_point;
    int n_atoms, n_sites;
    
    if (frame_source->dynamic_types == 1) report_invalid_setting("dynamic_types", "XTC");
    if (frame_source->molecule_flag == 1) report_invalid_setting("molecule_flag", "XTC");
	if (frame_source->dynamic_state_sampling == 1) report_invalid_setting("dynamic_state_sampling", "XTC"); 
    
    if (frame_source->position_dimension != 3) {
    	printf("XTC frames only support 3 dimensional positions!\n");
    	exit(EXIT_FAILURE);
    };
    
    
    // Get the number of sites in this initial frame and allocate memory to store their forces and positions.
    read_xtc_natoms(frame_source->trajectory_filename, &n_sites);
    read_xtc_natoms(frame_source->gromacs_data->extra_trajectory_filename, &n_atoms);
    if (n_sites != n_atoms) {
        printf("Atom numbers are not consistent between two xtc files!\n");
        exit(EXIT_FAILURE);
    }
    frame_source->frame_config = new FrameConfig(n_sites);
    frame_source->gromacs_data = new XRDData(n_sites);
    
    // Check that the trajectory is consistent with the desired CG model.
    check_molecule_sites(n_cg_sites, frame_source->frame_config->current_n_sites);
    
    // Use Gromacs xtc library routines to read all the data stored in the .xtc files for a single frame.
    frame_source->gromacs_data->trajectory_filepointer = xdrfile_open(frame_source->trajectory_filename, "r");
    frame_source->gromacs_data->extra_trajectory_filepointer = xdrfile_open(frame_source->gromacs_data->extra_trajectory_filename, "r");
    if ((read_xtc(frame_source->gromacs_data->extra_trajectory_filepointer, frame_source->frame_config->current_n_sites, &frame_source->current_timestep, &frame_source->time, frame_source->simulation_box_limits, frame_source->gromacs_data->f, &junk_floating_point)
         != exdrOK) ||
        (read_xtc(frame_source->gromacs_data->trajectory_filepointer, frame_source->frame_config->current_n_sites, &junk_integer, &junk_floating_point, frame_source->simulation_box_limits, frame_source->gromacs_data->x, &junk_floating_point)
         != exdrOK)) {
            printf("Can not read the first frame!\n");
            exit(EXIT_FAILURE);
    }
    frame_source->current_frame_n = 1;

	if (frame_source->bootstrapping_flag == 1) {
		frame_source->mt_rand_gen = std::mt19937(frame_source->random_num_seed);
	}
    
    // Finish up by changing information simply determined by the data just read.
    frame_source->gromacs_data->convert_rvec_to_vector(frame_source->frame_config->x, frame_source->frame_config->f, frame_source->frame_config->current_n_sites);
    for (int i = 0; i < DIMENSION; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
    return;
    #endif
}

// Read the initial frame of a lammps-dump-format trajectory.

void read_initial_lammps_frame(FrameSource* const frame_source, const int n_cg_sites, int* cg_site_types, int* mol_ids)
{
	assert(n_cg_sites > 0);	
	frame_source->lammps_data = new LammpsData;
	frame_source->lammps_data->header_size = 0;
	int n_sites = 0;

	// Set read_lammps_body function_pointer based on scalar_matching_flag.
	if (frame_source->scalar_matching_flag == 1) {
		frame_source->lammps_data->read_lammps_body = read_scalar_lammps_body;
	} else {
		frame_source->lammps_data->read_lammps_body = read_dimension_lammps_body;
	}
	
    // Get the number of sites in this initial frame and allocate memory to store their forces and positions.
	frame_source->lammps_data->trajectory_stream.open(frame_source->trajectory_filename, std::ifstream::in);
	if (frame_source->lammps_data->trajectory_stream.fail()) {
		printf("Problem opening lammps trajcetory %s\n", frame_source->trajectory_filename);
		exit(EXIT_FAILURE);
	}
	
	//read header for first frame 
	read_lammps_header(frame_source->lammps_data, &n_sites, &frame_source->current_timestep, &frame_source->time, frame_source->simulation_box_limits, frame_source->dynamic_types, frame_source->molecule_flag, frame_source->dynamic_state_sampling);
	if(n_sites <= 0) {
		exit(EXIT_FAILURE);
    }
    
    //allocate position and force vectors
    frame_source->frame_config = new FrameConfig(n_sites);
    frame_source->lammps_data->elements = new std::string[frame_source->lammps_data->header_size];
    if (frame_source->dynamic_state_sampling == 1) frame_source->lammps_data->cg_site_state_probabilities = new double[n_sites];
	else frame_source->lammps_data->state_pos = -1;
    if ( (frame_source->dynamic_types == 1) || (frame_source->dynamic_state_sampling == 1) ) {
    	frame_source->frame_config->cg_site_types = cg_site_types;
    } else {
    	frame_source->lammps_data->type_pos = -1;
    }
    if (frame_source->molecule_flag == 1) {
    	frame_source->frame_config->molecule_ids = mol_ids;
    }
        
    // Check that the trajectory is consistent with the desired CG model.
    check_molecule_sites(n_cg_sites, frame_source->frame_config->current_n_sites);
    
    if ( (frame_source->dynamic_types == 1) && (frame_source->dynamic_state_sampling == 1) ) {
		printf("Warning: Dynamic_state_sampling will override dynamic_types!\n");
	}
    
    //read the body of the frame into memory
    if ( frame_source->lammps_data->read_lammps_body(frame_source->lammps_data, frame_source->frame_config, frame_source->dynamic_types, frame_source->molecule_flag, frame_source->dynamic_state_sampling) != 1 ) {
    	printf("Cannot read the first frame!\n");			
    	if ( (frame_source->dynamic_types == 1) || (frame_source->dynamic_state_sampling == 1) ) frame_source->frame_config->cg_site_types = NULL; //undo aliasing to cg.topo_data.cg_site_types
		if (frame_source->molecule_flag == 1) {
			frame_source->frame_config->molecule_ids = NULL;
		}
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
    for (int i = 0; i < DIMENSION; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
    return;
}

void initial_nothing(FrameSource* const frame_source, const int n_cg_sites, int* cg_site_types, int* mol_ids)
{
}

// Read a frame of a .trr-format trajectory after the first has been read.

int read_next_trr_frame(FrameSource* const frame_source)
{
    int return_val = 0;
    
    #if _exclude_gromacs == 1
	#else
	real junk_floating_point;
    // Use Gromacs xtc library routines to read all the data stored in the .trr file for a single frame.
    if (read_trr(frame_source->gromacs_data->trajectory_filepointer, frame_source->frame_config->current_n_sites, &frame_source->current_timestep, &frame_source->time, &junk_floating_point, frame_source->simulation_box_limits, frame_source->gromacs_data->x, NULL, frame_source->gromacs_data->f) == exdrOK) return_val = 1;
    else return_val = 0;
    
    // Finish up by changing information simply determined by the data just read.
    frame_source->gromacs_data->convert_rvec_to_vector(frame_source->frame_config->x, frame_source->frame_config->f, frame_source->frame_config->current_n_sites);
    for (int i = 0; i < DIMENSION; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
    frame_source->current_frame_n += 1;
	#endif
    
    // Return any relevant error codes.
    return return_val;
}

// Read a frame of a .xtc-format trajectory after the first has been read

int read_next_xtc_frame(FrameSource* const frame_source)
{
    int return_val = 0;
    
    #if _exclude_gromacs == 1
    #else
    int junk_integer;
    real junk_floating_point;
    // Use Gromacs xtc library routines to read all the data stored in the .xtc files for a single frame.
    if ((read_xtc(frame_source->gromacs_data->extra_trajectory_filepointer, frame_source->frame_config->current_n_sites, &frame_source->current_timestep, &frame_source->time, frame_source->simulation_box_limits, frame_source->gromacs_data->f, &junk_floating_point)
         == exdrOK) &&
        (read_xtc(frame_source->gromacs_data->trajectory_filepointer, frame_source->frame_config->current_n_sites, &junk_integer, &junk_floating_point, frame_source->simulation_box_limits, frame_source->gromacs_data->x, &junk_floating_point)
         == exdrOK)) return_val = 1;
    else return_val = 0;
    
    // Finish up by changing information simply determined by the data just read.
    frame_source->gromacs_data->convert_rvec_to_vector(frame_source->frame_config->x, frame_source->frame_config->f, frame_source->frame_config->current_n_sites);
    for (int i = 0; i < DIMENSION; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
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

	read_lammps_header(frame_source->lammps_data, &frame_source->frame_config->current_n_sites, &frame_source->current_timestep, &frame_source->time, frame_source->simulation_box_limits, frame_source->dynamic_types, frame_source->molecule_flag, frame_source->dynamic_state_sampling);    

 	if (reference_atoms != frame_source->frame_config->current_n_sites) {
 		printf("Warning: Number of CG sites defined in top.in is not consistent with trajectory!\n");
 		return_value = 0;
 	} else if ( frame_source->lammps_data->read_lammps_body(frame_source->lammps_data, frame_source->frame_config, frame_source->dynamic_types, frame_source->molecule_flag, frame_source->dynamic_state_sampling) != 1) {
    	printf("Cannot read the frame at time %lf!\n", frame_source->time);
    	return_value = 0;
    }
 
    // Finish up by changing information simply determined by the data just read.
	for (int i = 0; i < DIMENSION; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
    frame_source->current_frame_n += 1;

 	// Return 1 if successful, 0 otherwise.
 	return return_value;
}

int read_junk_lammps_frame(FrameSource* const frame_source)
{
	int return_value = 1;  
	int reference_atoms  = frame_source->frame_config->current_n_sites;
	std::string line;

	read_lammps_header(frame_source->lammps_data, &frame_source->frame_config->current_n_sites, &frame_source->current_timestep, &frame_source->time, frame_source->simulation_box_limits, frame_source->dynamic_types, frame_source->molecule_flag, frame_source->dynamic_state_sampling);    

 	if (reference_atoms != frame_source->frame_config->current_n_sites) {
 		printf("Warning: Number of CG sites defined in top.in is not consistent with trajectory!\n");
 		return_value = 0;
 	} else {
 		// Skip through expected number of lines in frame body without parsing.
		for(int i=0; i < frame_source->frame_config->current_n_sites; i++) {
			check_and_read_next_line(frame_source->lammps_data->trajectory_stream, line);
    	}
	}
	 
    // Finish up by changing information simply determined by the data just read.
	for (int i = 0; i < DIMENSION; i++) frame_source->frame_config->simulation_box_half_lengths[i] = frame_source->simulation_box_limits[i][i] * 0.5;
    frame_source->current_frame_n += 1;

 	// Return 1 if successful, 0 otherwise.
 	return return_value;
}

int next_nothing(FrameSource* const frame_source)
{
	return 1;
}

void default_move_to_starting_frame(FrameSource* const frame_source) {
    for (int i = 0; i < frame_source->starting_frame - 1; i++) {
        if ((*frame_source->get_junk_frame)(frame_source) == 0) {
            printf("Failure attempting to skip frame %d. Check the trajectory file for errors.\n", i);
            exit(EXIT_FAILURE);
        }
    }
}

//-------------------------------------------------------------
// Helper functions for reading LAMMPS header and body
//-------------------------------------------------------------

void read_lammps_header(LammpsData *const lammps_data, int* const current_n_sites, int *const timestep, real *const time, matrix box, const int dynamic_types, const int molecule_flag, const int dynamic_state_sampling)
{
	double low = 0.0;
	double high = 0.0;
	std::string line;
	int flag = 1; 
	
	while(flag == 1) {
		//read next line of header (and wrap-up if end-of-file)
		check_and_read_next_line(lammps_data->trajectory_stream, line);
		
		//test if it is a labeled line (all LAMMPS labels start with "ITEM:")
		if( line.compare(0, 5, "ITEM:") == 0 ) {
			//find out which label matched (skip space after ITEM:)
			if( line.compare(6, 15, "NUMBER OF ATOMS") == 0) {
				
				//read number of atoms
				lammps_data->trajectory_stream >> *current_n_sites;
				
			} else if( line.compare(6, 10, "BOX BOUNDS") == 0) {
					
				//read in bounds (low high) for each dimensions
				for(int pos=0; pos <  DIMENSION; pos++) {
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
				int set_mol = 0;
				lammps_data->header_size = 0;
				
				//check for xpos and fpos as we tokenize string to determine number of columns in body
				while ((next = line.find_first_of(" ", prev)) != std::string::npos) {
					if( (next - prev) == 0 ) { //check if empty
						prev++;
						continue;
					} else if( (line.compare(prev, 1, "x") == 0) ||
							   (line.compare(prev, 1, "xu") == 0) ) {
						lammps_data->x_pos = lammps_data->header_size;
						set_x = 1;
					} else if( line.compare(prev, 2, "fx") == 0 ) {
						lammps_data->f_pos = lammps_data->header_size;
						set_f = 1;
					} else if( line.compare(prev, 3, "mol") == 0 ) {
						lammps_data->mol_pos = lammps_data->header_size;
						set_mol = 1;
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
                	if( (line.compare(prev, 1, "x") == 0 )  ||
					    (line.compare(prev, 1, "xu") == 0) ) {
                        lammps_data->x_pos = lammps_data->header_size;
                    	set_x = 1;
                    } else if( line.compare(prev, 2, "fx") == 0 ) {
                        lammps_data->f_pos = lammps_data->header_size;
                    	set_f = 1;
                    } else if( line.compare(prev, 3, "mol") == 0) {
						lammps_data->mol_pos = lammps_data->header_size;
						set_mol = 1;
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
				if ( (molecule_flag == 1) && (set_mol == 0) ) {
					printf("Warning: Molecule information not detected in header when parsing LAMMPS frame header!\n");
					exit(EXIT_FAILURE);
				}
			} else {
				printf("Unrecognized line in frame header: %s", line.c_str() );
				
			}	//close inner if/else if structure
		}		//close ITEM match
	}			//close while	
	return;
}

int read_dimension_lammps_body(LammpsData *const lammps_data, FrameConfig *const frame_config, const int dynamic_types, const int molecule_flag, const int dynamic_state_sampling)
{
	//read in current_n_sites lines to extract position and force information
	int j = 0;
	int return_value = 1;
	std::string line;
	
	//allocate space for array of strings based on lammps_data->header_size
	for(int i=0; i < frame_config->current_n_sites; i++)
	{
		//read in next line and tokenize 
		check_and_read_next_line(lammps_data->trajectory_stream, line);
		if(  (j = StringSplit(line, " \t", lammps_data->elements)) != lammps_data->header_size ) {	//allow for trailing white space
			printf("Warning: Number of fields detected in frame body");
			printf(" (%d) does not agree with number expected from frame header (%d)!\n", j, lammps_data->header_size);
			return_value = -1;
			break;
			}
			
		//extract position information
		for(j = 0; j < DIMENSION; j++) {
			frame_config->x[i][j] = atof( lammps_data->elements[j + lammps_data->x_pos].c_str() );
		}
		
		//extract force information
		for(j = 0; j < DIMENSION; j++) {
			frame_config->f[i][j] = atof( lammps_data->elements[j + lammps_data->f_pos].c_str() );
		}

		//extract type information
		if(dynamic_types == 1) { //check if dynamic_type is set
			frame_config->cg_site_types[i] = atoi( lammps_data->elements[lammps_data->type_pos].c_str() );
		}
		if(molecule_flag == 1) { //check if molecule_flag is set
			frame_config->molecule_ids[i] = atoi( lammps_data->elements[lammps_data->mol_pos].c_str() ) - 1;
		}
		if(dynamic_state_sampling == 1) { // check if dynamic_state_sampling is set
			lammps_data->cg_site_state_probabilities[i] = atof( lammps_data->elements[lammps_data->state_pos].c_str() );
		}
	}
	return return_value;
}

int read_scalar_lammps_body(LammpsData *const lammps_data, FrameConfig *const frame_config, const int dynamic_types, const int molecule_flag, const int dynamic_state_sampling)
{
	//read in current_n_sites lines to extract position and force information
	int j = 0;
	int return_value = 1;
	std::string line;
	
	//allocate space for array of strings based on lammps_data->header_size
	for(int i=0; i < frame_config->current_n_sites; i++)
	{
		//read in next line and tokenize 
		check_and_read_next_line(lammps_data->trajectory_stream, line);	
		if(  (j = StringSplit(line, " \t", lammps_data->elements)) != lammps_data->header_size ) {	//allow for trailing white space
			printf("Warning: Number of fields detected in frame body");
			printf(" (%d) does not agree with number expected from frame header (%d)!\n", j, lammps_data->header_size);
			return_value = -1;
			}
			
		//extract position information
		for(j = 0; j < DIMENSION; j++) {
			frame_config->x[i][j] = atof( lammps_data->elements[j + lammps_data->x_pos].c_str() );
		}
		
		//extract scalar "force" information
		frame_config->f[i][0] = atof( lammps_data->elements[lammps_data->f_pos].c_str() );

		//extract type information
		if(dynamic_types == 1) { //check if dynamic_type is set
			frame_config->cg_site_types[i] = atoi( lammps_data->elements[lammps_data->type_pos].c_str() );
		}
		if(molecule_flag == 1) { //check if molecule_flag is set
			frame_config->molecule_ids[i] = atoi( lammps_data->elements[lammps_data->mol_pos].c_str() ) - 1;
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

void read_frame_weights(FrameSource* const frame_source, const int start_frame, const int n_frames, const std::string &extension)
{
	std::string filename = "frame_weights." + extension;
    read_frame_values(filename.c_str(), start_frame, n_frames, frame_source->frame_weights);
	
    double total = 0.0;
    for (int i = 0; i < n_frames; i++) {
 		total += frame_source->frame_weights[i];
 	}
    frame_source->total_frame_weights = 1.0 / total;
}

// Read framewise information from an auxiliary file.
// Note: This supports:
// 1) the virial constraint for all frames in file 'p_con.in'.
// The pressure_constraint_rhs_vector is the RHS of eq. (12) in JCP,123,134105,2005
// 2) the framewise relative entropy observable method in file 'observable.in'. 
// 3) frame weights for statistical reweighting through 'frame_weights.in', 'frame_weights.ref', and 'frame_weights.cg'
// 4) framewise observable values for newobs.x through 'observables.ref' and 'observables.cg'

void read_frame_values(const char* filename, const int start_frame, const int n_frames, double* &values)
{
    std::ifstream vals_in;
    values = new double[n_frames];
    check_and_open_in_stream(vals_in, filename);
    read_stream_into_array(vals_in, start_frame, n_frames, values);
    vals_in.close();
}

inline void read_stream_into_array(std::ifstream &in_file, const int start_frame, const int n_frames, double* &values)
{
	double junk;
    for (int i = 0; i < start_frame - 1; i++) {
        in_file >> junk;
    }
    for (int i = 0; i < n_frames; i++) {
    	in_file >> values[i];
    }
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

// Local function prototypes for this section.
void recursive_cell_index_loop(const std::vector<int>& cell_number, std::vector<int> &indices, std::vector<int> &stencil, const std::vector<int> &hash_offset, const int current_dimension, dimension_neighbor_action action_todo);
void setup_pair_cell_stencil(const std::vector<int>& cell_number, std::vector<int> &cell_indices, std::vector<int> &stencil, const std::vector<int> &hash_offset);
void setup_3B_cell_stencil(const std::vector<int>& cell_number, std::vector<int> &cell_indices, std::vector<int> &stencil, const std::vector<int> &hash_offset);
int recursive_template_loop(const std::vector<int>& cell_number, const std::vector<int> &cell_indices, std::vector<int> &shift_indices, std::vector<int> &stencil, const std::vector<int> &hash_offset, int stencil_counter, const int current_dimension, add_stencil_element action_todo);
int add_pair_stencil_element(const std::vector<int>& cell_number,  const std::vector<int> &cell_indices, std::vector<int> &shift_indices, std::vector<int> &stencil, const std::vector<int> &hash_offset, int stencil_counter);
int add_3B_stencil_element(const std::vector<int>& cell_number, const std::vector<int> &cell_indices, std::vector<int> &shift_indices, std::vector<int> &stencil, const std::vector<int> &hash_offset, int stencil_counter);

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
 	for (int i = 0; i < DIMENSION; i++) {
 		if (cutoff > simulation_box_half_lengths[i]) {
	        printf("Cutoff is larger than half of the simulation box size!\n");
    	    exit(EXIT_FAILURE);
    	}
    }
	
	cell_number = std::vector<int>(DIMENSION);
	cell_size = std::vector<double>(DIMENSION);
	int too_small = 0;
	
	// Determine the number of cells in the box first by calculateng the number of cells needed to span each dimension.
    for (int i = 0; i < DIMENSION; i++) {
    	cell_number[i] = (int)(2.0 * simulation_box_half_lengths[i] / cutoff);
    }

    // Check that there are enough cells to make a cell list worthwhile.
    for (int i = 0; i < DIMENSION; i++) {
    	// NOTE: I am not sure that throwing away the cell list entirely is the smartest thing to do.
    	// It would be possible to shrink the cell list dimension by the number of dimensions that are too small.
    	
    	// This is a special case where the number of cells is smaller than the stencil in at least one dimension.
    	if (cell_number[i] < 3) {
    		std::fill(cell_number.begin(), cell_number.end(), 3);
    		std::fill(cell_size.begin(), cell_size.end(), -1.0);
    		too_small = 1;
    		break;
    	}
    }
    // Otherwise, continue with cell list setup.
    if (too_small == 0) {
    	for (int i = 0; i < DIMENSION; i++) {
	        cell_size[i] = 2.0 * simulation_box_half_lengths[i] / (double)(cell_number[i]);
	    }
	}
	
	// Allocate arrays based on the total number of cells needed to cover the entire simulation box.
    size = 1;
    for (int i = 0; i < DIMENSION; i++) {
    	size *= cell_number[i];
    }
    
    head = std::vector<int>(size);
    list = std::vector<int>(current_n_sites);
}

// Populate the cell lists.

void BaseCellList::populateList(const int n_particles, std::array<double, DIMENSION>* const &particle_positions)
{
    assert(n_particles > 0);
    int icell; // The index for the cell that a particle is in;
    
    // Pre-compute offsets in the cell array for steps in each dimension.
    hash_offset.resize(DIMENSION);
    hash_offset[0] = 1;
    for (int i = 1; i < DIMENSION; i++) {
    	hash_offset[i] = hash_offset[i - 1] * cell_number[i - 1];
    }
    
    // Calculate the inverse of the size of a cell in each dimension.
	std::vector<double> cell_inv(DIMENSION);
	for (int i = 0; i < DIMENSION; i++) {
		cell_inv[i] = 1.0 / cell_size[i];
	}
	
	// If we are actually using cell_lists.
	// // At the moment this is checked by only looking at the first dimension,
	// // but if cell list use is NOT all-or-none then this check would be insufficient.
    if (cell_size[0] > 0.0) {
		// Assign the particles to cells and build the neighbor list for each cell (backwards).
        // Initialize the cell heads.
		for (int i = 0; i < size; i++) {
            head[i] = -1;
        }
        for (int i = 0; i < n_particles; i++) {
			// Determine each particles cell. 
			// This "hash" for each cell refers to the cell's index since that data is stored in a flat array (x + y * x_offset + z * x_offset * y_offsets + ...).
            icell = 0;
            for (int j = 0; j < DIMENSION; j++) {
            	icell += (int)( particle_positions[i][j] * cell_inv[j] ) * hash_offset[j];
            }
			// Make this particle the head of its cell's list.
			// This particle now "points" to the index that was previously the head of this cell's list (-1 if it is the first particle).
            list[i] = head[icell];
            head[icell] = i;
        }
    } else {
		// In this special case, it does not make sense to use actual cells.
		// So, this is simply a list of all particles with each particle connected to the particles adjacent to it (by index) in the list.
        for (int i = 1; i < size; i++) {
            head[i] = -1;
        }
        head[0] = 0;
        for (int i = 0; i < n_particles - 1; i++) {
            list[i] = i + 1;
        }
        list[n_particles - 1] = -1;
    }
}

// Set up a pair list stencil.

void PairCellList::setUpCellListStencil()
{
    std::vector<int> indices(DIMENSION);

    // Pre-compute offsets in the cell array for steps in each dimension.
    hash_offset.resize(DIMENSION);
    hash_offset[0] = 1;
    for (int i = 1; i < DIMENSION; i++) {
    	hash_offset[i] = hash_offset[i - 1] * cell_number[i - 1];
    }

	// Determine how many neighboring cells each cell has, 
	// but only half need to be looked at using Newton's third law.
	int neighbor_cells = (int)( pow( 3.0, (double)(DIMENSION) ) ) - 1;
	// Get the total number of cells.
	int number_cells = head.size();
	// The stencil vector is a flat vector that includes
	// the neighboring cells that need to be looked at for all cells.
	stencil_size = neighbor_cells / 2;
    stencil = std::vector<int>(number_cells * stencil_size);
        
    // Note: The acceptable numbers could be modeled as
    // Run number from 0 or 1 to 2^(dimension) - 1.
    // Treat 1 as +/-1 for all digits other than MSB
    // Exclude -1 for numbers with (dimension - 1) or more 0's

	recursive_cell_index_loop(cell_number, indices, stencil, hash_offset, 0, setup_pair_cell_stencil);
}

// Set up a three body list stencil.
void ThreeBCellList::setUpCellListStencil()
{
    std::vector<int> indices(DIMENSION);

    // Pre-compute offsets in the cell array for steps in each dimension.
    hash_offset.resize(DIMENSION);
    hash_offset[0] = 1;
    for (int i = 1; i < DIMENSION; i++) {
    	hash_offset[i] = hash_offset[i - 1] * cell_number[i - 1];
    }

	// Determine how many neighboring cells each cell has.
	// For three_body_interactions all neighboring cells need to be looked at.
	int neighbor_cells = (int)( pow( 3.0, (double)(DIMENSION) ) ) - 1;
	// Get the total number of cells.
	int number_cells = head.size();
	// The stencil vector is a flat vector that includes
	// the neighboring cells that need to be looked at for all cells.
	stencil_size = neighbor_cells;
    stencil = std::vector<int>(number_cells * stencil_size);
        
    // Note: The acceptable numbers could be modeled as
    // Run number from 0 or 1 to 2^(dimension) - 1.
    // Treat 1 as +/-1 for all digits other than MSB
    // Exclude -1 for numbers with (dimension - 1) or more 0's

	recursive_cell_index_loop(cell_number, indices, stencil, hash_offset, 0, setup_3B_cell_stencil);
}

void recursive_cell_index_loop(const std::vector<int>& cell_number, std::vector<int> &indices, std::vector<int> &stencil, const std::vector<int> &hash_offset, const int current_dimension, dimension_neighbor_action action_todo) 
{
	if (current_dimension == (DIMENSION - 1) ) {
		// base case: loop over last dimension and call action_todo.
		for (int i = 0; i < cell_number[current_dimension]; i++) {
			indices[current_dimension] = i;
			action_todo(cell_number, indices, stencil, hash_offset);
		}
	} else {
		// otherwise, loop over this dimension (setting its index) and call the next layer.
		for (int i = 0; i < cell_number[current_dimension]; i++) {
			indices[current_dimension] = i;
			recursive_cell_index_loop(cell_number, indices, stencil, hash_offset, current_dimension + 1, action_todo);
		}
	}
}

void setup_pair_cell_stencil(const std::vector<int> &cell_number, std::vector<int> &cell_indices, std::vector<int> &stencil, const std::vector<int> &hash_offset) 
{
	// Determine the index for the current cell.
	std::vector<int> shift_indices(DIMENSION);
	int stencil_counter = 0;
	
	// Try combinations of 0, +/- 1 offsets for all but the first dimension (only, 0 and 1).
	for (int i = 0; i <= 1; i++) {
		shift_indices[0] = i;
		stencil_counter = recursive_template_loop(cell_number, cell_indices, shift_indices, stencil, hash_offset, stencil_counter, 1, add_pair_stencil_element);
	}
}

void setup_3B_cell_stencil(const std::vector<int> &cell_number, std::vector<int> &cell_indices, std::vector<int> &stencil, const std::vector<int> &hash_offset) 
{
	// Determine the index for the current cell.
	std::vector<int> shift_indices(DIMENSION);
	int stencil_counter = 0;
	
	// Try combinations of 0, +1, and -1 offsets for all dimensions.
	stencil_counter = recursive_template_loop(cell_number, cell_indices, shift_indices, stencil, hash_offset, stencil_counter, 0, add_3B_stencil_element);
}

int recursive_template_loop(const std::vector<int> &cell_number, const std::vector<int> &cell_indices, std::vector<int> &shift_indices, std::vector<int> &stencil, const std::vector<int> &hash_offset, int stencil_counter, const int current_dimension, add_stencil_element action_todo) 
{
	if (current_dimension == (DIMENSION - 1) ) {
	// base case: loop over last dimension and call action_todo.
		for (int i = -1; i <= 1; i++) {
			shift_indices[current_dimension] = i;
			stencil_counter = action_todo(cell_number, cell_indices, shift_indices, stencil, hash_offset, stencil_counter);
		}
	} else if (current_dimension >= DIMENSION) {
		// Check for 1D edge case.
		stencil_counter = action_todo(cell_number, cell_indices, shift_indices, stencil, hash_offset, stencil_counter);
	} else {
		// otherwise, loop over this dimension (setting its index) and call the next layer.
		for (int i = -1; i <= 1; i++) {
			shift_indices[current_dimension] = i;
			stencil_counter = recursive_template_loop(cell_number,cell_indices, shift_indices, stencil, hash_offset, stencil_counter, current_dimension + 1, action_todo);
		}
	}
	return stencil_counter;
}

int add_pair_stencil_element(const std::vector<int> &cell_number, const std::vector<int> &cell_indices, std::vector<int> &shift_indices, std::vector<int> &stencil, const std::vector<int> &hash_offset, int stencil_counter)
{
	// Check that this set of offsets is acceptable.
	// The first non-zero offset must be positive.
	int non_zero = 0;
	for (int i = 0; i < DIMENSION; i++) {
		if (shift_indices[i] == -1) {
			return stencil_counter;
		} else if (shift_indices[i] == 1) {
			non_zero = 1;
			break;
		}
	}
	if (non_zero == 0) {
		return stencil_counter;
	}
	// Otherwise, this set of offsets is acceptable.
	
	// Determine the cell's hash index.
	int cell_index = 0;
	for (int i = 0; i < DIMENSION; i++) {
		cell_index += cell_indices[i] * hash_offset[i];
	}
	
	// Now, determine this cell's beginning index in the stencil vector.
	int neighbor_cells = (int)( pow( 3.0, (double)(DIMENSION) ) ) - 1;
	int cell_stencil = cell_index * neighbor_cells / 2;
	
	// Determine the shifted cell's hash index.
	int shifted_cell_index = 0;
	for (int i = 0; i < DIMENSION; i++) {
		shifted_cell_index += ( (cell_indices[i] + shift_indices[i] + cell_number[i]) % cell_number[i] ) * hash_offset[i];
	}
	
	// Add this cell to the stencil.
	stencil[cell_stencil + stencil_counter] = shifted_cell_index;
	stencil_counter++;
	return stencil_counter;
}

int add_3B_stencil_element(const std::vector<int> &cell_number, const std::vector<int> &cell_indices, std::vector<int> &shift_indices, std::vector<int> &stencil, const std::vector<int> &hash_offset, int stencil_counter)
{
	// Check that this set of offsets is acceptable.
	int non_zero = 0;
	for (int i = 0; i < DIMENSION; i++) {
		if (shift_indices[i] != 0) {
			non_zero = 1;
			break;
		}
	}
	if (non_zero == 0) return stencil_counter;
	// Otherwise this set of shift_indices is acceptable.
	
	// Determine the cell's hash index.
	int cell_index = 0;
	for (int i = 0; i < DIMENSION; i++) {
		cell_index += cell_indices[i] * hash_offset[i];
	}
	
	// Now, determine this cell's beginning index in the stencil vector.
	int neighbor_cells = (int)( pow( 3.0, (double)(DIMENSION) ) ) - 1;
	int cell_stencil = cell_index * neighbor_cells;
	
	// Determine the shifted cell's hash index.
	int shifted_cell_index = 0;
	for (int i = 0; i < DIMENSION; i++) {
		shifted_cell_index += ( (cell_indices[i] + shift_indices[i] + cell_number[i]) % cell_number[i] ) * hash_offset[i];
	}
	
	// Add this cell to the stencil.
	stencil[cell_stencil + stencil_counter] = shifted_cell_index;
	stencil_counter++;
	return stencil_counter;
}
