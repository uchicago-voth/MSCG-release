#ifndef _read_lammps_h_
#define _read_lammps_h_

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <vector>
#include <stdint.h>
#include <string>

//-------------------------------------------------------------
// struct for keeping track of LAMMPS frame data
//-------------------------------------------------------------

struct rLammpsData {
	std::ifstream trajectory_stream;
	int type_pos;			// Index for type element in frame body
	int x_pos;				// Starting index for position elements in frame body
	int f_pos;				// Starting index for force elements in frame body
	int state_pos;			// Starting index for state probabilities in frame_body
	int header_size;		// Number of columns for header/body of frame
	std::string* elements; 	// Array to store tokenized header elements 
	double* cg_site_state_probabilities;   // A list of the probabilities for all states of all CG particles (used if dynamic_state_sampling = 1) (currently only for 2 states)
};

// Function prototypes
rLammpsData* rread_initial_lammps_frame(rLammpsData* lammps_data, const int n_cg_sites, int* cg_site_types, double* simulation_box_half_lengths, double* x, double* f, char* filename);
rLammpsData* rread_next_lammps_frame(rLammpsData* lammps_data, int const reference_atoms, double* simulation_box_half_lengths, double* x, double* f);
void rfinish_lammps_reading(rLammpsData* lammps_data);

#endif