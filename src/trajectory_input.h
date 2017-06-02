//
//  trajectory_input.h
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#ifndef _trajectory_input_h
#define _trajectory_input_h

# if _exclude_gromacs == 1
# else
#  undef _exclude_gromacs
#  define _exclude_gromacs 0
# endif

#include "misc.h"

#ifndef DIMENSION
#define DIMENSION 3
#endif

#include <array>
#include <cstdint>
#include <random>
#include <fstream>
#include <string>
#include <vector>

struct ControlInputs;
struct LammpsData;
struct XRDData;

typedef real matrix[3][3];

enum TrajectoryType {kGromacsTRR = 0, kGromacsXTC = 1, kLAMMPSDump = 2};

typedef void (*dimension_neighbor_action)(const std::vector<int> &cell_number, std::vector<int> &indices, std::vector<int> &stencil, const std::vector<int> &hash_offset);
typedef int (*add_stencil_element)(const std::vector<int> &cell_number, const std::vector<int> &cell_indices, std::vector<int> &shift_indices, std::vector<int> &stencil, const std::vector<int> &hash_offset, int stencil_counter);

//-------------------------------------------------------------
// High-level struct for passing one frame's configuration data
//-------------------------------------------------------------

struct FrameConfig {
    int current_n_sites;                // Total number of sites in this frame
    real* simulation_box_half_lengths;  // A list of half the box length in each dimension
    std::array<double, DIMENSION>* x;   // A list of all CG particle positions for a single frame stored in a flat array, x,y,z components contiguous 
    std::array<double, DIMENSION>* f;   // A list of all CG particle positions for a single frame stored in a flat array, x,y,z components contiguous    
	int* cg_site_types;				   	   // A list of all CG particle types (used if dynamic_types = 1)
	int* molecule_ids;					// A list of which molecule each particle is in (used if molecule_flag = 1)
	
	inline FrameConfig(const int n_sites) {
		current_n_sites = n_sites;
		x = new std::array<double, DIMENSION>[current_n_sites + 1];
		f = new std::array<double, DIMENSION>[current_n_sites + 1];
		simulation_box_half_lengths = new real[DIMENSION];
	};
	
	inline FrameConfig(const int n_sites, int* site_types) {
		current_n_sites = n_sites;
		x = new std::array<double, DIMENSION>[current_n_sites + 1];
		f = new std::array<double, DIMENSION>[current_n_sites + 1];

		simulation_box_half_lengths = new real[DIMENSION];
		cg_site_types = site_types;		
	};
	
	inline FrameConfig(const int n_sites, int* site_types, int* mol_ids) {
		current_n_sites = n_sites;
		x = new std::array<double, DIMENSION>[current_n_sites + 1];
		f = new std::array<double, DIMENSION>[current_n_sites + 1];

		simulation_box_half_lengths = new real[DIMENSION];
		cg_site_types = site_types;		
		molecule_ids = mol_ids;
	};
	
	inline ~FrameConfig() {
		delete [] simulation_box_half_lengths;
		delete [] x;
		delete [] f;
	};
};

//-------------------------------------------------------------
// High-level struct for keeping track of all frame reading data
//-------------------------------------------------------------

struct FrameSource {
    // Type-independent source specifications set in control.in.
    int use_statistical_reweighting;        // 1 to use per-frame statistical reweighting; 0 otherwise
    int pressure_constraint_flag;           // 1 to use the virial constraint; 0 otherwise
    int dynamic_types;						// 1 to use dynamic type tracking; 0 otherwise
    int molecule_flag;						// 2 to use molecule tracking from LAMMPS trajectory; 1 to set molecule information from top.in; 0 otherwise
    int dynamic_state_sampling;				// 1 to use dynamic state sampling; 0 otherwise
    int dynamic_state_samples_per_frame;	// Number of times each frame is resampled if dynamic_state_sampling is 1
    int bootstrapping_flag;					// 1 to use bootstrapping; 0 otherwise
	int bootstrapping_num_subsamples;		// Number of subsamples per estimate (amount of discrete frame weight to distribute) if bootstrapping_flag is 1
	int bootstrapping_num_estimates;			// Number of estimates (separate bootstrap estimates constructed) if bootstrapping_flag is 1
	uint_fast32_t random_num_seed;			// Random number seed only used if dynamic_state_sampling or bootstrapping_flag is 1
    int starting_frame;                     // Trajectory frame number to start from
    int n_frames;                           // Total number of frames to read for this force matching
    char trajectory_filename[1000];         // Trajectory file name (positions for .xtc, forces and positions for .trr)
    std::mt19937 mt_rand_gen;    // A Mersenne Twister random number generator for dynamic state sampling.
	int position_dimension;					// The number of elements in each particle's position vector.
	int scalar_matching_flag;				// Whether to match DIMENSION sized forces (0) or match scalar sized forces (1)
	
    // Type-dependent source data and functions
    TrajectoryType trajectory_type;         // 0 to use .trr format trajectories; 1 to use .xtc format trajectories; 2 to use LAMMPS trajectories
	XRDData* gromacs_data;
	LammpsData* lammps_data;

    // Type-dependent function to read the first frame of a given source
    // Performs initial sanity checks to make sure the frame is consistent 
    // with input specifications.
    void (*get_first_frame)(FrameSource * const frame_source, const int n_cg_sites, int* cg_site_types, int* mol_ids);
    // An optionally type-dependent function to skip frames of a given source
    void (*move_to_start_frame)(FrameSource * const frame_source);
    // Type-dependent function to read but not process the next frame of a given source
    int (*get_junk_frame)(FrameSource * const frame_source);
    // Type-dependent function to provide the next frame of a given source
    int (*get_next_frame)(FrameSource * const frame_source);
    // Type-dependent function to clean up after reading all desired frames
    void (*cleanup)(FrameSource * const frame_source);

    // Data for a single frame.
    int current_timestep;                           // The timestep of the current frame
    int current_frame_n;
    real time;                                      // The time value of the current frame
    matrix simulation_box_limits;                   // A matrix describing the limits of an orthogonal simulation box along each dimension

    // Data read for all frames at once, if at all.
    double* frame_weights;                          // A list of weights for statistical reweighting, one per trajectory frame
    double total_frame_weights;                     // The sum of the frame weights.
    double* pressure_constraint_rhs_vector;         // pressure_constraint_rhs_vector is the RHS of eq. (12) in JCP,123,134105,2005
	double* cg_observables;							// A list of frame-wise observables for relative entropy, one per trajectory frame for CG observable
	double* ref_observables;						// A list of frame-wise observables for relative entropy, one per trajectory frame for reference observable
	
	// Generate data for all frame at once, it at all.
	double** bootstrapping_weights;
	
	FrameConfig* frame_config;
    inline FrameConfig* getFrameConfig() { return frame_config;}
    void sampleTypesFromProbs();
};

//-------------------------------------------------------------
// Construct a frame source object.
//-------------------------------------------------------------

void parse_command_line_arguments(const int num_arg, char** arg, FrameSource* const frame_source);
void parse_entropy_command_line_arguments(const int num_arg, char** arg, FrameSource* const frame_source_cg, FrameSource* const frame_source_ref);

// Copy trajectory-reading specifications from ControlInputs to FRAME_DATA.
void parse_command_line_set(const char* arg1, const char* arg2, FrameSource* const frame_source_cg, FrameSource* const frame_source_ref, int& checker_cg, int& checker_ref);
void copy_control_inputs_to_frd(struct ControlInputs* const control_input, FrameSource* const frame_source);
void copy_control_inputs_to_frd(ControlInputs * const control_input, FrameSource* fs_ref, FrameSource* fs_cg);

//-------------------------------------------------------------
// Auxiliary-trajectory reading functions.
//-------------------------------------------------------------

// Read statistical weights for all needed frames.
void read_frame_weights(FrameSource* const frame_source, const int start_frame, const int n_frames, const std::string &extension);

// Read information relating to the virial constraint or frame-wise observable for all needed frames.
void read_frame_values(const char* filename, const int start_frame, const int n_frames, double* &vals);

//-------------------------------------------------------------
// Auxiliary-trajectory generating functions.
//-------------------------------------------------------------

// Generate bootstrapping weights for all estimates.
void generate_bootstrapping_weights(FrameSource* const frame_source, const int num_frames);

// Combine bootstrapping and reweighting weights.
void combine_reweighting_and_boostrapping_weights(FrameSource* const frame_source);

// Clean-up bootstrapping weights at end.
void free_bootstrapping_weights(FrameSource* const frame_source);

//--------------------------------------------------------------------
// Cell list routines for two- or three-body nonbonded interactions.
//--------------------------------------------------------------------

class BaseCellList {

	// This list is traversed for a given cell by looking up the head for a given cell (the index of head).
	// The value gives the index of the first particle in that cells list.
	// The next particle index is determined by looking up the value at the index of the previous particle's index in the "list".
	// This is repeated until the value is -1.

public:
    void init(const double cutoff, const FrameSource* const fr);
    void populateList(const int n_particles, std::array<double, DIMENSION>* const &particle_positions);
    inline int get_stencil_size() const { return stencil_size; };
    inline double get_cell_size(int i) const {return cell_size[i]; };
    int size;					// The total number of cells to cover the simulation box.
    std::vector<int> list;		// "Linked list" for neighbor list. 
								// The value at each particle's index is the next particle index in that cell's list. 
								// If it is the last particle in the list, its value is -1.
    std::vector<int> head;		// List of the first particle in each cell for all cells.
    std::vector<int> stencil;	// List of neighboring cells to look through during force computation.
    std::vector<int> hash_offset;
	
protected:
    // The number of cells in each dimension.
    std::vector<int> cell_number;
	// The size (in each dimension) that a given cell spans.
    std::vector<double> cell_size;
	int stencil_size;			// The number of neighboring cells surrounding a given cell that need to be searched through during force computation.

    void setUpCellListCells(const double cutoff, const real* simulation_box_half_lengths, const int current_n_sites);
    virtual void setUpCellListStencil() = 0;
};

class PairCellList: public BaseCellList {
protected:
    virtual void setUpCellListStencil();
};

class ThreeBCellList: public BaseCellList {
protected:
    virtual void setUpCellListStencil();
};

#endif
