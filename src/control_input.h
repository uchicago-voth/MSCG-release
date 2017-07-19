//
// control_input.h
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#ifndef _control_input_h
#define _control_input_h
#include <cstdint>

//-------------------------------------------------------------
// Control input struct definition
//-------------------------------------------------------------

typedef struct ControlInputs {
	// Data settings
    int starting_frame;
    int n_frames;
    int frames_per_traj_block;
    
    // Input specifications
    int use_statistical_reweighting;
	int pressure_constraint_flag;
	
	// Additional features
	int dynamic_types;
    int dynamic_state_sampling;
    int dynamic_state_samples_per_frame;
    int bootstrapping_flag;
    int bootstrapping_full_output_flag;
	int bootstrapping_num_estimates;
	int bootstrapping_num_subsamples;
    uint_fast32_t random_num_seed;					// Only used when dynamic_state_sampling or bootstrapping_flag is 1

    // Interaction style specifications.
    int angle_interaction_style;                // 1 to use distance-based angular interactions; 0 for angle-based angle interactions.
    int dihedral_interaction_style;             // 1 to use distance-based dihedral interactions; 0 for angle-based dihedral interactions.
    int three_body_flag;
    int three_body_nonbonded_exclusion_flag;
    int excluded_style;						// 0 no exclusions; 2 exclude 1-2 bonded; 3 exclude 1-2 and 1-3 bonded; 4 exclude 1-2, 1-3 and 1-4 bonded interactions
    double gamma;
    double pair_nonbonded_cutoff;
    int max_pair_bonds_per_site;
    int max_angles_per_site;
    int max_dihedrals_per_site;

    // Basis set specifications.
    double pair_nonbonded_fm_binwidth;
    double pair_bond_fm_binwidth;
    double angle_fm_binwidth;
    double dihedral_fm_binwidth;
    double three_body_fm_binwidth;
    int nonbonded_bspline_k;                // B-spline k value for nonbonded pair interactions
    int pair_bond_bspline_k;                // B-spline k value for bonded pair interactions
    int angle_bspline_k;                    // B-spline k value for bonded angular interactions
    int dihedral_bspline_k;                 // B-spline k value for bonded dihedral interactions
    int three_body_bspline_k;               // B-spline k value for nonbonded three body interactions
    int basis_set_type;
    
    // Output specifications. 
    int output_style;
    int output_solution_flag;    
    int output_residual;
    int output_spline_coeffs_flag;
    int output_normal_equations_rhs_flag;
    double pair_nonbonded_output_binwidth;
    double pair_bond_output_binwidth;
    double angle_output_binwidth;
    double dihedral_output_binwidth;
    double three_body_nonbonded_output_binwidth;

	// Rangefinder only output specifications
	int output_pair_nonbonded_parameter_distribution;
    int output_pair_bond_parameter_distribution;
    int output_angle_parameter_distribution;
    int output_dihedral_parameter_distribution;
    
    //REM and IBI specification
  	double iteration_step_size;
  	double temperature;
  	double boltzmann;

    // Matrix specifications
    int matrix_type;
    int itnlim;
    int iterative_calculation_flag;
    double tikhonov_regularization_param;
    int regularization_style;
    int bayesian_flag;
	int bayesian_max_iter;
    double rcond;
	double sparse_safety_factor; 
	int num_sparse_threads;
	
	ControlInputs(void);
	~ControlInputs(void);
} ControlInputs;

#endif
