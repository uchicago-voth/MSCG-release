//
//  interaction_model.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#ifndef _interaction_model_h
#define _interaction_model_h

#include <list>
#include <string>
#include <vector>

#include "interaction_hashing.h"
#include "splines.h"
#include "topology.h"
#include "misc.h"
#include "control_input.h"

struct MATRIX_DATA;
class PairCellList;
class ThreeBCellList;
struct InteractionClassComputer;
struct ThreeBodyNonbondedClassComputer;

// Function called externally
void free_interaction_data(CG_MODEL_DATA* cg);

//-------------------------------------------------------------
// Enumerated type definitions
//-------------------------------------------------------------

enum InteractionClassType {kPairNonbonded, kPairBonded, kAngularBonded, kDihedralBonded, kThreeBodyNonbonded};
// function pointer "type" used for polymorphism of matrix element calculation (for order parameter and pair nonbonded types)
typedef void (*calc_pair_matrix_elements)(InteractionClassComputer* const, const rvec*, const real*, MATRIX_DATA* const);

//-------------------------------------------------------------
// Interaction-model-related type definitions
//-------------------------------------------------------------

// This stores parameters that define an interaction class.

struct InteractionClassSpec {
	protected:
	BasisType basis_type;
    int bspline_k;
    double fm_binwidth;

	public:
    InteractionClassType class_type;
    int class_subtype;
    int basis_funcs_per_evaluation;
    
    // output only parameters
    int output_spline_coeffs_flag;
    double output_binwidth;
    int output_parameter_distribution;
    FILE** output_range_file_handles;

    int n_cg_types;    
    int n_defined;
    std::vector<unsigned> defined_to_possible_intrxn_index_map;
    int n_to_force_match;
    std::vector<unsigned> defined_to_matched_intrxn_index_map;
    std::vector<unsigned> interaction_column_indices;
    int n_tabulated;
    std::vector<unsigned> defined_to_tabulated_intrxn_index_map;
    
    double *lower_cutoffs;
    double *upper_cutoffs;
    double cutoff;
    
    double external_table_spline_binwidth;
    double **external_table_spline_coefficients;

	public:
	// These virtual functions need to be implemented for every new interaction type.
	virtual void determine_defined_intrxns(TopologyData*) = 0;
	virtual int get_n_body (void) const = 0;
    virtual std::string get_full_name(void) const = 0;
    virtual std::string get_short_name(void) const = 0;
    virtual std::string get_table_name(void) const = 0;
    virtual char get_char_id(void) const = 0;

	// Helper and implementation functions.
	void adjust_cutoffs_for_basis(int i);
    void adjust_cutoffs_for_type(int i);
    void setup_for_defined_interactions(TopologyData* topo_data); 
	void setup_indices_in_fm_matrix(void);
	void read_interaction_class_ranges(FILE *range_in); 
    int read_table(FILE* external_spline_table, int line, int offset);
	int read_bspline_table(FILE* external_spline_table, int line, int offset);
	void free_tabulated_interaction_data(FILE* spline_output_filep);
	void free_force_tabulated_interaction_data(void);
	
    std::string get_interaction_name(char **type_names, const int intrxn_index_among_defined) const;
    std::vector<int> get_interaction_types(const int intrxn_index_among_defined) const;
    inline int get_index_from_hash(const int hash_val) const {if (defined_to_possible_intrxn_index_map.size() == 0) return hash_val; else return SearchIntTable(defined_to_possible_intrxn_index_map, hash_val);}
    inline int get_hash_from_index(const int index) const {if (defined_to_possible_intrxn_index_map.size() > 0) return defined_to_possible_intrxn_index_map[index]; else return index;}

	// Functions meant to be eliminated.	
	// set_n_defined is only used in topology for three body interaction it should be eliminated eventually
	inline void set_n_defined(int n) {
		n_defined = n;
	};
	
	// Accessor (getter) functions.
	inline int get_num_basis_func(void) const {
		return interaction_column_indices[n_to_force_match];
    };

	inline BasisType get_basis_type(void) {
		return basis_type;
	};
	inline int get_n_defined(void) {
		return n_defined;
	};
	inline int get_bspline_k(void) {
		return bspline_k;
	};
	inline double get_fm_binwidth(void) {
		return fm_binwidth;
	};
	
	~InteractionClassSpec() {
		delete [] lower_cutoffs;
		delete [] upper_cutoffs;
	}
};

// Info needed for FM calculation of each interaction class, very closely
// related to the below struct. (Will be rebuilt from the below struct later.)

struct InteractionClassComputer {
	
    // Raw interaction class specifications
    InteractionClassSpec *ispec;
    double cutoff2;                             // Squared cutoff; used only for nonbonded interactions

    // Matrix-locations for storing results of computation
    int trajectory_block_frame_index;           // Index of the current frame in the current block of frames
    int current_frame_starting_row;             // Starting row number for the block of the FM matrix determined by the current frame
    int interaction_class_column_index;         // The starting column of the FM matrix block corresponding to the current class of interactions
    int basis_function_column_index;            // Starting column index for the matrix block corresponding to the current active interaction in the current class

    // Interacting particle indices: 
    // pair interactions: k-l; 
    // three-body interactions: k-j-l; 
    // four-body interactions: k-i-j-l.
    int k;
    int l;
    int i;
    int j;
    
    // Temps for determining which interaction the particles interact with.
    int index_among_defined_intrxns;
    int index_among_matched_interactions;
    int index_among_tabulated_interactions;

    // Calculation intermediates for the interaction
    double intrxn_param;                       // The interaction parameter for any single-parameter interaction (ie distance, angle, dihedral angle)
    double intrxn_param_less_lower_cutoff;     // Pair distance from the pair_nonbonded_interaction_lower_cutoffs_XOR_lower_cutoffs
    double stillinger_weber_angle_parameter;   // Current interaction's SW angle param

    // Function called to calculate matrix elements corresponding to an interaction in the current class of interactions.
    void (*calculate_fm_matrix_elements)(InteractionClassComputer* const self, const rvec *x, const real *simulation_box_half_lengths, MATRIX_DATA* const mat);
    // Function called to evaluate the values of functions in an interaction's basis set
    void (*set_up_fm_bases)(void);

	virtual void class_set_up_computer() = 0;  
    // Function to calculate index of an actual interaction among all possible interactions for the current class
	virtual int calculate_hash_number(int* const cg_site_types, const int n_cg_types) = 0;
	
	void set_up_computer(InteractionClassSpec* const ispec_pt, int *curr_iclass_col_index);	

	void calc_external_spline_interaction(void);
	void calc_grid_of_force_vals(const std::vector<double> &spline_coeffs, const int index_among_defined_intrxns, const double binwidth, std::vector<double> &axis_vals, std::vector<double> &force_vals);
	void calc_grid_of_force_and_deriv_vals(const std::vector<double> &spline_coeffs, const int index_among_defined_intrxns, const double binwidth, std::vector<double> &axis_vals, std::vector<double> &force_vals, std::vector<double> &deriv_vals);
	
	void walk_neighbor_list(MATRIX_DATA* const mat, calc_pair_matrix_elements calc_matrix_elements, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, const rvec* x, const real* simulation_box_half_lengths);

    // Spline computation objects for force matched and
    // tabulated interactions.
    SplineComputer* fm_s_comp;
    SplineComputer* table_s_comp;

    // Preallocating this temporary is worth ~20% of runtime in serial_fm.
    std::vector<double> fm_basis_fn_vals;
    std::vector<double> table_basis_fn_vals;

	protected:
	void calculate_bspline_matrix_elements(void);
	void calculate_linear_spline_matrix_elements(void);
};

struct PairNonbondedClassSpec: InteractionClassSpec {
	inline PairNonbondedClassSpec(ControlInputs* control_input) {
		class_type = kPairNonbonded;
		cutoff = control_input->pair_nonbonded_cutoff;
		basis_type = (BasisType) control_input->basis_set_type;
		output_spline_coeffs_flag = control_input->output_spline_coeffs_flag;
		fm_binwidth = control_input->pair_nonbonded_fm_binwidth;
		bspline_k = control_input->nonbonded_bspline_k;
		output_binwidth = control_input->pair_nonbonded_output_binwidth;
		output_parameter_distribution = control_input->output_pair_nonbonded_parameter_distribution;
	}
	
	void determine_defined_intrxns(TopologyData *topo_data) {
        n_defined = calc_n_distinct_pairs(topo_data->n_cg_types);
	}
	
	inline int get_n_body () const { return 2;}
    inline std::string get_full_name(void) const {return "pair nonbonded";}
    inline std::string get_short_name(void) const {return "";}
    inline std::string get_table_name(void) const {return "short_range";}
    inline char get_char_id(void) const {return 'n';}
};

struct PairBondedClassSpec: InteractionClassSpec {
	inline PairBondedClassSpec(ControlInputs* control_input) {
		class_type = kPairBonded;
		cutoff = control_input->pair_nonbonded_cutoff;
		basis_type = (BasisType) control_input->basis_set_type;
		output_spline_coeffs_flag = control_input->output_spline_coeffs_flag;
		fm_binwidth = control_input->pair_bond_fm_binwidth;
		bspline_k = control_input->pair_bond_bspline_k;
		output_binwidth = control_input->pair_bond_output_binwidth;
		output_parameter_distribution = control_input->output_pair_bond_parameter_distribution;
	}
	
	inline ~PairBondedClassSpec() {}
	
	void determine_defined_intrxns(TopologyData *topo_data) {
		int n_possible_interactions = calc_n_distinct_pairs(topo_data->n_cg_types);
        n_defined = calc_n_active_interactions(topo_data->bond_type_activation_flags, n_possible_interactions);
        defined_to_possible_intrxn_index_map = std::vector<unsigned>(n_defined, 0);
        set_up_interaction_type_hash_array(topo_data->bond_type_activation_flags, n_possible_interactions, defined_to_possible_intrxn_index_map);
	}
	
	int get_n_body () const { return 2;}
    inline std::string get_full_name(void) const {return "pair bonded";}
    inline std::string get_short_name(void) const {return "bon";}
    inline std::string get_table_name(void) const {return "bond";}
    inline char get_char_id(void) const {return 'b';}
};

struct AngularClassSpec: InteractionClassSpec {
	inline AngularClassSpec(ControlInputs* control_input) {
		class_type = kAngularBonded;
    	basis_type = (BasisType) control_input->basis_set_type;
    	output_spline_coeffs_flag = control_input->output_spline_coeffs_flag;
    	class_subtype = control_input->angle_interaction_style;
    	fm_binwidth = control_input->angle_fm_binwidth;
    	bspline_k = control_input->angle_bspline_k;
    	output_binwidth = control_input->angle_output_binwidth;
		output_parameter_distribution = control_input->output_angle_parameter_distribution;
    }
	
	inline ~AngularClassSpec() {}
	
	void determine_defined_intrxns(TopologyData *topo_data) {
		int n_possible_interactions = calc_n_distinct_triples(topo_data->n_cg_types);
        n_defined = calc_n_active_interactions(topo_data->angle_type_activation_flags, n_possible_interactions);
        defined_to_possible_intrxn_index_map = std::vector<unsigned>(n_defined, 0);
        set_up_interaction_type_hash_array(topo_data->angle_type_activation_flags, n_possible_interactions, defined_to_possible_intrxn_index_map);
	}
	
	int get_n_body () const { return 3;}
    inline std::string get_full_name(void) const {return "angular bonded";}
    inline std::string get_short_name(void) const {return "ang";}
    inline std::string get_table_name(void) const {return "angle";}
    inline char get_char_id(void) const {return 'a';}
};

struct DihedralClassSpec: InteractionClassSpec {
	inline DihedralClassSpec(ControlInputs* control_input) {
		class_type = kDihedralBonded;
    	basis_type = (BasisType) control_input->basis_set_type;
    	output_spline_coeffs_flag = control_input->output_spline_coeffs_flag;
    	class_subtype = control_input->dihedral_interaction_style;
    	fm_binwidth = control_input->dihedral_fm_binwidth;
		bspline_k = control_input->dihedral_bspline_k;
		output_binwidth = control_input->dihedral_output_binwidth;
		output_parameter_distribution = control_input->output_dihedral_parameter_distribution;
	}
	
	inline ~DihedralClassSpec() {}
	
	void determine_defined_intrxns(TopologyData *topo_data) {
		int n_possible_interactions = calc_n_distinct_quadruples(topo_data->n_cg_types);
        n_defined = calc_n_active_interactions(topo_data->dihedral_type_activation_flags, n_possible_interactions);
        defined_to_possible_intrxn_index_map = std::vector<unsigned>(n_defined, 0);
        set_up_interaction_type_hash_array(topo_data->dihedral_type_activation_flags, n_possible_interactions, defined_to_possible_intrxn_index_map);
	}
	
	int get_n_body () const { return 4;}
    inline std::string get_full_name(void) const {return "dihedral bonded";}
    inline std::string get_short_name(void) const {return "dih";}
    inline std::string get_table_name(void) const {return "dihedral";}
    inline char get_char_id(void) const {return 'd';}
};

struct ThreeBodyNonbondedClassSpec: InteractionClassSpec {

	double three_body_gamma;
    double* three_body_nonbonded_cutoffs;
    double* stillinger_weber_angle_parameters_by_type;
    double stillinger_weber_angle_parameter;
 	
	inline ThreeBodyNonbondedClassSpec(ControlInputs* control_input) {
		class_type = kThreeBodyNonbonded;
		basis_type = (BasisType) control_input->basis_set_type;
		output_spline_coeffs_flag = control_input->output_spline_coeffs_flag;
		class_subtype = control_input->three_body_flag;
		fm_binwidth = control_input->three_body_fm_binwidth;
		bspline_k = control_input->three_body_bspline_k;
		output_binwidth = control_input->three_body_nonbonded_output_binwidth;
		output_parameter_distribution = 0;
		three_body_gamma = control_input->gamma;
		n_defined = 0;
	}
	
	inline ~ThreeBodyNonbondedClassSpec() {
	    if (class_subtype > 0) {
			delete [] three_body_nonbonded_cutoffs;
			delete [] stillinger_weber_angle_parameters_by_type;
		}
	}
	
	void determine_defined_intrxns(TopologyData *topo_data) {n_defined = 0;}
	int get_n_body () const { return 3;}
    inline std::string get_full_name(void) const {return "three body nonbonded";}
    inline std::string get_short_name(void) const {return "";}
    inline std::string get_table_name(void) const {return "three_body";}
    inline char get_char_id(void) const {return '3';}
};

struct PairNonbondedClassComputer : InteractionClassComputer {
	void class_set_up_computer(void);
	void calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const PairCellList& pair_cell_list, const rvec* x, const real* simulation_box_half_lengths);
    int calculate_hash_number(int* const cg_site_types, const int n_cg_types) {
	    return calc_two_body_interaction_hash(cg_site_types[k], cg_site_types[l], n_cg_types);
	}
};

struct PairBondedClassComputer : InteractionClassComputer {
	void class_set_up_computer(void);
	void calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const rvec* x, const real* simulation_box_half_lengths); 
    int calculate_hash_number(int* const cg_site_types, const int n_cg_types) {
	    return calc_two_body_interaction_hash(cg_site_types[k], cg_site_types[l], n_cg_types);
	}
};

struct AngularClassComputer : InteractionClassComputer {
	void class_set_up_computer(void);
	void calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const rvec* x, const real* simulation_box_half_lengths);
    int calculate_hash_number(int* const cg_site_types, const int n_cg_types) {
	    return calc_three_body_interaction_hash(cg_site_types[j], cg_site_types[k], cg_site_types[l], n_cg_types);
	}
};

struct DihedralClassComputer : InteractionClassComputer {
	void class_set_up_computer(void);
	void calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const rvec* x, const real* simulation_box_half_lengths);
    int calculate_hash_number(int* const cg_site_types, const int n_cg_types) {
		return calc_four_body_interaction_hash(cg_site_types[i], cg_site_types[j], cg_site_types[k], cg_site_types[l], n_cg_types);
	}
};

struct ThreeBodyNonbondedClassComputer : InteractionClassComputer {
	double coef1[100];
 	
	void special_set_up_computer(InteractionClassSpec* const ispec_pt, int *curr_iclass_col_index);
	void class_set_up_computer(void) {} ;
	void calculate_interactions(MATRIX_DATA* const mat, int traj_block_frame_index, int curr_frame_starting_row, const int n_cg_types, const TopologyData& topo_data, const ThreeBCellList& three_body_cell_list, const rvec* x, const real* simulation_box_half_lengths);
    void calculate_bspline_elements_and_deriv_elements(double* coef1);
	void calculate_bspline_deriv_elements(double* coef1);
    int calculate_hash_number(int* const cg_site_types, const int n_cg_types) {
	    return calc_three_body_interaction_hash(cg_site_types[j], cg_site_types[k], cg_site_types[l], n_cg_types);
	}
};

// Major struct responsible for keeping track of all cg model parameters, interaction definitions,
// and basis set specifications.

struct CG_MODEL_DATA {
   
    // Cutoff specifications.
    double pair_nonbonded_cutoff;           // Nonbonded pair interaction cutoff
    double pair_nonbonded_cutoff2;          // Squared cutoff distance for pair nonbonded interactions
    double three_body_nonbonded_cutoff2;    // Squared cutoff distance for three body nonbonded interactions

    // Topology specifications.
    TopologyData topo_data;
    
    // CG site number and type specifications.
    int n_cg_types;
    int n_cg_sites;
    char **name;

    // Interaction class specification structs.
    PairNonbondedClassSpec pair_nonbonded_interactions;
    PairBondedClassSpec pair_bonded_interactions;
    AngularClassSpec angular_interactions;
    DihedralClassSpec dihedral_interactions;
    ThreeBodyNonbondedClassSpec three_body_nonbonded_interactions;

    // Interaction class computation structs.
    PairNonbondedClassComputer pair_nonbonded_computer;
    PairBondedClassComputer pair_bonded_computer;
    AngularClassComputer angular_computer;
    DihedralClassComputer dihedral_computer;
    ThreeBodyNonbondedClassComputer three_body_nonbonded_computer;

	// List for interactions and computer used by iterators.
	std::list<InteractionClassSpec*> iclass_list;
	std::list<InteractionClassComputer*> icomp_list;
	
    // Three body "dummy topology" hack for carrying info from top.in into 
    // the range input reading function (where it is freed).
    int* tb_n;
    int** tb_list;

    // Non-matrix-associated output flags.
    int output_spline_coeffs_flag;          // 1 to output spline coefficients as well as force tables; 0 otherwise

	inline CG_MODEL_DATA(ControlInputs* control_input) :
		pair_nonbonded_cutoff(control_input->pair_nonbonded_cutoff),
		topo_data(control_input->max_pair_bonds_per_site, control_input->max_angles_per_site, control_input->max_dihedrals_per_site),
		pair_nonbonded_interactions(control_input), pair_bonded_interactions(control_input),
		angular_interactions(control_input), dihedral_interactions(control_input),
		three_body_nonbonded_interactions(control_input),
        output_spline_coeffs_flag(control_input->output_spline_coeffs_flag)
	{
    	topo_data.excluded_style = control_input->excluded_style;
		pair_nonbonded_cutoff2 = pair_nonbonded_cutoff * pair_nonbonded_cutoff;
		
		iclass_list.push_back(&pair_nonbonded_interactions);
		iclass_list.push_back(&pair_bonded_interactions);
		iclass_list.push_back(&angular_interactions);
		iclass_list.push_back(&dihedral_interactions);
		
		icomp_list.push_back(&pair_nonbonded_computer);
		icomp_list.push_back(&pair_bonded_computer);
		icomp_list.push_back(&angular_computer);
		icomp_list.push_back(&dihedral_computer);	
	}
		
	~CG_MODEL_DATA() {};
};

// Uniform-interface wrappers for interaction hashing
// In header for inlining
inline int calculate_two_body_interaction_hash(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types) {
    return calc_two_body_interaction_hash(cg_site_types[info->k], cg_site_types[info->l], n_cg_types);
}

inline int calculate_three_body_interaction_hash(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types) {
    return calc_three_body_interaction_hash(cg_site_types[info->j], cg_site_types[info->k], cg_site_types[info->l], n_cg_types);
}

inline int calculate_four_body_interaction_hash(InteractionClassComputer* const info, int* const cg_site_types, const int n_cg_types) {
    return calc_four_body_interaction_hash(cg_site_types[info->i], cg_site_types[info->j], cg_site_types[info->k], cg_site_types[info->l], n_cg_types);
}

//--------------------------------------------------------------------
// Functions for setting up the potential model that will be used
// in the CG model from a range.in file.
//--------------------------------------------------------------------

// Read interaction ranges and assign the interactions to be force matched, tabulated, or null.
void read_all_interaction_ranges(CG_MODEL_DATA* const cg);

// Read tabulated interaction data from file
void read_tabulated_interaction_file(CG_MODEL_DATA* const cg, int n_cg_types);

#endif
