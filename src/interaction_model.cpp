//
//  interaction_model.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include "interaction_model.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <sstream>

//---------------------------------------------------------------
// Prototypes for internal implementations.
//---------------------------------------------------------------

// Input extraction functions.
inline void read_types(const int n_body, std::vector<int> &types, std::string* elements, const int n_types, char** name);

// Non-standard setup functions.
void three_body_setup_for_defined_interactions(InteractionClassSpec* ispec, TopologyData* topo_data);
void three_body_setup_indices_in_fm_matrix(InteractionClassSpec* ispec);

// Input checking and error reporting functions.
void check_nonbonded_interaction_range_cutoffs(PairNonbondedClassSpec *const ispec, double const cutoff);
void report_tabulated_interaction_format_error(const int line);
void report_tabulated_interaction_data_consistency_error(const int line);
void report_fields_error(const std::string &full_name, const int n_expected, const int n_fields); 
inline void check_mode(char* mode);

// Get the name of a single defined interaction via its index among
// defined interactions.

std::string InteractionClassSpec::get_interaction_name(char **type_names, const int intrxn_index_among_defined, const std::string &delimiter) const
{
    // Name first by the types involved in the interaction.
    // Name as type1_type2..._typeN.
    std::vector<int> types = get_interaction_types(intrxn_index_among_defined);
    
    // Reorder for special angles and dihedrals.
    if (format == 1) {
    	if (types.size() == 3) {
    		swap_pair(types[0], types[1]); // (start) B A C => A B C (end)
    	} else { //size is 4
    		swap_pair(types[2], types[1]); // (start) B C A D => B A C D
    		swap_pair(types[0], types[1]); //         B A C D => A B C D (end)
    	}
    }
    
    std::string namestring = std::string(type_names[types[0] - 1]);
    for(unsigned i = 1; i < types.size(); i++) {
        namestring += delimiter + type_names[types[i] - 1];
    }
    return namestring;
}

std::vector<int> InteractionClassSpec::get_interaction_types(const int index_among_defined_intrxns) const 
{
    std::vector<int> types(get_n_body(), 0);
	invert_interaction_hash(get_hash_from_index(index_among_defined_intrxns), n_cg_types, types);
    return types;
}

// Check that specified nonbonded interactions do not extend past the nonbonded cutoff
void check_nonbonded_interaction_range_cutoffs(PairNonbondedClassSpec *const ispec, double const cutoff)
{
	for (int i = 0; i < ispec->n_defined; i++) {
		if (ispec->defined_to_matched_intrxn_index_map[i] != 0) {
			if (ispec->upper_cutoffs[i] > (cutoff + ispec->output_binwidth + VERYSMALL) ) {
				fprintf(stderr, "An upper cutoff (%lf) specified in the range file is larger than the pair nonbonded cutoff speicified in the control file (%lf).\n", ispec->upper_cutoffs[i], cutoff);
				fprintf(stderr, "This can lead to unphysical looking interactions.\n");
				fprintf(stderr, "Please adjust and run again.\n");
				fflush(stderr);
				exit(EXIT_FAILURE);
			}
		}
	}
}

// Functions for reading a range.in file and assigning the FM matrix column indices for each basis function.

void InteractionClassSpec::read_interaction_class_ranges(std::ifstream &range_in)
{
    printf("Reading interaction ranges for %d %s interactions.\n", get_n_defined(), get_full_name().c_str());    

    int total_tabulated = 0;
    int total_to_fm = 0;
    
    int n_fields;
    int n_expected = 3 + get_n_body();
	std::vector<int> types(get_n_body());
	std::string* elements = new std::string[n_expected + 1];
 	std::string line;
    char mode[10] = "";
	
    for (int i = 0; i < n_defined; i++) {

    	check_and_read_next_line(range_in, line);		
		// Check that this line has enough fields.
		if ( (n_fields = StringSplit(line, " \t\n", elements)) < n_expected ) report_fields_error(get_full_name(), n_expected, n_fields);
		read_rmin_class(elements, get_n_body(), i, mode);

        // If the mode is none, the interaction is not actually in the model.
        // If the mode is fm or fm+tab/tab+fm, the interaction should be force matched.
        // If the mode is tab or fm+tab/tab+fm, the interaction should be tabulated.
        
        check_mode(mode);

        if (strcmp(mode,"fm") == 0 || strcmp(mode,"fm+tab") == 0 || strcmp(mode,"tab+fm") == 0) {
			// This interaction is to be force matched.
            // Increment running total and set the new index.
            total_to_fm++;
            defined_to_matched_intrxn_index_map[i] = total_to_fm;
            // Adjust for a basis by rounding the cutoffs to even
            // numbers of bins.
			adjust_cutoffs_for_basis(i);
            // Adjust nonbonded interactions to match the global cutoff.
            adjust_cutoffs_for_type(i);
		}
        if (strcmp(mode,"tab") == 0 || strcmp(mode,"fm+tab") == 0 || strcmp(mode,"tab+fm") == 0) {
			// This interaction is tabulated.
			// Increment the running total of the tabulated interactions.
			total_tabulated++;
			defined_to_tabulated_intrxn_index_map[i] = total_tabulated;
		}			
    }
    n_to_force_match = total_to_fm;
    n_tabulated = total_tabulated;
    printf("Will force match %d %s interactions and %d are interactions tabulated", n_to_force_match, get_full_name().c_str(), n_tabulated);
    printf(".\n"); 
    delete [] elements; 
}

void InteractionClassSpec::smart_read_interaction_class_ranges(std::ifstream &range_in, char** name)
{
    printf("Reading interaction ranges for %d %s interactions.\n", get_n_defined(), get_full_name().c_str());    

    int total_tabulated = 0;
    int total_to_fm = 0;
    int total_intrxns = 0;
    int n_fields;
    int n_expected = 3 + get_n_body();
    
    std::vector<int> types(get_n_body());
	std::string* elements = new std::string[n_expected + 1];
    std::string line;
    char mode[10];
	
	while( std::getline(range_in, line) != NULL) {
	
		// Check that this line has enough fields.
		if ( (n_fields = StringSplit(line, " \t\n", elements)) < n_expected ) {	//allow for trailing white space
			if(total_intrxns == 0) report_fields_error(get_full_name(), n_expected, n_fields);
			// This allows the reading to terminate if 
			// a line of whitespace is encountered before the end of file.
			break;
		}
		
		// Extract the useful information.
		read_types(get_n_body(), types, &elements[0], n_cg_types, name);
		int index_among_defined = calc_interaction_hash(types, n_cg_types);
	
		// Read the low and high parameter values;
		read_rmin_class(elements, get_n_body(), index_among_defined, mode);
		check_mode(mode);
		total_intrxns++;

        if (strcmp(mode,"fm") == 0 || strcmp(mode,"fm+tab") == 0 || strcmp(mode,"tab+fm") == 0) {
			// This interaction is to be force matched.
            // Increment running total and set the new index.
            total_to_fm++;
            defined_to_matched_intrxn_index_map[index_among_defined] = total_to_fm;
            // Adjust for a basis by rounding the cutoffs to even
            // numbers of bins.
			adjust_cutoffs_for_basis(index_among_defined);
            // Adjust nonbonded interactions to match the global cutoff.
            adjust_cutoffs_for_type(index_among_defined);
		}
        if (strcmp(mode,"tab") == 0 || strcmp(mode,"fm+tab") == 0 || strcmp(mode,"tab+fm") == 0) {
			// This interaction is tabulated.
			// Increment the running total of the tabulated interactions.
			total_tabulated++;
			defined_to_tabulated_intrxn_index_map[index_among_defined] = total_tabulated;
		}			
    }

    n_to_force_match = total_to_fm;
    n_tabulated = total_tabulated;
    printf("Will force match %d %s interactions or %d are interactions tabulated", n_to_force_match, get_full_name().c_str(), n_tabulated);
    printf(".\n");
	delete [] elements;
}

// For use with smart_read_interaction_class

void InteractionClassSpec::read_rmin_class(std::string* &elements, const int position, const int index_among_defined, char* mode) 
{
	lower_cutoffs[index_among_defined] = atof(elements[position].c_str());
	upper_cutoffs[index_among_defined] = atof(elements[position + 1].c_str());
	sprintf(mode, "%s", elements[position + 2].c_str());	
}
		
inline void read_types(const int n_body, std::vector<int> &types, std::string* elements, const int n_types, char** name)
{
	for (int j = 0; j < n_body; j++) {
    	types[j] = match_type(elements[j], name, n_types);
        if( types[j] == -1) {
        	fprintf(stderr, "Unrecognized type %s!\n", elements[j].c_str());
        	fflush(stderr);
    	}
    }
}

inline void check_mode(char* mode)
{
    if (strcmp(mode,"none") != 0 && strcmp(mode,"fm") != 0 && strcmp(mode,"tab") != 0 && strcmp(mode,"fm+tab") != 0 && strcmp(mode,"tab+fm") != 0 ){
        fprintf(stderr, "Interaction mode %s is not recognized\n", mode);
        fflush(stderr);
    	exit(EXIT_FAILURE);
	}
}

void InteractionClassSpec::adjust_cutoffs_for_basis(int i)
{
    if (basis_type == kLinearSpline) {
        lower_cutoffs[i] = floor(lower_cutoffs[i] / output_binwidth + 0.5) * output_binwidth;
        if (lower_cutoffs[i] < 0.0) lower_cutoffs[i] = 0.0;
        upper_cutoffs[i] = lower_cutoffs[i] + floor((upper_cutoffs[i] - lower_cutoffs[i]) / fm_binwidth + 0.5) * fm_binwidth;
    } else if ((basis_type == kBSpline) ||  (basis_type == kBSplineAndDeriv)) {
        upper_cutoffs[i] = (((int)(upper_cutoffs[i] / output_binwidth)) + 1) * output_binwidth;
        lower_cutoffs[i] = upper_cutoffs[i] - ((int)((upper_cutoffs[i] - lower_cutoffs[i]) / fm_binwidth) + 1) * fm_binwidth;
    }
}

void InteractionClassSpec::adjust_cutoffs_for_type(int i)
{
    if (basis_type == kLinearSpline) {
        if (class_type == kPairNonbonded && fabs(upper_cutoffs[i] - cutoff - fm_binwidth) < VERYSMALL_F) upper_cutoffs[i] -= fm_binwidth;
    }
}

// Determine number of columns for each interaction to be force matched.

void InteractionClassSpec::setup_indices_in_fm_matrix(void)
{
	int counter = 0;
	int grid_i;
	interaction_column_indices = std::vector<unsigned>(n_to_force_match + 1, 0);

	for (int i = 0; i < n_defined; i++) {
		if (defined_to_matched_intrxn_index_map[i] != 0) {
			grid_i = floor((upper_cutoffs[i] - lower_cutoffs[i]) / fm_binwidth + 0.5) + 1;
			if (grid_i > 1000) {
				fprintf(stderr, "\nWarning: An individual interaction has more than 1000 bins associated with it!\n");
				fprintf(stderr, "Please check that this is intentional.\n");
				fprintf(stderr, "This may be a sign that the wrong angle_style and dihedral_style is selected.\n\n");
				fflush(stderr);
			}
			
			interaction_column_indices[counter + 1] = interaction_column_indices[counter] + grid_i;
			// BSplines include an extra bspline_k - 2 knots.
			if ((basis_type == kBSpline) || (basis_type == kBSplineAndDeriv)) interaction_column_indices[counter + 1] += get_bspline_k() - 2;
			counter++;
		}
	}
}

void ThreeBodyNonbondedClassSpec::setup_indices_in_fm_matrix(void)
{ 
    if (class_subtype > 0) {
		interaction_column_indices = std::vector<unsigned>(get_n_defined(), 0);
		interaction_column_indices[0] = 0;
		
		n_tabulated = 0;
		n_to_force_match  = n_force = n_defined;
		
        if (class_subtype == 3) {
            // For this style, the whole interaction contributes only one single basis function.
            for (int i = 1; i < get_n_defined() + 1; i++) interaction_column_indices[i] = i;
        } else {  
            // For other styles, the potential contributes more basis functions.
            if ((get_basis_type() == kBSpline) || (get_basis_type() == kBSplineAndDeriv)) { // Set up a B-spline basis for this interaction.
                for (int i = 1; i < get_n_defined() + 1; i++) {
                	interaction_column_indices[i] = interaction_column_indices[i - 1] 
                		+ i * (get_bspline_k() - 2 + floor(180.0 / get_fm_binwidth() + 0.5) + 1);
                }
            } else if (get_basis_type() == kLinearSpline) { // Set up a linear spline basis for this interaction.
                for (int i = 1; i < get_n_defined() + 1; i++) {
                	interaction_column_indices[i] = interaction_column_indices[i - 1] 
                		+ i * (floor(180.0 / get_fm_binwidth() + 0.5) + 1);
                }
			}
		}
	}
}

// Allocate space for interactions that will be used.

void InteractionClassSpec::setup_for_defined_interactions(TopologyData* topo_data)
{
	n_cg_types = (int)(topo_data->n_cg_types);
    determine_defined_intrxns(topo_data);
	defined_to_matched_intrxn_index_map = std::vector<unsigned>(n_defined, 0);
	defined_to_tabulated_intrxn_index_map = std::vector<unsigned>(n_defined, 0);
	lower_cutoffs = new double[n_defined]();
	upper_cutoffs = new double[n_defined]();
	n_to_force_match = 0;
	n_tabulated = 0;
}

void InteractionClassSpec::dummy_setup_for_defined_interactions(TopologyData* topo_data)
{
	n_defined = 0;
	n_to_force_match = 0;
	n_tabulated = 0;
	lower_cutoffs = new double[1]();
	upper_cutoffs = new double[1]();
}

void three_body_setup_for_defined_interactions(InteractionClassSpec* ispec, TopologyData* topo_data)
{
    // This is equivalent to determine_defined_intrxns functions inside of setup_for_defined_interactions
    // for class_subtype > 0.
    ThreeBodyNonbondedClassSpec* tb_spec = static_cast<ThreeBodyNonbondedClassSpec*>(ispec);
	tb_spec->n_cg_types = topo_data->n_cg_types;
	
	// This is equivelant to the rest of setup_for_defined_interactions.
    if (tb_spec->class_subtype > 0) {
    
    	tb_spec->determine_defined_intrxns(topo_data);
	
	    // Allocate space for the three body nonbonded hash tables analogously to the bonded interactions.
        tb_spec->defined_to_matched_intrxn_index_map = std::vector<unsigned>(tb_spec->get_n_defined(), 0);   
		tb_spec->defined_to_tabulated_intrxn_index_map = std::vector<unsigned>(tb_spec->get_n_defined(), 0);   
        tb_spec->lower_cutoffs = new double[tb_spec->get_n_defined()];
        tb_spec->upper_cutoffs = new double[tb_spec->get_n_defined()];

        // The three body interaction basis functions depend only 
        // on a single angle by default.
        for (int i = 0; i < tb_spec->get_n_defined(); i++) {
            tb_spec->defined_to_matched_intrxn_index_map[i] = i + 1;
			tb_spec->defined_to_tabulated_intrxn_index_map[i] = 0;
            tb_spec->lower_cutoffs[i] = 0.0;
            tb_spec->upper_cutoffs[i] = 180.0;
        }

	} else {
        tb_spec->defined_to_matched_intrxn_index_map = std::vector<unsigned>(1, 0);
        tb_spec->defined_to_tabulated_intrxn_index_map = std::vector<unsigned>(1, 0);  
        tb_spec->lower_cutoffs = new double[1]();
		tb_spec->upper_cutoffs = new double[1]();
		tb_spec->interaction_column_indices = std::vector<unsigned>(1, 0);
    	
    }
}

//-------------------------------------------------------------------------------
// Functions for setting up the potential model that will be used in the CG model
//-------------------------------------------------------------------------------

// Tabulated potential reading error reporting functions

void report_tabulated_interaction_format_error(const int line)
{
    fprintf(stderr, "Wrong format in table.in:line %d!\n", line);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

void report_tabulated_interaction_data_consistency_error(const int line)
{
    fprintf(stderr, "Numbers of tabulated interactions from lower_cutoffs.in/pair_bond_interaction_lower_cutoffs.in and table.in are not consistent:line %d!\n", line);
    fflush(stderr);
    exit(EXIT_FAILURE);
}

void report_fields_error(const std::string &full_name, const int n_expected, const int n_fields)
{
	fprintf(stderr, "This %s interaction requires at least %d entries, but only %d were detected!\n", full_name.c_str(), n_expected, n_fields);
	fflush(stderr);
	exit(EXIT_FAILURE);	
}

void read_all_interaction_ranges(CG_MODEL_DATA* const cg)
{
    // Determine the number of interactions that are actually present in the model for each class of interactions, 
    // allocate a hash array and an index array, then set up the hash array.
    // The index array must be filled in from the range specifications in rmin.in and rmin_b.in.
	std::list<InteractionClassSpec*>::iterator iclass_iterator;
	for(iclass_iterator=cg->iclass_list.begin(); iclass_iterator != cg->iclass_list.end(); iclass_iterator++) {
		(*iclass_iterator)->setup_for_defined_interactions(&cg->topo_data);
	}

    // Now deal with interactions that do not fit the normal scheme.	
	three_body_setup_for_defined_interactions(&cg->three_body_nonbonded_interactions, &cg->topo_data);

    // Read normal range specifications.
    // Open the range files.
    std::ifstream nonbonded_range_in, bonded_range_in;    
   	check_and_open_in_stream(nonbonded_range_in, "rmin.in"); 
	check_and_open_in_stream(bonded_range_in, "rmin_b.in"); 

	// Read the ranges.
	for(iclass_iterator=cg->iclass_list.begin(); iclass_iterator != cg->iclass_list.end(); iclass_iterator++) {
        if ((*iclass_iterator)->n_defined == 0) continue;
        if ((*iclass_iterator)->class_type == kPairNonbonded) {
            (*iclass_iterator)->smart_read_interaction_class_ranges(nonbonded_range_in, cg->name);
		} else {
            (*iclass_iterator)->read_interaction_class_ranges(bonded_range_in);
        }
    }	
	// Close the range files.
	nonbonded_range_in.close();
    bonded_range_in.close();
	
	// Check that specified nonbonded interactions do not extend past the nonbonded cutoff
	check_nonbonded_interaction_range_cutoffs(&cg->pair_nonbonded_interactions, cg->pair_nonbonded_cutoff);
	
    // Allocate space for the column index of each block of basis functions associated with each class of interactions active
    // in the model and meant for force matching, then fill them in class by class.
	for(iclass_iterator=cg->iclass_list.begin(); iclass_iterator != cg->iclass_list.end(); iclass_iterator++) {
		(*iclass_iterator)->setup_indices_in_fm_matrix();
	}
	
	// Now, handle similar actions for non-standard interactions.
	cg->three_body_nonbonded_interactions.setup_indices_in_fm_matrix();	
}

void read_tabulated_interaction_file(CG_MODEL_DATA* const cg, int n_cg_types)
{
    int line = 0;
    std::ifstream external_spline_table;
    check_and_open_in_stream(external_spline_table, "table.in");
	
	std::list<InteractionClassSpec*>::iterator iclass_iterator;
	for(iclass_iterator=cg->iclass_list.begin(); iclass_iterator != cg->iclass_list.end(); iclass_iterator++) {
		line = (*iclass_iterator)->read_table(external_spline_table, line, cg->n_cg_types);
	}
	
	external_spline_table.close();
}

int InteractionClassSpec::read_table(std::ifstream &external_spline_table, int line, int num_cg_types) 
{
    std::string buff;
    char parameter_name[50];
    int hash_val, index_among_defined;
    int n_external_splines_to_read = 0;
    
    // Read the file header.
    check_and_read_next_line(external_spline_table, buff, line);
    sscanf(buff.c_str(), "%s %d %lf", parameter_name, &n_external_splines_to_read, &external_table_spline_binwidth);
    if (strcmp(parameter_name, get_table_name().c_str()) != 0) report_tabulated_interaction_format_error(line);
    if (n_external_splines_to_read != n_tabulated) report_tabulated_interaction_data_consistency_error(line);
    if (n_external_splines_to_read <= 0) return line;
    
    // Read each of the tabulated interactions.
    external_table_spline_coefficients = new double*[n_external_splines_to_read];
    for (int i = 0; i < n_external_splines_to_read; i++) {
        // Read the types of the interaction.
        check_and_read_next_line(external_spline_table, buff, line);
        std::istringstream buffss(buff);
        std::vector<int> types(get_n_body());
        for (unsigned j = 0; j < types.size(); j++) buffss >> types[j];
        // Find it in the defined interactions.
        hash_val = calc_interaction_hash(types, num_cg_types);
        index_among_defined = get_index_from_hash(hash_val);
        // Read the values.
        line = read_bspline_table(external_spline_table, line, index_among_defined);
    }
    return line;
}

int InteractionClassSpec::read_bspline_table(std::ifstream &external_spline_table, int line, int index_among_defined)
{
    std::string buff;
    int index_among_tabulated, n_external_spline_control_points;

    check_and_read_next_line(external_spline_table, buff, line);
    sscanf(buff.c_str(), "%lf%lf", &lower_cutoffs[index_among_defined], &upper_cutoffs[index_among_defined]);
    n_external_spline_control_points = floor((upper_cutoffs[index_among_defined] - lower_cutoffs[index_among_defined]) / external_table_spline_binwidth + 0.5) + 1;
    index_among_tabulated = defined_to_tabulated_intrxn_index_map[index_among_defined] - 1;
    external_table_spline_coefficients[index_among_tabulated] = new double[n_external_spline_control_points];
    for (int j = 0; j < n_external_spline_control_points; j++) {
        check_and_read_next_line(external_spline_table, buff, line);
		sscanf(buff.c_str(), "%lf", &external_table_spline_coefficients[index_among_tabulated][j]);
     }
     return line;
}

void InteractionClassComputer::calc_grid_of_force_vals(const std::vector<double> &spline_coeffs, const int index_among_defined, const double binwidth, std::vector<double> &axis_vals, std::vector<double> &force_vals) 
{
    // Size the output vectors of positions and forces conservatively.
    unsigned num_entries = int( (ispec->upper_cutoffs[index_among_defined] - ispec->lower_cutoffs[index_among_defined]) / binwidth  + 1.0);
    if (num_entries == 0) num_entries = 1;
    axis_vals = std::vector<double>(num_entries);
    force_vals = std::vector<double>(num_entries);    
    // Calculate forces by iterating over the grid points from low to high.
    double min = ((int)(ispec->lower_cutoffs[index_among_defined] / binwidth) + 1) * binwidth;
    double max = ispec->upper_cutoffs[index_among_defined];
    if (min >= max) {
    	fprintf(stderr, "No output will be generated for this interaction since the rounded lower cutoff is greater than or equal to the upper cutoff!\n");
    }
    
    unsigned counter = 0;
    for (double axis = min; axis < max; axis += binwidth) {
        force_vals[counter] = fm_s_comp->evaluate_spline(index_among_defined, interaction_class_column_index, spline_coeffs, axis);
        axis_vals[counter] = axis;
        counter++;
    }
    if (counter == 0) {
    	double axis = (ispec->upper_cutoffs[index_among_defined] + ispec->lower_cutoffs[index_among_defined]) * 0.5;
    	force_vals[counter] = fm_s_comp->evaluate_spline(index_among_defined, interaction_class_column_index, spline_coeffs, axis);
        axis_vals[counter] = axis;
        counter++;
    }
    // Set the correct size for the output vectors of positions and forces.
    axis_vals.resize(counter);
    force_vals.resize(counter);
}

void InteractionClassComputer::calc_grid_of_force_and_deriv_vals(const std::vector<double> &spline_coeffs, const int index_among_defined, const double binwidth, std::vector<double> &axis_vals, std::vector<double> &force_vals, std::vector<double> &deriv_vals)
{
    BSplineAndDerivComputer* s_comp_ptr = static_cast<BSplineAndDerivComputer*>(fm_s_comp);
    unsigned num_entries = int( (ispec->upper_cutoffs[index_among_defined] - ispec->lower_cutoffs[index_among_defined]) / binwidth  + 1.0);
    axis_vals = std::vector<double>(num_entries);
    force_vals = std::vector<double>(num_entries);
	deriv_vals = std::vector<double>(num_entries);
	
    // Calculate forces by iterating over the grid points from low to high.
    double min = ((int)(ispec->lower_cutoffs[index_among_defined] / binwidth) + 1) * binwidth;
    double max = ispec->upper_cutoffs[index_among_defined];
    if (min >= max) {
    	fprintf(stderr, "No output will be generated for this interaction since the rounded lower cutoff is greater than or equal to the upper cutoff!\n");
    }
    
    unsigned counter = 0;
    for (double axis = min; axis < max; axis += binwidth) {
    	axis_vals[counter] = axis;
        force_vals[counter] = s_comp_ptr->evaluate_spline(index_among_defined, interaction_class_column_index, spline_coeffs, axis);
        deriv_vals[counter] = s_comp_ptr->evaluate_spline_deriv(index_among_defined, interaction_class_column_index, spline_coeffs, axis);   
    	counter++;
    }
    
    // Set the correct size for the output vectors of positions and forces.
    axis_vals.resize(counter);
    force_vals.resize(counter);
}
