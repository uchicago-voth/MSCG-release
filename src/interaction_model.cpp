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
int extract_index_among_defined_from_line(InteractionClassSpec* ispec, std::string &buff, const int num_cg_types);
inline void extract_types_and_range(std::string* elements, const InteractionClassType type, const int n_body, std::vector<int> &types, double &low, double &high, char** name, const int n_cg_types, int &index_among_defined);
inline void read_types(const int n_body, std::vector<int> &types, std::string* elements, const int n_types, char** name);

// Non-standard setup functions.
void setup_periodic_index(InteractionClassSpec* iclass);
void three_body_setup_for_defined_interactions(InteractionClassSpec* ispec, TopologyData* topo_data);
void three_body_setup_indices_in_fm_matrix(InteractionClassSpec* ispec);

// Input checking and error reporting functions.
void check_nonbonded_interaction_range_cutoffs(PairNonbondedClassSpec *const ispec, double const cutoff);
void check_angular_interaction_range_cutoffs(AngularClassSpec *const ispec);
void report_tabulated_interaction_format_error(const int line, const std::string& string);
void report_tabulated_interaction_data_consistency_error(const int line);
void report_fields_error(const std::string &full_name, const int n_expected, const int n_fields); 
inline void check_mode(char* mode);


// Get the name of a single defined interaction via its index among
// defined interactions.  Also, add the short name at the end, if it exists.

std::string InteractionClassSpec::get_basename(char **name, const int intrxn_index_among_defined, const std::string &delimiter) 
{
  std::string basename = get_interaction_name(name, intrxn_index_among_defined, "_");
  if (!get_short_name().empty())
  {
    basename += "_" + get_short_name();
  }
  return basename;
}

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

// Select the correct type name array for the interaction.
char** select_name(InteractionClassSpec* const ispec, char ** const cg_name)
{
	return cg_name;
}

void check_input_values(CG_MODEL_DATA* cg)
{
	std::list<InteractionClassSpec*>::iterator iclass_iterator;
	for(iclass_iterator=cg->iclass_list.begin(); iclass_iterator != cg->iclass_list.end(); iclass_iterator++) {
		if ( (*iclass_iterator)->class_subtype < 0 || (*iclass_iterator)->class_subtype > 4) {
			printf("Invalid class_subtype (%d) for %s!\n", (*iclass_iterator)->class_subtype, (*iclass_iterator)->get_full_name().c_str());
			(*iclass_iterator)->class_subtype = 0;
		}
		if ( (*iclass_iterator)->cutoff <= 0.0 ) {
			printf("Invalid cutoff (%lf) for %s!\n", (*iclass_iterator)->cutoff, (*iclass_iterator)->get_full_name().c_str());
			(*iclass_iterator)->cutoff = cg->pair_nonbonded_cutoff;
		}
		if ( (*iclass_iterator)->output_spline_coeffs_flag < 0 || (*iclass_iterator)->output_spline_coeffs_flag > 1) {
			printf("Invalid output_spline_coeffs_flag (%d) for %s!\n", (*iclass_iterator)->output_spline_coeffs_flag, (*iclass_iterator)->get_full_name().c_str());
			(*iclass_iterator)->output_spline_coeffs_flag = 0;
		}
		if ( (*iclass_iterator)->get_fm_binwidth() < 0.0 || (*iclass_iterator)->get_fm_binwidth() > (*iclass_iterator)->cutoff )  {
			printf("Invalid fm_binwidth (%lf) for %s!\n", (*iclass_iterator)->get_fm_binwidth(), (*iclass_iterator)->get_full_name().c_str());
			exit(EXIT_FAILURE);
		}
		if ( (*iclass_iterator)->get_bspline_k() < 2 ) {
			printf("Invalid bspline_k (%d) for %s!\n (Must be at least 3)\n", (*iclass_iterator)->get_bspline_k(), (*iclass_iterator)->get_full_name().c_str());
			exit(EXIT_FAILURE);
		}
		if ( (*iclass_iterator)->output_parameter_distribution < 0 || (*iclass_iterator)->output_parameter_distribution > 2 ) {
			 printf("Invalid output_parameter_distribution (%d) for %s!\n", (*iclass_iterator)->output_parameter_distribution, (*iclass_iterator)->get_full_name().c_str());
			 (*iclass_iterator)->output_parameter_distribution = 0;
		}
	}
	
	if (cg->three_body_nonbonded_interactions.class_subtype < 0 || cg->three_body_nonbonded_interactions.class_subtype > 3) {
		printf("Invalid class_subtype (%d) for %s!\n", cg->three_body_nonbonded_interactions.class_subtype, cg->three_body_nonbonded_interactions.get_full_name().c_str());
		cg->three_body_nonbonded_interactions.class_subtype = 0;
	}
		if ( cg->three_body_nonbonded_interactions.get_fm_binwidth() < 0.0 ) {
		printf("Invalid fm_binwidth (%lf) for %s!\n", cg->three_body_nonbonded_interactions.get_fm_binwidth(), cg->three_body_nonbonded_interactions.get_full_name().c_str());
		exit(EXIT_FAILURE);
	}
	if ( cg->three_body_nonbonded_interactions.get_bspline_k() < 2 ) {
		printf("Invalid bspline_k (%d) for %s!\n (Must be at least 3)\n", cg->three_body_nonbonded_interactions.get_bspline_k(), cg->three_body_nonbonded_interactions.get_full_name().c_str());
		exit(EXIT_FAILURE);
	}
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

// Check that specified angular interactions do not extend past 0 or 180.
void check_angular_interaction_range_cutoffs(AngularClassSpec *const ispec)
{
	for (int i = 0; i < ispec->n_defined; i++) {
		if (ispec->defined_to_matched_intrxn_index_map[i] != 0) {
			if (ispec->upper_cutoffs[i] > 180.0 + VERYSMALL_F || ispec->lower_cutoffs[i] < -VERYSMALL_F) {
				fprintf(stderr, "The interaction range (%lf to %lf) specified in the range file is goes outside of the allowed range (0 to 180).\n", ispec->lower_cutoffs[i], ispec->upper_cutoffs[i]);
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
        }
        if (strcmp(mode,"tab") == 0 || strcmp(mode,"fm+tab") == 0 || strcmp(mode,"tab+fm") == 0) {
			// This interaction is tabulated.
			// Increment the running total of the tabulated interactions.
			total_tabulated++;
			defined_to_tabulated_intrxn_index_map[i] = total_tabulated;
		}	
    }
    n_to_force_match = total_to_fm;
    n_force = total_to_fm;
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
	
	std::getline(range_in, line);
	while(range_in.good() == 1) {
	
		// Check that this line has enough fields.
		if ( (n_fields = StringSplit(line, " \t\n", elements)) < n_expected ) {	//allow for trailing white space
			if(total_intrxns == 0) report_fields_error(get_full_name(), n_expected, n_fields);
			// This allows the reading to terminate if 
			// a line of whitespace is encountered before the end of file.
			break;
		}
		
		// Extract the useful information.
		int index_among_defined;
		read_types(get_n_body(), types, &elements[0], n_cg_types, name);
		index_among_defined = calc_interaction_hash(types, n_cg_types);
		
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
        }
        if (strcmp(mode,"tab") == 0 || strcmp(mode,"fm+tab") == 0 || strcmp(mode,"sym+tab") == 0
         || strcmp(mode,"tab+fm") == 0 || strcmp(mode,"tab+sym") == 0) {
			// This interaction is tabulated.
			// Increment the running total of the tabulated interactions.
			total_tabulated++;
			defined_to_tabulated_intrxn_index_map[index_among_defined] = total_tabulated;
		}		
    	std::getline(range_in, line);
    }
    n_to_force_match = total_to_fm;
    n_force = total_to_fm;
    n_tabulated = total_tabulated;
    printf("Will force match %d %s interactions and %d are interactions tabulated", n_to_force_match, get_full_name().c_str(), n_tabulated);
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

int extract_index_among_defined_from_line(InteractionClassSpec* ispec, std::string &buff, const int num_cg_types)
{
	int hash_val;
	std::vector<int> types(ispec->get_n_body());
	// Find it in the defined interactions.
	std::istringstream buffss(buff);
	for (unsigned j = 0; j < types.size(); j++) buffss >> types[j];
	hash_val = calc_interaction_hash(types, num_cg_types);

	return ispec->get_index_from_hash(hash_val);
}

inline void extract_types_and_range(std::string* elements, const InteractionClassType type, const int n_body, std::vector<int> &types, double &low, double &high, char** name, const int n_cg_types, int &index_among_defined)
{
	// Read the base types and 
	// look-up the index for this sub-interaction
	read_types(n_body, types, &elements[0], n_cg_types, name);
	index_among_defined = calc_interaction_hash(types, n_cg_types);
	
	// Read the low and high parameter values;
	low  = atof(elements[n_body].c_str());
	high = atof(elements[n_body + 1].c_str());
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
    if (strcmp(mode,"none") != 0 && strcmp(mode,"fm") != 0 && strcmp(mode,"tab") != 0 && strcmp(mode,"fm+tab") != 0 &&
        strcmp(mode,"tab+fm") != 0){
        fprintf(stderr, "Interaction mode %s is not recognized\n", mode);
        fflush(stderr);
    	exit(EXIT_FAILURE);
	}
}

void setup_periodic_index(InteractionClassSpec* iclass) 
{
	// Setup periodic flags for dihedral interactions based on upper and lower cutoff values.
	// Only do this for angle-based dihedrals (not distance-based dihedrals).
	if (iclass->class_subtype == 1) return;
	for (int i = 0; i < iclass->n_defined; i++) {
		// Only for matched interactions
		if (iclass->defined_to_matched_intrxn_index_map[i] == 0) continue;
		
		// Check if this interaction meets the criteria for periodic interactions.
		if (iclass->upper_cutoffs[i] > 179.0 && iclass->lower_cutoffs[i] < -179.0) {
			// This only looks for an interaction that covers the whole range.
			// It allows some flexibility in values near the edges since this will also get called for BI calculation directly from rangerinder.
			printf("Treating dihedral interaction range (%lf to %lf) as fully periodic (-180 to 180).\n", iclass->lower_cutoffs[i], iclass->upper_cutoffs[i]);
			iclass->defined_to_periodic_intrxn_index_map[i] = 1;
			iclass->upper_cutoffs[i] = 180.0;
			iclass->lower_cutoffs[i] = -180.0;
		} else if (iclass->upper_cutoffs[i] < iclass->lower_cutoffs[i] || iclass->upper_cutoffs[i] > 180.0 + VERYSMALL_F || iclass->lower_cutoffs[i] + VERYSMALL_F < -180.0) {
			// This looks for a range that passes through the periodic boundary.
			// For example, this can be specified as either +140 to -140 or +140 to +220.
			printf("Treating dihedral interaction range (%lf to %lf) as going through periodic wrapping.\n", iclass->lower_cutoffs[i], iclass->upper_cutoffs[i]);
			iclass->defined_to_periodic_intrxn_index_map[i] = 2;	
			
			if (iclass->lower_cutoffs[i] < -180.0) {
				iclass->lower_cutoffs[i] += 360.0;
				iclass->upper_cutoffs[i] += 360.0;
			}
			
			if (iclass->upper_cutoffs[i] < 180.0) {
				iclass->upper_cutoffs[i] += 360.0;
			}
			if (iclass->upper_cutoffs[i] - iclass->lower_cutoffs[i] > 350.0 - VERYSMALL_F) {
				printf("This dihedral interaction is nearly the entire interaction range!\n");
				printf("Please change the interaction range to -180 to 180 to take advantage of periodic features.\n");
				exit(EXIT_FAILURE);
			}
		}
		// Make sure that these cutoff values are within the right range.
		if (iclass->lower_cutoffs[i] < -180.0 - VERYSMALL_F || iclass->lower_cutoffs[i] > 540.0 ||
				   iclass->upper_cutoffs[i] > 540.0 + VERYSMALL_F) {
			printf("Invalid interaction range for dihedral angle (%lf to %lf)\n", iclass->lower_cutoffs[i], iclass->upper_cutoffs[i]);
			if (iclass->lower_cutoffs[i] < -180.0 - VERYSMALL_F ) {
				printf("Adjusting lower cutoff to -180.\n");
				iclass->lower_cutoffs[i] = -180.0;
			} else if (iclass->upper_cutoffs[i] > 540.0 && iclass->defined_to_periodic_intrxn_index_map[i] == 2) {
				printf("Adjusting unwrapped upper cutoff to 540.0\n");
				iclass->upper_cutoffs[i] = 540.0;
			} else {			
				exit(EXIT_FAILURE);
			}
		}
	}
}

void InteractionClassSpec::adjust_cutoffs_for_basis(int i)
{
   if ((basis_type == kLinearSpline) || (basis_type == kBSpline) ||  (basis_type == kBSplineAndDeriv)) {
    	// Do not adjust upper cutoff if the pair nonbonded interaction is already at the user-defined cutoff in control.in
    	// Otherwise, adjust so that the upper cutoff is divisible by the output binwidth.
        if (!(class_type == kPairNonbonded && fabs(upper_cutoffs[i] - cutoff)) < VERYSMALL_F) {
         upper_cutoffs[i] = floor( (upper_cutoffs[i] / output_binwidth) + 0.5 ) * output_binwidth;
        }
        // Now, round down the lower cutoff so that there is an integer number of a bin.
        lower_cutoffs[i] = upper_cutoffs[i] - floor( ((upper_cutoffs[i] - lower_cutoffs[i]) / fm_binwidth) + 0.5 ) * fm_binwidth;
        if ((class_type != kDihedralBonded) && (lower_cutoffs[i] < 0.0)) lower_cutoffs[i] = 0.0;
    } else {
    	// Nothing to be done here.
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
			if ((basis_type == kBSpline) || (basis_type == kBSplineAndDeriv)) {
				interaction_column_indices[counter + 1] -= 2;
				if(defined_to_periodic_intrxn_index_map[i] != 1) interaction_column_indices[counter + 1] += get_bspline_k();
			}
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
                		+ i * (get_bspline_k() - 2 + floor(180.0 / get_fm_binwidth() + 0.5));
                }
            } else if (get_basis_type() == kLinearSpline) { // Set up a linear spline basis for this interaction.
                for (int i = 1; i < get_n_defined() + 1; i++) {
                	interaction_column_indices[i] = interaction_column_indices[i - 1] 
                		+ i * (floor(180.0 / get_fm_binwidth() + 0.5));
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
	defined_to_periodic_intrxn_index_map = std::vector<unsigned>(n_defined, 0);
	lower_cutoffs = new double[n_defined]();
	upper_cutoffs = new double[n_defined]();
	n_to_force_match = 0;
	n_force = 0;
	n_tabulated = 0;
}

void InteractionClassSpec::dummy_setup_for_defined_interactions(TopologyData* topo_data)
{
	n_defined = 0;
	n_to_force_match = 0;
	n_force = 0;
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
        tb_spec->defined_to_periodic_intrxn_index_map = std::vector<unsigned>(tb_spec->get_n_defined(), 0);
        tb_spec->lower_cutoffs = new double[tb_spec->get_n_defined()]();
        tb_spec->upper_cutoffs = new double[tb_spec->get_n_defined()];

        // The three body interaction basis functions depend only 
        // on a single angle by default.
        for (int i = 0; i < tb_spec->get_n_defined(); i++) {
            tb_spec->defined_to_matched_intrxn_index_map[i] = i + 1;
            tb_spec->lower_cutoffs[i] = 0.0;
            tb_spec->upper_cutoffs[i] = 180.0;
        }

	} else {
        tb_spec->defined_to_matched_intrxn_index_map = std::vector<unsigned>(1, 0);
        tb_spec->defined_to_tabulated_intrxn_index_map = std::vector<unsigned>(1, 0);
        tb_spec->defined_to_periodic_intrxn_index_map = std::vector<unsigned>(1, 0);
        tb_spec->lower_cutoffs = new double[1]();
		tb_spec->upper_cutoffs = new double[1]();
		tb_spec->interaction_column_indices = std::vector<unsigned>(1, 0);
    	
    }
}

//--------------------------------------------------------------------
// Functions for setting up the potential model that will be used
// in the CG model from a range.in file using lots of hard-to-maintain
// tricks like 'if the range starts from -1, don't include this
// interaction in the model'
//--------------------------------------------------------------------

// Tabulated potential reading error reporting functions

void report_tabulated_interaction_format_error(const int line, const std::string& string)
{
    fprintf(stderr, "Wrong format in table.in:line %d! Expected keyword %s\n", line, string.c_str());
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

void reset_interaction_cutoff_arrays(CG_MODEL_DATA* const cg)
{
	// Reset upper and lower cutoffs for all interactions before starting BI calculation
	std::list<InteractionClassSpec*>::iterator iclass_iterator;
	for(iclass_iterator=cg->iclass_list.begin(); iclass_iterator != cg->iclass_list.end(); iclass_iterator++) {
		delete [] (*iclass_iterator)->upper_cutoffs;
		delete [] (*iclass_iterator)->lower_cutoffs;
	}
	delete [] cg->three_body_nonbonded_interactions.upper_cutoffs;
	delete [] cg->three_body_nonbonded_interactions.lower_cutoffs;
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
        char** name = select_name((*iclass_iterator), cg->name);
	
        if ((*iclass_iterator)->class_type == kPairNonbonded) {
            (*iclass_iterator)->smart_read_interaction_class_ranges(nonbonded_range_in, name);
		} else {
            (*iclass_iterator)->read_interaction_class_ranges(bonded_range_in);
        }
    }	
	// Close the range files.
	nonbonded_range_in.close();
    bonded_range_in.close();
		
	// Check that specified nonbonded interactions do not extend past the nonbonded cutoff
	check_nonbonded_interaction_range_cutoffs(&cg->pair_nonbonded_interactions, cg->pair_nonbonded_cutoff);
	// Check that specified nonbonded interactions do not extend past the nonbonded cutoff
	check_angular_interaction_range_cutoffs(&cg->angular_interactions);
	// Also, setup periodic flags for dihedral interactions based on upper and lower cutoff values.
	setup_periodic_index(&cg->dihedral_interactions);
	
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
    int index_among_defined;
    int n_external_splines_to_read = 0;
    
    // Read the file header.
    check_and_read_next_line(external_spline_table, buff, line);
    sscanf(buff.c_str(), "%s %d %lf", parameter_name, &n_external_splines_to_read, &external_table_spline_binwidth);
    if (strcmp(parameter_name, get_table_name().c_str()) != 0) report_tabulated_interaction_format_error(line, get_table_name());
    if (n_external_splines_to_read != n_tabulated) report_tabulated_interaction_data_consistency_error(line);
    if (n_external_splines_to_read <= 0) return line;
    
    // Read each of the tabulated interactions.
    external_table_spline_coefficients = new double*[n_external_splines_to_read];
    for (int i = 0; i < n_external_splines_to_read; i++) {
        // Read the types of the interaction.
        check_and_read_next_line(external_spline_table, buff, line);
        std::size_t pos = buff.find("types");
        if (pos != std::string::npos) { // this uses the "new" multiple types for same table format
    		int offset, prev_line;
    		std::vector<int> index_list(1);
    		std::string remainder = buff.substr(pos + 5);
    		index_list[0] = extract_index_among_defined_from_line(this, remainder, num_cg_types);
			while (1) {
				// save this line position 		
				offset = external_spline_table.tellg();
				check_and_read_next_line(external_spline_table, buff, line);
    			pos = buff.find("types"); // each line of types must begin with the "types" tag
    			if (pos != std::string::npos) { 
    				// add this to the vector and keep going
    				remainder = buff.substr(pos + 5);
    				index_list.push_back(extract_index_among_defined_from_line(this, remainder, num_cg_types));
    			} else {
    				// reset the file stream back one line and exit the loop
    				external_spline_table.clear();
    				external_spline_table.seekg(offset);
    				line--;
    				break;
    			}
    		}
    		
    		// read the table
    		prev_line = line;
    		line = read_bspline_table(external_spline_table, line, index_list[0]);
			
			// copy the table to all of the other defined interactions in index_list
			// also, increment the number of "tables read" (i)
			for (unsigned j = 1; j < index_list.size(); j++) {
				copy_table(index_list[0], index_list[j], line - prev_line);
				i++;
			}
			
    	} else { // process as single type per table ("old" format)
    		// determine the hash and index for this set of types
			index_among_defined = extract_index_among_defined_from_line(this, buff, num_cg_types);
			// Read the values.
			line = read_bspline_table(external_spline_table, line, index_among_defined);
		}
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

void InteractionClassSpec::copy_table(const int base_defined, const int target_defined, const int num_lines)
{
	lower_cutoffs[target_defined] = lower_cutoffs[base_defined];
	upper_cutoffs[target_defined] = upper_cutoffs[base_defined];
	
	int base_tab    = defined_to_tabulated_intrxn_index_map[base_defined  ] - 1;
	int target_tab  = defined_to_tabulated_intrxn_index_map[target_defined] - 1;
	
	int n_external_spline_control_points = num_lines - 1;
	external_table_spline_coefficients[target_tab] = new double[n_external_spline_control_points];
	for (int j = 0; j < n_external_spline_control_points; j++) {
		external_table_spline_coefficients[target_tab][j] = external_table_spline_coefficients[base_tab][j];
	}
}

void InteractionClassComputer::calc_grid_of_table_force_vals(const int index_among_defined, const double binwidth, std::vector<double> &axis_vals, std::vector<double> &force_vals) 
{
    // Calculate forces by iterating over the grid points from low to high.
    // The lower value is adjusted so that the difference between the output upper cutoff and output lower cutoff is divisible by the output binwidth and that the lower cutoff is always greater than basis lower cutoff.
    double max = ispec->upper_cutoffs[index_among_defined];
    double min = max - ((int)((max - ispec->lower_cutoffs[index_among_defined]) / binwidth)) * binwidth;

    // Size the output vectors of positions and forces conservatively.
    unsigned num_entries = int((max - min)/binwidth) + 2;
    if (num_entries <= 0) { 
    	num_entries = 1;
    	fprintf(stderr, "No output will be generated for this interaction since the rounded lower cutoff is greater than or equal to the upper cutoff!\n");
	}
	
    axis_vals = std::vector<double>(num_entries);
    force_vals = std::vector<double>(num_entries);    
    
    unsigned counter = 0;
    std::vector<double> junk_vector; // The external spline coefficients are grabed internally through evaluate_spline and calculate_basis_fn_vals.
    								 // So, they don't need to be converted from double* to std::vector<double> to be passed in here. (junk instead).
    for (double axis = min; axis <= max+VERYSMALL_F; axis += binwidth) {
        force_vals[counter] = table_s_comp->evaluate_spline(index_among_defined, interaction_class_column_index, junk_vector, axis);
        axis_vals[counter] = axis;
        counter++;
    }
    if (counter == 0) {
    	double axis = (ispec->upper_cutoffs[index_among_defined] + ispec->lower_cutoffs[index_among_defined]) * 0.5;
    	force_vals[counter] = table_s_comp->evaluate_spline(index_among_defined, interaction_class_column_index, junk_vector, axis);
        axis_vals[counter] = axis;
        counter++;
    }
    // Set the correct size for the output vectors of positions and forces.
    if (counter > num_entries) printf("Grid evaluation vector was undersized\n"); 
    else {
    	axis_vals.resize(counter);
    	force_vals.resize(counter);
    }
}

void InteractionClassComputer::calc_grid_of_force_vals(const std::vector<double> &spline_coeffs, const int index_among_defined, const double binwidth, std::vector<double> &axis_vals, std::vector<double> &force_vals) 
{
    // Calculate forces by iterating over the grid points from low to high.
    // The lower value is adjusted so that the difference between the output upper cutoff and output lower cutoff is divisible by the output binwidth and that the lower cutoff is always greater than basis lower cutoff.
    double max = ispec->upper_cutoffs[index_among_defined];
    double min = max - ((int)((max - ispec->lower_cutoffs[index_among_defined]) / binwidth)) * binwidth;
    
    // Size the output vectors of positions and forces conservatively.
    unsigned num_entries = int((max - min)/binwidth) + 2;
    if (num_entries <= 0) { 
    	num_entries = 1;
    	fprintf(stderr, "No output will be generated for this interaction since the rounded lower cutoff is greater than or equal to the upper cutoff!\n");
	}
    
    axis_vals = std::vector<double>(num_entries);
    force_vals = std::vector<double>(num_entries);    
    unsigned counter = 0;
    for (double axis = min; axis <= max + VERYSMALL_F; axis += binwidth) {
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
    if (counter > num_entries) printf("Grid evaluation vector was undersized\n"); 
    else {
    	axis_vals.resize(counter);
    	force_vals.resize(counter);
	}
}

void InteractionClassComputer::calc_grid_of_force_and_deriv_vals(const std::vector<double> &spline_coeffs, const int index_among_defined, const double binwidth, std::vector<double> &axis_vals, std::vector<double> &force_vals, std::vector<double> &deriv_vals)
{
    BSplineAndDerivComputer* s_comp_ptr = static_cast<BSplineAndDerivComputer*>(fm_s_comp);	
    
    // Calculate forces by iterating over the grid points from low to high.
    // The lower value is adjusted so that the difference between the output upper cutoff and output lower cutoff is divisible by the output binwidth and that the lower cutoff is always greater than basis lower cutoff.
	double max = ispec->upper_cutoffs[index_among_defined];
    double min = max - ((int)((max - ispec->lower_cutoffs[index_among_defined]) / binwidth)) * binwidth;
    
    // Size the output vectors of positions and forces conservatively.
    unsigned num_entries = int((max - min)/binwidth) + 2;
    if (num_entries <= 0) { 
    	num_entries = 1;
    	fprintf(stderr, "No output will be generated for this interaction since the rounded lower cutoff is greater than or equal to the upper cutoff!\n");
	}
	
	axis_vals = std::vector<double>(num_entries);
    force_vals = std::vector<double>(num_entries);
    deriv_vals = std::vector<double>(num_entries);

    unsigned counter = 0;
    for (double axis = min; axis <= max+VERYSMALL_F; axis += binwidth) {
    	axis_vals[counter] = axis;
   		force_vals[counter] = s_comp_ptr->evaluate_spline(index_among_defined, interaction_class_column_index, spline_coeffs, axis);
        deriv_vals[counter] = s_comp_ptr->evaluate_spline_deriv(index_among_defined, interaction_class_column_index, spline_coeffs, axis);   
    	counter++;
    }
    
    // Set the correct size for the output vectors of positions and forces.
    if (counter > num_entries) printf("Grid evaluation vector was undersized\n"); 
    else {
    	axis_vals.resize(counter);
    	force_vals.resize(counter);
    }
}