//
//  fm_output.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include "interaction_model.h"
#include "interaction_hashing.h"
#include "matrix.h"
#include "misc.h"

#include "fm_output.h"

//----------------------------------------------------------------------------
// Prototypes for private implementation routines.
//----------------------------------------------------------------------------

void write_interaction_data_to_file(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat);
void write_three_body_interaction_data(ThreeBodyNonbondedClassComputer* const icomp, MATRIX_DATA* const mat, char ** const name);

void pad_and_print_table_files(const char char_id, const std::string& basename, std::vector<double>& axis_vals, std::vector<double>& force_vals, std::vector<double>& potential_vals, const double cutoff);
void pad_and_print_single_table(const char char_id, const std::string& basename, std::vector<double>& axis_vals, std::vector<double>& force_vals, const double cutoff);
void print_table_files(const char char_id, std::string& basename, std::vector<double>& axis_vals, std::vector<double>& force_vals, std::vector<double>& potential_vals);

void write_one_param_table_files(InteractionClassComputer* const icomp, char ** const name, const std::vector<double> &spline_coeffs, const int index_among_defined_intrxns, const double cutoff);
void write_one_param_table_files_energy(InteractionClassComputer* const icomp, char ** const name, const std::vector<double> &spline_coeffs, const int index_among_defined_intrxns, const double cutoff);
void write_two_param_bspline_table_file(InteractionClassComputer* const icomp, char ** const name, MATRIX_DATA* const mat, const int index_among_defined);

void write_one_param_linear_spline_file(InteractionClassComputer* const icomp, char ** const name, MATRIX_DATA* const mat, const int index_among_defined_intrxns);
void write_one_param_bspline_file(InteractionClassComputer* const icomp, char ** const name, MATRIX_DATA* const mat, const int index_among_defined);
void write_output_solution(MATRIX_DATA* const mat);

void write_MSCGFM_table_output_file(const std::string& filename_base, const std::vector<double>& axis, const std::vector<double>& force);
void write_LAMMPS_table_output_file(const char i_type, const std::string& interaction_name, const int num_entries, double* axis_values, double* potential_vals, double* force_vals, int min_index, double min_val);
void write_LAMMPS_table_output_file(const char i_type, const std::string& interaction_name, std::vector<double>& axis_vals, std::vector<double>& potential_vals, std::vector<double>& force_vals);
void write_full_bootstrapping_MSCGFM_table_output_file(const std::string& filename_base, const std::vector<double>& axis, std::vector<double> const master_force, std::vector<double>* const force, int const bootstrapping_num_estimates);
void write_bootstrapping_MSCGFM_table_output_file(const std::string& filename_base, const std::vector<double>& axis, std::vector<double> const master_force, std::vector<double>* const force, int const bootstrapping_num_estimates);

void write_bootstrapping_one_param_table_files(InteractionClassComputer* const icomp, char **name, std::vector<double> const master_coeffs, std::vector<double>* const spline_coeffs, const int index_among_defined_intrxns, const int bootstrapping_num_estimates, const int bootstrapping_full_output_flag);
void write_bootstrapping_one_param_table_files_energy(InteractionClassComputer* const icomp, char **name, std::vector<double> const master_coeffs, std::vector<double>* const spline_coeffs, const int index_among_defined_intrxns, const int bootstrapping_num_estimates, const int bootstrapping_full_output_flag, const double cutoff);
void write_bootstrapping_one_param_bspline_file(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, const int index_among_defined_intrxns);
void write_bootstrapping_one_param_linear_spline_file(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, const int index_among_defined_intrxns);
std::vector<double> calculate_bootstrapping_standard_error(const std::vector<double> master_vals, const std::vector<double>* vals, const int bootstrapping_num_estimates);

//------------------------------------------------------------------------
//    Implementation
//------------------------------------------------------------------------

void write_fm_interaction_output_files(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat)
{

    // Write a binary copy of the solution vector if desired.
    if (mat->output_solution_flag == 1) {
    	write_output_solution(mat);
    }

    // Write all interaction-by-interaction output files.
    write_interaction_data_to_file(cg, mat);
}

void write_interaction_data_to_file(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat)
{
    // Erase the "b-spline.out" file if it exists, create it empty if not.
    // Move this to a b-spline output constructor.
    FILE *spline_output_filep = open_file("b-spline.out", "w");
    fclose(spline_output_filep);
    

    // For each class of interactions, perform output for the active 
    // interactions in that class. For one-parameter interactions right now.
	std::list<InteractionClassComputer*>::iterator icomp_iterator;
	for(icomp_iterator = cg->icomp_list.begin(); icomp_iterator != cg->icomp_list.end(); icomp_iterator++) {
        // For every defined interaction,
        for (unsigned i = 0; i < (*icomp_iterator)->ispec->defined_to_matched_intrxn_index_map.size(); i++) {
            // If that interaction is being matched (includes forces and symmetric/DOOM),
            if ((*icomp_iterator)->ispec->defined_to_matched_intrxn_index_map[i] != 0) {
            	    
            	// Select the correct type name array for the interaction.
				char** name = select_name((*icomp_iterator)->ispec, cg->name);

				// Write output based on energy splines
	        	if(mat->matrix_type == kREM || mat->matrix_type == kDummy){
	      	 		if(mat->matrix_type == kDummy && (*icomp_iterator)->ispec->output_parameter_distribution == 0) continue;
		         	
		         	if (mat->bootstrapping_flag == 1) {
	                	// Write tabular output, regardless of spline type.
 	   	            	write_bootstrapping_one_param_table_files_energy(*icomp_iterator, name, mat->fm_solution, mat->bootstrap_solutions, i, mat->bootstrapping_num_estimates, mat->bootstrapping_full_output_flag, cg->pair_nonbonded_cutoff);
        	        	// Write special output files for the specific spline types.
            	    	if ((*icomp_iterator)->ispec->get_basis_type() == kBSpline ||
            	    	    (*icomp_iterator)->ispec->get_basis_type() == kBSplineAndDeriv ) {
                		    write_bootstrapping_one_param_bspline_file(*icomp_iterator, name, mat, i);
	                	} else if ((*icomp_iterator)->ispec->get_basis_type() == kLinearSpline ||
	                			   (*icomp_iterator)->ispec->get_basis_type() == kDelta) {
    	            	    write_bootstrapping_one_param_linear_spline_file(*icomp_iterator, name, mat, i);
        	        	} else {
            	    	    printf("Unrecognized basis type.\n");
                		    exit(EXIT_FAILURE);
                		}
   					} else {
						// Write tabular output, regardless of spline type.
						write_one_param_table_files_energy(*icomp_iterator, name, mat->fm_solution, i, cg->pair_nonbonded_cutoff);
					
						// Write special output files for the specific spline types.	      
						if ((*icomp_iterator)->ispec->get_basis_type() == kBSpline ||
							(*icomp_iterator)->ispec->get_basis_type() == kBSplineAndDeriv ) {
							write_one_param_bspline_file(*icomp_iterator, name, mat, i);
						} else if ((*icomp_iterator)->ispec->get_basis_type() == kLinearSpline ||
						   (*icomp_iterator)->ispec->get_basis_type() == kDelta) {
							write_one_param_linear_spline_file(*icomp_iterator, name, mat, i);
						} else {
							printf("Unrecognized basis type.\n");
							exit(EXIT_FAILURE);
						}
					}
	      		} else {		
	      		// This is for force matching, MS-CODE (sym), and relative entropy of framewise observables (kObs)
	      		// Write output based on force splines		
   					if (mat->bootstrapping_flag == 1) {
	                	// Write tabular output, regardless of spline type.
 	   	            	write_bootstrapping_one_param_table_files(*icomp_iterator, name, mat->fm_solution, mat->bootstrap_solutions, i, mat->bootstrapping_num_estimates, mat->bootstrapping_full_output_flag);
        	        	// Write special output files for the specific spline types.
            	    	if ((*icomp_iterator)->ispec->get_basis_type() == kBSpline ||
            	    	    (*icomp_iterator)->ispec->get_basis_type() == kBSplineAndDeriv ) {
                		    write_bootstrapping_one_param_bspline_file(*icomp_iterator, name, mat, i);
	                	} else if ((*icomp_iterator)->ispec->get_basis_type() == kLinearSpline ||
	                			   (*icomp_iterator)->ispec->get_basis_type() == kDelta) {
    	            	    write_bootstrapping_one_param_linear_spline_file(*icomp_iterator, name, mat, i);
        	        	} else {
            	    	    printf("Unrecognized basis type.\n");
                		    exit(EXIT_FAILURE);
                		}
   					} else {
	                	// Write tabular output, regardless of spline type.
				  		write_one_param_table_files(*icomp_iterator, name, mat->fm_solution, i, cg->pair_nonbonded_cutoff);
        	        	// Write special output files for the specific spline types.
            	    	if ((*icomp_iterator)->ispec->get_basis_type() == kBSpline ||
            	    	    (*icomp_iterator)->ispec->get_basis_type() == kBSplineAndDeriv ) {
                		    write_one_param_bspline_file(*icomp_iterator, name, mat, i);
		    			} else if ((*icomp_iterator)->ispec->get_basis_type() == kLinearSpline ||
	                			   (*icomp_iterator)->ispec->get_basis_type() == kDelta) {
    	                	write_one_param_linear_spline_file(*icomp_iterator, name, mat, i);
    	            	} else {
            	        	printf("Unrecognized basis type.\n");
                	    	exit(EXIT_FAILURE);
                		}
					}
	      		}
	    	}
        }
	}
      
    // Write three body nonbonded interaction data.
	write_three_body_interaction_data(&cg->three_body_nonbonded_computer, mat, cg->name);
	
	printf("Done with output.\n"); fflush(stdout);

    for (int i = 0; i < cg->n_cg_types; i++) delete [] cg->name[i];
    delete [] cg->name;
}

// Write output for three-body non-bonded interactions.

void write_three_body_interaction_data(ThreeBodyNonbondedClassComputer* const icomp, MATRIX_DATA* const mat, char ** const name)
{
	if (icomp->ispec->class_subtype <= 0) return;
	
	ThreeBodyNonbondedClassSpec* iclass = static_cast<ThreeBodyNonbondedClassSpec*>(icomp->ispec);	
	FILE* tb_out = NULL;	
	if (iclass->class_subtype == 3) tb_out = open_file("3b.dat", "w");
	
	// For every defined interaction,
	for (unsigned i = 0; i < iclass->defined_to_matched_intrxn_index_map.size(); i++) {
		// If that interaction is being matched,
		if (iclass->defined_to_matched_intrxn_index_map[i] != 0) {
			if (iclass->class_subtype == 3) {
				int index_among_matched = iclass->defined_to_matched_intrxn_index_map[i];
				fprintf(tb_out, "%.15le\n", mat->fm_solution[icomp->interaction_class_column_index + iclass->interaction_column_indices[index_among_matched - 1]]);
			} else {
				if(iclass->get_basis_type() == kBSpline || iclass->get_basis_type() == kBSplineAndDeriv) {
					write_one_param_bspline_file(icomp, name, mat, i);
					write_two_param_bspline_table_file(icomp, name, mat, i);
				} else {
					write_one_param_linear_spline_file(icomp, name, mat, i);
				}
			}
		}
	}
		
	if (iclass->class_subtype == 3) fclose(tb_out);
}

void write_MSCGFM_table_output_file(const std::string& filename_base, const std::vector<double>& axis, const std::vector<double>& force) 
{
	std::string filename_tmp = filename_base + ".dat";
    FILE *curr_table_output_file = open_file(filename_tmp.c_str(), "w");
    for (unsigned i = 0; i < axis.size(); i++) {
        fprintf(curr_table_output_file, "%lf %.15le\n", axis[i], force[i]);
    }
    fclose(curr_table_output_file);
}

void write_LAMMPS_table_output_file(const char i_type, const std::string& interaction_name, std::vector<double>& axis_vals, std::vector<double>& potential_vals, std::vector<double>& force_vals) 
{
	// Set-up LAMMPS table file
	std::string filename = interaction_name + ".table";
	FILE* curr_table_output_file = open_file(filename.c_str(), "w");

	// Write header
	fprintf(curr_table_output_file, "# Header information on force file\n");
	fprintf(curr_table_output_file, "\n");
	fprintf(curr_table_output_file, "%s\n", interaction_name.c_str());
	
    // Adjust bonded interactions so that the minimum potential is at 0.0.
    if( (i_type == 'b') || (i_type == 'a') || (i_type == 'd') ) {
        standardize_potential(potential_vals);
    }

    // Write special header lines for specific interaction types.
	if ( (i_type == 'b') || (i_type == 'a') ) {
		// bond and angle
        int min_index = get_min_index(potential_vals);
		fprintf(curr_table_output_file, "N %lu FP 0.0 0.0 EQ %lf\n", axis_vals.size(), axis_vals[min_index]);
	} else if (i_type == 'd') { // dihedral
		fprintf(curr_table_output_file, "N %lu DEGREES\n", axis_vals.size());
	} else { // pair non-bonded, density, RG
		fprintf(curr_table_output_file, "N %lu R %lf %lf\n", axis_vals.size(), axis_vals[0], axis_vals[axis_vals.size() - 1]);
	}

	fprintf(curr_table_output_file, "\n");

	// Write body of table.
	for(unsigned k = 0; k < axis_vals.size(); k++)
	{
		fprintf(curr_table_output_file, "%d %lf %lf %lf\n", (k+1), axis_vals[k], potential_vals[k], force_vals[k]);
	}
	fclose(curr_table_output_file);
}

// Write the tabular output for a single interaction.
void write_one_param_table_files_energy(InteractionClassComputer* const icomp, char ** const name, const std::vector<double> &spline_coeffs, const int index_among_defined_intrxns, const double cutoff) 
{
  std::vector<double> axis_vals, force_vals, potential_vals;
  if (dynamic_cast<OneBodyClassComputer*>(icomp) != NULL)
    {
      icomp->calc_one_force_val(spline_coeffs, index_among_defined_intrxns, icomp->ispec->output_binwidth, axis_vals, force_vals, potential_vals);
    } 
  else
    {
      icomp->calc_grid_of_force_and_deriv_vals(spline_coeffs, index_among_defined_intrxns, icomp->ispec->output_binwidth, axis_vals, potential_vals, force_vals);
      make_negative(force_vals);
    }
    // Determine base for output filenames.
    std::string basename;
    if(icomp->ispec->class_type == kDensity)
      {
        DensityClassSpec* ispec = static_cast<DensityClassSpec*>(icomp->ispec);
        basename = ispec->get_basename(ispec->density_group_names, index_among_defined_intrxns, "_");
      }
    else
      {
        basename = icomp->ispec->get_basename(name, index_among_defined_intrxns, "_");
      }

    // Select symmetric basename modifier if appropriate (i.e. DOOM interactions)
   if(icomp->ispec->defined_to_symmetric_intrxn_index_map[index_among_defined_intrxns] != 0)
     {
        basename += "_sym";	
     }
	
	// Print out tabulated output files in MSCGFM style and LAMMPS style.
    write_MSCGFM_table_output_file(basename, axis_vals, potential_vals);
	pad_and_print_table_files(icomp->ispec->get_char_id(), basename, axis_vals, force_vals, potential_vals, cutoff);
}
				 
void write_one_param_table_files(InteractionClassComputer* const icomp, char ** const name, const std::vector<double> &spline_coeffs, const int index_among_defined_intrxns, const double cutoff) 
{	
    // Compute forces over a grid of parameter values.
    std::vector<double> axis_vals, force_vals, potential_vals;
    int index_among_tabulated = icomp->ispec->defined_to_tabulated_intrxn_index_map[index_among_defined_intrxns];
    if (dynamic_cast<OneBodyClassComputer*>(icomp) != NULL)
      {
    	icomp->calc_one_force_val(spline_coeffs, index_among_defined_intrxns, icomp->ispec->output_binwidth, axis_vals, force_vals, potential_vals);
      } 
    else
      {
    	icomp->calc_grid_of_force_vals(spline_coeffs, index_among_defined_intrxns, icomp->ispec->output_binwidth, axis_vals, force_vals);
    	// Integrate force starting from cutoff = 0.0 potential.
    	integrate_force(axis_vals, force_vals, potential_vals);
      }
    
    // Determine base for output filenames.
    std::string basename;
    if(icomp->ispec->class_type == kDensity)
      {
        DensityClassSpec* ispec = static_cast<DensityClassSpec*>(icomp->ispec);
        basename = ispec->get_basename(ispec->density_group_names, index_among_defined_intrxns, "_");
      }
    else
      {
        basename = icomp->ispec->get_basename(name, index_among_defined_intrxns, "_");
      }
    // Select symmetric basename modifier if appropriate (i.e. DOOM interactions)
   if(icomp->ispec->defined_to_symmetric_intrxn_index_map[index_among_defined_intrxns] != 0)
     {
        basename += "_sym";	
     }
	
	// Print out tabulated output files in MSCGFM style and LAMMPS style.
    write_MSCGFM_table_output_file(basename, axis_vals, force_vals);
    if( index_among_tabulated == 0) {
    	pad_and_print_table_files(icomp->ispec->get_char_id(), basename, axis_vals, force_vals, potential_vals, cutoff);
    } else {
    	print_table_files(icomp->ispec->get_char_id(), basename, axis_vals, force_vals, potential_vals);
    	if (icomp->ispec->class_type != kOneBody) {
		    std::vector<double> tab_axis_vals, tab_force_vals;
		    icomp->calc_grid_of_table_force_vals(index_among_defined_intrxns, icomp->ispec->output_binwidth, tab_axis_vals, tab_force_vals);
    		// sum the forces between the fm and tab.
    		add_force_vals(axis_vals, force_vals, tab_axis_vals, tab_force_vals);
    		// Integrate force starting from cutoff
    		integrate_force(axis_vals, force_vals, potential_vals);
    		
    		write_MSCGFM_table_output_file(basename + "_sum", axis_vals, force_vals);
			pad_and_print_table_files(icomp->ispec->get_char_id(), basename + "_sum", axis_vals, force_vals, potential_vals, cutoff);
    	}
    }
}

void pad_and_print_table_files(const char char_id, const std::string& basename, std::vector<double>& axis_vals, std::vector<double>& force_vals, std::vector<double>& potential_vals, const double cutoff)
{	
	if (char_id == 'n') {
    	pad_and_print_single_table(char_id, basename, axis_vals, force_vals, 0.0);
	} else if (char_id == 'b') {
		pad_and_print_single_table(char_id, basename, axis_vals, force_vals, cutoff);
	} else if (char_id == 'a') {
    	pad_and_print_single_table(char_id, basename, axis_vals, force_vals, 180.0);
    } else if (char_id == 'g') {
    	int size = axis_vals.size();
    	std::vector<double> rg_potential_vals;
    	std::vector<double> sqrt_axis_vals(size);
    	for (int i = 0; i < size; i++) {
    		sqrt_axis_vals[i] = sqrt(axis_vals[i]);
	    }
    	integrate_force(sqrt_axis_vals, force_vals, rg_potential_vals);
	    write_LAMMPS_table_output_file(char_id, basename, axis_vals, rg_potential_vals, force_vals);
    } else if (char_id == 'd') {
   	 trim_excess_axis(-180.0, 180.0, axis_vals, force_vals);
   	 std::vector<double> corrected_potential_vals;
   	 integrate_force(axis_vals, force_vals, corrected_potential_vals);
   	 write_LAMMPS_table_output_file(char_id, basename, axis_vals, potential_vals, force_vals); 
    } else {		
    	write_LAMMPS_table_output_file(char_id, basename, axis_vals, potential_vals, force_vals);   
    }
}

void print_table_files(const char char_id, std::string& basename, std::vector<double>& axis_vals, std::vector<double>& force_vals, std::vector<double>& potential_vals)
{	
	if (char_id == 'g') {
    	int size = axis_vals.size();
    	std::vector<double> rg_potential_vals;
    	std::vector<double> sqrt_axis_vals(size);
    	for (int i = 0; i < size; i++) {
    		sqrt_axis_vals[i] = sqrt(axis_vals[i]);
	    }
	    integrate_force(sqrt_axis_vals, force_vals, rg_potential_vals);
	    write_LAMMPS_table_output_file(char_id, basename, axis_vals, rg_potential_vals, force_vals);
    } else {		
    	write_LAMMPS_table_output_file(char_id, basename, axis_vals, potential_vals, force_vals);   
    }
}

void pad_and_print_single_table(const char char_id, const std::string& basename, std::vector<double>& axis_vals, std::vector<double>& force_vals, const double cutoff)
{
   std::vector<double> padded_potential_vals;
   
   // pad front
   int status = pad_values_front_with_fix(axis_vals,force_vals);
   if(char_id == 'a'){
     pad_values_front(0.0,axis_vals,force_vals,force_vals[0]);
   }
   if (status == -1) {
	printf("Error encountered when padding lower end of %s! Please check the output tables carefully before using!\n", basename.c_str());
   }
   
   // pad back
   status = pad_values_back_with_fix(cutoff,axis_vals,force_vals);
   if(char_id == 'a') {
     pad_values_back(180.0,axis_vals,force_vals,force_vals[force_vals.size()-1]);
     trim_excess_axis(0.0, 180.0, axis_vals, force_vals);
   }
   if (status == -1) {
	printf("Error encountered when padding upper end of %s! Please check the output tables carefully before using!\n", basename.c_str());
   }
   
   // Generate padded_potential_vals
   integrate_force(axis_vals, force_vals, padded_potential_vals);
   
   // write LAMMPS table using padded forces and potentials
   write_LAMMPS_table_output_file(char_id, basename, axis_vals, padded_potential_vals, force_vals);
}

void write_two_param_bspline_table_file(InteractionClassComputer* const icomp, char ** const name, MATRIX_DATA* const mat, const int index_among_defined)
{
    // Print out a table of the interaction forces.
    std::string filename = icomp->ispec->get_basename(name, index_among_defined, "_") + ".dat";
    FILE* curr_spline_output_file = open_file(filename.c_str(), "w");
	std::vector<double> axis_vals, force_vals, deriv_vals;
	icomp->calc_grid_of_force_and_deriv_vals(mat->fm_solution, index_among_defined, icomp->ispec->output_binwidth, axis_vals, force_vals, deriv_vals);
	for (unsigned i = 0; i < axis_vals.size(); i++) {
        fprintf(curr_spline_output_file, "%lf %.15le %.15le\n", axis_vals[i], force_vals[i], deriv_vals[i]);
    }
    fclose(curr_spline_output_file);
}

// Write the linear spline coefficient output for a single interaction.

void write_one_param_linear_spline_file(InteractionClassComputer* const icomp, char ** const name, MATRIX_DATA* const mat, const int index_among_defined_intrxns)
{
    if (icomp->ispec->output_spline_coeffs_flag == 1) {
        char name_tmp1[100], name_tmp2[100];
        FILE *raw_rhs_output_file, *normal_form_rhs_output_file;
        std::string basename;
		
		if (icomp->ispec->class_type == kDensity) {
			DensityClassSpec* ispec = static_cast<DensityClassSpec*>(icomp->ispec);
			basename = ispec->get_basename(name, index_among_defined_intrxns, "_");
		} else {
			basename = icomp->ispec->get_basename(name, index_among_defined_intrxns, "_");
		}

        int index_among_matched_interactions = icomp->ispec->defined_to_matched_intrxn_index_map[index_among_defined_intrxns];
    	// Select symmetric basename modifier if appropriate (i.e. DOOM interactions)
		if(icomp->ispec->defined_to_symmetric_intrxn_index_map[index_among_defined_intrxns] != 0) {
			basename += "_sym";	
		}

        // Output the linear spline bin points and coefficients.
        sprintf(name_tmp1, "%s.b", basename.c_str());
        raw_rhs_output_file = open_file(name_tmp1, "w");
        for (unsigned i = icomp->interaction_class_column_index + icomp->ispec->interaction_column_indices[index_among_matched_interactions - 1]; i < icomp->interaction_class_column_index + icomp->ispec->interaction_column_indices[index_among_matched_interactions]; i++) {
            fprintf(raw_rhs_output_file, "%lf %.15le\n", icomp->ispec->lower_cutoffs[index_among_defined_intrxns] + icomp->ispec->get_fm_binwidth() * (i - icomp->interaction_class_column_index - icomp->ispec->interaction_column_indices[index_among_matched_interactions]), mat->fm_solution[i]);
        }
        fclose(raw_rhs_output_file);
        
        // Output normal equation right hand side for this interaction.
        if (mat->output_normal_equations_rhs_flag == 1) {
            sprintf(name_tmp2, "%s.dense_fm_normal_rhs_vector", basename.c_str());
            normal_form_rhs_output_file = open_file(name_tmp2, "w");
            for (unsigned i = icomp->interaction_class_column_index + icomp->ispec->interaction_column_indices[index_among_matched_interactions - 1]; i < icomp->interaction_class_column_index + icomp->ispec->interaction_column_indices[index_among_matched_interactions]; i++) {
            	fprintf(normal_form_rhs_output_file, "%lf %.15le\n", icomp->ispec->lower_cutoffs[index_among_defined_intrxns] + icomp->ispec->get_fm_binwidth() * (i - icomp->interaction_class_column_index - icomp->ispec->interaction_column_indices[index_among_matched_interactions - 1]), mat->dense_fm_normal_rhs_vector[i]);
            }
            fclose(normal_form_rhs_output_file);
        }
    }
}

// Write the bspline coefficient output for a single interaction.

void write_one_param_bspline_file(InteractionClassComputer* const icomp, char ** const name, MATRIX_DATA* const mat, const int index_among_defined)
{
    std::string type_names;
    if (icomp->ispec->class_type == kDensity) {
   		DensityClassSpec* dspec = static_cast<DensityClassSpec*>(icomp->ispec);
   	    type_names = dspec->get_interaction_name(name, index_among_defined, " ");
	} else {
   	    type_names = icomp->ispec->get_interaction_name(name, index_among_defined, " ");
	}

	FILE* spline_output_filep = open_file("b-spline.out", "a");
	
	fprintf(spline_output_filep, "%c: ", icomp->ispec->get_char_id());
	fprintf(spline_output_filep, "%s ", type_names.c_str());

    // Print number of splines and cutoffs for the interactions.
    int index_among_matched = icomp->ispec->defined_to_matched_intrxn_index_map[index_among_defined];
    int n_basis_funcs = icomp->ispec->interaction_column_indices[index_among_matched] - icomp->ispec->interaction_column_indices[index_among_matched - 1];
    int n_to_print_minus_bspline_k = n_basis_funcs - icomp->ispec->get_bspline_k() + 2;
   
   fprintf(spline_output_filep, "%d %d ", icomp->ispec->get_bspline_k(), n_to_print_minus_bspline_k);
   fprintf(spline_output_filep, "%.15le %.15le\n", icomp->ispec->lower_cutoffs[index_among_defined], icomp->ispec->upper_cutoffs[index_among_defined]);

    // Print the spline coefficients.
    for (int k = 0; k < n_basis_funcs; k++) {
        fprintf(spline_output_filep, "%.15le ", mat->fm_solution[icomp->ispec->interaction_column_indices[index_among_matched - 1] + k]);
    }
    // Complete the line.
    fprintf(spline_output_filep, "\n");
	fclose(spline_output_filep);
}

void write_output_solution(MATRIX_DATA* const mat)
{
	FILE* xout = open_file("x.out", "wb");
	if (mat->bootstrapping_flag == 1) {
		for (int i = 0; i < mat->bootstrapping_num_estimates; i++) {
			fwrite(&(mat->bootstrap_solutions[i][0]), sizeof(double), mat->fm_matrix_columns, xout);
		}
	} else {
		fwrite(&(mat->fm_solution[0]), sizeof(double), mat->fm_matrix_columns, xout);
	}
	fclose(xout);
}

void write_full_bootstrapping_MSCGFM_table_output_file(const std::string& filename_base, const std::vector<double>& axis, std::vector<double> const master_force, std::vector<double>* const force, int const bootstrapping_num_estimates) 
{
	std::string filename_tmp = filename_base + ".dat";
    FILE *curr_table_output_file = open_file(filename_tmp.c_str(), "w");
    for (unsigned i = 0; i < axis.size(); i++) {
        fprintf(curr_table_output_file, "%lf\t", axis[i]);
        fprintf(curr_table_output_file, "%lf\t", master_force[i]);
        for (int j = 0; j < bootstrapping_num_estimates; j++) {
        	fprintf(curr_table_output_file, " %.15le", force[j][i]);
    	}
    	fprintf(curr_table_output_file, "\n");
    }
    fclose(curr_table_output_file);
}

std::vector<double> calculate_bootstrapping_standard_error(const std::vector<double> master_vals, const std::vector<double>* vals, const int bootstrapping_num_estimates)
{
	double sum, squared;
	std::vector<double> standard_error(master_vals.size());
	for (unsigned i = 0; i < master_vals.size(); i++) {
		sum = master_vals[i];
		squared = master_vals[i] * master_vals[i];
		for (int j = 0; j < bootstrapping_num_estimates; j++) {
			sum += vals[j][i];
			squared += vals[j][i] * vals[j][i];
		}
		standard_error[i] = sqrt(squared - (sum * sum / (double)(bootstrapping_num_estimates + 1)) ) / (double)(bootstrapping_num_estimates + 1);
	}
	return standard_error;
}

void write_bootstrapping_MSCGFM_table_output_file(const std::string& filename_base, const std::vector<double>& axis, std::vector<double> const master_force, std::vector<double>* const force, int const bootstrapping_num_estimates) 
{
	std::string filename_tmp = filename_base + ".dat";
    FILE *curr_table_output_file = open_file(filename_tmp.c_str(), "w");

    // Calculate standard error for all samples
    std::vector<double> standard_error = calculate_bootstrapping_standard_error(master_force, force, bootstrapping_num_estimates);
    for (unsigned i = 0; i < axis.size(); i++) {
        fprintf(curr_table_output_file, "%lf\t", axis[i]);
        fprintf(curr_table_output_file, "%lf\t", master_force[i]);
        fprintf(curr_table_output_file, "%lf", standard_error[i]);
    	fprintf(curr_table_output_file, "\n");
    }
    fclose(curr_table_output_file);
}

void write_bootstrapping_one_param_table_files(InteractionClassComputer* const icomp, char **name, std::vector<double> const master_coeffs, std::vector<double>* const spline_coeffs, const int index_among_defined_intrxns, const int bootstrapping_num_estimates, const int bootstrapping_full_output_flag) 
{
    // Compute forces over a grid of parameter values.
    std::vector<double> axis_vals;
    std::vector<double> master_force_vals;
    std::vector<double> master_potential_vals;	
    std::vector<double>* force_vals = new std::vector<double>[bootstrapping_num_estimates];
    std::vector<double>* potential_vals = new std::vector<double>[bootstrapping_num_estimates];	
    int i;
    
    // Master
    if (dynamic_cast<OneBodyClassComputer*>(icomp) != NULL) {
    	icomp->calc_one_force_val(master_coeffs, index_among_defined_intrxns, icomp->ispec->output_binwidth, axis_vals, master_force_vals, master_potential_vals);
 	} else {
   		icomp->calc_grid_of_force_vals(master_coeffs, index_among_defined_intrxns, icomp->ispec->output_binwidth, axis_vals, master_force_vals);
    	integrate_force(axis_vals, master_force_vals, master_potential_vals);
    }
    
    // Bootstrap estimates
    for (i = 0; i < bootstrapping_num_estimates; i++) {
	    if (dynamic_cast<OneBodyClassComputer*>(icomp) != NULL) {
    		icomp->calc_one_force_val(spline_coeffs[i], index_among_defined_intrxns, icomp->ispec->output_binwidth, axis_vals, force_vals[i], potential_vals[i]);
 		} else {
 			icomp->calc_grid_of_force_vals(spline_coeffs[i], index_among_defined_intrxns, icomp->ispec->output_binwidth, axis_vals, force_vals[i]);
			// Integrate force starting from cutoff = 0.0 potential.
    		integrate_force(axis_vals, force_vals[i], potential_vals[i]);
		}
		
	}
    
    // Print out tabulated output files in MSCGFM style and LAMMPS style.
    std::string basename = icomp->ispec->get_basename(name, index_among_defined_intrxns, "_");

    // Select symmetric basename modifier if appropriate (i.e. DOOM interactions)
	if(icomp->ispec->defined_to_symmetric_intrxn_index_map[index_among_defined_intrxns] != 0) {
		basename += "_sym";	
	}

  	if (bootstrapping_full_output_flag == 1) {
  		// Write master followed by all estimates.
   		write_full_bootstrapping_MSCGFM_table_output_file(basename, axis_vals, master_force_vals, force_vals, bootstrapping_num_estimates);
    } else {
    	// Write master and standard error only.
    	write_bootstrapping_MSCGFM_table_output_file(basename, axis_vals, master_force_vals, force_vals, bootstrapping_num_estimates);
    }
    // Only write master copy for LAMMPS output.
    write_LAMMPS_table_output_file(icomp->ispec->get_char_id(), basename, axis_vals, potential_vals[0], force_vals[0]); 
    delete [] force_vals;
    delete [] potential_vals;  
}

void write_bootstrapping_one_param_table_files_energy(InteractionClassComputer* const icomp, char **name, std::vector<double> const master_coeffs, std::vector<double>* const spline_coeffs, const int index_among_defined_intrxns, const int bootstrapping_num_estimates, const int bootstrapping_full_output_flag, const double cutoff) 
{
    // Compute forces over a grid of parameter values.
    std::vector<double> axis_vals;
    std::vector<double> master_force_vals;
    std::vector<double> master_potential_vals;	
    std::vector<double>* force_vals = new std::vector<double>[bootstrapping_num_estimates];
    std::vector<double>* potential_vals = new std::vector<double>[bootstrapping_num_estimates];	
    int i;
    
    // Master
    if (dynamic_cast<OneBodyClassComputer*>(icomp) != NULL) {
    	icomp->calc_one_force_val(master_coeffs, index_among_defined_intrxns, icomp->ispec->output_binwidth, axis_vals, master_force_vals, master_potential_vals);
 	} else {
        icomp->calc_grid_of_force_and_deriv_vals(master_coeffs, index_among_defined_intrxns, icomp->ispec->output_binwidth, axis_vals, master_potential_vals, master_force_vals);
        make_negative(master_force_vals);
    }
    
    // Bootstrap estimates
    for (i = 0; i < bootstrapping_num_estimates; i++) {
	    if (dynamic_cast<OneBodyClassComputer*>(icomp) != NULL) {
    		icomp->calc_one_force_val(spline_coeffs[i], index_among_defined_intrxns, icomp->ispec->output_binwidth, axis_vals, force_vals[i], potential_vals[i]);
 		} else {
 			icomp->calc_grid_of_force_and_deriv_vals(spline_coeffs[i], index_among_defined_intrxns, icomp->ispec->output_binwidth, axis_vals, potential_vals[i], force_vals[i]);
			make_negative(force_vals[i]);
		}
	}

    // Determine base for output filenames.
    std::string basename;
    if(icomp->ispec->class_type == kDensity)
      {
        DensityClassSpec* ispec = static_cast<DensityClassSpec*>(icomp->ispec);
        basename = ispec->get_basename(ispec->density_group_names, index_among_defined_intrxns, "_");
      }
    else
      {
        basename = icomp->ispec->get_basename(name, index_among_defined_intrxns, "_");
      }

    // Select symmetric basename modifier if appropriate (i.e. DOOM interactions)
	if(icomp->ispec->defined_to_symmetric_intrxn_index_map[index_among_defined_intrxns] != 0) {
		basename += "_sym";	
	}

  	if (bootstrapping_full_output_flag == 1) {
  		// Write master followed by all estimates.
   		write_full_bootstrapping_MSCGFM_table_output_file(basename, axis_vals, master_force_vals, force_vals, bootstrapping_num_estimates);
    } else {
    	// Write master and standard error only.
    	write_bootstrapping_MSCGFM_table_output_file(basename, axis_vals, master_force_vals, force_vals, bootstrapping_num_estimates);
    }
    // Only write master copy for LAMMPS output. 
	pad_and_print_table_files(icomp->ispec->get_char_id(), basename, axis_vals, master_force_vals, master_potential_vals, cutoff);
    delete [] force_vals;
    delete [] potential_vals;  
}

void write_bootstrapping_one_param_bspline_file(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, const int index_among_defined_intrxns)
{
    FILE *spline_output_filep = open_file("b-spline.out", "a");

    // Print out class character & the types involved.
    fprintf(spline_output_filep, "%c: ", icomp->ispec->get_char_id());
    std::string type_names = icomp->ispec->get_interaction_name(name, index_among_defined_intrxns, " ");
	fprintf(spline_output_filep, "%s ", type_names.c_str());
    
    // Print number of splines and cutoffs for the interactions.
    int index_among_matched_interactions = icomp->ispec->defined_to_matched_intrxn_index_map[index_among_defined_intrxns];
    unsigned interaction_column_indices = icomp->ispec->interaction_column_indices[index_among_matched_interactions] - icomp->ispec->interaction_column_indices[index_among_matched_interactions - 1];
    unsigned n_to_print_minus_bspline_k = interaction_column_indices - icomp->ispec->get_bspline_k() + 2;
    fprintf(spline_output_filep, "%d %d %.15le %.15le\n", icomp->ispec->get_bspline_k(), n_to_print_minus_bspline_k, icomp->ispec->lower_cutoffs[index_among_defined_intrxns], icomp->ispec->upper_cutoffs[index_among_defined_intrxns]);
    
    for (int i = 0; i < mat->bootstrapping_num_estimates; i++) {
    	// Print the spline coefficients.
    	for (unsigned k = 0; k < interaction_column_indices; k++) {
        	fprintf(spline_output_filep, "%.15le ", mat->bootstrap_solutions[i][icomp->ispec->interaction_column_indices[index_among_matched_interactions - 1] + k]);
    	}
    
	    // Complete the line and flush.
    	fprintf(spline_output_filep, "\n"); fflush(spline_output_filep); 
    }
    
    // Close the file.
    fclose(spline_output_filep);
}

void write_bootstrapping_one_param_linear_spline_file(InteractionClassComputer* const icomp, char **name, MATRIX_DATA* const mat, const int index_among_defined_intrxns)
{
    if (icomp->ispec->output_spline_coeffs_flag == 1) {
        char name_tmp1[100], name_tmp2[100];
        FILE *raw_rhs_output_file, *normal_form_rhs_output_file;
        std::string basename = icomp->ispec->get_basename(name, index_among_defined_intrxns, "_");
        
        int index_among_matched_interactions = icomp->ispec->defined_to_matched_intrxn_index_map[index_among_defined_intrxns];
    	// Select symmetric basename modifier if appropriate (i.e. DOOM interactions)
		if(icomp->ispec->defined_to_symmetric_intrxn_index_map[index_among_defined_intrxns] != 0) {
			basename += "_sym";	
		}

        // Output the linear spline bin points and coefficients.
        sprintf(name_tmp1, "%s.b", basename.c_str());
        raw_rhs_output_file = open_file(name_tmp1, "w");
        for (unsigned i = icomp->interaction_class_column_index + icomp->ispec->interaction_column_indices[index_among_matched_interactions - 1]; i < icomp->interaction_class_column_index + icomp->ispec->interaction_column_indices[index_among_matched_interactions]; i++) {
               fprintf(raw_rhs_output_file, "%lf", icomp->ispec->lower_cutoffs[index_among_defined_intrxns] + icomp->ispec->get_fm_binwidth() * (i - icomp->interaction_class_column_index - icomp->ispec->interaction_column_indices[index_among_matched_interactions]));
               for (int j = 0; j < mat->bootstrapping_num_estimates; j++) {               
                  fprintf(raw_rhs_output_file, " %.15le", mat->bootstrap_solutions[j][i]);
               }
               fprintf(raw_rhs_output_file, "\n");
        }
        fclose(raw_rhs_output_file);
        
        // Output normal equation right hand side for this interaction.
        if (mat->output_normal_equations_rhs_flag == 1) {
            sprintf(name_tmp2, "%s.dense_fm_normal_rhs_vector", basename.c_str());
            normal_form_rhs_output_file = open_file(name_tmp2, "w");
            for (unsigned i = icomp->interaction_class_column_index + icomp->ispec->interaction_column_indices[index_among_matched_interactions - 1]; i < icomp->interaction_class_column_index + icomp->ispec->interaction_column_indices[index_among_matched_interactions]; i++) {
                fprintf(normal_form_rhs_output_file, "%lf", icomp->ispec->lower_cutoffs[index_among_defined_intrxns] + icomp->ispec->get_fm_binwidth() * (i - icomp->interaction_class_column_index - icomp->ispec->interaction_column_indices[index_among_matched_interactions - 1]));
            	for (int j = 0; j < mat->bootstrapping_num_estimates; j++) {
	            	fprintf(normal_form_rhs_output_file, "%.15le",  mat->bootstrapping_dense_fm_normal_rhs_vectors[j][i]);
    			}
            	fprintf(normal_form_rhs_output_file, "\n");
            }
            fclose(normal_form_rhs_output_file);
        }
    }
}
