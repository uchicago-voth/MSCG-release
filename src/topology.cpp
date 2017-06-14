//
//  topology.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cstdlib>
#include <cstring>
#include <cassert>
#include <fstream>
#include <vector>

#include "interaction_model.h"
#include "misc.h"

#include "topology.h"

// Helper functions for reading a topology file for the CG model

// Determine the maximum number of items possible for a given exclusion.
inline unsigned get_max_exclusion_number(TopologyData const* topo_data, const int excluded_style);
// Read the molecule groups for radius of gyration types defined in the topology file.
int read_molecule_groups(TopologyData* topo_data, CG_MODEL_DATA* const cg, std::ifstream &top_in, int line_num, std::string* &rg_lines);
void process_molecule_groups(TopologyData* topo_data, std::string* &rg_lines);
// Read the specification for three-body interactions.
int read_three_body_topology(TopologyData* topo_data, CG_MODEL_DATA* const cg, std::ifstream &top_in, int line);
// Read the density groups defined in the topology file.
int read_density_groups(TopologyData* topo_data, CG_MODEL_DATA* const cg, std::ifstream &top_in, int line_num);
int read_density_weights(TopologyData* topo_data, CG_MODEL_DATA* const cg, std::ifstream &top_in, int line_num);
// Read a single molecule's topology specification.
void read_molecule_definition(TopologyData* const mol, TopologyData* const cg, std::ifstream &top_in, int &line);
// Determine appropriate non-bonded exclusions based on bonded topology and excluded_style setting
bool new_excluded_partner(unsigned** const excluded_partners, const int location, const unsigned current_size, const unsigned new_site);
// Report an error in the topology input file.
void report_topology_input_format_error(const int line, char *parameter_name);
// Search function for molecular exclusion.
void recursive_exclusion_search(TopologyData const* topo_data, TopoList* &exclusion_list, std::vector<int> &path_list);

//---------------------------------------------------------------
// Functions for managing TopoList structs
//---------------------------------------------------------------

TopoList::TopoList(unsigned n_sites, unsigned partners_per, unsigned max_partners) : 
    n_sites_(n_sites), partners_per_(partners_per), max_partners_(max_partners) {
    modified = 0;
    partner_numbers_ = new unsigned[n_sites_]();
    partners_ = new unsigned*[n_sites_];
    for (unsigned i = 0; i < n_sites_; i++) {
        partners_[i] = new unsigned[partners_per_ * max_partners_]();
    }
}

TopoList::~TopoList() {
	if (modified == 0) {
    	for (unsigned i = 0; i < n_sites_; i++) {
        	delete [] partners_[i];
    	}
    	delete [] partners_;
    	delete [] partner_numbers_;
	}
}

//---------------------------------------------------------------
// Functions for managing TopologyData structs
//---------------------------------------------------------------

inline unsigned get_max_exclusion_number(TopologyData const* topo_data, const int excluded_style) 
{
	assert( (excluded_style == 0) || (excluded_style == 2) || (excluded_style == 3) || (excluded_style == 4) || (excluded_style == 5) );
	switch (excluded_style) {
		case 2:
			return topo_data->max_pair_bonds_per_site;
			break;
			
		case 3:
			return topo_data->max_pair_bonds_per_site + topo_data->max_angles_per_site;
			break;
			
		case 4:
			return topo_data->max_pair_bonds_per_site + topo_data->max_angles_per_site + topo_data->max_dihedrals_per_site;
			break;
			
		case 5:
			return topo_data->max_pair_bonds_per_site + topo_data->max_angles_per_site + topo_data->max_dihedrals_per_site;
			break;
		
		default:
			return 1;
			break;
	}
}

// Initialize a TopologyData struct given that the number of sites and maximum
// numbers of bonds, angles, and dihedrals for each site are fixed.
void initialize_topology_data(TopologyData* const topo_data)
{
	unsigned nsites = topo_data->n_cg_sites;
    topo_data->cg_site_types = new int[nsites]();
    topo_data->molecule_ids  = new int[nsites]();

    topo_data->name = new char*[topo_data->n_cg_types]();
    for (unsigned i = 0; i < topo_data->n_cg_types; i++) {
        topo_data->name[i] = new char[MAX_CG_TYPE_NAME_LENGTH + 1]();
    }
    topo_data->max_quints_per_site = topo_data->max_dihedrals_per_site + 2;
    // Ready pair bond data structures.
    topo_data->bond_type_activation_flags = new int[calc_n_distinct_pairs(topo_data->n_cg_types)]();
    topo_data->bond_list = new TopoList(nsites, 1, topo_data->max_pair_bonds_per_site);
    // Ready angle bond data structures.
    topo_data->angle_type_activation_flags = new int[calc_n_distinct_triples(topo_data->n_cg_types)]();
    topo_data->angle_list = new TopoList(nsites, 2, topo_data->max_angles_per_site);
    // Ready dihedral bond data structures.
    topo_data->dihedral_type_activation_flags = new int[calc_n_distinct_quadruples(topo_data->n_cg_types)]();
    topo_data->dihedral_list = new TopoList(nsites, 3, topo_data->max_dihedrals_per_site);
    // Ready dihedral bond data structures.
    topo_data->quint_type_activation_flags = new int[calc_n_distinct_quints(topo_data->n_cg_types)]();
    topo_data->quint_list = new TopoList(nsites, 4, topo_data->max_quints_per_site);
    // Ready exclusion list data structures.
	topo_data->exclusion_list = new TopoList(nsites, 1, get_max_exclusion_number(topo_data, topo_data->excluded_style));
	// Read density exclusion list data structures.
	topo_data->density_exclusion_list = new TopoList(nsites, 1, get_max_exclusion_number(topo_data, topo_data->density_excluded_style));
}

// Free an initialized TopologyData struct.
void TopologyData::free_topology_data(void) {
  printf("Freeing topology information.\n");
  
  if (cg_site_types != NULL) delete [] cg_site_types;
  if (molecule_ids != NULL) delete [] molecule_ids;
    
    if (name != NULL) {
        for (unsigned i = 0; i < n_cg_types; i++) {
            if (name[i] != NULL) delete [] name[i];
        }
        delete [] name;
    }

    // Delete the topology lists.
    if (bond_list != NULL) delete bond_list;
    if (angle_list != NULL) delete angle_list;
    if (dihedral_list != NULL) delete dihedral_list;
    if (quint_list != NULL) delete quint_list;
    if (molecule_list != NULL) delete molecule_list;
    delete exclusion_list;
	delete density_exclusion_list;
	
    // Delete the topology interaction type activation flags.
    if (bond_type_activation_flags != NULL) delete [] bond_type_activation_flags;
    if (angle_type_activation_flags != NULL) delete [] angle_type_activation_flags;
    if (dihedral_type_activation_flags != NULL) delete [] dihedral_type_activation_flags;
    if (quint_type_activation_flags != NULL) delete [] quint_type_activation_flags;
}

//---------------------------------------------------------------
// Functions for reading a topology file for the CG model
//---------------------------------------------------------------

void read_topology_file(TopologyData* topo_data, CG_MODEL_DATA* const cg) 
{
    int line = 0;
    unsigned i;
    char parameter_name[50];
    std::string* rg_lines = NULL;
    std::string buff;
    std::ifstream top_in;
    check_and_open_in_stream(top_in, "top.in");
    
    // Read the number of sites and number of types of sites.
    
    check_and_read_next_line(top_in, buff, line);
    sscanf(buff.c_str(), "%s%u", parameter_name, &topo_data->n_cg_sites);
    if (strcmp(parameter_name, "cgsites") != 0) report_topology_input_format_error(line,parameter_name);
    cg->n_cg_sites = int(topo_data->n_cg_sites);
    check_and_read_next_line(top_in, buff, line);
    sscanf(buff.c_str(), "%s%u", parameter_name, &topo_data->n_cg_types);
    if (strcmp(parameter_name, "cgtypes") != 0) report_topology_input_format_error(line,parameter_name);
    cg->n_cg_types = int(topo_data->n_cg_types);

    // Allocate bond, angle, and dihedral interaction
    // class arrays indexed by all possible such interactions
    // using that information.
    initialize_topology_data(topo_data);
    cg->name = new char*[cg->n_cg_types]();
    for (i = 0; i < unsigned(cg->n_cg_types); i++) {
        cg->name[i] = new char[MAX_CG_TYPE_NAME_LENGTH + 1];
    }
    
    // Read the string names of each site type.
    for (i = 0; i < topo_data->n_cg_types; i++) {
        check_and_read_next_line(top_in, buff, line);
        sscanf(buff.c_str(), "%s", topo_data->name[i]);
    }
    
    // Read the section of top.in defining radius of gyration interactions, if present.
    if (cg->radius_of_gyration_interactions.class_subtype > 0 || cg->helical_interactions.class_subtype > 0) {
    	line = read_molecule_groups(topo_data, cg, top_in, line, rg_lines);
    }
    
    // Read the section of top.in defining three body nonbonded interactions, if present.
    if (cg->three_body_nonbonded_interactions.class_subtype > 0) {
		line = read_three_body_topology(topo_data, cg, top_in, line);
	}
    
	// Read the section of top.in defining density dependent nonbonded interactions, if present.
	topo_data->n_density_groups = 0;
	if(cg->density_interactions.class_subtype > 0) {
		line = read_density_groups(topo_data, cg, top_in, line);
		line = read_density_weights(topo_data, cg, top_in, line);
	}
	
    // Read the molecules section of top.in to define all topology lists.
    
    // Read the number of distinct types of molecule
    check_and_read_next_line(top_in, buff, line);
    unsigned n_molecule_types;
    sscanf(buff.c_str(), "%s%d", parameter_name, &n_molecule_types);
    if (strcmp(parameter_name, "moltypes") != 0) report_topology_input_format_error(line,parameter_name);
    
    // Allocate space for an array containing information pertaining
    // to each individual molecule type treated as an entire system
    // and initialize the topology data sections of them.
    TopologyData* mol = new TopologyData[n_molecule_types];
    for (i = 0; i < n_molecule_types; i++) {
        mol[i] = TopologyData(topo_data->max_pair_bonds_per_site, topo_data->max_angles_per_site, topo_data->max_dihedrals_per_site);
	}

    // Read the definitions for each molecule to get the topological
    // building blocks of the full system.
    // Assume that one should input/output angles and dihedrals in bonding order (ABC, ABCD)
    int max_mol_size = 0;
    topo_data->angle_format = 1;
	topo_data->dihedral_format = 1;
    for (i = 0; i < n_molecule_types; i++) {
        printf("Reading definition of molecule type %d.\n", i);
        read_molecule_definition(&mol[i], topo_data, top_in, line);
        if ((int)(mol[i].n_cg_sites) > max_mol_size) {
        	max_mol_size = (int)(mol[i].n_cg_sites);
        }
        // Make sure that mol_style is consistent.
        if (i == 0) {
        	topo_data->angle_format = mol[i].angle_format;
        	topo_data->dihedral_format = mol[i].dihedral_format;
        } else {
        	if (topo_data->angle_format != mol[i].angle_format) {
        		fprintf(stderr, "Cannot have two different angle formats!\n");
        		exit(EXIT_FAILURE);
        	}
        	if (topo_data->dihedral_format != mol[i].dihedral_format) {
        	    fprintf(stderr, "Cannot have two different dihedral formats!\n");
        		exit(EXIT_FAILURE);
        	}
        }
    }
    
    // Read the system section of top.in to define the full system's
    // topology.
    int total_cg_molecules = 0;
    unsigned j, k, l;
    unsigned block, molecule_type, n_molecules_in_system, first_site_in_molecule = 0;
    check_and_read_next_line(top_in, buff, line);
    sscanf(buff.c_str(), "%s%d", parameter_name, &block);
    if (strcmp(parameter_name, "system") != 0) report_topology_input_format_error(line, parameter_name);
    
    // For each line in this section, add a given number of molecules
    // of a given type to the system in order.
    for (i = 0; i < block; i++) {
        // Get a number and type of molecules.
        check_and_read_next_line(top_in, buff, line);
        sscanf(buff.c_str(), "%d%d", &molecule_type, &n_molecules_in_system);
        molecule_type--;
        total_cg_molecules += n_molecules_in_system;
        assert(total_cg_molecules <= cg->n_cg_sites);
        
        // For each molecules of that type to add to the growing
        // total system definition, add an appropriate number
        // of CG sites to the topology lists using the molecule's
        // topology as a guide.
        for (j = 0; j < n_molecules_in_system; j++) {
            for (k = 0; k < mol[molecule_type].n_cg_sites; k++) {
            	if (mol[molecule_type].bond_list->partner_numbers_[k] > topo_data->max_pair_bonds_per_site) {
                	printf("Warning: Too many bonds (%d) to handle based on max_pair_bonds_per_site (%d) !\n", mol[molecule_type].bond_list->partner_numbers_[k], topo_data->max_pair_bonds_per_site);
    				exit(EXIT_FAILURE);
                }
                if (mol[molecule_type].angle_list->partner_numbers_[k] > topo_data->max_angles_per_site) {
                	printf("Warning: Too many bonds (%d) to handle based on max_angles_per_site (%d) !\n", mol[molecule_type].angle_list->partner_numbers_[k], topo_data->max_angles_per_site);
    				exit(EXIT_FAILURE);
                }
                if (mol[molecule_type].dihedral_list->partner_numbers_[k] > topo_data->max_dihedrals_per_site) {
                	printf("Warning: Too many bonds (%d) to handle based on max_dihedrals_per_site (%d) !\n", mol[molecule_type].dihedral_list->partner_numbers_[k], topo_data->max_dihedrals_per_site);
    				exit(EXIT_FAILURE);
                }
                topo_data->cg_site_types[k + first_site_in_molecule] = mol[molecule_type].cg_site_types[k];
                topo_data->molecule_ids[k + first_site_in_molecule] = j;
                topo_data->bond_list->partner_numbers_[k + first_site_in_molecule] = mol[molecule_type].bond_list->partner_numbers_[k];
                for (l = 0; l < mol[molecule_type].bond_list->partner_numbers_[k]; l++) {
                    topo_data->bond_list->partners_[k + first_site_in_molecule][l] = mol[molecule_type].bond_list->partners_[k][l] + first_site_in_molecule;
                }
                topo_data->angle_list->partner_numbers_[k + first_site_in_molecule] = mol[molecule_type].angle_list->partner_numbers_[k];
                for (l = 0; l < 2 * mol[molecule_type].angle_list->partner_numbers_[k]; l++) {
                    topo_data->angle_list->partners_[k + first_site_in_molecule][l] = mol[molecule_type].angle_list->partners_[k][l] + first_site_in_molecule;
                }
                topo_data->dihedral_list->partner_numbers_[k + first_site_in_molecule] = mol[molecule_type].dihedral_list->partner_numbers_[k];
                for (l = 0; l < 3 * mol[molecule_type].dihedral_list->partner_numbers_[k]; l++) {
                    topo_data->dihedral_list->partners_[k + first_site_in_molecule][l] = mol[molecule_type].dihedral_list->partners_[k][l] + first_site_in_molecule;
                }
                topo_data->quint_list->partner_numbers_[k + first_site_in_molecule] = mol[molecule_type].quint_list->partner_numbers_[k];
                for (l = 0; l < 4 * mol[molecule_type].quint_list->partner_numbers_[k]; l++) {
                    topo_data->quint_list->partners_[k + first_site_in_molecule][l] = mol[molecule_type].quint_list->partners_[k][l] + first_site_in_molecule;
                }
            }
            first_site_in_molecule += mol[molecule_type].n_cg_sites;
        }
    }
    
    // Give the CG_MODEL_DATA struct its own copy of the type names.
    for (i = 0; i < topo_data->n_cg_types; i++) {
        sscanf(topo_data->name[i], "%s", cg->name[i]);
    }
    
	// Set-up appropriate bonded exclusions list from non-bonded interactions based on excluded_style
	setup_excluded_list(topo_data, topo_data->exclusion_list, topo_data->excluded_style);
	setup_excluded_list(topo_data, topo_data->density_exclusion_list, topo_data->density_excluded_style);
	setup_molecule_list(topo_data, topo_data->molecule_list, total_cg_molecules, max_mol_size);
	
	// Process the molecule group information now, if necessary.
	if (cg->radius_of_gyration_interactions.class_subtype > 0 || cg->helical_interactions.class_subtype > 0) {
    	process_molecule_groups(topo_data, rg_lines);
    }
	
    // Close the file and free the memory used to store the temporary
    // single-molecule topologies.
    top_in.close();
	for (i = 0; i < n_molecule_types; i++) {
        mol[i].free_topology_data();
    }
    delete [] mol;
}

int read_molecule_groups(TopologyData* topo_data, CG_MODEL_DATA* const cg, std::ifstream &top_in, int line_num, std::string* &rg_lines)
{
	std::string buff;
	char parameter_name[50] = "";
	
	// Check label.
	check_and_read_next_line(top_in, buff);
    line_num++;	
	sscanf(buff.c_str(), "%s%d", parameter_name, &(topo_data->n_molecule_groups));
  	if (strcmp(parameter_name, "molecule_groups") != 0) report_topology_input_format_error(line_num, parameter_name);
	if (topo_data->n_molecule_groups > 5) {
		printf("Cannot define more than 5 molecule groups!\n");
		exit(EXIT_FAILURE);
	}
	
	// Allocate data, read the lines, and store them for later processing
	rg_lines = new std::string[topo_data->n_molecule_groups];
	printf("Read definitions for %d molecule groups.\n", topo_data->n_molecule_groups);
	for(int i=0; i < topo_data->n_molecule_groups; i++) {
		check_and_read_next_line(top_in, rg_lines[i]);
		line_num++;
	}
	
	return line_num;
}

void process_molecule_groups(TopologyData* topo_data, std::string* &rg_lines)
{
	char type_num[5] = "";
	int mol_id, num_entries;
	unsigned n_sites = topo_data->molecule_list->n_sites_;
	int n_cg_molecules = (int)(n_sites);
	std::string line;
	
	std::string* type_list = new std::string[n_cg_molecules + 1];
    topo_data->molecule_groups = new bool[topo_data->n_molecule_groups * n_cg_molecules]();
	topo_data->molecule_group_names = new char*[topo_data->n_molecule_groups];
    for (int i = 0; i < topo_data->n_molecule_groups; i++) {
        topo_data->molecule_group_names[i] = new char[MAX_CG_TYPE_NAME_LENGTH + 1];
    }
	  
	for(int i=0; i < topo_data->n_molecule_groups; i++) {
		// Split line using string tools. 
		num_entries = StringSplit(rg_lines[i], " \t", type_list);
		// First entry is the name of the CG group.
		sprintf(topo_data->molecule_group_names[i], "%s", type_list[0].c_str());
		// Add cg_type to molecule_group definition.
		for (int j=1; j < num_entries; j++) {
			sprintf(type_num, "%s", type_list[j].c_str());
			if( isalpha(type_num[0]) ) {
				printf("CG types in molecule group need to be stated in term of numbers based on their order in the systems defined in top.in.\n");
				printf("Problem with label: %s\n", type_num);
				fflush(stdout);
				exit(EXIT_FAILURE);
			}
			mol_id = atoi( type_num );
			// Check that read was successful and then correct mol_id from 1-base to 0-base.
			if( (mol_id > 0) && (mol_id <= n_cg_molecules) ){
				topo_data->molecule_groups[i*(n_cg_molecules) + mol_id - 1] = true;
			} else {
				printf("invalid molecule number: %s\n", type_num); fflush(stdout);
			}
		}
	}
	delete [] type_list;
	delete [] rg_lines;
}

int read_three_body_topology(TopologyData* topo_data, CG_MODEL_DATA* const cg, std::ifstream &top_in, int line)
{
	std::string buff;
	char parameter_name[50];
	unsigned tbtype, i;
	int* tb_i;
	int* tb_j;
	int* tb_k;
	
	check_and_read_next_line(top_in, buff, line);
	sscanf(buff.c_str(), "%s%d", parameter_name, &tbtype);
	if (strcmp(parameter_name, "threebody") != 0) report_topology_input_format_error(line, parameter_name);

	cg->three_body_nonbonded_interactions.stillinger_weber_angle_parameters_by_type = new double[tbtype];
	cg->three_body_nonbonded_interactions.three_body_nonbonded_cutoffs = new double[tbtype];

	tb_i = new int[tbtype];
	tb_j = new int[tbtype];
	tb_k = new int[tbtype];
	
	for (i = 0; i < tbtype; i++) {
		check_and_read_next_line(top_in, buff, line);
		sscanf(buff.c_str(), "%d%d%d%lf%lf", tb_i + i, tb_j + i, tb_k + i, cg->three_body_nonbonded_interactions.stillinger_weber_angle_parameters_by_type + i, cg->three_body_nonbonded_interactions.three_body_nonbonded_cutoffs + i);
		tb_i[i]--;
		tb_j[i]--;
		tb_k[i]--;
	}
	
	cg->three_body_nonbonded_interactions.tb_n = new int[topo_data->n_cg_types]();
	cg->three_body_nonbonded_interactions.tb_list = new int*[topo_data->n_cg_types];
	
	for (i = 0; i < tbtype; i++) {
		cg->three_body_nonbonded_interactions.tb_n[tb_i[i]]++;
	}
	
	for (i = 0; i < topo_data->n_cg_types; i++) {
		cg->three_body_nonbonded_interactions.tb_list[i] = new int[cg->three_body_nonbonded_interactions.tb_n[i] * 2];
		cg->three_body_nonbonded_interactions.tb_n[i] = 0;
	}
	
	for (i = 0; i < tbtype; i++) {
		cg->three_body_nonbonded_interactions.tb_list[tb_i[i]][cg->three_body_nonbonded_interactions.tb_n[tb_i[i]] * 2] = tb_j[i] + 1;
		cg->three_body_nonbonded_interactions.tb_list[tb_i[i]][cg->three_body_nonbonded_interactions.tb_n[tb_i[i]] * 2 + 1] = tb_k[i] + 1;
		cg->three_body_nonbonded_interactions.tb_n[tb_i[i]]++;
	}

	cg->three_body_nonbonded_interactions.set_n_defined(tbtype);
	delete [] tb_i;
	delete [] tb_j;
	delete [] tb_k;
	
	return line;
}

int read_density_groups(TopologyData* topo_data, CG_MODEL_DATA* const cg, std::ifstream &top_in, int line_num)
{
	// Allocate data.
	std::string buff;
	char parameter_name[50] = "";
	int cg_type, num_entries;
	std::string* type_list = new std::string[cg->n_cg_types + 1];
    
	// Check label.
	check_and_read_next_line(top_in, buff);
	line_num++;	
	sscanf(buff.c_str(), "%s%d", parameter_name, &(topo_data->n_density_groups));
  	if (strcmp(parameter_name, "density_nonbonded") != 0) report_topology_input_format_error(line_num, parameter_name);
	if (topo_data->n_density_groups > 5) {
		printf("Cannot define more than 5 density groups!\n");
		delete [] type_list;
		exit(EXIT_FAILURE);
	}
	
	// Allocate data requiring n_density_groups
	topo_data->density_groups = new bool[topo_data->n_density_groups * cg->n_cg_types]();
	topo_data->density_group_names = new char*[topo_data->n_density_groups];
    for (int i = 0; i < topo_data->n_density_groups; i++) {
        topo_data->density_group_names[i] = new char[MAX_CG_TYPE_NAME_LENGTH + 1];
    }
	  
	printf("Read definitions for %d density groups.\n", topo_data->n_density_groups);
	for(int i=0; i < topo_data->n_density_groups; i++) {
		check_and_read_next_line(top_in, buff);
    	line_num++;
      	
		// Split line using string tools. 
		num_entries = StringSplit(buff, " \t", type_list);
		// First entry is the name of the CG group.
		sprintf(topo_data->density_group_names[i], "%s", type_list[0].c_str());
		// Add cg_type to density_group definition.
		for (int j=1; j < num_entries; j++) {
			cg_type = match_type(type_list[j], topo_data->name, cg->n_cg_types);
	        if( cg_type == -1) {
    	    	printf("Unrecognized type %s!\n", type_list[j].c_str());
    		} else {
    			// Convert this number for 1-based to 0-based.
    			cg_type--;
    		}
			
			// Check that read was successful for cg_type.
			// This assumes that cg_type is 1-based (not 0-based).
			if ( (cg_type >= 0) && (cg_type < cg->n_cg_types) ) {
			//if ( (cg_type > 0) && (cg_type <= cg->n_cg_types) ) {
				topo_data->density_groups[i*(cg->n_cg_types) + cg_type] = true;
			} else {
				printf("Out-of-range type number in density group definition: %s\n", type_list[j].c_str()); fflush(stdout);
			}
		}
	}
	delete [] type_list;
	return line_num;
}

int read_density_weights(TopologyData* topo_data, CG_MODEL_DATA* const cg, std::ifstream &top_in, int line_num)
{
	// Need to rewrite this so that it does not use iclass
	// Then need to pass information from topology into density
	// Allocate density_weights to match density_groups
	if(topo_data->n_density_groups <= 0) return line_num;

	topo_data->density_weights = new double[topo_data->n_density_groups * cg->n_cg_types];
	
	// Default initialization is for number density (all weights are 1.0)
	for (int i = 0; i < topo_data->n_density_groups * cg->n_cg_types; i++) {
		topo_data->density_weights[i] = 1.0;
	}
	
	// Check if that is all that needs to be done.
	if (cg->density_interactions.density_weights_flag == 0) return line_num;
	
	// Otherwise read speicific weights
	std::string buff;
	std::string* elements = new std::string[4];
	char parameter_name[50] = "";
	int num_elements, num_entries;
	double weight;
	
	// Check label
	check_and_read_next_line(top_in, buff);
    line_num++;	
	sscanf(buff.c_str(), "%s%d", parameter_name, &num_entries);
  	if (strcmp(parameter_name, "density_weights") != 0) report_topology_input_format_error(line_num, parameter_name);
	
	// There should be an entry for each type listed in each density group.
	// So, go through density_group array to find specified elements.
	printf("Read density weights.\n");
	for (int i = 0; i < topo_data->n_density_groups; i++) {
		for(int j = 0; j < cg->n_cg_types; j++) {
	
			// Read in value if this type (j) is part of density_group (i).
			// The first element in the line should be the density group.
			// The second element in the line should be the CG type.
			// The third element of the line should be the weight.
			if (topo_data->density_groups[i * cg->n_cg_types + j] == true) {
				check_and_read_next_line(top_in, buff);
    			line_num++;
      	
				// Split line using string tools. 
				if( (num_elements = StringSplit(buff, " \t\n", elements)) < 3 ) {
					printf("Not enough elements (%d) in lines of den_weights.in. Expected at least 3 elements.\n", num_elements);
					exit(EXIT_FAILURE);
				}
				
				// Only use weight if it is non-zero.
				// Weight could also be 0 if this element is non-numeric
				if ( (weight = atof(elements[2].c_str())) != 0.0) {
					topo_data->density_weights[i * cg->n_cg_types + j] = weight;
				} else {
					printf("Warning: Detected a density weight of 0 for type %d in density_group %s!\n", j, topo_data->density_group_names[i]);
					printf("Ignoring this weight and setting weight for this type to 1.0!\n");
				}
			}
		}
	}
	delete [] elements;
	return line_num;
}

void read_molecule_definition(TopologyData* const mol, TopologyData *topo_data, std::ifstream &top_in, int &line)
{
    unsigned i, j, hash_val;
    int automatic_topology_style = -5;
    std::string buff;
    char parameter_name[50];
    
    // Read in the first line of the single-molecule topology input section
    // Find the number of sites in the molecule and the style of topology
    // data included here. 'automatic_topology_style' is -1 if angles and
    // dihedrals are manually specified, 1 if only pair bonding is used,
    // 2 if angles are assigned automatically between bonded triples, and
    // 3 if dihedrals are assigned automatically between bonded quadruples.
    check_and_read_next_line(top_in, buff, line);
    sscanf(buff.c_str(), "%s%d%d", parameter_name, &mol->n_cg_sites, &automatic_topology_style);
    
    if (strcmp(parameter_name, "mol") != 0) report_topology_input_format_error(line, parameter_name);
    if (automatic_topology_style == -5) {
        printf("Wrong format in top.in: line %d, provide a valid topology inference style.\n", line);
        exit(EXIT_FAILURE);
    }
    
    mol->n_cg_types = 0;
    mol->excluded_style = 0;
	mol->density_excluded_style = 0;
    // Allocate space for bond, angle, and dihedral lists.
    initialize_topology_data(mol);
    
    // Read the the CG site type of each site in the CG molecule
    check_and_read_next_line(top_in, buff, line);
    sscanf(buff.c_str(), "%s", parameter_name);
    if (strcmp(parameter_name, "sitetypes") != 0) report_topology_input_format_error(line, parameter_name);
    for (i = 0; i < mol->n_cg_sites; i++) {
        check_and_read_next_line(top_in, buff, line);
        sscanf(buff.c_str(), "%d", &mol->cg_site_types[i]);
    }
    
    // Define all topology lists.        
    // Read in the number of pair bonds in the molecule
    unsigned n_pair_bonded_interactions;
    int cg_site1, cg_site2, cg_site3, cg_site4;
    int cg_type1, cg_type2, cg_type3, cg_type4;
    check_and_read_next_line(top_in, buff, line);
    sscanf(buff.c_str(), "%s%d", parameter_name, &n_pair_bonded_interactions);
    if (strcmp(parameter_name, "bonds") != 0) report_topology_input_format_error(line,  parameter_name);
    
    // Read in the pair bond topology.
    // The input file indexes particles from 1, this program
    // indexes them from zero.
    for (i = 0; i < n_pair_bonded_interactions; i++) {
        
        // Read the next line in the bond table.
        check_and_read_next_line(top_in, buff, line);
        sscanf(buff.c_str(), "%d%d", &cg_site1, &cg_site2);
        
        // Correct to index from zero.
        cg_site1--;
        cg_site2--;
        
        // Determine the bond type and activate it.
        cg_type1 = mol->cg_site_types[cg_site1];
        cg_type2 = mol->cg_site_types[cg_site2];
        hash_val = calc_two_body_interaction_hash(cg_type1, cg_type2, topo_data->n_cg_types);
        topo_data->bond_type_activation_flags[hash_val] = 1;
        
        // Add the bond to the bond table.
        mol->bond_list->partners_[cg_site1][mol->bond_list->partner_numbers_[cg_site1]] = cg_site2;
        mol->bond_list->partner_numbers_[cg_site1]++;
        mol->bond_list->partners_[cg_site2][mol->bond_list->partner_numbers_[cg_site2]] = cg_site1;
        mol->bond_list->partner_numbers_[cg_site2]++;
    }
    printf("Read pair bond topology; %d bonds of %d bond types.\n", n_pair_bonded_interactions, calc_n_active_interactions(topo_data->bond_type_activation_flags, calc_n_distinct_pairs(topo_data->n_cg_types)));
    
    if (automatic_topology_style == -1) {
        
        // If all topology info is set manually in the table, read in
        // the angle and dihedral tables.
        
        // Read the angle tables.
        unsigned n_angle_interactions;
        mol->angle_format = -1;
        check_and_read_next_line(top_in, buff, line);
        sscanf(buff.c_str(), "%s%d%d", parameter_name, &n_angle_interactions, &mol->angle_format);
        if (strcmp(parameter_name, "angles") != 0) report_topology_input_format_error(line, parameter_name);
        if (mol->angle_format == 0) {
			printf("Reading angle interactions for A-B and B-C as B A C.\n");
		} else if(mol->angle_format == 1) {
			printf("Reading angle interactions for A-B and B-C as A B C.\n");
			
		} else {
        	fprintf(stderr, "WARNING: 2 numbers are now required after the angles tag in top.in.\n");
        	fprintf(stderr, "To specify the angle of A-B and B-C as B A C choose option 0 or as A B C choose option 1.\n");
        	fprintf(stderr, "Setting option to 0 by default.\n");
        	mol->angle_format = 0;
		}
		
        for (i = 0; i < n_angle_interactions; i++) {
            
            // Read an angle definition.
            check_and_read_next_line(top_in, buff, line);
            sscanf(buff.c_str(), "%d%d%d", &cg_site1, &cg_site2, &cg_site3);
            
            // Re-index from zero.
            cg_site1--;
            cg_site2--;
            cg_site3--;
            
            // Determine bond type.
            cg_type1 = mol->cg_site_types[cg_site1];
            cg_type2 = mol->cg_site_types[cg_site2];
            cg_type3 = mol->cg_site_types[cg_site3];
            
            if (mol->angle_format == 1) { // Swap cg_type1 and cg_type2 to look like angle_format 0.
            	swap_pair(cg_type1, cg_type2); // (start) A B C => B A C (end)
            	swap_pair(cg_site1, cg_site2); // mirror for sites (as well as types)
            }
            hash_val = calc_three_body_interaction_hash(cg_type1, cg_type2, cg_type3, topo_data->n_cg_types);
            topo_data->angle_type_activation_flags[hash_val] = 1;
            
            // Add to the list.
            mol->angle_list->partners_[cg_site2][2 * mol->angle_list->partner_numbers_[cg_site2]] = cg_site1;
            mol->angle_list->partners_[cg_site2][2 * mol->angle_list->partner_numbers_[cg_site2] + 1] = cg_site3;
            mol->angle_list->partner_numbers_[cg_site2]++;
            mol->angle_list->partners_[cg_site3][2 * mol->angle_list->partner_numbers_[cg_site3]] = cg_site1;
            mol->angle_list->partners_[cg_site3][2 * mol->angle_list->partner_numbers_[cg_site3] + 1] = cg_site2;
            mol->angle_list->partner_numbers_[cg_site3]++;
        }
        printf("Read manually specified angle topology; %d angles of %d angle types.\n", n_angle_interactions, calc_n_active_interactions(topo_data->angle_type_activation_flags, calc_n_distinct_triples(topo_data->n_cg_types)));

        
        // Read the dihedral tables.
        unsigned n_dihedral_interactions;
        mol->dihedral_format = -1;
        check_and_read_next_line(top_in, buff, line);
        sscanf(buff.c_str(), "%s%d%d", parameter_name, &n_dihedral_interactions, &mol->dihedral_format);
        if (strcmp(parameter_name, "dihedrals") != 0) report_topology_input_format_error(line, parameter_name);
        if (mol->dihedral_format == 0) {
			printf("Reading angle interactions for A-B, B-C, and C-D as B C A D.\n");
		} else if (mol->dihedral_format == 1) {
			printf("Reading angle interactions for A-B, B-C, and C-D as A B C D.\n");
		} else {
        	fprintf(stderr, "WARNING: 2 numbers are now required after the dihedral tag in top.in.\n");
        	fprintf(stderr, "To specify the dihedral of A-B, B-C, and C-D as B C A D choose option 0 or as A B C D choose option 1.\n");
			fprintf(stderr, "Setting option to 0.\n");
			mol->dihedral_format = 0;
		}
		
        for (i = 0; i < n_dihedral_interactions; i++) {
            
            // Read a single dihedral definition.
            check_and_read_next_line(top_in, buff, line);
            sscanf(buff.c_str(), "%d%d%d%d", &cg_site1, &cg_site2, &cg_site3, &cg_site4);
            
            // Re-index from zero.
            cg_site1--;
            cg_site2--;
            cg_site3--;
            cg_site4--;
            
            // Determine bond type.
            cg_type1 = mol->cg_site_types[cg_site1];
            cg_type2 = mol->cg_site_types[cg_site2];
            cg_type3 = mol->cg_site_types[cg_site3];
            cg_type4 = mol->cg_site_types[cg_site4];
            
            if (mol->dihedral_format == 1) { // rearrange indices so that it looks like dihedral_style 0.
            	swap_pair(cg_type1, cg_type2); // (start) A B C D => B A C D
            	swap_pair(cg_type2, cg_type3); // B A C D => B C A D (done)
            	
            	swap_pair(cg_site1, cg_site2); // mirror for sites as well
            	swap_pair(cg_site2, cg_site3);
            }
            
            hash_val = calc_four_body_interaction_hash(cg_type1, cg_type2, cg_type3, cg_type4, topo_data->n_cg_types);
            topo_data->dihedral_type_activation_flags[hash_val] = 1;
            
            // Add the dihedral to the table.
            mol->dihedral_list->partners_[cg_site3][3 * mol->dihedral_list->partner_numbers_[cg_site3]] = cg_site1;
            mol->dihedral_list->partners_[cg_site3][3 * mol->dihedral_list->partner_numbers_[cg_site3] + 1] = cg_site2;
            mol->dihedral_list->partners_[cg_site3][3 * mol->dihedral_list->partner_numbers_[cg_site3] + 2] = cg_site4;
            mol->dihedral_list->partner_numbers_[cg_site3]++;
            mol->dihedral_list->partners_[cg_site4][3 * mol->dihedral_list->partner_numbers_[cg_site4]] = cg_site2;
            mol->dihedral_list->partners_[cg_site4][3 * mol->dihedral_list->partner_numbers_[cg_site4] + 1] = cg_site1;
            mol->dihedral_list->partners_[cg_site4][3 * mol->dihedral_list->partner_numbers_[cg_site4] + 2] = cg_site3;
            mol->dihedral_list->partner_numbers_[cg_site4]++;
        }
        printf("Read manually specified dihedral topology; %d dihedrals of %d dihedral types.\n", n_dihedral_interactions, calc_n_active_interactions(topo_data->dihedral_type_activation_flags, calc_n_distinct_quadruples(topo_data->n_cg_types)));
        
    } else {
        
        // If the automatic_topology_style specifies that this topology
        // reader is responsible for inferring angular and dihedral
        // interactions from pair bond interactions, then infer those
        // interactions as necessary.
        unsigned l, n_angles;
        
        // Loop over CG sites in the molecule
        for (cg_site1 = 0; cg_site1 < (int)(mol->n_cg_sites); cg_site1++) {
        	cg_type1 = mol->cg_site_types[cg_site1];
        	n_angles = 0;
            // Loop over sites bonded to that molecule.
            for (j = 0; j < mol->bond_list->partner_numbers_[cg_site1]; j++) {
                
                // Overdetermine this bond style's activation flag.
                cg_site2 = mol->bond_list->partners_[cg_site1][j];
                cg_type2 = mol->cg_site_types[cg_site2];
                if (cg_type1 <= cg_type2) {
                    hash_val = calc_two_body_interaction_hash(cg_type1, cg_type2, topo_data->n_cg_types);
                    topo_data->bond_type_activation_flags[hash_val] = 1;
                }
                
                // That's enough for some styles.
                if (automatic_topology_style == 1) continue;
                
                // Otherwise, check the sites bonded to that one
                // and add consider adding the angle to the angle table.
                for (l = 0; l < mol->bond_list->partner_numbers_[cg_site2]; l++) {
                    
                    cg_site3 = mol->bond_list->partners_[cg_site2][l];
                    cg_type3 = mol->cg_site_types[cg_site3];
                    if (unsigned(cg_site3) == cg_site1) continue;
                    
                    mol->angle_list->partners_[cg_site1][2 * n_angles]     = cg_site2;
                    mol->angle_list->partners_[cg_site1][2 * n_angles + 1] = cg_site3;
                    n_angles++;
                    
                    hash_val = calc_three_body_interaction_hash(cg_type2, cg_type1, cg_type3, topo_data->n_cg_types);
                    topo_data->angle_type_activation_flags[hash_val] = 1;
                }
            }
            mol->angle_list->partner_numbers_[cg_site1] = n_angles;
        }
        if (automatic_topology_style != 1) {
	        int total_angles=0;
	        for (i = 0; i < mol->n_cg_sites; i++) {
	            total_angles += mol->angle_list->partner_numbers_[i];
	        }
            printf("Automatically generated angle topology; %d angles of %d angle types.\n", total_angles/2, calc_n_active_interactions(topo_data->angle_type_activation_flags, calc_n_distinct_triples(topo_data->n_cg_types)));
        }
        
        unsigned k, n_dihedrals;
        
        // Loop over CG sites in the molecule 
        for (cg_site1 = 0; cg_site1 < (int)(mol->n_cg_sites); cg_site1++) {
            cg_type1 = mol->cg_site_types[cg_site1];
            
            // For some styles of automatic topology generation, every
            // iteration of this loop should be skipped.
            if (automatic_topology_style <= 2) continue;
            
            // Infer dihedrals for this site using a combination of the 
            // angle topology and the pair bond topology.
            n_dihedrals = 0;
            for (j = 0; j < mol->angle_list->partner_numbers_[cg_site1]; j++) {
            	cg_site2 = mol->angle_list->partners_[cg_site1][2 * j];
                cg_site3 = mol->angle_list->partners_[cg_site1][2 * j + 1];  
                cg_type2 = mol->cg_site_types[cg_site2];
                cg_type3 = mol->cg_site_types[cg_site3];
                for (k = 0; k < mol->bond_list->partner_numbers_[cg_site3]; k++) {
                    cg_site4 = mol->bond_list->partners_[cg_site3][k];
                    cg_type4 = mol->cg_site_types[cg_site4];
                    if (cg_site4 == cg_site2) continue;
                    if (cg_site4 == cg_site1) continue;
                    mol->dihedral_list->partners_[cg_site1][3 * n_dihedrals] = cg_site2;
                    mol->dihedral_list->partners_[cg_site1][3 * n_dihedrals + 1] = cg_site3;
                    mol->dihedral_list->partners_[cg_site1][3 * n_dihedrals + 2] = cg_site4;
                    hash_val = calc_four_body_interaction_hash(cg_type2, cg_type3, cg_type1, cg_type4, topo_data->n_cg_types);
                    topo_data->dihedral_type_activation_flags[hash_val] = 1;
                    n_dihedrals++;
                }
            }
            mol->dihedral_list->partner_numbers_[cg_site1] = n_dihedrals;
        }
        if (automatic_topology_style > 2) {
            int total_dihedrals=0;
	        for (i = 0; i < mol->n_cg_sites; i++) {
	            total_dihedrals += mol->dihedral_list->partner_numbers_[i];
	        }
            printf("Automatically generated dihedral topology; %d dihedrals of %d dihedral types.\n", total_dihedrals/2, calc_n_active_interactions(topo_data->dihedral_type_activation_flags, calc_n_distinct_quadruples(topo_data->n_cg_types)));
        }
    }
    
    // Set-up quint list by connecting angles with angles that do not go "back" on themselves
    // Loop over CG sites in the molecule
    int cg_site5, cg_type5;
    for (cg_site1 = 0; cg_site1 < (int)(mol->n_cg_sites); cg_site1++) {
    	cg_type1 = mol->cg_site_types[cg_site1];
    	int n_quints = 0;
        // Loop over angles involving that CG site.
        for (j = 0; j < mol->angle_list->partner_numbers_[cg_site1]; j++) {
        	// Then look at angles involving that partner (j).
        	cg_site2 = mol->angle_list->partners_[cg_site1][2 * j    ]; // The sites in this angle are mol->angle_list->partner_[i][2 * j] 
        	cg_site3 = mol->angle_list->partners_[cg_site1][2 * j + 1]; //                         and mol->angle_list->partner_[i][2 * j + 1].
			cg_type2 = mol->cg_site_types[cg_site2];
			cg_type3 = mol->cg_site_types[cg_site3];
			for (unsigned k = 0; k < mol->angle_list->partner_numbers_[mol->angle_list->partners_[cg_site1][2 * j + 1]]; k++) {
				// Make sure that this angle is not folding back on the first angle
				cg_site4 = mol->angle_list->partners_[cg_site3][2 * k];
				cg_site5 = mol->angle_list->partners_[cg_site3][2 * k + 1]; 
				cg_type4 = mol->cg_site_types[cg_site4];
				cg_type5 = mol->cg_site_types[cg_site5];
				// // The end of the 2nd angle can not match a site in the 1st angle.
				if (cg_site5 == cg_site1) continue;
				if (cg_site5 == cg_site2) continue;
				// // The middle of the 2nd angle can not match a site in the 1st angle.
				if (cg_site4 == cg_site1) continue;
				if (cg_site4 == cg_site2) continue;
				
				// Add this quint to this partner list
				mol->quint_list->partners_[cg_site1][4 * n_quints    ] = cg_site2;
				mol->quint_list->partners_[cg_site1][4 * n_quints + 1] = cg_site3;
				mol->quint_list->partners_[cg_site1][4 * n_quints + 2] = cg_site4;
				mol->quint_list->partners_[cg_site1][4 * n_quints + 3] = cg_site5;
				
				// Set the activation flag for this site. 
				hash_val = calc_five_body_interaction_hash(cg_type2, cg_type3, cg_type4, cg_type1, cg_type5, topo_data->n_cg_types);
				topo_data->quint_type_activation_flags[hash_val] = 1;
				n_quints++;
			}
		}
		// set actual number of quints for this site
		mol->quint_list->partner_numbers_[cg_site1] = n_quints;
    }
}

void recursive_exclusion_search(TopologyData const* topo_data, TopoList* &exclusion_list, std::vector<int> &path_list)
{
	int path_size = path_list.size();
	int base_id = path_list[0];
	int last_id = path_list[path_size - 1];
	int bond_partner_size = topo_data->bond_list->partner_numbers_[ path_list[path_size - 1] ];
	int partner_id;
	
	for (int i = 0; i < bond_partner_size; i++) {
		partner_id = topo_data->bond_list->partners_[last_id][i];
		
		// Check if this site is already excluded
		if ((base_id != partner_id) &&
		    (new_excluded_partner(exclusion_list->partners_, base_id, exclusion_list->partner_numbers_[base_id], partner_id) == true)) {
	        exclusion_list->partners_[base_id][ exclusion_list->partner_numbers_[base_id] ] = partner_id;
    	    exclusion_list->partner_numbers_[base_id]++;
 	        
            path_list.push_back(partner_id);
            recursive_exclusion_search(topo_data, exclusion_list, path_list);
            path_list.pop_back();
		}
	}
}

void setup_excluded_list( TopologyData const* topo_data, TopoList* &exclusion_list, const int excluded_style) 
{
	// Automatically determine topology to set appropriate
    // bond, angle, and/or dihedral exclusion as appropriate.
	printf("Setting up exclusion list for excluded_style %d.\n", excluded_style);
    if (excluded_style == 0) return;
	
	unsigned max_excluded_number = get_max_exclusion_number(topo_data, excluded_style);
       
    if (excluded_style == 5) {
    	for (unsigned i = 0; i < topo_data->n_cg_sites; i++) {
    		std::vector<int> path_list = {(int)(i)};
    		recursive_exclusion_search(topo_data, exclusion_list, path_list);
    		
    		if ( exclusion_list->partner_numbers_[i] > max_excluded_number ) {
    			printf("Warning: Too many excluded interactions (%d) to handle based on max_pair_bonds_per_site (%d)!\n", exclusion_list->partner_numbers_[i], max_excluded_number);
    			exit(EXIT_FAILURE);
    		}
		}
    	return;	
    }
            
    // Loop over each CG site and its bonds to copy into exclusion list
    for (unsigned i = 0; i < topo_data->n_cg_sites; i++) {
        for (unsigned j = 0; j < topo_data->bond_list->partner_numbers_[i]; j++) {            
            exclusion_list->partners_[i][ exclusion_list->partner_numbers_[i] ] = topo_data->bond_list->partners_[i][j];
            exclusion_list->partner_numbers_[i]++;
        }
    
    	if ( exclusion_list->partner_numbers_[i] > max_excluded_number ) {
    		printf("Warning: Too many excluded interactions (%d) to handle based on max_pair_bonds_per_site (%d)!\n", exclusion_list->partner_numbers_[i], max_excluded_number);
    		exit(EXIT_FAILURE);
    	}
    }
    
    // If we also need angles, check the sites bonded to known bonds and
    // if we need dihedrals, check the sites bonded to know angles to add to exclusion list.
    if (excluded_style == 2) return;
    
    // Loop over CG sites in the molecule
    for (unsigned i = 0; i < topo_data->n_cg_sites; i++) {
        
        // Loop over sites (cg_site1) bonded to that CG site (i).
        for (unsigned j = 0; j < topo_data->bond_list->partner_numbers_[i]; j++) {
        
            unsigned cg_site1 = topo_data->bond_list->partners_[i][j];
    
    		// Loop over potential angles (cg_site2: bonded to cg_site1 which is bonded to i)
            for (unsigned k = 0; k < topo_data->bond_list->partner_numbers_[cg_site1]; k++) {
                    
                unsigned cg_site2 = topo_data->bond_list->partners_[cg_site1][k];
                
                //add this as an angle interaction as long as it does not loop back to the starting CG site
                if (cg_site2 == i) continue; 
                if ( new_excluded_partner(exclusion_list->partners_, i, exclusion_list->partner_numbers_[i], cg_site2) == true ) {
	                exclusion_list->partners_[i][ exclusion_list->partner_numbers_[i] ] = cg_site2;
    	            exclusion_list->partner_numbers_[i]++;
        		}
        		        
                //if we need dihedrals, loop over all potential dihedrals (cg_site3: bonded to angle)
                if (excluded_style == 4) {
                	for(unsigned l = 0; l < topo_data->bond_list->partner_numbers_[cg_site2]; l++) {
                		unsigned cg_site3 = topo_data->bond_list->partners_[cg_site2][l];
                		
                		//add this as a dihedral interaction as long as it does not loop back to a CG site in the existing angle
                		if( (cg_site3 == i) || (cg_site3 == cg_site1) || (cg_site3 == cg_site2) ) continue;
                		if ( new_excluded_partner(exclusion_list->partners_, i, exclusion_list->partner_numbers_[i], cg_site3) == true ) {
                			exclusion_list->partners_[i][ exclusion_list->partner_numbers_[i] ] = cg_site3;
               		 		exclusion_list->partner_numbers_[i]++;
               		 	}
               		}
            	}
            
        	}
        }
        
    	if ( exclusion_list->partner_numbers_[i] > max_excluded_number ) {
    		printf("Warning: Too many excluded interactions (%d) to handle based on max_pair_bonds_per_site, max_angles_per_site, and max_dihedrals_per_site (total %d)!\n", exclusion_list->partner_numbers_[i], max_excluded_number);
    		exit(EXIT_FAILURE);
    	}
    }
}

void setup_molecule_list(TopologyData const* topo_data, TopoList* &molecule_list, const int n_molecules, const int max_molecule_size)
{
	if (molecule_list != NULL) {
		delete molecule_list;
	}
	molecule_list = new TopoList(n_molecules, 1, max_molecule_size);
	int my_id;
	for (unsigned i = 0; i < topo_data->n_cg_sites; i++) {
		my_id = topo_data->molecule_ids[i];
		molecule_list->partners_[my_id][molecule_list->partner_numbers_[my_id]] = i;
		molecule_list->partner_numbers_[my_id]++;
	}
}

void update_molecule_list(TopologyData const* topo_data, TopoList* &molecule_list)
{
	// Check number of molecules
	int n_molecules = molecule_list->n_sites_;
	for (unsigned i = 0; i < topo_data->n_cg_sites; i++) {
		if (topo_data->molecule_ids[i] > n_molecules) {
			n_molecules = topo_data->molecule_ids[i];
		}
	}
	
	// Determine max molecule size
	int max_molecule_size = 1;
	int* molecule_size = new int[n_molecules]();
	for (int i = 0; i < n_molecules; i++) {
		molecule_size[i] = 0;
	}
	for (unsigned i = 0; i < topo_data->n_cg_sites; i++) {
		molecule_size[topo_data->molecule_ids[i]]++;
	}
	for (int i = 0; i  < n_molecules; i++) {
		if(molecule_size[i] > max_molecule_size) {
			max_molecule_size = molecule_size[i];
		}
	}
	delete [] molecule_size;
	
	// Setup new molecule list
	setup_molecule_list(topo_data, molecule_list, n_molecules, max_molecule_size);
}

bool new_excluded_partner(unsigned** const excluded_partners, const int location, const unsigned current_size, const unsigned new_site) {
	for(unsigned i = 0; i < current_size; i++) {
		if (excluded_partners[location][i] == new_site)	return false;
	}
	return true;
}

void report_topology_input_format_error(const int line, char *parameter_name)
{
    printf("Wrong format in top.in: line %d, unexpected parameter %s.\n", line, parameter_name);
    exit(EXIT_FAILURE);
}