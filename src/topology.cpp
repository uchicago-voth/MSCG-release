//
//  topology.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cstdlib>
#include <cstring>
#include <cassert>

#include "interaction_model.h"
#include "misc.h"

#include "topology.h"

// Helper functions for reading a topology file for the CG model
unsigned get_max_exclusion_number(TopologyData const* const topo_data);

// Read the specification for three-body interactions.
int read_three_body_topology(TopologyData* topo_data, CG_MODEL_DATA* const cg, FILE* top_in, int line);
// Read a single molecule's topology specification.
void read_molecule_definition(TopologyData* const mol, TopologyData* const cg, FILE* const top_in, int* line);
// Determine appropriate non-bonded exclusions based on bonded topology and exclusion_style setting
bool new_excluded_partner(unsigned** const excluded_partners, const int location, const unsigned current_size, const unsigned new_site);
// Report an error in the topology input file.
void report_topology_input_format_error(const int line, char *parameter_name);

//---------------------------------------------------------------
// Functions for managing TopoList structs
//---------------------------------------------------------------

TopoList::TopoList(unsigned n_sites, unsigned partners_per, unsigned max_partners) : 
    n_sites_(n_sites), partners_per_(partners_per), max_partners_(max_partners) {
    partner_numbers_ = new unsigned[n_sites_]();
    partners_ = new unsigned*[n_sites_];
    modified = 0;
    for (unsigned i = 0; i < n_sites_; i++) {
        partners_[i] = new unsigned[partners_per_ * max_partners_];
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

unsigned get_max_exclusion_number(TopologyData const* const topo_data) {
	assert( (topo_data->excluded_style == 0) || (topo_data->excluded_style == 2) || (topo_data->excluded_style == 3) || (topo_data->excluded_style == 4) );
	switch (topo_data->excluded_style) {
		case 2:
			return topo_data->max_pair_bonds_per_site;
			break;
			
		case 3:
			return topo_data->max_pair_bonds_per_site + topo_data->max_angles_per_site;
			break;
			
		case 4:
			return topo_data->max_pair_bonds_per_site + topo_data->max_angles_per_site + topo_data->max_dihedrals_per_site;
			break;
			
		default:
			return 0;
			break;
	}
}

// Initialize a TopologyData struct given that the number of sites and maximum
// numbers of bonds, angles, and dihedrals for each site are fixed.
void initialize_topology_data(TopologyData* const topo_data)
{
	unsigned nsites = topo_data->n_cg_sites;
    topo_data->cg_site_types = new int[nsites];

    topo_data->name = new char*[topo_data->n_cg_types];
    for (unsigned i = 0; i < topo_data->n_cg_types; i++) {
        topo_data->name[i] = new char[MAX_CG_TYPE_NAME_LENGTH + 1];
    }
    // Ready pair bond data structures.
    topo_data->bond_type_activation_flags = new int[calc_n_distinct_pairs(topo_data->n_cg_types)]();
    topo_data->bond_list = new TopoList(nsites, 1, topo_data->max_pair_bonds_per_site);
    // Ready angle bond data structures.
    topo_data->angle_type_activation_flags = new int[calc_n_distinct_triples(topo_data->n_cg_types)]();
    topo_data->angle_list = new TopoList(nsites, 2, topo_data->max_angles_per_site);
    // Ready dihedral bond data structures.
    topo_data->dihedral_type_activation_flags = new int[calc_n_distinct_quadruples(topo_data->n_cg_types)]();
    topo_data->dihedral_list = new TopoList(nsites, 3, topo_data->max_dihedrals_per_site);
    // Ready exclusion list data structures.
    if( topo_data->excluded_style != 0) topo_data->exclusion_list = new TopoList(nsites, 1, get_max_exclusion_number(topo_data));
}

// Free an initialized TopologyData struct.
void free_topology_data(TopologyData *topo_data)
{
    delete [] topo_data->cg_site_types;
    if (topo_data->name != NULL) {
        for (unsigned i = 0; i < topo_data->n_cg_types; i++) {
            if (topo_data->name[i] != NULL) delete [] topo_data->name[i];
        }
        delete [] topo_data->name;
    }

    // Delete the topology lists.
    if (topo_data->bond_list != NULL) delete topo_data->bond_list;
    if (topo_data->angle_list != NULL) delete topo_data->angle_list;
    if (topo_data->dihedral_list != NULL) delete topo_data->dihedral_list;
    if (topo_data->excluded_style != 0) delete topo_data->exclusion_list;

    // Delete the topology interaction type activation flags.
    if (topo_data->bond_type_activation_flags != NULL) delete [] topo_data->bond_type_activation_flags;
    if (topo_data->angle_type_activation_flags != NULL) delete [] topo_data->angle_type_activation_flags;
    if (topo_data->dihedral_type_activation_flags != NULL) delete [] topo_data->dihedral_type_activation_flags;
}

//---------------------------------------------------------------
// Functions for reading a topology file for the CG model
//---------------------------------------------------------------

void read_topology_file(TopologyData* topo_data, CG_MODEL_DATA* const cg)
{
    unsigned i;
    char buff[100], parameter_name[50];
    FILE* top_in;
    top_in = open_file("top.in", "r");
    int line = 0;
    
    // Read the number of sites and number of types of sites.
    fgets(buff, 100, top_in);
    line++;
    sscanf(buff, "%s%d", parameter_name, &topo_data->n_cg_sites);
    if (strcmp(parameter_name, "cgsites") != 0) report_topology_input_format_error(line,parameter_name);
    cg->n_cg_sites = int(topo_data->n_cg_sites);
    fgets(buff, 100, top_in);
    line++;
    sscanf(buff, "%s%d", parameter_name, &topo_data->n_cg_types);
    if (strcmp(parameter_name, "cgtypes") != 0) report_topology_input_format_error(line,parameter_name);
    cg->n_cg_types = int(topo_data->n_cg_types);

    // Allocate bond, angle, and dihedral interaction
    // class arrays indexed by all possible such interactions
    // using that information.
    initialize_topology_data(topo_data);
    cg->name = new char*[cg->n_cg_types];
    for (i = 0; i < unsigned(cg->n_cg_types); i++) {
        cg->name[i] = new char[MAX_CG_TYPE_NAME_LENGTH + 1];
    }
    
    // Read the string names of each site type.
    for (i = 0; i < topo_data->n_cg_types; i++) {
        fgets(buff, 100, top_in);
        line++;
        sscanf(buff, "%s", topo_data->name[i]);
    }
    
    // Read the section of top.in defining three body nonbonded interactions.
    if (cg->three_body_nonbonded_interactions.class_subtype > 0) {
        line = read_three_body_topology(topo_data, cg, top_in, line);
    }
    
    // Read the molecules section of top.in to define all topology lists.
    
    // Read the number of distinct types of molecule
    fgets(buff, 100, top_in);
    line++;
    unsigned n_molecule_types;
    sscanf(buff, "%s%d", parameter_name, &n_molecule_types);
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
    for (i = 0; i < n_molecule_types; i++) {
        printf("Reading definition of molecule type %d.\n", i);
        read_molecule_definition(&mol[i], topo_data, top_in, &line);
    }
    
    // Read the system section of top.in to define the full system's
    // topology.
    
    unsigned j, k, l;
    unsigned block, molecule_type, n_molecules_in_system, first_site_in_molecule = 0;
    fgets(buff, 100, top_in);
    line++;
    sscanf(buff, "%s%d", parameter_name, &block);
    if (strcmp(parameter_name, "system") != 0) report_topology_input_format_error(line, parameter_name);
    
    // For each line in this section, add a given number of molecules
    // of a given type to the system in order.
    for (i = 0; i < block; i++) {
        // Get a number and type of molecules.
        fgets(buff, 100, top_in);
        line++;
        sscanf(buff, "%d%d", &molecule_type, &n_molecules_in_system);
        molecule_type--;
        
        assert(n_molecules_in_system <= unsigned(cg->n_cg_sites));
        
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
            }
            first_site_in_molecule += mol[molecule_type].n_cg_sites;
        }
    }
    
    // Give the CG_MODEL_DATA struct its own copy of the type names.
    for (i = 0; i < topo_data->n_cg_types; i++) {
        sscanf(topo_data->name[i], "%s", cg->name[i]);
    }

	// Set-up appropriate bonded exclusions list from non-bonded interactions based on exclusion_style
	setup_excluded_list(topo_data);

    // Close the file and free the memory used to store the temporary
    // single-molecule topologies.
    fclose(top_in);
    for (i = 0; i < n_molecule_types; i++) {
        free_topology_data(&mol[i]);
    }
    delete [] mol;
}

int read_three_body_topology(TopologyData* topo_data, CG_MODEL_DATA* const cg, FILE* top_in, int line)
{
	char buff[100], parameter_name[50];
	unsigned tbtype, i;
	int* tb_i;
	int* tb_j;
	int* tb_k;
	
	fgets(buff, 100, top_in);
	line++;
	sscanf(buff, "%s%d", parameter_name, &tbtype);
	if (strcmp(parameter_name, "threebody") != 0) report_topology_input_format_error(line, parameter_name);

	cg->three_body_nonbonded_interactions.stillinger_weber_angle_parameters_by_type = new double[tbtype];
	cg->three_body_nonbonded_interactions.three_body_nonbonded_cutoffs = new double[tbtype];

	tb_i = new int[tbtype];
	tb_j = new int[tbtype];
	tb_k = new int[tbtype];
	
	for (i = 0; i < tbtype; i++) {
		fgets(buff, 100, top_in);
		line++;
		sscanf(buff, "%d%d%d%lf%lf", tb_i + i, tb_j + i, tb_k + i, cg->three_body_nonbonded_interactions.stillinger_weber_angle_parameters_by_type + i, cg->three_body_nonbonded_interactions.three_body_nonbonded_cutoffs + i);
		tb_i[i]--;
		tb_j[i]--;
		tb_k[i]--;
	}
	
	cg->tb_n = new int[topo_data->n_cg_types]();
	cg->tb_list = new int*[topo_data->n_cg_types];
	
	for (i = 0; i < tbtype; i++) {
		cg->tb_n[tb_i[i]]++;
	}
	
	for (i = 0; i < topo_data->n_cg_types; i++) {
		cg->tb_list[i] = new int[cg->tb_n[i] * 2];
		cg->tb_n[i] = 0;
	}
	
	for (i = 0; i < tbtype; i++) {
		cg->tb_list[tb_i[i]][cg->tb_n[tb_i[i]] * 2] = tb_j[i] + 1;
		cg->tb_list[tb_i[i]][cg->tb_n[tb_i[i]] * 2 + 1] = tb_k[i] + 1;
		cg->tb_n[tb_i[i]]++;
	}

	cg->three_body_nonbonded_interactions.set_n_defined(tbtype);
	delete [] tb_i;
	delete [] tb_j;
	delete [] tb_k;
	
	return line;
}

void read_molecule_definition(TopologyData* const mol, TopologyData *topo_data, FILE* const top_in, int* line)
{
    unsigned i, j;
    int automatic_topology_style = -5;
    char buff[100], parameter_name[50];
    
    // Read in the first line of the single-molecule topology input section
    // Find the number of sites in the molecule and the style of topology
    // data included here. 'automatic_topology_style' is -1 if angles and
    // dihedrals are manually specified, 1 if only pair bonding is used,
    // 2 if angles are assigned automatically between bonded triples, and
    // 3 if dihedrals are assigned automatically between bonded quadruples.
    fgets(buff, 100, top_in);
    (*line)++;
    sscanf(buff, "%s%d%d", parameter_name, &mol->n_cg_sites, &automatic_topology_style);
    
    if (strcmp(parameter_name, "mol") != 0) report_topology_input_format_error(*line, parameter_name);
    if (automatic_topology_style == -5) {
        printf("Wrong format in top.in: line %d, provide a valid topology inference style.\n", *line);
        exit(EXIT_FAILURE);
    }
    
    mol->n_cg_types = 0;
    mol->excluded_style = 0;
    // Allocate space for bond, angle, and dihedral lists.
    initialize_topology_data(mol);
    
    // Read the the CG site type of each site in the CG molecule
    fgets(buff, 100, top_in);
    (*line)++;
    sscanf(buff, "%s", parameter_name);
    if (strcmp(parameter_name, "sitetypes") != 0) report_topology_input_format_error(*line, parameter_name);
    for (i = 0; i < mol->n_cg_sites; i++) {
        fgets(buff, 100, top_in);
        (*line)++;
        sscanf(buff, "%d", &mol->cg_site_types[i]);
    }
    
    // Define all topology lists.        
    // Read in the number of pair bonds in the molecule
    unsigned n_pair_bonded_interactions;
    int cg_site1, cg_site2, cg_site3, cg_site4;
    int cg_type1, cg_type2, cg_type3, cg_type4;
    int hash_val;
    fgets(buff, 100, top_in);
    (*line)++;
    sscanf(buff, "%s%d", parameter_name, &n_pair_bonded_interactions);
    if (strcmp(parameter_name, "bonds") != 0) report_topology_input_format_error(*line,  parameter_name);
    
    // Read in the pair bond topology.
    // The input file indexes particles from 1, this program
    // indexes them from zero.
    for (i = 0; i < n_pair_bonded_interactions; i++) {
        
        // Read the next line in the bond table.
        fgets(buff, 100, top_in);
        (*line)++;
        sscanf(buff, "%d%d", &cg_site1, &cg_site2);
        
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
        fgets(buff, 100, top_in);
        (*line)++;
        sscanf(buff, "%s%d", parameter_name, &n_angle_interactions);
        if (strcmp(parameter_name, "angles") != 0) report_topology_input_format_error(*line, parameter_name);
        
        for (i = 0; i < n_angle_interactions; i++) {
            
            // Read an angle definition.
            fgets(buff, 100, top_in);
            (*line)++;
            sscanf(buff, "%d%d%d", &cg_site1, &cg_site2, &cg_site3);
            
            // Re-index from zero.
            cg_site1--;
            cg_site2--;
            cg_site3--;
            
            // Determine bond type.
            cg_type1 = mol->cg_site_types[cg_site1];
            cg_type2 = mol->cg_site_types[cg_site2];
            cg_type3 = mol->cg_site_types[cg_site3];
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
        fgets(buff, 100, top_in);
        (*line)++;
        sscanf(buff, "%s%d", parameter_name, &n_dihedral_interactions);
        if (strcmp(parameter_name, "dihedrals") != 0) report_topology_input_format_error(*line, parameter_name);
        for (i = 0; i < n_dihedral_interactions; i++) {
            
            // Read a single dihedral definition.
            fgets(buff, 100, top_in);
            (*line)++;
            sscanf(buff, "%d%d%d%d", &cg_site1, &cg_site2, &cg_site3, &cg_site4);
            
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
        for (i = 0; i < mol->n_cg_sites; i++) {
            n_angles = 0;
            // Loop over sites bonded to that molecule.
            for (j = 0; j < mol->bond_list->partner_numbers_[i]; j++) {
                
                // Overdetermine this bond style's activation flag.
                cg_site1 = mol->bond_list->partners_[i][j];
                if (mol->cg_site_types[i] <= mol->cg_site_types[cg_site1]) {
                    hash_val = calc_two_body_interaction_hash(mol->cg_site_types[i], mol->cg_site_types[cg_site1], topo_data->n_cg_types);
                    topo_data->bond_type_activation_flags[hash_val] = 1;
                }
                
                // That's enough for some styles.
                if (automatic_topology_style == 1) continue;
                
                // Otherwise, check the sites bonded to that one
                // and add consider adding the angle to the angle table.
                for (l = 0; l < mol->bond_list->partner_numbers_[cg_site1]; l++) {
                    
                    cg_site2 = mol->bond_list->partners_[cg_site1][l];
                    if (unsigned(cg_site2) == i) continue;
                    
                    mol->angle_list->partners_[i][2 * n_angles] = cg_site1;
                    mol->angle_list->partners_[i][2 * n_angles + 1] = cg_site2;
                    n_angles++;
                    
                    cg_type1 = mol->cg_site_types[cg_site1];
                    cg_type2 = mol->cg_site_types[i];
                    cg_type3 = mol->cg_site_types[cg_site2];
                    hash_val = calc_three_body_interaction_hash(cg_type1, cg_type2, cg_type3, topo_data->n_cg_types);
                    topo_data->angle_type_activation_flags[hash_val] = 1;
                }
            }
            mol->angle_list->partner_numbers_[i] = n_angles;
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
        for (i = 0; i < mol->n_cg_sites; i++) {
            
            // For some styles of automatic topology generation, every
            // iteration of this loop should be skipped.
            if (automatic_topology_style <= 2) continue;
            
            // Infer dihedrals for this site using a combination of the 
            // angle topology and the pair bond topology.
            n_dihedrals = 0;
            for (j = 0; j < mol->angle_list->partner_numbers_[i]; j++) {
                for (k = 0; k < mol->bond_list->partner_numbers_[mol->angle_list->partners_[i][2 * j + 1]]; k++) {
                    if (mol->bond_list->partners_[mol->angle_list->partners_[i][2 * j + 1]][k] == mol->angle_list->partners_[i][2 * j]) continue;
                    if (mol->bond_list->partners_[mol->angle_list->partners_[i][2 * j + 1]][k] == i) continue;
                    mol->dihedral_list->partners_[i][3 * n_dihedrals] = mol->angle_list->partners_[i][2 * j];
                    mol->dihedral_list->partners_[i][3 * n_dihedrals + 1] = mol->angle_list->partners_[i][2 * j + 1];
                    mol->dihedral_list->partners_[i][3 * n_dihedrals + 2] = mol->bond_list->partners_[mol->angle_list->partners_[i][2 * j + 1]][k];
                    cg_site2 = calc_four_body_interaction_hash(mol->cg_site_types[mol->dihedral_list->partners_[i][3 * n_dihedrals]], mol->cg_site_types[mol->dihedral_list->partners_[i][3 * n_dihedrals + 1]], mol->cg_site_types[i], mol->cg_site_types[mol->dihedral_list->partners_[i][3 * n_dihedrals + 2]], topo_data->n_cg_types);
                    topo_data->dihedral_type_activation_flags[cg_site2] = 1;
                    n_dihedrals++;
                }
            }
            mol->dihedral_list->partner_numbers_[i] = n_dihedrals;
        }
        if (automatic_topology_style > 2) {
            int total_dihedrals=0;
	        for (i = 0; i < mol->n_cg_sites; i++) {
	            total_dihedrals += mol->dihedral_list->partner_numbers_[i];
	        }
            printf("Automatically generated dihedral topology; %d dihedrals of %d dihedral types.\n", total_dihedrals/2, calc_n_active_interactions(topo_data->dihedral_type_activation_flags, calc_n_distinct_quadruples(topo_data->n_cg_types)));
        }
    }
}

void setup_excluded_list( TopologyData const* topo_data ) 
{
	// Automatically determine topology to set appropriate
    // bond, angle, and/or dihedral exclusion as appropriate.
    if (topo_data->excluded_style == 0) return;
	printf("Setting up exclusion list for excluded_style %d.\n", topo_data->excluded_style);
    unsigned max_excluded_number = get_max_exclusion_number(topo_data);
            
    // Loop over each CG site and its bonds to copy into exclusion list
    for (unsigned i = 0; i < topo_data->n_cg_sites; i++) {
        for (unsigned j = 0; j < topo_data->bond_list->partner_numbers_[i]; j++) {            
            topo_data->exclusion_list->partners_[i][ topo_data->exclusion_list->partner_numbers_[i] ] = topo_data->bond_list->partners_[i][j];
            topo_data->exclusion_list->partner_numbers_[i]++;
        }
    
    	if ( topo_data->exclusion_list->partner_numbers_[i] > max_excluded_number ) {
    		printf("Warning: Too many excluded interactions (%d) to handle based on max_pair_bonds_per_site (%d)!\n", topo_data->exclusion_list->partner_numbers_[i], max_excluded_number);
    		exit(EXIT_FAILURE);
    	}
    }
    
    // If we also need angles, check the sites bonded to known bonds and
    // if we need dihedrals, check the sites bonded to know angles to add to exclusion list.
    if (topo_data->excluded_style == 2) return;
    
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
                if ( new_excluded_partner(topo_data->exclusion_list->partners_, i, topo_data->exclusion_list->partner_numbers_[i], cg_site2) == true ) {
	                topo_data->exclusion_list->partners_[i][ topo_data->exclusion_list->partner_numbers_[i] ] = cg_site2;
    	            topo_data->exclusion_list->partner_numbers_[i]++;
        		}
        		        
                //if we need dihedrals, loop over all potential dihedrals (cg_site3: bonded to angle)
                if (topo_data->excluded_style == 4) {
                	for(unsigned l = 0; l < topo_data->bond_list->partner_numbers_[cg_site2]; l++) {
                		unsigned cg_site3 = topo_data->bond_list->partners_[cg_site2][l];
                		
                		//add this as a dihedral interaction as long as it does not loop back to a CG site in the existing angle
                		if( (cg_site3 == i) || (cg_site3 == cg_site1) || (cg_site3 == cg_site2) ) continue;
                		if ( new_excluded_partner(topo_data->exclusion_list->partners_, i, topo_data->exclusion_list->partner_numbers_[i], cg_site3) == true ) {
                			topo_data->exclusion_list->partners_[i][ topo_data->exclusion_list->partner_numbers_[i] ] = cg_site3;
               		 		topo_data->exclusion_list->partner_numbers_[i]++;
               		 	}
               		}
            	}
            
        	}
        }
        
    	if ( topo_data->exclusion_list->partner_numbers_[i] > max_excluded_number ) {
    		printf("Warning: Too many excluded interactions (%d) to handle based on max_pair_bonds_per_site, max_angles_per_site, and max_dihedrals_per_site (total %d)!\n", topo_data->exclusion_list->partner_numbers_[i], max_excluded_number);
    		exit(EXIT_FAILURE);
    	}
    }
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
