//
//  topology.h
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#ifndef _topology_h
#define _topology_h

struct CG_MODEL_DATA;

struct TopoList {
    const unsigned n_sites_;		// The size of the patner_numbers_ and the primary index of the partners_ arrays.
									// This is equal to the current number of CG sites.
    const unsigned partners_per_;	// The number of partners each interaction has.
									// This is determined by the interaction type (1 for bond, 2 for angle, 3 for dihedral).
    const unsigned max_partners_;	// The maximum number of partners (set in control.in).
    unsigned* partner_numbers_;		// The number of partners each CG site has for this topological attribute.
    unsigned** partners_;			// For a given CG site, it lists the CG site indices of all partnered particles (for this topological attribute).

	int modified;					// A flag indicating if the pointers for partners_ and partner_numbers_ arrays are shared (1 for yes, 0 for no).
									// This is primarily useful in the LAMMPS fix when these arrays are allocated, freed, and owned by LAMMPS.
									
    inline TopoList() : TopoList(0, 0, 0) {}
    TopoList(unsigned n_sites, unsigned partners_per, unsigned max_partners);
    ~TopoList();
};

// Struct responsible for keeping track of all cg site numbers, types, bonds,
// angles, and dihedrals.

struct TopologyData {
    unsigned n_cg_sites;                    // Total number of CG sites
    unsigned n_cg_types;                    // Total number of CG site types
    char** name;                            // Human-readable names for each CG site type
    int* cg_site_types;                     // List of the CG types for each site
	int* molecule_ids;						// List of the molecule id for each site
    
    unsigned max_pair_bonds_per_site;            // Max pair bonds defined by connection to a single site
    unsigned max_angles_per_site;                // Max angles defined by connection to a single site
    unsigned max_dihedrals_per_site;             // Max dihedrals defined by connection to a single site

	// Pointers to TopoLists containing the sites connected to a given site through a given topological feature.
    TopoList* bond_list;
    TopoList* angle_list;			// This list only has "partners" for the end sites of an angle (not the middle).
    								// The first partner in this list for each angle is the "center" followed by the remaining site.
    TopoList* dihedral_list;		// This list only has "partners" for the end sites of an angle (not the central "bond").
    								// The first two partners in this list for each dihedral are the "center bond sites" followed by the remaining site.
    TopoList* exclusion_list;
    TopoList* molecule_list;
    
    int* bond_type_activation_flags;        // 0 if a given type of pair bonded interaction is active in the model; 1 otherwise
    int* angle_type_activation_flags;       // 0 if a given type of angular bonded interaction is active in the model; 1 otherwise
    int* dihedral_type_activation_flags;    // 0 if a given type of dihedral bonded interaction is active in the model; 1 otherwise
	
	int excluded_style;					// 0 no exclusions; 2 exclude 1-2 bonded; 3 exclude 1-2 and 1-3 bonded; 4 exclude 1-2, 1-3 and 1-4 bonded interactions
	int angle_format;
	int dihedral_format;
	
	// Radius of gyration specific variables
	int n_molecule_groups;
	bool* molecule_groups;
	char** molecule_group_names;
	
	// Density basis specific variables
	int n_density_groups;
	bool* density_groups;
	double* density_weights;
	char** density_group_names;
	int density_excluded_style;
	TopoList* density_exclusion_list;
	
	inline TopologyData () {
		TopologyData(4,12,36);
	}
	
	inline TopologyData (int max_bonds, int max_angles, int max_dihedrals) :
		max_pair_bonds_per_site(max_bonds),
		max_angles_per_site(max_angles),
		max_dihedrals_per_site(max_dihedrals) {
		bond_list = angle_list = dihedral_list  = NULL;
		exclusion_list = density_exclusion_list = NULL;
		molecule_list = NULL;
		cg_site_types = NULL;
  		molecule_ids = NULL;
		bond_type_activation_flags = angle_type_activation_flags = dihedral_type_activation_flags = NULL;
		};
	
	void free_topology_data(void);
};

// Helper functions for managing TopologyData structs.
void free_topology_data(TopologyData *topo_data);

// Primary function for reading a topology file for the CG model
void read_topology_file(TopologyData* topo_data, CG_MODEL_DATA* const cg);

// Initialize topology data structure (used for LAMMPS fix).
void initialize_topology_data(TopologyData* const topo_data);

// Determine appropriate non-bonded exclusions based on bonded topology and exclusion_style setting (used for LAMMPS fix).
void setup_excluded_list(TopologyData const* topo_data,  TopoList* &exclusion_list, const int excluded_style);

// Determine molecule lists based on molecule id information.
void setup_molecule_list(TopologyData const* topo_data, TopoList* &molecule_list, const int n_molecules, const int max_size_per_molecule);
void update_molecule_list(TopologyData const* topo_data, TopoList* &molecule_list);
#endif
