// mscg_driver.cpp
#include <cstdio>
#include "read_lammps.h"

// "Core" setup functions for running MSCG as a library.
void* mscg_startup_part1(void* mscg_struct);
//void* setup_frame_config(void* void_in, const int n_cg_sites, int * cg_site_types, double* box_half_lengths);
void* setup_topology_and_frame(void* void_in, const int n_cg_sites, const int n_cg_types, char ** type_names, int * cg_site_types, double* box_half_lengths);
void* set_bond_topology(void* void_in, unsigned** bond_partners, unsigned* bond_partner_numbers);
void* set_angle_topology(void* void_in, unsigned** angle_partners, unsigned* angle_partner_numbers);
void* set_dihedral_topology(void* void_in, unsigned** dihedral_partners, unsigned* dihedral_partner_numbers);
void* set_exclusion_topology(void* void_in, unsigned** exclusion_partners, unsigned* exclusion_partner_numbers);
void* mscg_startup_part2(void* mscg_struct);

// This function processes each frame for MSCG.
void* mscg_process_frame(void* void_in, double* const x, double* const f);

// This function finalizes MSCG by solving the matrix and generating outputs.
void* mscg_solve_and_output(void* void_in);


// Other library functions
void* update_frame_config(void* void_in, const int n_cg_sites, int * cg_site_types, double* box_half_lengths);
void* generate_exclusion_topology(void* void_in);
//void* rebuild_neighbor_cells(void* void_in, double* simulation_box_half_lengths);

int main(void) 
{
	void* mscg_struct;

	// Begin setup for mscg
	mscg_struct = mscg_startup_part1(mscg_struct);

	// Setup inputs for mscg
	int n_cg_sites = 1000;
	int n_cg_types = 1;
	int* cg_site_types = new int[n_cg_sites]();
	for (int i = 0; i < n_cg_sites; i++) cg_site_types[i] = 1;
	
	double box_half_lengths[3] = {20.6917, 20.6917, 20.6917};
	char ** type_names = new char*[n_cg_types];
	type_names[0] = new char[24];
	type_names[0] = "MeOH\0";
	double* x = new double[n_cg_sites*3]();
	double* f = new double[n_cg_sites*3]();
	char* filename = "MeOH_example.dat";
	int n_frames = 19;
	int status;
	
	int max_bonds = 4;
	int max_angles = 12;
	int max_dihedrals = 36;
	int max_exclusions = 16;
	unsigned** bond_partners = new unsigned*[n_cg_sites];
	unsigned** angle_partners = new unsigned*[n_cg_sites];
	unsigned** dihedral_partners = new unsigned*[n_cg_sites];
	unsigned** exclusion_partners = new unsigned*[n_cg_sites];
	for(int i = 0; i < n_cg_sites; i++) {
		bond_partners[i] = new unsigned[1 * max_bonds]();
		angle_partners[i] = new unsigned[2 * max_angles]();
		dihedral_partners[i] = new unsigned[3 * max_dihedrals]();
		exclusion_partners[i] = new unsigned[1 * max_exclusions]();
	}
	unsigned* bond_partner_numbers = new unsigned[n_cg_sites]();
	unsigned* angle_partner_numbers = new unsigned[n_cg_sites]();
	unsigned* dihedral_partner_numbers = new unsigned[n_cg_sites]();
	unsigned* exclusion_partner_numbers = new unsigned[n_cg_sites]();
	
	// Setup frame_config
	//printf("call setup_frame_config\n"); fflush(stdout);
	//mscg_struct = setup_frame_config(mscg_struct, n_cg_sites, cg_site_types, box_half_lengths);
	
	// Setup general topology and bond/angle/dihedral information
	printf("call setup_topology_and_frame\n\n"); fflush(stdout);
	mscg_struct = setup_topology_and_frame(mscg_struct, n_cg_sites, n_cg_types, type_names, cg_site_types, box_half_lengths);
		
	// These topologies also need to generate activation flags!
	printf("call set_bond_topology\n"); fflush(stdout);
	mscg_struct = set_bond_topology(mscg_struct, bond_partners, bond_partner_numbers);
	printf("call set_angle_topology\n"); fflush(stdout);
	mscg_struct = set_angle_topology(mscg_struct, angle_partners, angle_partner_numbers);
	printf("call set_dihedral_topology\n"); fflush(stdout);
	mscg_struct = set_dihedral_topology(mscg_struct, dihedral_partners, dihedral_partner_numbers);
	printf("call set_exclusion_topology\n"); fflush(stdout);
	mscg_struct = set_exclusion_topology(mscg_struct, exclusion_partners, exclusion_partner_numbers);
	
	// Finish setup for mscg
	printf("call mscg_startup_part2\n"); fflush(stdout);
	mscg_struct = mscg_startup_part2(mscg_struct);
	
	//Open trajectory file and open first frame
	rLammpsData* lammps_data;
	printf("read_initial_lammps_frame\n"); fflush(stdout);
	lammps_data = rread_initial_lammps_frame(lammps_data, n_cg_sites, cg_site_types, box_half_lengths, x, f, filename);
	
	// Process this frame
	printf("\nCall mscg_process_frame\n"); fflush(stdout);
	mscg_struct = mscg_process_frame(mscg_struct, x, f);

	// Continue until the specified number of frames have been processed.
	for(int i = 1; i < n_frames; i++) {
		printf("\nFRAME: %d\n", i); fflush(stdout);
		printf("call rread_next_lammps_frame\n"); fflush(stdout);
		lammps_data = rread_next_lammps_frame(lammps_data, n_cg_sites, box_half_lengths, x, f);
		printf("call mscg_process_frame\n"); fflush(stdout);
		mscg_struct = mscg_process_frame(mscg_struct, x, f);
	}

	// Finish mscg and output the determined interactions.
	printf("call mscg_solve_and_output\n"); fflush(stdout);
	mscg_solve_and_output(mscg_struct);

	// Close input trajectory
	printf("call finish_lammps_reading\n"); fflush(stdout);
	rfinish_lammps_reading(lammps_data);

	return 0;
}