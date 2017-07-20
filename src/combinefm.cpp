//
//  combinefm.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//


#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "batch_fm_combination.h"
#include "control_input.h"
#include "force_computation.h"
#include "fm_output.h"
#include "interaction_model.h"
#include "matrix.h"
#include "misc.h"

int main(int argc, char* argv[])
{
    double start_cputime = clock();

    
    printf("Reading high level control parameters.\n");
    ControlInputs control_input;
    CG_MODEL_DATA cg(&control_input); // CG model parameters and data; put here to initialize without default constructor

    printf("Reading topology file.\n");
    read_topology_file(&cg.topo_data, &cg);

    printf("Reading interaction ranges.\n");
    read_all_interaction_ranges(&cg);

    if (cg.pair_nonbonded_interactions.n_tabulated > 0 ||
        cg.pair_bonded_interactions.n_tabulated > 0 ||
        cg.angular_interactions.n_tabulated > 0 ||
        cg.dihedral_interactions.n_tabulated > 0 ||
		cg.density_interactions.n_tabulated > 0) {
        printf("Reading tabulated reference potentials.\n");
        read_tabulated_interaction_file(&cg, cg.topo_data.n_cg_types);
    } 

    MATRIX_DATA mat(&control_input, &cg);

    set_up_force_computers(&cg);
 
    read_binary_matrix(&mat);

    mat.finish_fm(&mat);

    write_fm_interaction_output_files(&cg, &mat);

    //print cpu time used
    double end_cputime = clock();
    double elapsed_cputime = ((double)(end_cputime - start_cputime)) / CLOCKS_PER_SEC;
    printf("\n%f seconds used.\n", elapsed_cputime);

    return 0;
}
