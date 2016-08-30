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

    ControlInputs control_input;

    reset_control_defaults_and_read_control_input(&control_input);
    CG_MODEL_DATA cg(&control_input); // CG model parameters and data; put here to initialize without default constructor
    read_topology_file(&cg.topo_data, &cg);
    read_all_interaction_ranges(&cg);

    MATRIX_DATA* mat = make_matrix(&control_input, &cg);

    if (cg.pair_nonbonded_interactions.get_basis_type() == kBSpline) {
        set_up_force_computers(&cg);
    }

    read_binary_matrix(mat);

    free_AV_fm_setup_temp_variables(&cg, mat);

    mat->finish_fm(mat);

    write_fm_interaction_output_files(&cg, mat);

    delete mat;
    //print cpu time used
    double end_cputime = clock();
    double elapsed_cputime = ((double)(end_cputime - start_cputime)) / CLOCKS_PER_SEC;
    printf("\n%f seconds used.\n", elapsed_cputime);

    return 0;
}
