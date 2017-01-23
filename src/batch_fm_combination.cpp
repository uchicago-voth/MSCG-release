//
//  batch_fm_combination.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include "batch_fm_combination.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "external_matrix_routines.h"
#include "topology.h"
#include "interaction_hashing.h"
#include "interaction_model.h"
#include "matrix.h"
#include "misc.h"

void solve_regularized_accumulation_form_fm_equations(MATRIX_DATA* const mat)
{
    int i, j, k;
    double mscg_residual;

    FILE* lambda_table;
    int n_l, l_mod;
    double* lam, *lam2, tx, tx1;
    double *xnorm_l, *resid_l, *xnorm_l1;
    FILE* lout;
    
    FILE* solution_file;
    double* singular_values;
    int lapack_setup_flag = -1;
    double* lapack_temp_workspace;
    int info_in;
    char schar = 'S';
    double* uu;
    double* vt;

    // Initialize the arrays and other temps needed for SVD
    mscg_residual = fabs(mat->dense_fm_normal_rhs_vector[mat->fm_matrix_columns]);
    singular_values = new double[mat->fm_matrix_columns];
    lapack_temp_workspace = new double[1];
    mat->fm_solution = std::vector<double>(mat->fm_matrix_columns);
    uu = new double[mat->accumulation_matrix_columns * mat->fm_matrix_columns];
    vt = new double[mat->fm_matrix_columns * mat->fm_matrix_columns];
    
    // The SVD routine is first called to determine the size of the needed
    // workspace, then that workspace is allocated, and then the SVD routine
    // is called again to perform the actual singular value decomposition.
    dgesvd_(&schar, &schar, &mat->accumulation_matrix_columns, &mat->fm_matrix_columns, mat->dense_fm_matrix->values, &mat->accumulation_matrix_columns, singular_values, uu, &mat->accumulation_matrix_columns, vt, &mat->fm_matrix_columns, lapack_temp_workspace, &lapack_setup_flag, &info_in);
    lapack_setup_flag = lapack_temp_workspace[0];
    free(lapack_temp_workspace);
    lapack_temp_workspace = new double[lapack_setup_flag];
    dgesvd_(&schar, &schar, &mat->accumulation_matrix_columns, &mat->fm_matrix_columns, mat->dense_fm_matrix->values, &mat->accumulation_matrix_columns, singular_values, uu, &mat->accumulation_matrix_columns, vt, &mat->fm_matrix_columns, lapack_temp_workspace, &lapack_setup_flag, &info_in);

    
    // Use the SVD as output and to calculate the final FM solution;
    // depends on output style choices and regularization style choices.
    
    double s, s1, s2, squared_regularization_parameter;
    double eps;
    double* tmp = new double[mat->fm_matrix_columns];
    
    if (mat->rcond < 0) eps = DBL_EPSILON;
    else eps = mat->rcond;

    if (mat->regularization_style == 0 || mat->regularization_style == 1) {
        
        // Print the singular values.
        solution_file = open_file("sol_info.out", "w");
        fprintf(solution_file, "Singular vector:\n");
        for (i = 0; i < mat->fm_matrix_columns; i++) {
            fprintf(solution_file, "%le\n", singular_values[i]);
        }

        if (mat->regularization_style == 0) {

            // Use the singular values to determine the FM solution
            // without introducing any regularization
            for (j = 0; j < mat->fm_matrix_columns; j++) {
                s = 0.0;
                if (singular_values[j] > eps) {
                    for (i = 0; i < mat->accumulation_matrix_columns; i++) {
                        s += uu[j * mat->accumulation_matrix_columns + i] * mat->dense_fm_normal_rhs_vector[i];
                    }
                    s /= singular_values[j];

                }
                tmp[j] = s;
            }

            for (j = 0; j < mat->fm_matrix_columns; j++) {
                s = 0.0;
                for (i = 0; i < mat->fm_matrix_columns; i++) {
                    s += vt[j * mat->fm_matrix_columns + i] * tmp[i];
                }
                mat->fm_solution[j] = s;
            }
        
        } else {
            
            // Use the singular values to determine the FM solution
            // using scalar Tikhonov regularization 
            squared_regularization_parameter = mat->tikhonov_regularization_param * mat->tikhonov_regularization_param;
            for (j = 0; j < mat->fm_matrix_columns; j++) {
                s = 0.0;
                for (i = 0; i < mat->accumulation_matrix_columns; i++) {
                    s += uu[j * mat->accumulation_matrix_columns + i] * mat->dense_fm_normal_rhs_vector[i];
                }
                s = s * singular_values[j] / (singular_values[j] * singular_values[j] + squared_regularization_parameter);
                tmp[j] = s;
            }

            for (j = 0; j < mat->fm_matrix_columns; j++) {
                s = 0.0;
                for (i = 0; i < mat->fm_matrix_columns; i++) {
                    s += vt[j * mat->fm_matrix_columns + i] * tmp[i];
                }
                mat->fm_solution[j] = s;
            }
        }

    } else {
        
        // Use the singular values to determine the FM solution
        // using many-parameter regularization 
        lambda_table = open_file("lambda.in", "r");
        fscanf(lambda_table, "%d%d", &n_l, &l_mod);
        lout = open_file("l.out", "w");

        lam = new double[n_l];
        lam2 = new double[n_l];
        xnorm_l = new double[n_l];
        resid_l = new double[n_l];
        xnorm_l1 = new double[n_l];

        for (k = 0; k < n_l; k++) {
            fscanf(lambda_table, "%lf", lam + k);
            if (l_mod == 1) lam[k] = pow(10.0, lam[k]);
            lam2[k] = lam[k] * lam[k];
            xnorm_l[k] = 0.0;
            resid_l[k] = 0.0;
            xnorm_l1[k] = 0.0;
        }


        for (j = 0; j < mat->fm_matrix_columns; j++) {
            s = 0.0;
            for (i = 0; i < mat->accumulation_matrix_columns; i++) {
                s += uu[j * mat->accumulation_matrix_columns + i] * mat->dense_fm_normal_rhs_vector[i];
            }
            tx1 = singular_values[j] * singular_values[j];
            for (k = 0; k < n_l; k++) {
                tx = tx1 + lam2[k];
                s2 = s / tx;
                s1 = s2 * singular_values[j];
                xnorm_l[k] += s1 * s1;
                s1 = s2 * lam2[k];
                resid_l[k] += s1 * s1;
                xnorm_l1[k] += s * s * lam2[k] * tx1 / (tx * tx * tx);
            }
        }

        for (i = 0; i < n_l; i++) {
            xnorm_l1[i] *= (-4.0 / lam[i]);
            tx1 = xnorm_l[i] * resid_l[i];
            tx = 2.0 * tx1 / xnorm_l1[i] * (lam2[i] * xnorm_l1[i] * resid_l[i] + 2.0 * lam[i] * tx1 + lam2[i] * lam2[i] * xnorm_l[i] * xnorm_l1[i]) / pow(lam2[i] * xnorm_l[i] * xnorm_l[i] + resid_l[i] * resid_l[i], 1.5);
            fprintf(lout, "%19.14le %19.14le %19.14le %19.14le\n", resid_l[i] + mscg_residual, xnorm_l[i], lam[i], -tx);

        }

        fclose(lout);
        delete [] lam;
        delete [] lam2;
        delete [] xnorm_l;
        delete [] xnorm_l1;
        delete [] resid_l;
    }


    delete [] singular_values;
    delete [] lapack_temp_workspace;
    delete [] mat->dense_fm_normal_rhs_vector;
    delete [] mat->dense_fm_matrix;
    if (mat->regularization_style == 2) exit(EXIT_SUCCESS);
}
